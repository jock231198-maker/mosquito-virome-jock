#!/usr/bin/env python3
"""
viral_db_tools.py — ayudantes de build_viral_db.sh (solo biblioteca estándar).

Subcomandos
  meta     data_report.jsonl (NCBI Datasets, virus) + taxdump -> TSV de metadatos
  protmap  annotation_report.jsonl + protein.faa + metadatos -> prot2taxid (DIAMOND)
  ntmap    FASTA + metadatos -> "accession taxid" (makeblastdb -taxid_map)
  select   IDs de la referencia de mapeo (completos o RefSeq)
  clstr    .clstr de CD-HIT -> TSV representante / miembro / identidad
  summary  conteos por fuente, completitud, familia
  hostscan metadatos por ramas del arbol viral -> accesiones con huesped bajo un taxon
  annot    mejor hit blastn/DIAMOND por query + linaje -> tabla por vOTU
  probe    muestra la estructura REAL de los ficheros de Datasets (antes de la descarga grande)

Por qué el JSON se RECORRE en vez de leerse por claves fijas: el esquema de
Datasets cambia entre versiones (camelCase, anidamiento). Recorrer el árbol y
cruzar accesiones contra conjuntos ya conocidos no depende de los nombres exactos
de los campos.
"""
import argparse, gzip, json, re, sys, os
from collections import Counter, defaultdict

ACC_RE = re.compile(r'^[A-Z]{1,3}_?[0-9]{5,}\.[0-9]+$')
RANKS = ["realm", "kingdom", "phylum", "class", "order", "family", "genus", "species"]


def opn(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def die(msg):
    sys.stderr.write(f"[ERROR] {msg}\n"); sys.exit(1)


def fasta_ids(path):
    with opn(path) as fh:
        for line in fh:
            if line.startswith(">"):
                yield line[1:].split(None, 1)[0], line[1:].rstrip("\n")


def walk_strings(obj):
    if isinstance(obj, dict):
        for v in obj.values():
            yield from walk_strings(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from walk_strings(v)
    elif isinstance(obj, str):
        yield obj


def dig(d, *paths, default=""):
    """Primer valor no vacío entre varias rutas candidatas ('a.b.c')."""
    for p in paths:
        cur = d
        for k in p.split("."):
            if isinstance(cur, dict) and k in cur:
                cur = cur[k]
            else:
                cur = None; break
        if cur not in (None, "", [], {}):
            return cur
    return default


# --------------------------------------------------------------------------- taxonomy
class Taxonomy:
    def __init__(self, taxdir):
        self.parent, self.rank, self.name, self.merged = {}, {}, {}, {}
        for f in ("nodes.dmp", "names.dmp"):
            if not os.path.exists(os.path.join(taxdir, f)):
                die(f"falta {f} en {taxdir}")
        with open(os.path.join(taxdir, "nodes.dmp")) as fh:
            for line in fh:
                p = line.split("\t|\t")
                self.parent[p[0]] = p[1]; self.rank[p[0]] = p[2]
        with open(os.path.join(taxdir, "names.dmp")) as fh:
            for line in fh:
                p = line.rstrip("\t|\n").split("\t|\t")
                if p[3] == "scientific name":
                    self.name[p[0]] = p[1]
        mf = os.path.join(taxdir, "merged.dmp")
        if os.path.exists(mf):
            with open(mf) as fh:
                for line in fh:
                    p = line.rstrip("\t|\n").split("\t|\t")
                    self.merged[p[0]] = p[1]
        self._cache = {}

    def current(self, tid):
        tid = str(tid)
        return self.merged.get(tid, tid)

    def lineage(self, tid):
        tid = self.current(tid)
        if tid in self._cache:
            return self._cache[tid]
        out = {r: "" for r in RANKS}
        seen, t = 0, tid
        while t in self.parent and seen < 60:
            r = self.rank.get(t, "")
            if r in out and not out[r]:
                out[r] = self.name.get(t, "")
            if self.parent[t] == t:
                break
            t = self.parent[t]; seen += 1
        out["_known"] = tid in self.parent
        self._cache[tid] = out
        return out


# --------------------------------------------------------------------------- meta
META_COLS = ["accession", "taxid", "virus_name", "source_db", "completeness", "length",
             "segment", "host_taxid", "host_name", "collection_date", "geo_location",
             "set"] + RANKS


def cmd_meta(a):
    tax = Taxonomy(a.taxdump)
    seen, n_in, n_dup, n_notax, n_unknown = set(), 0, 0, 0, 0
    out = open(a.out, "w")
    out.write("\t".join(META_COLS) + "\n")
    for spec in a.reports:
        label, path = spec.split("=", 1)
        with opn(path) as fh:
            for line in fh:
                if not line.strip():
                    continue
                n_in += 1
                r = json.loads(line)
                acc = dig(r, "accession")
                if not acc:
                    die(f"registro sin 'accession' en {path}. Claves: {sorted(r)}")
                if acc in seen:
                    n_dup += 1; continue
                seen.add(acc)
                tid = str(dig(r, "virus.taxId", "virus.taxid", "taxId", "taxid"))
                if not tid:
                    n_notax += 1
                lin = tax.lineage(tid) if tid else {k: "" for k in RANKS}
                if tid and not lin.get("_known"):
                    n_unknown += 1
                row = [acc, tax.current(tid) if tid else "",
                       dig(r, "virus.organismName", "virus.sciName", "organismName"),
                       dig(r, "sourceDatabase", "sourcedb"),
                       dig(r, "completeness"),
                       str(dig(r, "length")),
                       dig(r, "segment"),
                       str(dig(r, "host.taxId", "host.taxid")),
                       dig(r, "host.organismName", "host.sciName"),
                       dig(r, "isolate.collectionDate", "collectionDate"),
                       dig(r, "location.geographicLocation", "geoLocation", "location.geographicRegion"),
                       label] + [lin[k] for k in RANKS]
                out.write("\t".join(str(x).replace("\t", " ") for x in row) + "\n")
    out.close()
    kept = len(seen)
    print(f"registros leidos: {n_in} | unicos: {kept} | duplicados por accesion: {n_dup}")
    print(f"sin taxid: {n_notax} | taxid ausente del taxdump: {n_unknown}")
    if kept == 0:
        die("metadatos vacios")
    if n_notax > 0.01 * kept:
        die("mas del 1% sin taxid: cambio el esquema del reporte. Corre 'probe'.")


def read_meta(path):
    m = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            d = dict(zip(hdr, line.rstrip("\n").split("\t")))
            m[d["accession"]] = d
    return m


# --------------------------------------------------------------------------- protmap
def cmd_protmap(a):
    meta = read_meta(a.meta)
    prot_ids, prot_org = set(), {}
    for faa in a.faa:
        for pid, hdr in fasta_ids(faa):
            prot_ids.add(pid)
            # "[organism=X]" (formato actual de Datasets) o "[X]" (clasico)
            mm = re.search(r"\[organism=([^\]]+)\]", hdr) or re.search(r"\[([^\[\]=]+)\]\s*$", hdr)
            if mm:
                prot_org[pid] = mm.group(1)
    p2t, conflicts = {}, 0
    # 1) annotation_report: toda accesion de proteina que aparece en el mismo
    #    registro que un nucleotido conocido hereda el taxid de ese nucleotido
    for rep in a.annot:
        with opn(rep) as fh:
            for line in fh:
                if not line.strip():
                    continue
                strs = [s for s in walk_strings(json.loads(line)) if ACC_RE.match(s)]
                nucs = {s for s in strs if s in meta}
                tids = {meta[n]["taxid"] for n in nucs if meta[n]["taxid"]}
                if len(tids) != 1:
                    continue
                t = tids.pop()
                for s in strs:
                    if s in prot_ids and s not in nucs:
                        if s in p2t and p2t[s] != t:
                            conflicts += 1
                        p2t.setdefault(s, t)
    # 1b) CDS sin accesion propia: ID "NC_139268.1_123-456" -> taxid del nucleotido
    via_nuc = 0
    for pid in prot_ids:
        if pid in p2t: continue
        mm = re.match(r"^([A-Z]+_?[0-9]+\.[0-9]+)_", pid)
        if mm and mm.group(1) in meta and meta[mm.group(1)]["taxid"]:
            p2t[pid] = meta[mm.group(1)]["taxid"]; via_nuc += 1
    via_annot = len(p2t) - via_nuc
    # 2) respaldo: nombre del organismo entre corchetes en la cabecera -> taxid (solo si es unico)
    name2tid = defaultdict(set)
    for d in meta.values():
        if d["virus_name"] and d["taxid"]:
            name2tid[d["virus_name"]].add(d["taxid"])
    via_name = 0
    for pid in prot_ids:
        if pid in p2t:
            continue
        org = prot_org.get(pid)
        if org and len(name2tid.get(org, ())) == 1:
            p2t[pid] = next(iter(name2tid[org])); via_name += 1
    with gzip.open(a.out, "wt") as out:
        out.write("accession.version\ttaxid\n")
        for pid in sorted(p2t):
            out.write(f"{pid}\t{p2t[pid]}\n")
    missing = sorted(prot_ids - set(p2t))
    with open(a.missing, "w") as fh:
        fh.write("\n".join(missing) + ("\n" if missing else ""))
    n = len(prot_ids)
    pct = 100 * len(p2t) / n if n else 0
    print(f"proteinas: {n} | via annotation_report: {via_annot} | via CDS-nucleotido: {via_nuc} | via nombre: {via_name} | "
          f"sin taxid: {len(missing)} ({100-pct:.2f}%) | conflictos: {conflicts}")
    if n and pct < a.min_pct:
        die(f"solo {pct:.1f}% de proteinas con taxid (< {a.min_pct}%). Corre 'probe'.")


# --------------------------------------------------------------------------- ntmap / select
def cmd_ntmap(a):
    meta = read_meta(a.meta)
    n = ok = 0
    with open(a.out, "w") as out:
        for sid, _ in fasta_ids(a.fasta):
            n += 1
            t = meta.get(sid, {}).get("taxid")
            if t:
                out.write(f"{sid} {t}\n"); ok += 1
    print(f"secuencias: {n} | con taxid: {ok} | sin taxid: {n-ok}")
    if n and ok / n < a.min_frac:
        die("demasiadas secuencias sin taxid")


def cmd_select(a):
    meta = read_meta(a.meta)
    n = 0
    with open(a.out, "w") as out:
        for sid, _ in fasta_ids(a.fasta):
            d = meta.get(sid)
            if not d:
                continue
            is_refseq = d["source_db"].lower() == "refseq" or sid[:3] in ("NC_", "AC_")
            if d["completeness"].upper() == "COMPLETE" or is_refseq:
                out.write(sid + "\n"); n += 1
    print(f"seleccionadas para la referencia de mapeo: {n}")


# --------------------------------------------------------------------------- clstr
def cmd_clstr(a):
    rows, members, rep = [], [], None
    def flush():
        for m, ident in members:
            rows.append((rep, m, ident))
    with open(a.clstr) as fh:
        for line in fh:
            if line.startswith(">Cluster"):
                if rep: flush()
                members, rep = [], None
                continue
            mm = re.search(r">(\S+?)\.\.\.\s*(\*|at\s+\S+?/?([\d.]+)%)", line)
            if not mm:
                continue
            sid = mm.group(1)
            if mm.group(2) == "*":
                rep = sid; members.append((sid, "100.00"))
            else:
                members.append((sid, mm.group(3)))
    if rep: flush()
    with open(a.out, "w") as out:
        out.write("representative\tmember\tidentity_pct\n")
        for r in rows:
            out.write("\t".join(r) + "\n")
    print(f"clusters: {len({r[0] for r in rows})} | secuencias: {len(rows)}")


# --------------------------------------------------------------------------- summary
def cmd_summary(a):
    meta = read_meta(a.meta)
    ids = {sid for sid, _ in fasta_ids(a.fasta)} if a.fasta else set(meta)
    rows = [meta[i] for i in ids if i in meta]
    def show(title, cnt, top=None):
        print(f"\n## {title}")
        for k, v in cnt.most_common(top):
            print(f"{v:>9}  {k or '(vacio)'}")
    print(f"secuencias con metadatos: {len(rows)}")
    show("source_db", Counter(r["source_db"] for r in rows))
    show("completeness", Counter(r["completeness"] for r in rows))
    show("set", Counter(r["set"] for r in rows))
    show("familia (top 40)", Counter(r["family"] for r in rows), 40)
    show("huesped (top 20)", Counter(r["host_name"] for r in rows), 20)


# --------------------------------------------------------------------------- hostscan
# Datasets filtra --host por el taxon EXACTO (probado 26-sep: --host 6656 solo da
# registros cuyo huesped es literalmente "Arthropoda"). Aqui se recorre el arbol
# viral del taxdump, se piden SOLO metadatos (datasets summary) por ramas y se
# filtra el huesped por su LINAJE. Las ramas enormes y sin interes (SARS-CoV-2,
# VIH, gripe...) se saltan con --skip para no bajar millones de registros.
def _host_lineage_ids(r, tax):
    h = r.get("host") or {}
    ids = {str(x.get("tax_id", x.get("taxId", ""))) for x in (h.get("lineage") or [])}
    tid = str(h.get("tax_id", h.get("taxId", "")) or "")
    if tid:
        ids.add(tid)
        if len(ids) == 1:                      # sin linaje en el registro: taxdump
            t, n = tax.current(tid), 0
            while t in tax.parent and n < 80:
                ids.add(t)
                if tax.parent[t] == t: break
                t = tax.parent[t]; n += 1
    return ids, (h.get("organism_name") or h.get("organismName") or "")


def cmd_hostscan(a):
    import subprocess, time
    tax = Taxonomy(a.taxdump)
    children = defaultdict(list)
    for t, par in tax.parent.items():
        if t != par:
            children[par].append(t)
    root = tax.current(str(a.root))
    skip = {}
    for x in a.skip.split(","):
        x = x.strip()
        if not x: continue
        cur = tax.current(x)
        if cur not in tax.parent:
            print(f"[aviso] skip {x}: no esta en el taxdump, se ignora"); continue
        skip[cur] = tax.name.get(cur, "?")
    anc = set()                                 # nodos por encima de algun skip
    for s_ in skip:
        t = tax.parent.get(s_)
        while t and t != "1":
            anc.add(t)
            if t == root: break
            t = tax.parent.get(t)
    queries, stack = [], [root]
    while stack:
        t = stack.pop()
        if t in skip: continue
        if t in anc: stack.extend(children[t])
        else: queries.append(t)
    queries.sort(key=int)
    os.makedirs(a.outdir, exist_ok=True)
    print(f"ramas a consultar: {len(queries)} | saltadas: " +
          ", ".join(f"{k} ({v})" for k, v in skip.items()))
    with open(os.path.join(a.outdir, "plan.tsv"), "w") as fh:
        for q in queries: fh.write(f"{q}\t{tax.name.get(q, '?')}\n")
    host = str(a.host)
    base = ["datasets", "summary", "virus", "genome", "taxon"]
    extra = ["--as-json-lines"] + (["--api-key", a.api_key] if a.api_key else [])
    failed, t0 = [], time.time()
    for i, q in enumerate(queries, 1):
        done = os.path.join(a.outdir, f"q_{q}.done")
        if os.path.exists(done): continue
        for attempt in range(1, a.retries + 1):
            scanned, kept, hosts = 0, [], Counter()
            p = subprocess.Popen(base + [q] + extra, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE, text=True)
            for line in p.stdout:
                if not line.strip(): continue
                scanned += 1
                r = json.loads(line)
                ids, hname = _host_lineage_ids(r, tax)
                if host in ids:
                    kept.append(r.get("accession", "")); hosts[hname] += 1
            err = p.stderr.read(); rc = p.wait()
            empty = rc != 0 and scanned == 0 and re.search(r"no .*(found|match)|0 genomes", err, re.I)
            if rc == 0 or empty:
                with open(os.path.join(a.outdir, f"q_{q}.acc"), "w") as fh:
                    fh.write("".join(x + "\n" for x in kept if x))
                with open(os.path.join(a.outdir, f"q_{q}.hosts"), "w") as fh:
                    fh.write("".join(f"{v}\t{k}\n" for k, v in hosts.items()))
                open(done, "w").write(f"{q}\t{tax.name.get(q,'?')}\t{scanned}\t{len(kept)}\n")
                if scanned:
                    print(f"[{i}/{len(queries)}] {q} {tax.name.get(q,'?')[:40]:40s} "
                          f"leidos {scanned:>9} | artropodo {len(kept):>7} | {time.time()-t0:7.0f}s", flush=True)
                break
            print(f"[{i}/{len(queries)}] {q} intento {attempt} fallo (rc={rc}): {err.strip()[:200]}", flush=True)
            time.sleep(30 * attempt)
        else:
            failed.append(q)
    # consolidar
    accs, hosts, tot_s = set(), Counter(), 0
    for q in queries:
        d = os.path.join(a.outdir, f"q_{q}.done")
        if not os.path.exists(d): continue
        tot_s += int(open(d).read().split("\t")[2])
        accs.update(l.strip() for l in open(os.path.join(a.outdir, f"q_{q}.acc")) if l.strip())
        for l in open(os.path.join(a.outdir, f"q_{q}.hosts")):
            v, k = l.rstrip("\n").split("\t", 1); hosts[k] += int(v)
    with open(os.path.join(a.outdir, "accessions.txt"), "w") as fh:
        fh.write("".join(x + "\n" for x in sorted(accs)))
    with open(os.path.join(a.outdir, "hosts.tsv"), "w") as fh:
        for k, v in hosts.most_common(): fh.write(f"{v}\t{k}\n")
    print(f"\nregistros leidos: {tot_s} | con huesped bajo {host}: {len(accs)}")
    print("huespedes mas frecuentes:")
    for k, v in hosts.most_common(10): print(f"  {v:>8}  {k}")
    if failed:
        die(f"{len(failed)} ramas fallaron ({','.join(failed[:10])}...). Relanza: retoma donde se quedo.")


# --------------------------------------------------------------------------- annot
# Mejor hit por query (bitscore) de blastn y de DIAMOND + linaje del taxdump.
def _best_hits(path, cols):
    best = {}
    if not path or not os.path.exists(path):
        return best
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(cols):
                f += [""] * (len(cols) - len(f))
            d = dict(zip(cols, f))
            q = d["qseqid"]
            try:
                b = float(d["bitscore"])
            except ValueError:
                continue
            if q not in best or b > float(best[q]["bitscore"]):
                best[q] = d
    return best


def cmd_annot(a):
    tax = Taxonomy(a.taxdump)
    BN = ["qseqid", "sacc", "pident", "length", "qlen", "slen", "evalue", "bitscore", "staxids", "sscinames", "qcovs"]
    DM = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore", "staxids", "sscinames", "qcovhsp", "stitle"]
    bn = _best_hits(a.blastn, BN)
    dm = _best_hits(a.diamond, DM)
    lens = {}
    with open(a.lengths) as fh:
        for line in fh:
            sid, ln = line.rstrip("\n").split("\t")[:2]
            lens[sid] = ln
    def lin(t):
        t = (t or "").split(";")[0]
        return tax.lineage(t) if t and t != "0" else {k: "" for k in RANKS}
    cols = ["votu", "longitud", "categoria",
            "bn_acc", "bn_pident", "bn_qcov", "bn_evalue", "bn_taxid", "bn_nombre", "bn_familia", "bn_genero",
            "dm_acc", "dm_pident", "dm_qcov_hsp", "dm_evalue", "dm_taxid", "dm_nombre", "dm_familia", "dm_genero",
            "dm_titulo", "familia", "orden", "realm"]
    cat = Counter(); fam = Counter()
    with open(a.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for q in lens:
            b, d = bn.get(q), dm.get(q)
            lb = lin(b["staxids"]) if b else {k: "" for k in RANKS}
            ld = lin(d["staxids"]) if d else {k: "" for k in RANKS}
            if b and float(b["pident"]) >= a.id_known and float(b["qcovs"] or 0) >= a.cov_known:
                c = "conocido_nt"
            elif b:
                c = "pariente_nt"
            elif d:
                c = "divergente_solo_aa"
            else:
                c = "sin_hit_viral"
            L = lb if b else ld
            row = [q, lens[q], c]
            row += [b[k] for k in ("sacc", "pident", "qcovs", "evalue", "staxids", "sscinames")] if b else [""] * 6
            row += [lb["family"], lb["genus"]]
            row += [d[k] for k in ("sseqid", "pident", "qcovhsp", "evalue", "staxids", "sscinames")] if d else [""] * 6
            row += [ld["family"], ld["genus"], d["stitle"] if d else ""]
            row += [L["family"], L["order"], L["realm"]]
            out.write("\t".join(str(x) for x in row) + "\n")
            cat[c] += 1; fam[L["family"] or ("(sin familia asignada)" if c != "sin_hit_viral" else "(sin hit)")] += 1
    print(f"vOTUs: {len(lens)} | con hit blastn: {len(set(bn) & set(lens))} | con hit DIAMOND: {len(set(dm) & set(lens))}")
    print("\n## categoria")
    for k, v in cat.most_common(): print(f"{v:>6}  {k}")
    print("\n## familia (mejor hit; blastn si lo hay, si no DIAMOND)")
    for k, v in fam.most_common(25): print(f"{v:>6}  {k}")


# --------------------------------------------------------------------------- probe
def cmd_probe(a):
    d = a.dir
    def find(name):
        for root, _, files in os.walk(d):
            if name in files:
                return os.path.join(root, name)
    for f in ("genomic.fna", "protein.faa", "data_report.jsonl", "annotation_report.jsonl"):
        p = find(f)
        print(f"{f:26s} {'OK  ' + p if p else 'NO ESTA'}")
    p = find("data_report.jsonl")
    if p:
        r = json.loads(open(p).readline())
        print("\n-- data_report, primer registro (recortado) --")
        print(json.dumps(r, indent=1)[:2500])
        print("\nvirus.taxId ->", dig(r, "virus.taxId", "virus.taxid", "taxId") or "NO ENCONTRADO")
        print("completeness ->", dig(r, "completeness") or "NO ENCONTRADO")
        print("sourceDatabase ->", dig(r, "sourceDatabase") or "NO ENCONTRADO")
    p = find("protein.faa")
    if p:
        print("\n-- protein.faa, primeras 3 cabeceras --")
        print("\n".join(h for _, h in list(fasta_ids(p))[:3]))
    p = find("annotation_report.jsonl")
    if p:
        r = json.loads(open(p).readline())
        print("\n-- annotation_report, accesiones halladas en el primer registro --")
        print(sorted({s for s in walk_strings(r) if ACC_RE.match(s)})[:20])


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("meta"); p.add_argument("--taxdump", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("reports", nargs="+", help="etiqueta=ruta/data_report.jsonl (el orden = prioridad)")
    p = sp.add_parser("protmap"); p.add_argument("--meta", required=True)
    p.add_argument("--faa", nargs="+", required=True); p.add_argument("--annot", nargs="*", default=[])
    p.add_argument("--out", required=True); p.add_argument("--missing", required=True)
    p.add_argument("--min-pct", type=float, default=90.0)
    p = sp.add_parser("ntmap"); p.add_argument("--meta", required=True)
    p.add_argument("--fasta", required=True); p.add_argument("--out", required=True)
    p.add_argument("--min-frac", type=float, default=0.99)
    p = sp.add_parser("select"); p.add_argument("--meta", required=True)
    p.add_argument("--fasta", required=True); p.add_argument("--out", required=True)
    p = sp.add_parser("clstr"); p.add_argument("clstr"); p.add_argument("--out", required=True)
    p = sp.add_parser("summary"); p.add_argument("--meta", required=True); p.add_argument("--fasta")
    p = sp.add_parser("probe"); p.add_argument("dir")
    p = sp.add_parser("hostscan"); p.add_argument("--taxdump", required=True)
    p.add_argument("--outdir", required=True); p.add_argument("--host", default="6656")
    p.add_argument("--skip", default=""); p.add_argument("--api-key", default="")
    p.add_argument("--retries", type=int, default=4)
    p.add_argument("--root", default="10239", help="raiz del recorrido (10239 = Viruses; otra para pruebas)")
    p = sp.add_parser("annot"); p.add_argument("--taxdump", required=True)
    p.add_argument("--lengths", required=True); p.add_argument("--blastn", default="")
    p.add_argument("--diamond", default=""); p.add_argument("--out", required=True)
    p.add_argument("--id-known", type=float, default=95.0); p.add_argument("--cov-known", type=float, default=80.0)
    a = ap.parse_args()
    {"meta": cmd_meta, "protmap": cmd_protmap, "ntmap": cmd_ntmap, "select": cmd_select,
     "clstr": cmd_clstr, "summary": cmd_summary, "probe": cmd_probe,
     "hostscan": cmd_hostscan, "annot": cmd_annot}[a.cmd](a)


if __name__ == "__main__":
    main()
