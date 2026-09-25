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
    via_annot = len(p2t)
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
    print(f"proteinas: {n} | via annotation_report: {via_annot} | via nombre: {via_name} | "
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
    a = ap.parse_args()
    {"meta": cmd_meta, "protmap": cmd_protmap, "ntmap": cmd_ntmap, "select": cmd_select,
     "clstr": cmd_clstr, "summary": cmd_summary, "probe": cmd_probe}[a.cmd](a)


if __name__ == "__main__":
    main()
