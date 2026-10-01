#!/usr/bin/env python3
"""arbiter_nr.py — DIAMOND vs nr como arbitro: ¿el mejor hit GLOBAL de cada vOTU es viral o celular?
Uso: python3 arbiter_nr.py <taxdump_dir> <nr.tsv> <annot_por_votu.tsv> <salida.tsv> [margen=1.1]"""
import sys, os, csv
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from viral_db_tools import Taxonomy

GROUPS = [("10239", "Virus"), ("6656", "Artropodo"), ("33208", "Metazoo_otro"), ("4751", "Hongo"),
          ("33090", "Planta"), ("2759", "Eucariota_otro"), ("2", "Bacteria"), ("2157", "Arquea")]
NR = ["qseqid", "sseqid", "pident", "length", "qlen", "evalue", "bitscore", "staxids", "sscinames", "qcovhsp", "stitle"]

def group(tax, tid):
    tid = (tid or "").split(";")[0]
    if not tid or tid == "0": return "sin_taxid"
    t, n, anc = tax.current(tid), 0, set()
    if t not in tax.parent: return "taxid_desconocido"
    while t in tax.parent and n < 80:
        anc.add(t)
        if tax.parent[t] == t: break
        t = tax.parent[t]; n += 1
    for g, name in GROUPS:
        if g in anc: return name
    return "otro"

taxdir, nrf, annf, outf = sys.argv[1:5]
margin = float(sys.argv[5]) if len(sys.argv) > 5 else 1.1
tax = Taxonomy(taxdir)
best, bestv, bestc = {}, {}, {}
for line in open(nrf):
    f = line.rstrip("\n").split("\t"); f += [""] * (len(NR) - len(f)); d = dict(zip(NR, f))
    try: b = float(d["bitscore"])
    except ValueError: continue
    d["grupo"] = group(tax, d["staxids"]); q = d["qseqid"]
    if q not in best or b > float(best[q]["bitscore"]): best[q] = d
    tgt = bestv if d["grupo"] == "Virus" else bestc
    if d["grupo"] not in ("sin_taxid", "taxid_desconocido") and (q not in tgt or b > float(tgt[q]["bitscore"])): tgt[q] = d

ann = list(csv.DictReader(open(annf), delimiter="\t"))
cols = ["votu", "longitud", "categoria_viral", "familia_viral", "dm_titulo_viral", "nr_acc", "nr_pident", "nr_qcov",
        "nr_evalue", "nr_bits", "nr_taxid", "nr_nombre", "nr_grupo", "nr_titulo", "bits_mejor_viral",
        "bits_mejor_celular", "celular_nombre", "veredicto"]
ver, cross, famv = Counter(), Counter(), defaultdict(Counter)
with open(outf, "w") as out:
    out.write("\t".join(cols) + "\n")
    for r in ann:
        q = r["votu"]; b = best.get(q); bv = bestv.get(q); bc = bestc.get(q)
        vb = float(bv["bitscore"]) if bv else 0.0; cb = float(bc["bitscore"]) if bc else 0.0
        if not b: v = "sin_hit_nr"
        elif vb and (not cb or vb >= cb): v = "viral" if (not cb or vb >= margin * cb) else "ambiguo"
        elif vb and vb >= cb / margin: v = "ambiguo"
        else: v = "celular_" + (bc["grupo"] if bc else b["grupo"])
        row = [q, r["longitud"], r["categoria"], r["familia"], r["dm_titulo"]]
        row += [b[k] for k in ("sseqid", "pident", "qcovhsp", "evalue", "bitscore", "staxids", "sscinames", "grupo", "stitle")] if b else [""] * 9
        row += [f"{vb:.0f}" if vb else "", f"{cb:.0f}" if cb else "", bc["sscinames"] if bc else "", v]
        out.write("\t".join(str(x) for x in row) + "\n")
        ver[v] += 1; cross[(r["categoria"], v)] += 1; famv[r["familia"] or "(sin familia)"][v] += 1

print(f"vOTUs: {len(ann)} | con hit en nr: {sum(1 for r in ann if r['votu'] in best)}")
print(f"(viral si su bitscore >= {margin} x el mejor celular; 'ambiguo' si estan a menos de eso)")
print("\n## veredicto nr")
for k, n in ver.most_common(): print(f"{n:>6}  {k}")
vs = [k for k, _ in ver.most_common()]
print("\n## categoria (base viral) x veredicto (nr)")
print(f"{'':22s}" + "".join(f"{v[:14]:>15s}" for v in vs))
for c in ["conocido_nt", "pariente_nt", "divergente_solo_aa", "sin_hit_viral"]:
    print(f"{c:22s}" + "".join(f"{cross[(c, v)]:>15d}" for v in vs))
print("\n## familia (base viral) x veredicto (nr), familias con >=3 vOTUs")
for fam, cn in sorted(famv.items(), key=lambda x: -sum(x[1].values())):
    if sum(cn.values()) >= 3: print(f"{fam[:24]:24s} " + ", ".join(f"{k}={n}" for k, n in cn.most_common()))
