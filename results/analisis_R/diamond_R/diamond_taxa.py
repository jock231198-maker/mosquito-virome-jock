#!/usr/bin/env python3
"""
diamond_taxa.py — añade a una tabla de DIAMOND (outfmt 6, sin cabecera) la
cabecera y cuatro columnas de taxonomía por hit, sacadas del taxdump:

    grupo    Virus / Artropodo / Metazoo_otro / Hongo / Planta / Eucariota_otro /
             Bacteria / Arquea / otro / sin_taxid / taxid_desconocido
             (la MISMA función que usa arbiter_nr.py, así que cuadra con el veredicto)
    realm, familia, genero   (vacíos si no aplica)

Si un hit trae varios taxids ("a;b") se usa el primero, igual que el árbitro.
Solo biblioteca estándar + viral_db_tools.py del pipeline.

Uso:
  python3 diamond_taxa.py --scripts-dir DIR --taxdump DIR \
      --cols "qseqid sseqid ... stitle" --in diamond.tsv --out salida.tsv.gz
"""
import argparse, gzip, sys

ap = argparse.ArgumentParser()
ap.add_argument("--scripts-dir", required=True, help="carpeta con viral_db_tools.py")
ap.add_argument("--taxdump", required=True)
ap.add_argument("--cols", required=True, help="nombres de columna, separados por espacios")
ap.add_argument("--in", dest="inp", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()

sys.path.insert(0, a.scripts_dir)
from viral_db_tools import Taxonomy, _group  # noqa: E402

cols = a.cols.split()
if "staxids" not in cols:
    sys.exit("ERROR: --cols no tiene 'staxids'")
it = cols.index("staxids")
tax = Taxonomy(a.taxdump)
cache, n, malas = {}, 0, 0
with open(a.inp) as fh, gzip.open(a.out, "wt") as out:
    out.write("\t".join(cols + ["grupo", "realm", "familia", "genero"]) + "\n")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) != len(cols):
            malas += 1
            continue
        tid = (f[it] or "").split(";")[0]
        if tid not in cache:
            g = _group(tax, tid)
            lin = tax.lineage(tid) if g not in ("sin_taxid", "taxid_desconocido") else {}
            cache[tid] = [g, lin.get("realm", ""), lin.get("family", ""), lin.get("genus", "")]
        out.write("\t".join(f + cache[tid]) + "\n")
        n += 1
if malas:
    sys.exit(f"ERROR: {malas} lineas con un numero de columnas distinto de {len(cols)}")
print(f"  {n} hits, {len(cache)} taxids distintos")
