#!/usr/bin/env python3
"""
integrar_arbitro.py — une el veredicto del arbitro de nr a la tabla final de vOTUs.

No toca TABLA_FINAL_votus.tsv: escribe una tabla nueva con 6 columnas mas.

  clase (tabla_final.py)  +  veredicto (arbiter_nr.py)  ->  clase_final

  rRNA / ADN_descartado            se mantienen
  sin_hit + viral                  viral_solo_nr     (nr ve virus que la base viral no vio)
  sin_hit + ambiguo                ambiguo
  sin_hit + celular_Artropodo      huesped
  sin_hit + celular_<otro>         celular_otro      (bacteria, hongo, planta, ...)
  sin_hit + sin_hit_nr             oscura            (sin hit en ninguna base)
  identificado/divergente/pariente_debil
          + celular_*              <clase> y conflicto_nr = si
          + otro                   <clase>

Uso:
  python3 integrar_arbitro.py \
      $RESULTS_DIR/qc_control/TABLA_FINAL_votus.tsv \
      $RESULTS_DIR/votu_annot/votus_nr_jr46/veredicto_por_votu.tsv \
      $RESULTS_DIR/qc_control/TABLA_FINAL_votus_arb.tsv
"""
import csv
import sys
from collections import Counter

if len(sys.argv) != 4:
    sys.exit(__doc__)
f_tabla, f_arb, f_out = sys.argv[1:]

COLS_ARB = ["veredicto", "nr_nombre", "nr_grupo", "nr_pident", "celular_nombre"]
CONSERVAR = {"rRNA", "ADN_descartado"}
VIRALES = {"identificado", "divergente", "pariente_debil"}


def leer(path):
    with open(path, newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        return r.fieldnames, list(r)


cab, tabla = leer(f_tabla)
cab_arb, arb_filas = leer(f_arb)
falta = [c for c in ["votu"] + COLS_ARB if c not in cab_arb]
if falta:
    sys.exit(f"ERROR: faltan columnas en {f_arb}: {falta}")
arb = {r["votu"]: r for r in arb_filas}


def clase_final(clase, ver):
    if clase in CONSERVAR:
        return clase, ""
    celular = ver.startswith("celular_")
    if clase in VIRALES:
        return clase, ("si" if celular else "")
    # clase == sin_hit (o cualquier otra no prevista)
    if ver == "viral":
        return "viral_solo_nr", ""
    if ver == "ambiguo":
        return "ambiguo", ""
    if ver == "celular_Artropodo":
        return "huesped", ""
    if celular:
        return "celular_otro", ""
    if ver == "sin_hit_nr":
        return "oscura", ""
    return clase, ""


sin_arb = 0
cruce = Counter()
salida = []
for r in tabla:
    a = arb.get(r["votu"])
    if a is None:
        sin_arb += 1
        a = {c: "" for c in COLS_ARB}
    ver = a["veredicto"] or "sin_arbitro"
    cf, conflicto = clase_final(r["clase"], ver)
    cruce[(r["clase"], ver)] += 1
    r2 = dict(r)
    r2["nr_veredicto"] = ver
    for c in COLS_ARB[1:]:
        r2[c] = a[c]
    r2["clase_final"] = cf
    r2["conflicto_nr"] = conflicto
    salida.append(r2)

cab_out = cab + ["nr_veredicto"] + COLS_ARB[1:] + ["clase_final", "conflicto_nr"]
with open(f_out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cab_out, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(salida)

# ---------------------------------------------------------------- informe
print(f"vOTUs en tabla final: {len(tabla)}   con veredicto: {len(tabla) - sin_arb}   sin veredicto: {sin_arb}")

clases = sorted({c for c, _ in cruce}, key=lambda c: -sum(n for (k, _), n in cruce.items() if k == c))
vers = sorted({v for _, v in cruce}, key=lambda v: -sum(n for (_, k), n in cruce.items() if k == v))
print("\n== clase (tabla_final) x veredicto (nr) ==")
print(f"  {'':16s}" + "".join(f"{v[:13]:>14s}" for v in vers))
for c in clases:
    print(f"  {c:16s}" + "".join(f"{cruce.get((c, v), 0) or '.':>14}" for v in vers))

print("\n== clase_final ==")
for c, n in Counter(r["clase_final"] for r in salida).most_common():
    print(f"  {c:16s} {n:4d}")

conf = [r for r in salida if r["conflicto_nr"] == "si"]
print(f"\n== conflictos: clase viral pero nr dice celular ({len(conf)}) ==")
for r in sorted(conf, key=lambda r: r["clase"]):
    print(f"  {r['clase']:15s} {r['especie'][:28]:28s} -> {r['nr_veredicto']:20s} {r['celular_nombre'][:30]}")

print(f"\n  escrito: {f_out}")
