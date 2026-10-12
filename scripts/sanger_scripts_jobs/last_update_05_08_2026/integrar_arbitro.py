#!/usr/bin/env python3
"""
integrar_arbitro.py — une a la tabla final el veredicto del arbitro de nr y,
opcionalmente, el cribado contra el genoma del huesped (cribado_genoma.sh).

No toca TABLA_FINAL_votus.tsv: escribe una tabla nueva con columnas extra y
`clase_final` + `es_virus`, que son las que usan las figuras.

Reglas (en este orden; EN_GENOMA = frac_en_genoma >= 0.5):

  rRNA / ADN_descartado                         se mantienen
  identificado / divergente / pariente_debil:
      EN_GENOMA                                 EVE_genoma     (en AaegL5 a >= 98%)
      nr celular_Artropodo, clase divergente    huesped        (homologia espuria: fagos, gigantes, retro)
      nr celular_Artropodo, resto               EVE_candidato  (no esta en la referencia; small RNA lo decide)
      nr otro celular                           celular_otro   (bacteria, hongo, planta, metazoo)
      pariente_debil sin proteina (sin hit DIAMOND, nr no viral)
        y blastn con qcov < 30 %                evidencia_debil (solo un tramo corto de nt)
      resto                                     se mantiene
  sin_hit:
      nr viral                                  viral_solo_nr  (EVE_genoma si EN_GENOMA)
      nr ambiguo                                ambiguo
      nr celular_Artropodo                      huesped
      nr otro celular                           celular_otro
      nr sin_hit_nr                             oscura         (huesped si EN_GENOMA)

es_virus = clase_final in {identificado, divergente, pariente_debil, viral_solo_nr}

nombre_virus (solo vOTUs virales): en "identificado", la especie de blastn; en el resto,
"afin a <X>" con X el mejor hit de proteina en nr si es viral (la evidencia mas fuerte),
y si no, el de blastn. Un hit de blastn sobre el 5 % de la vOTU no da nombre.
grupo_virus: "fago" si cualquiera de los dos nombres es de fago; "artropodo" si no.

Uso:
  python3 integrar_arbitro.py TABLA_FINAL_votus.tsv veredicto_por_votu.tsv SALIDA.tsv [votus_en_genoma.tsv]
"""
import csv
import re
import sys
from collections import Counter

if len(sys.argv) not in (4, 5):
    sys.exit(__doc__)
f_tabla, f_arb, f_out = sys.argv[1:4]
f_gen = sys.argv[4] if len(sys.argv) == 5 else None

EN_GENOMA = 0.5
QCOV_MIN = 30.0          # cobertura minima de blastn para nombrar sin apoyo de proteina
COLS_ARB = ["veredicto", "nr_nombre", "nr_grupo", "nr_pident", "celular_nombre"]
CONSERVAR = {"rRNA", "ADN_descartado"}
VIRALES = {"identificado", "divergente", "pariente_debil"}
ES_VIRUS = VIRALES | {"viral_solo_nr"}


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

gen = {}
if f_gen:
    _, gen_filas = leer(f_gen)
    gen = {r["votu"]: r for r in gen_filas}


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def clase_final(r, ver, fg):
    clase = r["clase"]
    if clase in CONSERVAR:
        return clase, ""
    cel = ver.startswith("celular_")
    en_gen = fg is not None and fg >= EN_GENOMA
    if clase in VIRALES:
        if en_gen:
            return "EVE_genoma", "si"
        if ver == "celular_Artropodo":
            return ("huesped" if clase == "divergente" else "EVE_candidato"), "si"
        if cel:
            return "celular_otro", "si"
        sin_proteina = num(r.get("dm_pident")) is None and ver != "viral"
        qcov = num(r.get("bn_qcov"))
        if clase == "pariente_debil" and sin_proteina and (qcov is None or qcov < QCOV_MIN):
            return "evidencia_debil", ""
        return clase, ""
    if ver == "viral":
        return ("EVE_genoma" if en_gen else "viral_solo_nr"), ""
    if ver == "ambiguo":
        return "ambiguo", ""
    if ver == "celular_Artropodo":
        return "huesped", ""
    if cel:
        return "celular_otro", ""
    if ver == "sin_hit_nr":
        return ("huesped" if en_gen else "oscura"), ""
    return clase, ""


sin_arb = sin_gen = 0
cruce = Counter()
salida = []
for r in tabla:
    v = r["votu"]
    a = arb.get(v)
    if a is None:
        sin_arb += 1
        a = {c: "" for c in COLS_ARB}
    ver = a["veredicto"] or "sin_arbitro"

    g = gen.get(v)
    fg = None
    if g is not None:
        fg = float(g["frac_en_genoma"])
    elif f_gen:
        sin_gen += 1

    cf, conflicto = clase_final(r, ver, fg)
    cruce[(r["clase"], ver)] += 1

    r2 = dict(r)
    r2["nr_veredicto"] = ver
    for c in COLS_ARB[1:]:
        r2[c] = a[c]
    r2["frac_en_genoma"] = "" if fg is None else f"{fg:.3f}"
    r2["cromosomas_genoma"] = g["cromosomas"] if g else ""
    r2["clase_final"] = cf
    r2["es_virus"] = "si" if cf in ES_VIRUS else "no"
    r2["conflicto"] = conflicto
    if r2["es_virus"] == "si":
        especie = (r.get("especie") or "").strip()
        nr_viral = (a["nr_nombre"] or "").strip() if ver == "viral" else ""
        if cf == "identificado":
            nombre = especie or nr_viral or v
        else:
            nombre = "afín a " + (nr_viral or especie or v)
        es_fago = any(re.search(r"phage|Leviviricetes|Caudoviricetes", x, re.I) for x in (especie, nr_viral))
        r2["nombre_virus"] = nombre
        r2["grupo_virus"] = "fago" if es_fago else "artropodo"
    else:
        r2["nombre_virus"] = ""
        r2["grupo_virus"] = ""
    salida.append(r2)

cab_out = cab + ["nr_veredicto"] + COLS_ARB[1:] + ["frac_en_genoma", "cromosomas_genoma",
                                                   "clase_final", "es_virus", "conflicto",
                                                   "nombre_virus", "grupo_virus"]
with open(f_out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cab_out, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(salida)

# ---------------------------------------------------------------- informe
print(f"vOTUs: {len(tabla)}   sin veredicto nr: {sin_arb}"
      + (f"   sin cribado de genoma: {sin_gen}" if f_gen else "   (sin cribado de genoma)"))

clases = sorted({c for c, _ in cruce}, key=lambda c: -sum(n for (k, _), n in cruce.items() if k == c))
vers = sorted({v for _, v in cruce}, key=lambda v: -sum(n for (_, k), n in cruce.items() if k == v))
print("\n== clase (tabla_final) x veredicto (nr) ==")
print(f"  {'':16s}" + "".join(f"{v[:13]:>14s}" for v in vers))
for c in clases:
    print(f"  {c:16s}" + "".join(f"{cruce.get((c, v), 0) or '.':>14}" for v in vers))

print("\n== clase_final ==")
for c, n in Counter(r["clase_final"] for r in salida).most_common():
    print(f"  {c:16s} {n:4d}")
print(f"\n  es_virus = si: {sum(r['es_virus'] == 'si' for r in salida)} vOTUs")

if f_gen:
    eve = [r for r in salida if r["clase_final"] == "EVE_genoma"]
    print(f"\n== EVE en el genoma de referencia ({len(eve)}) ==")
    for r in sorted(eve, key=lambda r: -float(r["frac_en_genoma"])):
        print(f"  {r['clase']:14s} {r['especie'][:30]:30s} frac={r['frac_en_genoma']}  {r['cromosomas_genoma'][:40]}")

    virus = [r for r in salida if r["es_virus"] == "si" and r["frac_en_genoma"] and float(r["frac_en_genoma"]) > 0]
    print(f"\n== virus con algo de secuencia en el genoma, por debajo del umbral ({len(virus)}) ==")
    for r in sorted(virus, key=lambda r: -float(r["frac_en_genoma"])):
        print(f"  {r['clase_final']:14s} {r['especie'][:30]:30s} frac={r['frac_en_genoma']}")

for g in ("artropodo", "fago"):
    esp = Counter(r["nombre_virus"] for r in salida if r["es_virus"] == "si" and r["grupo_virus"] == g)
    print(f"\n== virus ({g}): {len(esp)} nombres, {sum(esp.values())} vOTUs ==")
    for e, n in esp.most_common():
        print(f"  {e[:45]:45s} {n:3d} vOTUs")

deb = [r for r in salida if r["clase_final"] == "evidencia_debil"]
print(f"\n== evidencia_debil: solo un tramo corto de nt, sin proteina ({len(deb)}) ==")
for r in deb:
    print(f"  {r['especie'][:32]:32s} bn={r['bn_pident']}%/{r['bn_qcov']}%  ensamblada en {r['n_muestras_asm']}")

print(f"\n  escrito: {f_out}")
