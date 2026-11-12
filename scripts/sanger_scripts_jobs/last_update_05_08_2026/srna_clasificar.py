#!/usr/bin/env python3
"""
srna_clasificar.py -- FASE 2 / NIVEL B, paso 3: replica activa frente a EVE.

Aplica a cada par vOTU x muestra los criterios de PLAN_smallRNA.md §4c, fijados
ANTES de ver los datos. Aqui solo se hacen explicitos los dos que el plan dejaba
en prosa; los umbrales son argumentos, de modo que mover uno queda registrado en
la linea de comandos y no escondido en el codigo:

  "lecturas confinadas a un tramo" (clase EVE)
      frac_cubierta <= --tramo-cubierto (0.30) Y span_frac <= --tramo-span (0.30)
      Un EVE es un fragmento integrado: las lecturas se apinan en el trozo
      homologo, no se reparten por el genoma viral entero.

  "ping-pong" (clase virus + piRNA)
      pp_z10 >= --pp-z (2.0), con al menos --pp-pares (100) pares de
      solapamiento. z de 2 es el corte habitual; por debajo de 100 pares el
      z-score de un vOTU con pocas lecturas no significa nada.

Clases (se evaluan en este orden; la primera que encaja gana):

  no_interpretable_rrna  el vOTU lleva rRNA (barrnap). En el catalogo de 500 bp
                         el 77-92 % de las bases mapeadas va a vOTUs de rRNA y
                         SortMeRNA nunca corrio: se marca y no se interpreta.
  bajo_umbral            rpm < --rpm-min (10)
  replica_activa         pico de 21 nt + siRNA en las dos hebras + cobertura
                         repartida + mayoria de lecturas NO en el genoma
                         (sufijo _con_piRNA si ademas hay ping-pong: Lau 2026
                         lo describe para ANPHV y sigue siendo replicacion)
  virus_mas_piRNA        siRNA simetricos y ping-pong, sin pico de 21 nt claro
  EVE_NIRVS              piRNA de una hebra, sin pico de 21 nt, confinados a un
                         tramo y presentes en AaegL5
  ambiguo                lo demas -> perfil de cobertura a mano

El motivo lleva los criterios que FALLARON, para poder discutir un ambiguo sin
volver a la tabla.

Consenso por vOTU: un virus que replica en un solo pool replica. Pero si un vOTU
sale replica en unos pools y EVE en otros se marca MIXTO, que es exactamente la
sospecha abierta de CFAV (infeccion real en M72/M61/M41 sobre un fondo de ~27 %
de cobertura en todas las librerias, RESUMEN_smallRNA.md §3).

Uso:
    srna_clasificar.py metricas.tsv --salida clases.tsv --consenso consenso.tsv \\
        [--rrna votus_rrna.txt] [--en-genoma votus_en_genoma.txt]
"""

import argparse
import csv
import sys
from collections import defaultdict


def num(x, d=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return d


def lista(path):
    if not path:
        return set()
    out = set()
    with open(path) as fh:
        for ln in fh:
            s = ln.split("\t")[0].strip()
            if s and not s.startswith("#"):
                out.add(s)
    return out


def clasificar(r, a, en_rrna):
    """Devuelve (clase, motivo)."""
    votu = r["votu"]
    if votu in en_rrna:
        return "no_interpretable_rrna", "vOTU con rRNA (barrnap)"

    rpm = num(r["rpm"])
    if rpm < a.rpm_min:
        return "bajo_umbral", f"rpm {rpm:.2f} < {a.rpm_min}"

    pico21 = r["pico21"] == "si"
    si_sent = num(r["si_frac_sentido"], -1.0)
    si_bihebra = 0.0 <= si_sent and a.hebra_min <= si_sent <= a.hebra_max
    frac_cub = num(r["frac_cubierta"])
    frac_gen = num(r["frac_en_genoma"])
    frac_pi = num(r["frac_pi"])
    pol = num(r["polaridad"])
    una_hebra = max(pol, 1 - pol) >= a.eve_hebra
    span = num(r["span_frac"])
    z = r["pp_z10"]
    pp = (z != "NA" and num(z) >= a.pp_z and num(r["pp_pares"]) >= a.pp_pares)

    # --- replica activa
    falla = []
    if not pico21:
        falla.append("sin pico de 21 nt")
    if not si_bihebra:
        falla.append(f"siRNA no bihebra (sentido {si_sent:.2f})" if si_sent >= 0
                     else "sin siRNA")
    if frac_cub < a.cob_min:
        falla.append(f"cubre {frac_cub:.2f} < {a.cob_min}")
    if frac_gen >= 0.5:
        falla.append(f"{frac_gen:.2f} de las lecturas en AaegL5")
    if not falla:
        return ("replica_activa_con_piRNA" if pp else "replica_activa"), "ok"

    # --- virus + piRNA: siRNA simetricos y ping-pong, aunque el 21 no domine
    if si_bihebra and pp and frac_gen < 0.5:
        return "virus_mas_piRNA", f"z10 {num(z):.1f}; " + "; ".join(falla)

    # --- EVE / NIRVS
    falla_eve = []
    if frac_pi < a.pi_min:
        falla_eve.append(f"piRNA {frac_pi:.2f} < {a.pi_min}")
    if not una_hebra:
        falla_eve.append(f"no es de una hebra (polaridad {pol:.2f})")
    if pico21:
        falla_eve.append("tiene pico de 21 nt")
    if not (frac_cub <= a.tramo_cubierto and span <= a.tramo_span):
        falla_eve.append(f"no confinado (cubre {frac_cub:.2f}, span {span:.2f})")
    if frac_gen < a.eve_genoma:
        falla_eve.append(f"solo {frac_gen:.2f} de las lecturas en AaegL5")
    if not falla_eve:
        return "EVE_NIRVS", "ok"

    return "ambiguo", "replica: " + "; ".join(falla) + " | EVE: " + "; ".join(falla_eve)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("metricas", help="TSV junto de srna_votu_metrics.py (todas las muestras)")
    ap.add_argument("--salida", required=True)
    ap.add_argument("--consenso", required=True)
    ap.add_argument("--rrna", help="lista de vOTUs con rRNA (barrnap): votus_rrna")
    ap.add_argument("--en-genoma", dest="en_genoma",
                    help="lista de vOTUs del cribado AaegL5 (solo informativa)")
    ap.add_argument("--rpm-min", type=float, default=10.0)
    ap.add_argument("--cob-min", type=float, default=0.50, help="replica: fraccion del vOTU cubierta")
    ap.add_argument("--hebra-min", type=float, default=0.10, help="replica: siRNA por hebra, minimo")
    ap.add_argument("--hebra-max", type=float, default=0.90, help="replica: siRNA por hebra, maximo")
    ap.add_argument("--pi-min", type=float, default=0.80, help="EVE: fraccion de piRNA")
    ap.add_argument("--eve-hebra", type=float, default=0.80, help="EVE: fraccion de una hebra")
    ap.add_argument("--eve-genoma", type=float, default=0.50, help="EVE: fraccion de lecturas en AaegL5")
    ap.add_argument("--tramo-cubierto", type=float, default=0.30)
    ap.add_argument("--tramo-span", type=float, default=0.30)
    ap.add_argument("--pp-z", type=float, default=2.0)
    ap.add_argument("--pp-pares", type=float, default=100)
    a = ap.parse_args()

    en_rrna = lista(a.rrna)
    en_gen = lista(a.en_genoma)

    por_votu = defaultdict(lambda: defaultdict(int))
    rpm_max = defaultdict(float)
    muestras = defaultdict(set)

    with open(a.metricas) as fh, open(a.salida, "w") as out:
        rd = csv.DictReader(fh, delimiter="\t")
        faltan = {"votu", "muestra", "rpm", "pico21", "si_frac_sentido",
                  "frac_cubierta", "frac_en_genoma", "frac_pi", "polaridad",
                  "span_frac", "pp_z10", "pp_pares"} - set(rd.fieldnames or [])
        if faltan:
            sys.exit(f"ERROR: a {a.metricas} le faltan columnas: {sorted(faltan)}")
        cols = list(rd.fieldnames) + ["en_cribado_genoma", "clase", "motivo"]
        out.write("\t".join(cols) + "\n")
        for r in rd:
            clase, motivo = clasificar(r, a, en_rrna)
            r["en_cribado_genoma"] = "si" if r["votu"] in en_gen else "no"
            r["clase"] = clase
            r["motivo"] = motivo
            out.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
            por_votu[r["votu"]][clase] += 1
            muestras[r["votu"]].add(r["muestra"])
            rpm_max[r["votu"]] = max(rpm_max[r["votu"]], num(r["rpm"]))

    with open(a.consenso, "w") as out:
        out.write("votu\tmuestras\treplica\tvirus_mas_piRNA\tEVE_NIRVS\tambiguo\t"
                  "bajo_umbral\trrna\trpm_max\tconsenso\n")
        for votu in sorted(por_votu):
            c = por_votu[votu]
            rep = c["replica_activa"] + c["replica_activa_con_piRNA"]
            vpi = c["virus_mas_piRNA"]
            eve = c["EVE_NIRVS"]
            amb = c["ambiguo"]
            bajo = c["bajo_umbral"]
            rrna = c["no_interpretable_rrna"]
            if rrna:
                cons = "no_interpretable_rrna"
            elif (rep + vpi) and eve:
                cons = "MIXTO_replica_y_EVE"
            elif rep + vpi:
                cons = "replica"
            elif eve:
                cons = "EVE_NIRVS"
            elif amb:
                cons = "ambiguo"
            else:
                cons = "bajo_umbral"
            out.write(f"{votu}\t{len(muestras[votu])}\t{rep}\t{vpi}\t{eve}\t{amb}\t"
                      f"{bajo}\t{rrna}\t{rpm_max[votu]:.2f}\t{cons}\n")


if __name__ == "__main__":
    main()
