#!/usr/bin/env python3
"""
fusionar_veredictos.py — une los veredictos de nr de las dos ramas.

POR QUE
  El catalogo de 500 pb tiene 5,899 vOTUs: 5,164 nuevas y 735 heredadas del
  catalogo de 1 kb. El arbitraje de nr se corrio solo sobre las nuevas, asi que
  integrar_arbitro.py deja esas 735 como 'sin_arbitro' — y son justamente las
  mejor caracterizadas (anfevirus, Verdadero, Chaq-like, Guadeloupe).

QUE HACE
  Lee los dos ficheros de veredictos, los proyecta a las columnas que
  integrar_arbitro.py realmente usa (COLS_ARB, linea 53 de ese script) y
  escribe un fichero unico. La proyeccion es POR NOMBRE de columna, no por
  posicion, asi que no importa que los dos ficheros tengan distinto ancho.

  Si una vOTU aparece en los dos, gana la de la rama nueva y se avisa.

USO
  python3 fusionar_veredictos.py SALIDA.tsv FUENTE1.tsv [FUENTE2.tsv ...]
  Las fuentes se leen en orden de PRIORIDAD DECRECIENTE: la primera gana.

  Ejemplo:
    python3 fusionar_veredictos.py \\
        $RESULTS_DIR/qc_control_min500/veredicto_todas.tsv \\
        $RESULTS_DIR/qc_control_min500/veredicto_votus_nuevas.tsv \\
        $RESULTS_DIR/votu_annot/votus_nr_jr46/veredicto_por_votu.tsv

  Con --esperadas puede comprobar la cobertura contra la tabla final:
    ... --esperadas $RESULTS_DIR/qc_control_min500/TABLA_FINAL_votus.tsv
"""
import csv
import sys

COLS_ARB = ["veredicto", "nr_nombre", "nr_grupo", "nr_pident", "celular_nombre"]
COL_ID = "votu"


def leer(path):
    """Devuelve (cabecera, filas como dict). Detecta el separador."""
    with open(path, newline="") as fh:
        muestra = fh.readline()
        fh.seek(0)
        sep = "\t" if "\t" in muestra else ","
        r = csv.DictReader(fh, delimiter=sep)
        if r.fieldnames is None:
            return [], []
        return [c.strip() for c in r.fieldnames], [dict(x) for x in r]


def main():
    args = [a for a in sys.argv[1:]]
    esperadas = None
    if "--esperadas" in args:
        i = args.index("--esperadas")
        esperadas = args[i + 1]
        del args[i:i + 2]

    if len(args) < 2:
        print(__doc__)
        sys.exit(1)

    salida, fuentes = args[0], args[1:]

    fusion = {}          # votu -> dict con COL_ID + COLS_ARB
    procedencia = {}     # votu -> fichero del que salio
    solapes = []

    for path in fuentes:
        cab, filas = leer(path)
        if not cab:
            print(f"  AVISO: {path} esta vacio o no tiene cabecera")
            continue
        if COL_ID not in cab:
            print(f"  ERROR: {path} no tiene columna '{COL_ID}'. Cabecera: {cab[:8]}")
            sys.exit(1)
        faltan = [c for c in COLS_ARB if c not in cab]
        if faltan:
            print(f"  ERROR: a {path} le faltan columnas que integrar_arbitro.py usa: {faltan}")
            print(f"         cabecera: {cab}")
            sys.exit(1)

        nuevas = 0
        for f in filas:
            vid = (f.get(COL_ID) or "").strip()
            if not vid:
                continue
            if vid in fusion:
                solapes.append(vid)
                continue           # gana la fuente de mayor prioridad
            fusion[vid] = {COL_ID: vid, **{c: (f.get(c) or "") for c in COLS_ARB}}
            procedencia[vid] = path
            nuevas += 1
        print(f"  {nuevas:>6} vOTUs de {path}")

    with open(salida, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[COL_ID] + COLS_ARB, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for vid in fusion:
            w.writerow(fusion[vid])

    print(f"\n  escrito: {salida}   ({len(fusion)} vOTUs)")
    if solapes:
        print(f"  {len(solapes)} vOTUs estaban en mas de una fuente; gano la de mayor prioridad")
        for v in solapes[:5]:
            print(f"      {v}")

    # reparto de veredictos, util para comprobar que la fusion tiene sentido
    reparto = {}
    for v in fusion.values():
        reparto[v["veredicto"] or "(vacio)"] = reparto.get(v["veredicto"] or "(vacio)", 0) + 1
    print("\n  reparto de veredictos:")
    for k in sorted(reparto, key=lambda x: -reparto[x]):
        print(f"      {reparto[k]:>6}  {k}")

    if esperadas:
        cab_e, filas_e = leer(esperadas)
        ids = {(f.get(COL_ID) or "").strip() for f in filas_e}
        ids.discard("")
        sin = sorted(ids - set(fusion))
        print(f"\n  comprobacion contra {esperadas}")
        print(f"      vOTUs en la tabla final: {len(ids)}")
        print(f"      sin veredicto tras la fusion: {len(sin)}")
        for v in sin[:10]:
            print(f"          {v}")
        if len(sin) > 10:
            print(f"          ... y {len(sin) - 10} mas")


if __name__ == "__main__":
    main()
