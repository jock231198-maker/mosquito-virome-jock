#!/usr/bin/env python3
"""
nivelB_objetivos.py — escribe la lista de vOTUs objetivo del nivel B.

POR QUE
  nivelB_clases.tsv va a tener una fila por cada vOTU con lecturas de ARN
  pequeno, de un catalogo de 5,899. Sin una lista de objetivos el RESUMEN.txt
  del colector es una lista plana de ids y hay que volver a cruzar a mano.
  Esta lista agrupa las vOTUs por LA PREGUNTA que cada una responde.

QUE ESCRIBE
  nivelB_objetivos.txt            un id por linea, lo que come el colector
  nivelB_objetivos_anotados.tsv   el mismo conjunto con motivo y contexto

LA COLUMNA QUE IMPORTA: prediccion_largo
  RESUMEN_min500.md 5.4b separo PCLV de sus propios EVE con los perfiles de
  cobertura del RNA largo: las 7 vOTUs de PCLV dan CV de profundidad 0.40-0.74
  y los 5 EVE de Phasivirus/HTV dan 0.85-2.37, sin solape. El nivel B decide lo
  mismo por una via independiente (siRNA de 21 nt frente a piRNA).

  Esta columna deja la prediccion ESCRITA ANTES de correr, que es la unica
  forma de que el acuerdo signifique algo. Al terminar se cruza
  prediccion_largo con la clase del nivel B: si coinciden, hay dos tecnicas
  independientes de acuerdo; si no, el criterio de CV se cae y hay que saberlo.

USO
  python3 nivelB_objetivos.py TABLA_FINAL_votus_arb.tsv DIR_SALIDA [PERFILES.tsv]

  TABLA   $RESULTS_DIR/qc_control_min500/TABLA_FINAL_votus_arb.tsv  (rama MERIDA)
  SALIDA  el results de MERIDA_SRNA, que es donde lo busca el colector
  PERFILES  qc_control_min500/perfil_cobertura/perfiles_cobertura.tsv
            opcional: sin el, prediccion_largo queda vacia
"""
import csv
import gzip
import os
import statistics as st
import sys
from collections import Counter, defaultdict

# El control positivo de perfil EVE (PLAN_nivelB_fase2.md 6). No es AEFE1:
# ese elemento no esta depositado en ninguna base.
CONTROL = "EVE_CFAV_AaegL5_chr1_469224864_469225891"

# Las cuatro especies nuevas de RESUMEN_min500.md 5.3, por su nombre en nr.
NUEVAS = ("Palmetto toti-like virus", "Grenada mosquito rhabdovirus 1",
          "Sinu virus", "High Island virus")

# Umbrales de la prediccion del RNA largo (RESUMEN_min500.md 5.4b).
CV_CORTE = 0.80
FRAC_CORTE = 0.50


def motivo_de(r):
    """Un motivo por vOTU, en orden de prioridad: el primero que encaje."""
    txt = " ".join((r.get(c) or "") for c in ("nombre_virus", "especie", "nr_nombre"))
    cf = (r.get("clase_final") or "").strip()
    es_virus = (r.get("es_virus") or "").strip().lower() in ("si", "sí", "1", "true")
    nombre = (r.get("nombre_virus") or "").strip()

    if "Phasivirus" in txt:
        return "PCLV"
    if "Humaita" in txt:
        return "HTV"
    if "Cell fusing agent" in txt:
        return "CFAV"
    if any(n in txt for n in NUEVAS):
        return "especie_nueva"
    if cf == "pariente":
        return "pariente_especie_nueva"
    if es_virus and nombre.startswith("af"):      # "afín a X" = los "-like"
        return "virus_like"
    if es_virus:
        return "virus_nombrado"
    if cf == "EVE_genoma":
        return "EVE_genoma"
    if cf == "EVE_candidato":
        return "EVE_candidato"
    return ""


def perfiles_cv(path, interesan):
    """CV de profundidad y cobertura de cada vOTU, en su mejor muestra."""
    if not path or not os.path.exists(path):
        return {}
    por = defaultdict(list)
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        cab = fh.readline()
        if "profundidad" not in cab:          # sin cabecera: rebobinar
            fh.seek(0)
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) < 5 or p[1] not in interesan:
                continue
            try:
                por[(p[1], p[0])].append(float(p[4]))
            except ValueError:
                continue
    out = {}
    for (v, s), xs in por.items():
        cub = sum(1 for x in xs if x > 0) / len(xs)
        m = st.mean(xs)
        cv = (st.pstdev(xs) / m) if m else 0.0
        # la mejor muestra es la de mas cobertura; a igualdad, la mas profunda
        clave = (cub, m)
        if v not in out or clave > out[v][0]:
            out[v] = (clave, s, cub, cv, m)
    return {v: (s, cub, cv, m) for v, (_, s, cub, cv, m) in out.items()}


def prediccion(frac, cv):
    """La prediccion preinscrita de 5.4b. Vacia si falta un dato."""
    if frac is None or cv is None:
        return ""
    if frac < FRAC_CORTE and cv < CV_CORTE:
        return "replica"
    if frac >= FRAC_CORTE and cv >= CV_CORTE:
        return "EVE"
    return "indeterminado"


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    tabla, salida = sys.argv[1], sys.argv[2]
    perfiles = sys.argv[3] if len(sys.argv) > 3 else None
    os.makedirs(salida, exist_ok=True)

    with open(tabla, newline="") as fh:
        filas = list(csv.DictReader(fh, delimiter="\t"))
    if not filas:
        sys.exit(f"ERROR: {tabla} esta vacia")

    obj = []
    for r in filas:
        m = motivo_de(r)
        if m:
            obj.append((r["votu"], m, r))

    ids = {v for v, _, _ in obj}
    cv = perfiles_cv(perfiles, ids)

    # el control va siempre, y no esta en la tabla del transcriptoma
    lista = [v for v, _, _ in obj] + [CONTROL]

    with open(os.path.join(salida, "nivelB_objetivos.txt"), "w") as fh:
        fh.write("\n".join(lista) + "\n")

    COLS = ["votu", "motivo", "clase_final", "nombre_virus", "longitud",
            "frac_en_genoma", "n_muestras_mapeo", "muestra_cob", "cobertura",
            "cv_profundidad", "prof_media", "prediccion_largo"]
    anot = os.path.join(salida, "nivelB_objetivos_anotados.tsv")
    with open(anot, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(COLS)
        w.writerow([CONTROL, "CONTROL", "EVE_genoma_conocido", "EVE de CFAV en AaegL5",
                    "1028", "1.000", "", "", "", "", "", "EVE"])
        for v, m, r in sorted(obj, key=lambda t: (t[1], t[0])):
            try:
                fg = float(r.get("frac_en_genoma") or "")
            except ValueError:
                fg = None
            s, cub, c, prof = cv.get(v, ("", None, None, None))
            w.writerow([
                v, m, r.get("clase_final", ""), r.get("nombre_virus", ""),
                r.get("longitud", ""), r.get("frac_en_genoma", ""),
                r.get("n_muestras_mapeo", ""), s,
                f"{cub:.3f}" if cub is not None else "",
                f"{c:.3f}" if c is not None else "",
                f"{prof:.1f}" if prof is not None else "",
                prediccion(fg, c),
            ])

    # --- resumen en pantalla ---
    print(f"objetivos: {len(lista)} vOTUs (incluido el control)\n")
    cnt = Counter(m for _, m, _ in obj)
    print(f"{'motivo':26} {'n':>5}")
    print("-" * 33)
    for m, n in sorted(cnt.items(), key=lambda t: -t[1]):
        print(f"{m:26} {n:5d}")
    print(f"{'CONTROL':26} {1:5d}")

    if cv:
        pred = Counter()
        for v, m, r in obj:
            try:
                fg = float(r.get("frac_en_genoma") or "")
            except ValueError:
                fg = None
            _, _, c, _ = cv.get(v, ("", None, None, None))
            pred[prediccion(fg, c) or "(sin perfil)"] += 1
        print(f"\nprediccion del RNA largo (5.4b), a comprobar con el nivel B:")
        for k, n in sorted(pred.items(), key=lambda t: -t[1]):
            print(f"  {k:24} {n:5d}")
    else:
        print("\n  AVISO: sin perfiles de cobertura, prediccion_largo va vacia")

    print(f"\nescrito: {os.path.join(salida, 'nivelB_objetivos.txt')}")
    print(f"escrito: {anot}")


if __name__ == "__main__":
    main()
