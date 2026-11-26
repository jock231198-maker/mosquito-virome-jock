#!/usr/bin/env python3
"""
Separa los contigs con hit a Wolbachia en rRNA y no-rRNA.

La hipotesis a matar: todos los hits son ARN ribosomal conservado, que pega
contra cualquier Wolbachia por igual, y no queda nada especifico de Wolbachia.
Si eso es cierto, no hay Wolbachia en las muestras y la senal de Kraken2 es un
artefacto de no haber depletado rRNA.

Si queda sequence no-rRNA concentrada en M20/M21, entonces si hay algo.

Modos:
  extraer   hits.tsv union.fasta salida.fasta        saca los contigs con hit
  resumir   hits.tsv barrnap.gff salida_prefijo      clasifica y resume

    python3 -I wolb_rrna.py extraer hits.tsv union.fasta wolb_contigs.fasta
    python3 -I wolb_rrna.py resumir hits.tsv barrnap.gff wolb
"""

import collections
import sys

IDENT_MIN = 90.0
LEN_MIN = 200


def leer_hits(p):
    """contig -> lista de (ref, pident, alen, sstart, send)"""
    h = collections.defaultdict(list)
    with open(p) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            try:
                pid, alen = float(f[2]), int(f[3])
                ss, se = int(f[6]), int(f[7])
            except ValueError:
                continue
            if pid < IDENT_MIN or alen < LEN_MIN:
                continue
            h[f[0]].append((f[1], pid, alen, min(ss, se), max(ss, se)))
    return h


def cmd_extraer(hits_p, fasta_p, out_p):
    quiero = set(leer_hits(hits_p))
    escritos, largo = 0, {}
    with open(fasta_p) as fh, open(out_p, "w") as out:
        guardar, nombre = False, None
        for ln in fh:
            if ln.startswith(">"):
                nombre = ln[1:].split()[0]
                guardar = nombre in quiero
                if guardar:
                    escritos += 1
                    largo[nombre] = 0
                    out.write(ln)
            elif guardar:
                out.write(ln)
                largo[nombre] += len(ln.strip())
    print(f"contigs con hit: {len(quiero)}")
    print(f"contigs escritos: {escritos}")
    faltan = quiero - set(largo)
    if faltan:
        print(f"AVISO: no encontre {len(faltan)} en el fasta, p.ej. "
              f"{sorted(faltan)[:3]}", file=sys.stderr)
    if largo:
        ls = sorted(largo.values())
        print(f"longitud: min {ls[0]}  mediana {ls[len(ls) // 2]}  max {ls[-1]}")


def cmd_resumir(hits_p, gff_p, pref):
    hits = leer_hits(hits_p)

    # barrnap: contig -> bp cubiertos por genes de rRNA, y que genes
    rrna_bp = collections.defaultdict(int)
    rrna_gen = collections.defaultdict(set)
    with open(gff_p) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "rRNA":
                continue
            try:
                a, b = int(f[3]), int(f[4])
            except ValueError:
                continue
            rrna_bp[f[0]] += abs(b - a) + 1
            for campo in f[8].split(";"):
                if campo.startswith("Name="):
                    rrna_gen[f[0]].add(campo[5:])

    # un contig es rRNA si el rRNA cubre >=50 % de lo que alineo contra Wolbachia
    filas = []
    for ctg, hs in hits.items():
        muestra = ctg.split("__")[0]
        alen_max = max(h[2] for h in hs)
        es_rrna = rrna_bp.get(ctg, 0) >= 0.5 * alen_max
        mejor = max(hs, key=lambda h: h[2])
        filas.append((muestra, ctg, alen_max, es_rrna,
                      ",".join(sorted(rrna_gen.get(ctg, []))) or "-",
                      mejor[0], mejor[1], mejor[3], mejor[4]))

    por_m = collections.defaultdict(lambda: [0, 0, 0, 0])
    for m, _c, alen, es_rrna, _g, _r, _p, _s, _e in filas:
        por_m[m][0] += 1
        if es_rrna:
            por_m[m][1] += 1
            por_m[m][2] += alen
        else:
            por_m[m][3] += alen

    print("== por muestra ==")
    print("muestra\tcontigs\tson_rRNA\tno_rRNA\tbp_rRNA\tbp_no_rRNA")
    for m, (tot, nr, bpr, bpn) in sorted(por_m.items(), key=lambda t: -t[1][3]):
        print(f"{m}\t{tot}\t{nr}\t{tot - nr}\t{bpr}\t{bpn}")

    print()
    print("== genes de rRNA encontrados ==")
    c = collections.Counter()
    for _m, ctg, _a, es, g, *_ in filas:
        if es:
            c[g] += 1
    for g, k in c.most_common():
        print(f"  {k:4d} contigs  {g}")

    print()
    nr = [f for f in filas if not f[3]]
    print(f"== contigs NO rRNA: {len(nr)} ==")
    if nr:
        print("muestra\tcontig\tbp_align\tmejor_ref\tident\tref_ini\tref_fin")
        for f in sorted(nr, key=lambda r: -r[2])[:40]:
            print(f"{f[0]}\t{f[1]}\t{f[2]}\t{f[5]}\t{f[6]}\t{f[7]}\t{f[8]}")
    else:
        print("  ninguno: todo lo que pega contra Wolbachia es rRNA conservado")

    with open(pref + "_clasificacion.tsv", "w") as out:
        out.write("muestra\tcontig\tbp_align\tes_rrna\tgenes_rrna"
                  "\tmejor_ref\tident\tref_ini\tref_fin\n")
        for f in filas:
            out.write("\t".join(str(x) for x in f) + "\n")
    print()
    print(f"escrito {pref}_clasificacion.tsv")


if len(sys.argv) < 2:
    sys.exit(__doc__)
if sys.argv[1] == "extraer" and len(sys.argv) == 5:
    cmd_extraer(sys.argv[2], sys.argv[3], sys.argv[4])
elif sys.argv[1] == "resumir" and len(sys.argv) == 5:
    cmd_resumir(sys.argv[2], sys.argv[3], sys.argv[4])
else:
    sys.exit(__doc__)
