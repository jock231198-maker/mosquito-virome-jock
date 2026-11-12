#!/usr/bin/env python3
"""
srna_votu_metrics.py -- FASE 2 / NIVEL B, paso 2: metricas por vOTU x muestra.

Lee por STDIN un flujo de alineaciones ORDENADO POR NOMBRE DE LECTURA, con las
columnas que produce srna_map_worker.sh (join del SAM viral con la tabla del
genoma):

    lectura   n_loci_genoma   votu   pos   hebra   seq

  lectura        nombre del FASTA, con el peso al final: r123_x45
  n_loci_genoma  loci en AaegL5 de esa lectura (0 = no esta en el genoma)
  pos            POS del SAM, 1-based, extremo izquierdo
  hebra          + / -
  seq            SEQ del SAM (tal como se alineo) o * en alineaciones secundarias

Por que ordenado por lectura y no el SAM crudo: (1) el join con el genoma se
hace en disco, no en un diccionario de decenas de millones de nombres; (2) las
alineaciones de una misma lectura quedan juntas POR CONSTRUCCION, asi que
"lectura unica" (una sola vOTU) es exacto y no depende de que bowtie las emita
seguidas. Las secuencias emparentadas del catalogo recogen las mismas lecturas
(RESUMEN_smallRNA.md §4.6: Hubei/Kwale/Wenzhou sobemo-like salen en los mismos
pools que Guadeloupe), asi que cada vOTU lleva su conteo total y su conteo
unico, y frac_unica deja ver la inflacion.

PESOS: cada lectura vale el _x<w> de su nombre (los crudos del proveedor vienen
colapsados). Todos los conteos son pesados.

Salida: TSV por STDOUT, una fila por vOTU con lecturas en la muestra.

Uso (lo llama el worker):
    ... | srna_votu_metrics.py --muestra M2 --denominador 31500000 \\
              --longitudes longitudes.tsv [--perfil perfil.tsv]
"""

import argparse
import math
import sys
from collections import defaultdict

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def peso(nombre):
    """r123_x45 -> 45. Sin sufijo, 1."""
    i = nombre.rfind("_x")
    if i == -1:
        return 1
    try:
        w = int(nombre[i + 2:])
    except ValueError:
        return 1
    return w if w > 0 else 1


class Votu:
    __slots__ = ("lon", "cub", "lecturas", "unicas", "sent", "anti", "len_sent",
                 "len_anti", "pi_anti_1u", "pi_anti_tot", "pi_sent_10a",
                 "pi_sent_tot", "pp_sent", "pp_anti", "en_genoma",
                 "loci", "ini", "fin")

    def __init__(self, lon):
        self.lon = lon
        self.cub = bytearray(lon) if lon else bytearray()
        self.lecturas = 0
        self.unicas = 0
        self.sent = 0
        self.anti = 0
        self.len_sent = defaultdict(int)
        self.len_anti = defaultdict(int)
        self.pi_anti_1u = 0
        self.pi_anti_tot = 0
        self.pi_sent_10a = 0
        self.pi_sent_tot = 0
        self.pp_sent = defaultdict(int)
        self.pp_anti = defaultdict(int)
        self.en_genoma = 0
        self.loci = []
        self.ini = None
        self.fin = None


def pingpong_z(sent, anti, kmax=25):
    """z-score del solapamiento 5'-5' de 10 nt (firma de ping-pong).

    sent[p]  = peso de piRNA sentido cuyo 5' esta en p
    anti[p]  = peso de piRNA antisentido cuyo 5' esta en p
    Un solapamiento de k nt entre un sentido que empieza en p y un antisentido
    cuyo 5' cae en p + k - 1. z sobre los solapamientos 1..kmax sin el 10.
    """
    if not sent or not anti:
        return None, 0
    c = []
    for k in range(1, kmax + 1):
        tot = 0
        for p, w in sent.items():
            w2 = anti.get(p + k - 1)
            if w2:
                tot += w * w2
        c.append(tot)
    n_pares = sum(c)
    if n_pares == 0:
        return None, 0
    otros = [c[k - 1] for k in range(1, kmax + 1) if k != 10]
    mu = sum(otros) / len(otros)
    var = sum((x - mu) ** 2 for x in otros) / (len(otros) - 1)
    sd = math.sqrt(var)
    if sd == 0:
        return None, n_pares
    return (c[9] - mu) / sd, n_pares


def mediana(xs):
    if not xs:
        return 0
    xs = sorted(xs)
    n = len(xs)
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2


COLS = ["muestra", "votu", "longitud", "lecturas", "lecturas_unicas", "frac_unica",
        "rpm", "rpm_unico", "r18_23", "r21", "r24_35", "sentido", "antisentido",
        "polaridad", "frac_si", "frac_pi", "pico21", "si_frac_sentido",
        "pos_cubiertas", "frac_cubierta", "span_frac", "u1_anti", "a10_sentido",
        "pp_z10", "pp_pares", "frac_en_genoma", "loci_genoma_mediana"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--muestra", required=True)
    ap.add_argument("--denominador", type=int, required=True,
                    help="lecturas 18-35 no rRNA de la libreria (denominador de RPM)")
    ap.add_argument("--longitudes", required=True,
                    help="TSV votu<TAB>longitud (del encabezado @SQ del SAM)")
    ap.add_argument("--perfil", help="TSV opcional: votu, longitud_nt, hebra, lecturas")
    ap.add_argument("--si", nargs=2, type=int, default=[18, 23], metavar=("MIN", "MAX"))
    ap.add_argument("--pi", nargs=2, type=int, default=[24, 35], metavar=("MIN", "MAX"))
    a = ap.parse_args()
    si_min, si_max = a.si
    pi_min, pi_max = a.pi

    lon = {}
    with open(a.longitudes) as fh:
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 2 and p[1].isdigit():
                lon[p[0]] = int(p[1])

    V = {}
    actual = None          # nombre de la lectura en curso
    grupo = []             # [(votu, pos, hebra)]
    seq = None
    w = 1
    nloci = 0

    def cerrar():
        if not grupo:
            return
        L = len(seq) if seq else 0
        es_si = si_min <= L <= si_max
        es_pi = pi_min <= L <= pi_max
        unica = len({g[0] for g in grupo}) == 1
        for votu, pos, hebra in grupo:
            v = V.get(votu)
            if v is None:
                v = V[votu] = Votu(lon.get(votu, 0))
            v.lecturas += w
            if unica:
                v.unicas += w
            if nloci:
                v.en_genoma += w
                v.loci.append(nloci)
            if hebra == "+":
                v.sent += w
                if L:
                    v.len_sent[L] += w
                if es_pi:
                    v.pi_sent_tot += w
                    if len(seq) >= 10 and seq[9] == "A":
                        v.pi_sent_10a += w
                    v.pp_sent[pos] += w
            else:
                v.anti += w
                if L:
                    v.len_anti[L] += w
                if es_pi:
                    v.pi_anti_tot += w
                    if seq and seq[0] == "T":
                        v.pi_anti_1u += w
                    v.pp_anti[pos + L - 1] += w
            if L and v.lon:
                ini = max(0, pos - 1)
                fin = min(v.lon, ini + L)
                if fin > ini:
                    v.cub[ini:fin] = b"\x01" * (fin - ini)
                    v.ini = ini if v.ini is None else min(v.ini, ini)
                    v.fin = fin if v.fin is None else max(v.fin, fin)

    raros = 0
    for linea in sys.stdin:
        c = linea.rstrip("\n").split("\t")
        if len(c) < 6:
            raros += 1
            continue
        nombre, g, votu, pos, hebra, s = c[0], c[1], c[2], c[3], c[4], c[5]
        if nombre != actual:
            cerrar()
            actual, grupo, seq, w = nombre, [], None, peso(nombre)
            nloci = int(g) if g.isdigit() else 0
        if s != "*" and seq is None:
            # SEQ viene en la orientacion de la referencia; se guarda la de la lectura
            seq = revcomp(s) if hebra == "-" else s
        grupo.append((votu, int(pos), hebra))
    cerrar()

    if raros:
        print(f"[metrics] {raros} lineas mal formadas ignoradas", file=sys.stderr)

    perfil = open(a.perfil, "w") if a.perfil else None
    if perfil:
        perfil.write("muestra\tvotu\tlongitud_nt\thebra\tlecturas\n")

    print("\t".join(COLS))
    M = 1e6 / a.denominador if a.denominador else 0.0
    for votu in sorted(V):
        v = V[votu]
        tot = v.lecturas
        if tot <= 0:
            continue
        r_si = sum(w2 for L, w2 in v.len_sent.items() if si_min <= L <= si_max) + \
               sum(w2 for L, w2 in v.len_anti.items() if si_min <= L <= si_max)
        r_pi = sum(w2 for L, w2 in v.len_sent.items() if pi_min <= L <= pi_max) + \
               sum(w2 for L, w2 in v.len_anti.items() if pi_min <= L <= pi_max)
        r21 = v.len_sent.get(21, 0) + v.len_anti.get(21, 0)
        otros_si = [v.len_sent.get(L, 0) + v.len_anti.get(L, 0)
                    for L in range(si_min, si_max + 1) if L != 21]
        pico21 = r21 > 0 and r21 >= (max(otros_si) if otros_si else 0)
        si_sent = sum(w2 for L, w2 in v.len_sent.items() if si_min <= L <= si_max)
        cub = sum(v.cub) if v.lon else 0
        frac_cub = cub / v.lon if v.lon else 0.0
        span = (v.fin - v.ini) / v.lon if (v.lon and v.ini is not None) else 0.0
        z, pares = pingpong_z(v.pp_sent, v.pp_anti)
        fila = [
            a.muestra, votu, v.lon, tot, v.unicas, f"{v.unicas / tot:.4f}",
            f"{tot * M:.3f}", f"{v.unicas * M:.3f}", r_si, r21, r_pi, v.sent, v.anti,
            f"{v.sent / tot:.4f}",
            f"{r_si / tot:.4f}", f"{r_pi / tot:.4f}",
            "si" if pico21 else "no",
            f"{si_sent / r_si:.4f}" if r_si else "NA",
            cub, f"{frac_cub:.4f}", f"{span:.4f}",
            f"{v.pi_anti_1u / v.pi_anti_tot:.4f}" if v.pi_anti_tot else "NA",
            f"{v.pi_sent_10a / v.pi_sent_tot:.4f}" if v.pi_sent_tot else "NA",
            f"{z:.3f}" if z is not None else "NA", pares,
            f"{v.en_genoma / tot:.4f}", f"{mediana(v.loci):.1f}",
        ]
        print("\t".join(str(x) for x in fila))
        if perfil:
            for L in sorted(set(v.len_sent) | set(v.len_anti)):
                if v.len_sent.get(L):
                    perfil.write(f"{a.muestra}\t{votu}\t{L}\t+\t{v.len_sent[L]}\n")
                if v.len_anti.get(L):
                    perfil.write(f"{a.muestra}\t{votu}\t{L}\t-\t{v.len_anti[L]}\n")
    if perfil:
        perfil.close()


if __name__ == "__main__":
    main()
