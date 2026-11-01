#!/usr/bin/env python3
"""
diamond_taxa.py — añade a una tabla de DIAMOND (outfmt 6, sin cabecera) la
cabecera y cuatro columnas de taxonomía por hit, sacadas del taxdump de NCBI:

    grupo    Virus / Artropodo / Metazoo_otro / Hongo / Planta / Eucariota_otro /
             Bacteria / Arquea / otro / sin_taxid / taxid_desconocido
    realm, familia, genero   (vacíos si no aplica)

Autocontenido (solo biblioteca estándar): la asignación de grupo es la MISMA
regla que usa arbiter_nr.py (primer taxid de "a;b", merged.dmp, ancestros).

Con --verificar VEREDICTO.tsv (tabla del árbitro) recalcula el veredicto de cada
contig a partir de los grupos de esta tabla y lo compara con el del árbitro.

Uso:
  python3 diamond_taxa.py --taxdump DIR --cols "qseqid ... stitle" \
      --in diamond.tsv --out salida.tsv.gz [--verificar veredicto_por_votu.tsv]
"""
import argparse, gzip, os, sys

GROUPS = [("10239", "Virus"), ("6656", "Artropodo"), ("33208", "Metazoo_otro"), ("4751", "Hongo"),
          ("33090", "Planta"), ("2759", "Eucariota_otro"), ("2", "Bacteria"), ("2157", "Arquea")]
RANKS = ("realm", "family", "genus")


class Taxonomy:
    def __init__(self, taxdir):
        self.parent, self.rank, self.name, self.merged = {}, {}, {}, {}
        for f in ("nodes.dmp", "names.dmp"):
            if not os.path.exists(os.path.join(taxdir, f)):
                sys.exit(f"ERROR: falta {f} en {taxdir}")
        with open(os.path.join(taxdir, "nodes.dmp")) as fh:
            for line in fh:
                p = line.split("\t|\t")
                self.parent[p[0]] = p[1]; self.rank[p[0]] = p[2]
        with open(os.path.join(taxdir, "names.dmp")) as fh:
            for line in fh:
                p = line.rstrip("\t|\n").split("\t|\t")
                if len(p) > 3 and p[3] == "scientific name":
                    self.name[p[0]] = p[1]
        mf = os.path.join(taxdir, "merged.dmp")
        if os.path.exists(mf):
            with open(mf) as fh:
                for line in fh:
                    p = line.rstrip("\t|\n").split("\t|\t")
                    self.merged[p[0]] = p[1]

    def current(self, tid):
        return self.merged.get(str(tid), str(tid))

    def annotate(self, tid):
        """(grupo, realm, familia, genero) — grupo con la regla de arbiter_nr.py"""
        tid = (tid or "").split(";")[0]
        if not tid or tid in ("0", "N/A"):
            return ["sin_taxid", "", "", ""]
        t = self.current(tid)
        if t not in self.parent:
            return ["taxid_desconocido", "", "", ""]
        anc, lin, n = set(), {}, 0
        while t in self.parent and n < 80:
            anc.add(t)
            r = self.rank.get(t, "")
            if r in RANKS and r not in lin:
                lin[r] = self.name.get(t, "")
            if self.parent[t] == t:
                break
            t = self.parent[t]; n += 1
        g = next((name for gid, name in GROUPS if gid in anc), "otro")
        return [g] + [lin.get(r, "") for r in RANKS]


def verificar(hits, vfile, margin=1.1):
    """Recalcula el veredicto del árbitro con estos grupos y lo compara."""
    best, bv, bc = {}, {}, {}
    for q, b, g in hits:
        if q not in best or b > best[q][0]: best[q] = (b, g)
        if g in ("sin_taxid", "taxid_desconocido"): continue
        t = bv if g == "Virus" else bc
        if q not in t or b > t[q][0]: t[q] = (b, g)
    with open(vfile) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        iv, iq = hdr.index("veredicto"), hdr.index("votu")
        n = dif = 0; ej = []
        for line in fh:
            f = line.rstrip("\n").split("\t"); q = f[iq]; n += 1
            vb = bv[q][0] if q in bv else 0.0; cb = bc[q][0] if q in bc else 0.0
            if q not in best: v = "sin_hit_nr"
            elif vb and (not cb or vb >= cb): v = "viral" if (not cb or vb >= margin * cb) else "ambiguo"
            elif vb and vb >= cb / margin: v = "ambiguo"
            else: v = "celular_" + (bc[q][1] if q in bc else best[q][1])
            if v != f[iv]:
                dif += 1
                if len(ej) < 5: ej.append(f"{q}: arbitro={f[iv]} recalculado={v}")
    print(f"  verificación contra el árbitro: {n - dif}/{n} contigs con el mismo veredicto")
    for e in ej: print("    " + e)
    return dif


ap = argparse.ArgumentParser()
ap.add_argument("--taxdump", required=True)
ap.add_argument("--cols", required=True, help="nombres de columna, separados por espacios")
ap.add_argument("--in", dest="inp", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--verificar", default="", help="veredicto_por_votu.tsv del árbitro (solo para nr)")
ap.add_argument("--scripts-dir", default="", help=argparse.SUPPRESS)  # compatibilidad, no se usa
a = ap.parse_args()

cols = a.cols.split()
for c in ("staxids", "qseqid", "bitscore"):
    if c not in cols: sys.exit(f"ERROR: --cols no tiene '{c}'")
it, iq, ib = cols.index("staxids"), cols.index("qseqid"), cols.index("bitscore")
tax = Taxonomy(a.taxdump)
cache, n, malas, hits = {}, 0, 0, []
with open(a.inp) as fh, gzip.open(a.out, "wt") as out:
    out.write("\t".join(cols + ["grupo", "realm", "familia", "genero"]) + "\n")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) != len(cols):
            malas += 1; continue
        tid = (f[it] or "").split(";")[0]
        if tid not in cache:
            cache[tid] = tax.annotate(tid)
        out.write("\t".join(f + cache[tid]) + "\n")
        if a.verificar:
            hits.append((f[iq], float(f[ib]), cache[tid][0]))
        n += 1
if malas:
    sys.exit(f"ERROR: {malas} lineas con un numero de columnas distinto de {len(cols)}")
print(f"  {n} hits, {len(cache)} taxids distintos, "
      f"{sum(1 for v in cache.values() if v[0] in ('sin_taxid', 'taxid_desconocido'))} sin taxonomía")
if a.verificar:
    dif = verificar(hits, a.verificar)
    if dif:
        print(f"  AVISO: {dif} contigs con veredicto distinto (¿otro taxdump en el árbitro?). "
              "Las figuras usan el veredicto del árbitro; los grupos por hit pueden no cuadrar en esos contigs.")
