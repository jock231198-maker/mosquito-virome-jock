#!/usr/bin/env python3
"""
Los EVE_candidato y los evidencia_debil que el ARN pequeno llama replicativos.

QUE PREGUNTA RESUELVE
  EVE_genoma (frac_en_genoma >= 0.5, confirmados en AaegL5) tiene CERO replica:
  las dos tecnicas concuerdan sin una sola excepcion. EVE_candidato es, por
  definicion, lo que NO esta en AaegL5 -- una etiqueta provisional. Los 5 con
  perfil de replica no son "EVE que replican": son candidatos a virus que el
  ARN largo no pudo confirmar y el ARN pequeno rescata.

  revisar_evidencia_debil es otra cosa y hay que mirarla de cerca: son las 21
  que el suelo de identidad degrado, diez con nombre de dengue. Si alguno de
  los que dan replica lleva nombre de dengue, eso hay que resolverlo antes de
  entregar nada.

QUE COMPRUEBA, ademas de listarlos
  Si los discordantes son virus de verdad, deberian distinguirse TAMBIEN en el
  ARN largo: menos frac_en_genoma que los EVE_candidato que si dieron perfil de
  EVE. Si no se distinguen, el rescate se apoya solo en una tecnica y hay que
  decirlo asi.

USO
  python3 -I nivelB_discordantes.py nivelB_consenso.tsv nivelB_objetivos_anotados.tsv
"""

import collections
import statistics
import sys

if len(sys.argv) < 3:
    sys.exit(__doc__)

CONSENSO, OBJETIVOS = sys.argv[1], sys.argv[2]

# nombres que obligan a parar y mirar: el punto mas sensible del documento
ALARMA = ("dengue", "denv", "zika", "chikung", "oropouche", "west nile",
          "yellow fever", "mayaro", "usutu", "japanese encephalitis")


def leer(path):
    with open(path) as fh:
        cab = fh.readline().rstrip("\n").split("\t")
        for ln in fh:
            v = ln.rstrip("\n").split("\t")
            if len(v) >= len(cab):
                yield dict(zip(cab, v))


def num(x, por_defecto=float("nan")):
    try:
        return float(x)
    except (TypeError, ValueError):
        return por_defecto


cons = {r["votu"]: r for r in leer(CONSENSO)}
obj = list(leer(OBJETIVOS))
print(f"consenso: {len(cons)} filas   objetivos: {len(obj)} filas")

falta = [r["votu"] for r in obj if r["votu"] not in cons]
if falta:
    print(f"AVISO: {len(falta)} objetivos sin fila en el consenso, p.ej. {falta[:3]}")

# --------------------------------------------------- 1. la matriz de control
print("\n== motivo x consenso ==")
mat = collections.defaultdict(collections.Counter)
clases = collections.Counter()
for r in obj:
    c = cons.get(r["votu"], {}).get("consenso", "SIN_FILA")
    mat[r["motivo"]][c] += 1
    clases[c] += 1
orden = [c for c, _ in clases.most_common()]
anchura = max(len(m) for m in mat) + 2
print("motivo".ljust(anchura) + "".join(c[:13].rjust(15) for c in orden) + "total".rjust(9))
for m in sorted(mat, key=lambda k: -sum(mat[k].values())):
    fila = "".join(str(mat[m][c] or "").rjust(15) for c in orden)
    print(m.ljust(anchura) + fila + str(sum(mat[m].values())).rjust(9))

# -------------------------------------------------- 2. los discordantes
MOTIVOS = ("EVE_candidato", "revisar_evidencia_debil")
REPLICA = ("replica", "MIXTO_replica_y_EVE")

disc = [r for r in obj
        if r["motivo"] in MOTIVOS
        and cons.get(r["votu"], {}).get("consenso", "") in REPLICA]

print(f"\n== los {len(disc)} discordantes, en detalle ==")
for r in sorted(disc, key=lambda r: (r["motivo"], -num(cons[r["votu"]]["rpm_max"], 0))):
    c = cons[r["votu"]]
    nombre = (r.get("nombre_virus") or c.get("nombre_virus") or "(sin nombre)").strip()
    bandera = ""
    bajo = nombre.lower()
    if any(a in bajo for a in ALARMA):
        bandera = "   <<<< NOMBRE DE ARBOVIRUS: RESOLVER ANTES DE ENTREGAR"
    print(f"\n  {r['votu']}{bandera}")
    print(f"    nombre            {nombre}")
    print(f"    motivo / clase    {r['motivo']} / {r['clase_final']}")
    print(f"    consenso srna     {c['consenso']}   rpm_max {c['rpm_max']}")
    print(f"    por libreria      replica {c['replica']} · virus+piRNA {c['virus_mas_piRNA']} · "
          f"EVE {c['EVE_NIRVS']} · ambiguo {c['ambiguo']} · bajo umbral {c['bajo_umbral']}"
          f"  (de {c['muestras']})")
    print(f"    ARN largo         longitud {r['longitud']} · frac_en_genoma {r['frac_en_genoma']} · "
          f"cobertura {r['cobertura']} en {r['muestra_cob']}")
    print(f"                      prof_media {r['prof_media']} · CV {r['cv_profundidad']} · "
          f"n_muestras_mapeo {r['n_muestras_mapeo']}")
    if r.get("prediccion_largo"):
        print(f"    prediccion_largo  {r['prediccion_largo']}")

alarmas = [r for r in disc
           if any(a in ((r.get("nombre_virus") or cons[r["votu"]].get("nombre_virus") or "").lower())
                  for a in ALARMA)]
print(f"\n  con nombre de arbovirus: {len(alarmas)}"
      + ("  -> RESOLVER" if alarmas else "  -> ninguno, no bloquea"))

# ------------------------------- 3. se distinguen tambien en el ARN largo?
print("\n== frac_en_genoma: los discordantes frente a los EVE_candidato que si dan EVE ==")
print("   (si los rescatados son virus, deberian tener MENOS frac_en_genoma)")


def resumen(rotulo, vals):
    vals = [v for v in vals if v == v]          # fuera los NaN
    if not vals:
        print(f"  {rotulo:<34} sin datos")
        return
    vals.sort()
    med = statistics.median(vals)
    print(f"  {rotulo:<34} n={len(vals):<4} min {vals[0]:.3f}  mediana {med:.3f}  max {vals[-1]:.3f}")


for motivo in MOTIVOS:
    rescatados = [num(r["frac_en_genoma"]) for r in obj
                  if r["motivo"] == motivo
                  and cons.get(r["votu"], {}).get("consenso", "") in REPLICA]
    comoeve = [num(r["frac_en_genoma"]) for r in obj
               if r["motivo"] == motivo
               and cons.get(r["votu"], {}).get("consenso", "") == "EVE_NIRVS"]
    print(f"\n  {motivo}")
    resumen("con perfil de replica", rescatados)
    resumen("con perfil de EVE", comoeve)

print("\n  referencia: el umbral de EVE_genoma es frac_en_genoma >= 0.50")
