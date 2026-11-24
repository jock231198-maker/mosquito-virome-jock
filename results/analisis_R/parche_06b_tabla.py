#!/usr/bin/env python3
"""
Dos cosas sobre 06b:

  1. El pie gana la aclaracion de que los dos ejes usan criterios distintos.
     Sin ella, los tres taxones con el naranja a la derecha del azul
     (Aedes aegypti totivirus 6 frente a 5, CFAV 3 frente a 2) parecen un
     error: ensamblar en mas pools de los que detecta el mapeo suena
     imposible. No lo es, porque el ensamblaje cuenta pools donde se armo
     cualquier contig de la especie y el mapeo cuenta pools donde la vOTU
     representante llega al 70 % de cobertura.

  2. Un chunk nuevo detras de la figura que CUENTA la discordancia en vez de
     dejar que se lea a ojo. Es la cifra que va al texto.

Tiene que estar en el mismo directorio que figuras.qmd, y despues de
parche_figuras_rev.py.
No escribe nada si cualquier aserto falla. Deja figuras.qmd.bak2.

    python3 parche_06b_tabla.py
"""

import os
import shutil
import sys

QMD = "figuras.qmd"

if not os.path.exists(QMD):
    sys.exit(f"ERROR: no encuentro {QMD} en {os.getcwd()}")

with open(QMD, "rb") as fh:
    crudo = fh.read()
crlf = b"\r\n" in crudo
nl = "\r\n" if crlf else "\n"
texto = crudo.decode("utf-8").replace("\r\n", "\n")

# ------------------------------------------------------------------ 1. pie
VIEJO_PIE = (
    '         "cases more directly, which is why this is a review signal and not a criterion of its own.")) +'
)
NUEVO_PIE = (
    '         "cases more directly, which is why this is a review signal and not a criterion of its own. The two\\n",\n'
    '         "axes use different criteria, so a few taxa sit the other way round: assembly counts a pool where any\\n",\n'
    '         "contig of that species was built, while mapping counts a pool where the representative vOTU reaches\\n",\n'
    '         "70 % of its length covered, and a partial contig does not.")) +'
)

n = texto.count(VIEJO_PIE)
assert n == 1, (f"pie de 06b: esperaba 1 coincidencia, encontre {n}. "
                f"Corriste antes parche_figuras_rev.py?")
texto = texto.replace(VIEJO_PIE, NUEVO_PIE)

# -------------------------------------------------------------- 2. chunk
ANCLA = 'p6b\n```\n\n```{r table-rate, eval = HAY_VOTUS}'

CHUNK = '''p6b
```

```{r tabla-discordancia, eval = HAY_VOTUS}
# La cifra que va al texto. RESUMEN_min500.md 5.6 fijo la regla -- mapeo alto
# con ensamblaje <= 3 casi siempre es un EVE o un acierto celular -- y aqui se
# cuenta en vez de leerse de la figura. El corte de mapeo >= 4 es deliberado y
# mas laxo que el 18/22 del 5.6: la pregunta no es cuantos casos extremos hay,
# sino si el patron es de unos pocos taxones o de una categoria entera.
disc <- comp |>
  mutate(
    grupo = if_else(
      str_detect(as.character(especie_fig), "-like|-related"),
      "Relative (-like / -related)", "Named species"),
    discordante = ensamblaje <= 1 & mapeo >= 4)

disc |>
  group_by(Group = grupo) |>
  summarise(Taxa = n(),
            `Assembled in ≤ 1 pool, mapped in ≥ 4` = sum(discordante),
            `% of the group` = paste0(round(100 * sum(discordante) / n()), " %"),
            .groups = "drop") |>
  knitr::kable()
```

The most extreme cases, assembled in a single pool:

```{r tabla-discordancia-extremos, eval = HAY_VOTUS}
disc |>
  filter(ensamblaje <= 1, mapeo >= 8) |>
  arrange(desc(mapeo)) |>
  transmute(Virus = as.character(especie_fig),
            `Assembled in` = ensamblaje,
            `Detected in` = mapeo) |>
  knitr::kable()
```

```{r table-rate, eval = HAY_VOTUS}'''

n = texto.count(ANCLA)
assert n == 1, f"ancla del chunk: esperaba 1 coincidencia, encontre {n}"
texto = texto.replace(ANCLA, CHUNK)

shutil.copy2(QMD, QMD + ".bak2")
with open(QMD, "w", encoding="utf-8", newline=nl) as fh:
    fh.write(texto)

print("  pie de 06b ampliado")
print("  2 chunks nuevos insertados tras la figura")
print(f"  finales de linea: {'CRLF' if crlf else 'LF'} (preservados)")
print(f"  copia de seguridad: {QMD}.bak2")
print("ok")
