# analisis_R — figuras del viroma MERIDA

Las figuras se hacen en local (Mac) a partir de tablas **congeladas** que vienen de la farm por git.
Nada aquí lee de la farm directamente.

```
analisis_R/
  congelar_datos.sh     se corre EN LA FARM: copia las tablas a datos_vN con MANIFEST (md5)
  metadata_MERIDA.tsv   metadatos de los 22 pools (editado a mano)
  datos_v1/             tablas congeladas; una versión no se modifica, se crea la siguiente
  R/leer.R              lectores (comprueban columnas) y cálculos: presencia, tasa de infección
  R/tema.R              paleta, tema y guardar() -> PDF + PNG
  figuras.qmd           una sección por paso del pipeline
  figuras/              salida: NN_nombre.pdf y .png
```

## Flujo

En la farm (lote MERIDA activo, tras `integrar_arbitro.py`):

```bash
cd ~/jkvirome/jm79/mosquito-virome-jock
bash results/analisis_R/congelar_datos.sh v1
git add results/analisis_R/datos_v1 && git commit -m "datos_v1 congelados" && git push
```

En el Mac:

```bash
git pull
cd results/analisis_R
quarto render figuras.qmd                     # usa datos_v3 (por defecto)
quarto render figuras.qmd -P datos:datos_v1   # otra versión
```

Sale `figuras.html` (informe con tablas) y los PDF/PNG en `figuras/`.

## Instalación (una vez)

```bash
brew install --cask quarto
Rscript -e 'install.packages(c("readr","dplyr","tidyr","ggplot2","scales","forcats","stringr","patchwork","knitr","rmarkdown"), repos = "https://cloud.r-project.org")'
```

## Criterios usados en las figuras

- **Presencia**: ≥ 70 % del genoma de la vOTU cubierto a ≥ 1×; una especie está presente si lo está alguna de sus vOTUs.
- **Abundancia**: profundidad media por millón de pares recortados (no por pares no huésped: M36 y M74 tienen ~60 % no huésped que no es virus).
- **Virus**: `es_virus == "si"` en `TABLA_FINAL_votus_arb.tsv` (identificado, divergente, pariente débil, viral solo en nr). Los EVE y lo celular quedan fuera.
- **Tasa de infección**: pools de 15; MIR y MLE (Chiang y Reeves) por 1,000 hembras, IC de Clopper-Pearson transformado.
