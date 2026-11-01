# diamond_R — DIAMOND sobre los 6,004 contigs de la unión

Tablas y figuras de las tres corridas de DIAMOND blastx sobre **los mismos 6,004 contigs**
de la unión MERIDA (>= 1,000 pb, 22 pools, 5 ensambladores, sin `spades_meta` ni `trinity_rescate`).

| Corrida | Base | Parámetros | Fichero en la farm |
|---|---|---|---|
| Base viral propia | `refs/viral_custom/20260925` (RefSeq viral + GenBank con huésped artrópodo, 680,459 proteínas; solo virus, con taxonomía) | `--very-sensitive -e 1e-3 -k 25` | `results/MERIDA/diamond_TODOS_vcustom.tsv` |
| RVDB v32.0 | `U-RVDBv32.0-prot_unique` (solo virus, con retroelementos y virus gigantes; sin taxid) | `--very-sensitive`, e <= 1e-3 | `results/MERIDA/diamond_TODOS_rvdb.tsv` |
| nr (jr46) | `nr_tax.dmnd` de jr46 (707 M secuencias; mismo contenido que el nr de feb-2024; virus + celulares) | `--sensitive -e 1e-5 -k 25` + árbitro | `results/MERIDA/contig_annot/contigs_union_MERIDA_nr_jr46/diamond.tsv` |

**El árbitro** (`arbiter_nr.py`) decide por contig con los hits de nr: *viral* si el mejor bitscore
viral es >= 1.1 x el mejor celular; *ambiguo* si están a menos del 10 %; si no, *celular_<grupo>*
según la taxonomía. Las bases solo virales no pueden decidir esto: su mejor hit siempre es un virus.

## Flujo

```
diamond_R/
  congelar_diamond.sh   EN LA FARM: valida y copia las tablas a datos_diamond_vN (con MANIFEST md5)
  diamond_taxa.py       lo llama congelar_diamond.sh: grupo/realm/familia/género por hit (taxdump)
  R/diamond.R           lectura, tabla por contig, clases y LISTAS EDITABLES de curación
  diamond_tablas.R      genera tablas/ y figuras/
  datos_diamond_vN/     tablas congeladas (no se tocan: si cambia algo, vN+1)
```

**1. En la farm**

```bash
cd ~/jkvirome/jm79/mosquito-virome-jock && git pull
bash results/analisis_R/diamond_R/congelar_diamond.sh check        # solo comprueba
bsub -Is -q normal -n 1 -M 8000 -R "select[mem>8000] rusage[mem=8000]" \
     bash results/analisis_R/diamond_R/congelar_diamond.sh v1     # congela (carga el taxdump)
git add results/analisis_R/diamond_R/datos_diamond_v1 && git commit -m "datos_diamond_v1 congelados" && git push
```

**2. En casa**

```bash
cd ~/JK/git/mosquito-virome-jock && git pull
cd results/analisis_R/diamond_R
Rscript diamond_tablas.R                  # la última datos_diamond_vN
Rscript diamond_tablas.R datos_diamond_v1 # una concreta
```

Paquetes: `readr dplyr tidyr stringr ggplot2 scales forcats` (+ `writexl` y `ragg`, opcionales). R >= 4.1.

## Salida

| Tabla | Contenido |
|---|---|
| `T1_runs` | las tres corridas: base, parámetros, hits, contigs con hit, e-value máximo observado |
| `T2_verdicts` | los 6,004 contigs por veredicto del árbitro |
| `T3_viraldb_vs_nr` | contigs con hit en cada base, repartidos por veredicto de nr |
| `T4_virus_catalogue` | un virus por fila: familia, pools, contigs, identidad aa, tipo (virus / candidato a EVE / fago) |
| `T5_per_pool` | por pool: contigs por veredicto, nº de virus y de candidatos a EVE |
| `T6_viral_contigs` | cada contig viral con su mejor hit en las tres bases |
| `T7_corrections` | correcciones manuales de veredicto aplicadas |

Todo junto en `tablas/diamond_tables.xlsx`, con una hoja README de definiciones.

| Figura | Qué muestra |
|---|---|
| `D1_databases_vs_arbiter` | contigs con hit en cada base, coloreados por el veredicto de nr |
| `D2_identity_by_verdict` | identidad del mejor hit en las bases virales, viral frente a celular |
| `D3_verdict_by_assembler` | veredicto por ensamblador |
| `D4_virus_by_pool_heatmap` | virus x pool, por clase (Known / Related / Divergent / EVE candidate) |
| `D5_identity_vs_length` | identidad aa frente a longitud de los contigs virales |
| `D6_nudivirus_contigs` | los contigs del nudivirus, con la proteína de su mejor hit |

## Curación (en `R/diamond.R`, arriba)

- `CORREGIR_VEREDICTO`: veredictos corregidos a mano (por defecto, el contig de M1_S1 cuyo mejor hit sin taxid es de *Aedes*: 204 -> 203 virales).
- `AGRUPAR`: mejores hits que son el mismo virus (los nudivirus de M71_S19; *Aedes aegypti anphevirus* = *Aedes anphevirus*).
- `FAMILIAS_EVE` y `UMBRAL_EVE`: qué se marca como candidato a EVE hasta pasar el cribado contra AaegL5.
- `UMBRAL_CONOCIDO` (95) y `UMBRAL_PARIENTE` (70): clases por identidad aa.

## Para no equivocarse al presentar

- Contigs no son virus: el mismo virus aparece en varios pools y ensambladores. La prevalencia es `pools`.
- Los e-values no son iguales: 1e-3 en las bases virales, 1e-5 en nr.
- Los candidatos a EVE no cuentan como virus hasta el cribado.
- Estas corridas son sobre contigs; las de los 748 vOTUs (`votu_annot/`) son otro subconjunto.
