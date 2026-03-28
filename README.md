# HIT_VALIDATION

Pipeline de validacion rigurosa de hits de screening virtual con DOCK6. Evalua cada hit en base a meritos absolutos (scoring, cobertura de zonas, interacciones PLIP) para generar evidencia de decision: seleccion como template Pharmit o compra. Datos pre-computados de referencia (ej. UDX de `reference_docking`) se importan como contexto informativo en tablas de salida, pero nunca se usan como input de ningun paso del pipeline.

 Parte del ciclo iterativo de descubrimiento:

```
reference_docking → molecular_docking → HIT_VALIDATION → Pharmit → reference_docking [iteracion 2] → ...
```

- **reference_docking**: Establece la linea base cristalografica (UDX en iteracion 1, mejor hit en iteracion 2+)
- **molecular_docking**: Screening virtual masivo contra librerias de moleculas
- **hit_validation** (ESTE PROYECTO): Validacion rigurosa de hits — evidencia per-molecula para template Pharmit o compra

## Diferencias clave vs reference_docking

| Aspecto | reference_docking | hit_validation |
|---------|-------------------|----------------|
| Sujeto | Ligando cristalografico unico (UDX) | Multiples hits de screening (10-20) |
| Coordenadas | Cristalograficas del PDB HETATM | Generadas por docking (sin pose cristal) |
| Prep. ligando | Sin antechamber (preserva geometria cristal) | Antechamber requerido (AM1-BCC + Sybyl) |
| Box margin | 1.0 A (tight, mantiene ligando cerca del cristal) | 5.0 A (espacio de exploracion para hits) |
| Metodo docking | Rigid (pose cristal se mantiene) | Flex (anchor-and-grow, busqueda conformacional) |
| Ligando ref. | UDX ES el sujeto | Ninguno (hits evaluados por merito propio; datos UDX importados como contexto) |
| Input PLIP | Complejo cristalino (receptor + UDX co-cristalizado) | Complejo dockeado (receptor + mejor pose) |
| GB/SA minimize | false (geometria cristal es sagrada) | true (poses dockeadas se benefician de relajacion) |
| Output | Perfil de referencia energetico | Reporte de decision per-molecula con argumentos de compra |
| 06a/06b | Farmacoforo del UDX cristal | Farmacoforo del mejor hit validado (iteracion 2) |

## Requisitos

| Dependencia | Instalacion | Notas |
|---|---|---|
| Python 3.9-3.12 | conda | |
| RDKit, OpenBabel | conda | Via `environment.yaml` |
| AmberTools >= 23.0 | conda | Provee antechamber para cargas AM1-BCC |
| PLIP | conda | `conda install -c conda-forge plip` |
| ionization | pip | `git+https://github.com/benjaminprieto/ionization.git` |
| **DOCK6** | Manual | Licencia academica gratuita de [UCSF](https://dock.compbio.ucsf.edu/DOCK_6/index.htm) |
| **ChimeraX** | Manual | [Descargar](https://www.cgl.ucsf.edu/chimerax/download.html) — necesario para 00b |

## Instalacion

```bash
# 1. Clonar
git clone https://github.com/benjaminprieto/hit_validation.git
cd hit_validation

# 2. Crear entorno conda
conda env create -f environment.yaml
conda activate hit_validation_env

# 3. Instalar paquete en modo editable
pip install -e ".[dev]"

# 4. Verificar dependencias externas
bash check_dependencies.sh

# 5. Correr tests
pytest
```

**Alternativa**: Compartir `molecular_docking_env` (ya tiene AmberTools + ionization):
```bash
conda activate molecular_docking_env
pip install -e ".[dev]"
```

## Estructura del proyecto

```
hit_validation/
├── 01_src/hit_validation/          Core modules (logica, sin CLI)
│   ├── m00_preparation/
│   │   ├── ligand_preparation.py       00a — antechamber batch (AM1-BCC + Sybyl)
│   │   ├── receptor_preparation.py     00b — ChimeraX DockPrep → rec_charged.mol2
│   │   └── binding_site_definition.py  00d — binding site por residuos/coordenadas
│   ├── m01_docking/
│   │   ├── grid_generation.py          01b — DMS → spheres → box 5.0A → grids
│   │   ├── dock6_runner.py             01c — DOCK6 docking (flex, grid_score_primary)
│   │   ├── footprint_rescore.py        01d — fps_primary per-residuo vdW+ES
│   │   ├── gbsa_rescore.py             01f — GB/SA Hawkins (con minimizacion)
│   │   └── score_collector.py          01e — scores → Excel/CSV
│   ├── m03_interaction_analysis/
│   │   └── plip_interaction_analysis.py 03a — PLIP en poses dockeadas (batch)
│   ├── m04_dock6_analysis/
│   │   ├── footprint_analysis.py       04b — consenso per-residuo + zonas + ref. context
│   │   └── footprint_rescoring.py      04b rescore
│   ├── m06_pharmit/
│   │   ├── pharmit_pharmacophore.py    06a — footprint+PLIP → queries Pharmit
│   │   └── pharmit_zone_selector.py    06b — mapeo de features por zonas
│   └── m07_decision_report/
│       └── decision_report.py          07a — reporte de decision (HTML + CSV)
├── 02_scripts/                     13 CLI wrappers (patron --config --campaigns)
├── 03_configs/                     13 YAML (uno por modulo)
├── 04_data/campaigns/              Inputs por campana
│   └── example_hit_validation/     Template para nuevas validaciones
├── 05_results/                     Outputs (gitignored)
├── tests/
│   └── test_pipeline.py            21 smoke tests
├── run_pipeline.sh                 Orquestador del pipeline
├── check_dependencies.sh
├── environment.yaml                Entorno conda (hit_validation_env)
├── pyproject.toml                  Build config (hit_validation v1.0.0)
├── CLAUDE.md                       Guia tecnica completa (ingles)
└── README.md                       Este archivo (espanol)
```

## Pipeline completo

```
00a → 00b → 00d → 01b → 01c → 01d → 01f → 01e → 03a → 04b → [06a → 06b] → 07a
```

### 00a — Preparacion de ligandos (REESCRITO)

Prepara N hits de screening usando antechamber con cargas AM1-BCC y atom types Sybyl.

**Fallback de 3 niveles:**
1. `antechamber -c bcc -at sybyl` (AM1-BCC completo + Sybyl)
2. `antechamber -c gas -at sybyl` (Gasteiger + Sybyl, rapido, sin QM)
3. `obabel` conversion (ultimo recurso, atom typing degradado)

**Pipeline por molecula:** SDF→mol2 → protonacion al pH (paquete ionization) → antechamber → validacion

**Restricciones criticas:**
- `atom_type` DEBE ser `sybyl`. GAFF2 causa que DOCK6 caiga silenciosamente a rigid docking con cero descriptores.
- El paquete `ionization` protona al `docking_ph` ANTES de antechamber.

**Parametros (03_configs/00a_ligand_preparation.yaml):**
- `atom_type: "sybyl"`
- `charge_method: "bcc"`
- `timeout_per_molecule: 300`

**Output:** `{name}/{name}.mol2`, `preparation_status.csv`, `preparation_summary.txt`

### 00b — Preparacion del receptor (IDENTICO)

ChimeraX DockPrep → `rec_charged.mol2` con cargas AMBER ff14SB y atom types Sybyl.

**Optimizacion:** Si el receptor ya fue preparado en reference_docking, se puede apuntar a `precomputed_mol2` en campaign_config.

### 00d — Definicion del binding site (se salta si hay grids precomputados)

Define el binding site por lista de residuos o coordenadas, y recorta el receptor para grid generation. **Se salta** si `precomputed_grids_dir` esta definido en campaign_config (flujo preferido — reusar grids de reference_docking o molecular_docking).

**Metodos:** residues (default) | coordinates | reference_ligand (legacy)

### 01b — Generacion de grids (se salta si hay grids precomputados)

DMS surface → esferas → seleccion → box → grids de energia. **Se salta** si `precomputed_grids_dir` esta definido.

**Config:** `box_margin: 5.0` (era 1.0 en reference_docking). Los hits necesitan libertad conformacional.

### 01c — DOCK6 docking (SOLO CONFIG)

DOCK6 con grid_score como funcion de scoring primaria.

**Cambios de config:**
- `search_method: "flex"` (era "rigid") — hits necesitan busqueda conformacional anchor-and-grow
- `pruning_clustering_cutoff: 2.0` — **CRITICO**: el default de DOCK6 (100 A) fusiona TODAS las poses en 1 cluster
- `num_preclustered_conformers: 1000` (era 500) — mas muestreo para hits
- `timeout_per_molecule: 600`

**Templates DOCK6:** `DOCK6_FLEX_TEMPLATE` y `DOCK6_RIGID_TEMPLATE`. Ninguno contiene parametros `_secondary`.

### 01d — Footprint rescore (auto-referencia)

Re-scoring con `footprint_similarity_score_primary` — descomposicion per-residuo vdW + ES. Cada hit usa su propia pose como referencia (auto-referencia) — no necesita UDX mol2.

### 01f — GB/SA Hawkins rescore (SOLO CONFIG)

Solvatacion implicita como rescore rigido separado con `gbsa_hawkins_score_primary=yes`.

**Cambio de config:** `minimize: true` (era false). Poses dockeadas se benefician de relajacion simplex antes del scoring GB/SA.

### 01e — Score collection (IDENTICO)

Parsea todos los scored mol2 y colecta scores en CSV/Excel.

**Campos capturados:** Grid_Score, Grid_vdw_energy, Grid_es_energy, Internal_energy_repulsive, DOCK_Rotatable_Bonds, Formal_Charge, Heavy_Atoms, Molecular_Weight, HBond_Acceptors, HBond_Donors

### 03a — Analisis PLIP (MODIFICADO para batch de poses dockeadas)

Corre PLIP en las mejores poses dockeadas (no en el complejo cristalino).

**Pipeline batch por molecula:** Extraer mejor pose → convertir a PDB → combinar con receptor → PLIP → JSON/CSV

**Parametros:**
- `pose_source: "01c_dock6_run"`
- `pose_selection: "best_score"` (Grid_Score mas bajo)

**Output por molecula:** `interactions.json`, `interactions_summary.csv`, `complex_for_plip.pdb`
**Consolidado:** `plip_summary_all.csv`

### 04b — Analisis de footprint (multi-molecula + contexto de referencia opcional)

Descomposicion per-residuo con remapeo de numeracion (mol2 secuencial → PDB original). Si `reference_context.footprint_csv` esta definido en campaign_config, los datos pre-computados se importan como fila de referencia en todas las tablas (nunca computados por este pipeline).

**ZONE_DEFINITIONS** (sistema XT1, propiedad del receptor, hardcoded):

| Zona | Residuos | Drug-like |
|------|----------|-----------|
| phosphate | ARG598, LYS599 | No |
| xylose | TRP392, TRP495, TYR565, SER575 | Si |
| ribose | HIS335, VAL333, THR390 | Si |
| uracil | ASP361, ARG363 | Si |
| catalytic | GLU529 | Si |

**Outputs:**
- `footprint_per_molecule.csv` — matriz residuo x molecula
- `residue_consensus.csv` — consenso per-residuo
- `hit_vs_udx_comparison.csv` — delta per-residuo (hit - referencia) por cada hit (solo si reference_context disponible)
- `subpocket_coverage.csv` — que zonas cubre cada hit (+ fila de referencia si disponible)
- `ranking_by_zone.csv` — hits rankeados por energia por sub-pocket (+ fila de referencia, no rankeada)
- `pharmacophore_residues.json` — residuos con >80% frecuencia de contacto
- `binding_site_zones.html` — reporte HTML interactivo de zonas

### 06a/06b — Pharmit (IDENTICOS, condicionales)

Generan queries de farmacoforo Pharmit. Solo corren si `pharmit_template_molecule` esta definido en campaign_config.

**Estrategias:** xylose, uracil, combined, analogues, druglike

### 07a — Reporte de decision (NUEVO)

Genera reporte HTML consolidado + CSV ranking con evidencia de validacion per-molecula. Si `reference_context` esta disponible, los datos de referencia aparecen como fila informativa en todas las tablas — pero NO participan en ranking ni clasificacion.

**Logica de clasificacion (umbrales absolutos — no relativos a UDX):**
- **"Strong Pharmit template candidate"**: cubre >= 2 zonas drug-like con energia significativa (< -0.5 kcal/mol por zona)
- **"Purchase candidate"**: Grid Score <= -30.0 Y cubre zona xylose O uracil
- **"Weak candidate"**: en otro caso

**Residuos clave** para comparacion detallada: TRP392, HIS335, GLU529, ASP494, ARG598

**Parametros (03_configs/07a_decision_report.yaml):**
- `zone_energy_cutoff: -0.5` (umbral absoluto de energia para cobertura de zona)
- `min_zones_for_template: 2`
- `grid_score_purchase_threshold: -30.0`

**Output:** `decision_report.html` (una tarjeta por molecula), `decision_summary.csv`

## Uso

### 1. Crear campana

```bash
cp -r 04_data/campaigns/example_hit_validation 04_data/campaigns/mi_campana
```

Colocar los archivos mol2/SDF de los hits en `mi_campana/ligands/`.
Editar `campaign_config.yaml`:

```yaml
campaign_id: "NAMIKI_top20_pH63"
target: {name: "XT1", pdb_id: "6EJ7"}
receptor: {pdb: "receptor/XT1_6EJ7.pdb", chain: "A", precomputed_mol2: null}
docking_ph: 6.3
ligands_dir: "ligands/"
grids:
  precomputed_grids_dir: null  # Definir para saltar 00d + 01b
  binding_site:
    method: "residues"
    residues: ["TRP392", "TRP495", "TYR565", "SER575", "HIS335", "ASP361", "ARG363", "GLU529", "ARG598", "LYS599"]
reference_context:  # Datos importados como fila de referencia en 04b/07a (informativo)
  label: "UDX (crystallographic)"
  scores_csv: null     # Path a dock6_scores.csv de reference_docking
  plip_json: null      # Path a interactions.json de reference_docking
  footprint_csv: null  # Path a residue_consensus.csv de reference_docking
pharmit_template_molecule: null  # Definir DESPUES de revisar 07a
```

### 2. Correr pipeline

```bash
# Pipeline completo
bash run_pipeline.sh 04_data/campaigns/mi_campana/campaign_config.yaml

# Desde un modulo especifico
bash run_pipeline.sh 04_data/campaigns/mi_campana/campaign_config.yaml 01c

# Modulo individual
python 02_scripts/01c_dock6_run.py \
    --config 03_configs/01c_dock6_run.yaml \
    --campaigns 04_data/campaigns/mi_campana/campaign_config.yaml
```

### 3. Revisar resultados

Despues de correr el pipeline:
- `07a_decision_report/decision_report.html` — reporte visual por molecula
- `07a_decision_report/decision_summary.csv` — ranking con recomendaciones
- `04b_footprint_analysis/subpocket_coverage.csv` — cobertura por sub-pocket (+ fila de referencia si disponible)
- `04b_footprint_analysis/hit_vs_udx_comparison.csv` — delta per-residuo vs referencia (si reference_context disponible)
- `04b_footprint_analysis/binding_site_zones.html` — mapa interactivo de zonas

### 4. Seleccionar template Pharmit (opcional)

Si un hit muestra buena cobertura de zonas, editarlo como template para iteracion 2:

```yaml
# En campaign_config.yaml:
pharmit_template_molecule: "NAMIKI_0042"
```

Y re-correr desde 06a:
```bash
bash run_pipeline.sh 04_data/campaigns/mi_campana/campaign_config.yaml 06a
```

## Estructura de outputs

```
05_results/{campaign_id}/
├── 00a_ligand_preparation/     {name}/{name}.mol2, preparation_status.csv
├── 00b_receptor_preparation/   rec_charged.mol2, rec_noH.pdb
├── 00d_binding_site/           rec_noH_site.pdb
├── 01b_grid_generation/        ligand.nrg, ligand.bmp, spheres_ligand.sph
├── 01c_dock6_run/              {name}/dock6.in, {name}_scored.mol2
├── 01d_footprint_rescore/      {name}/{name}_fps_scored.mol2
├── 01f_gbsa_rescore/           {name}/{name}_gbsa_scored.mol2
├── 01e_score_collection/       dock6_scores.csv, dock6_scores.xlsx, best_poses/
├── 03a_plip_analysis/          {name}/interactions.json, plip_summary_all.csv
├── 04_dock6_analysis/
│   └── 04b_footprint_analysis/ footprint_per_molecule.csv, residue_consensus.csv,
│                                hit_vs_udx_comparison.csv, subpocket_coverage.csv,
│                                ranking_by_zone.csv, binding_site_zones.html
├── 06a_pharmit/                pharmacophore_{strategy}.json (condicional)
├── 06b_pharmit_zones/          pharmit_{strategy}.json (condicional)
└── 07a_decision_report/        decision_report.html, decision_summary.csv
```

## Arquitectura DOCK6.13

DOCK6.13 ignora silenciosamente todas las funciones de scoring `_secondary` durante flex docking. Por lo tanto:
- **01c** usa solo `grid_score_primary` (paso de docking)
- **01d** usa `footprint_similarity_score_primary` (rescore separado)
- **01f** usa `gbsa_hawkins_score_primary` (rescore separado)

Cada funcion de scoring corre como un paso independiente de rescore rigido con `orient_ligand=no`.

Los template strings de DOCK6 (en dock6_runner.py, footprint_rescore.py, gbsa_rescore.py, footprint_rescoring.py) contienen **cero** parametros `_secondary`. Esto esta verificado por tests.

## Notas tecnicas

**DOCK6 y el limite de 80 caracteres.** Los programas Fortran de DOCK6 (sphgen, showbox) truncan paths a ~80 chars. El pipeline usa symlinks y filenames cortos automaticamente.

**`pruning_clustering_cutoff` DEBE ser 2.0 A.** El default de DOCK6 (100 A) fusiona todas las poses en 1 cluster, produciendo solo 1 pose de salida. Configurado en `03_configs/01c_dock6_run.yaml`.

**`atom_type` DEBE ser sybyl.** GAFF2 causa que DOCK6 caiga silenciosamente a rigid docking con cero descriptores moleculares. Aplicado en 00a.

**PLIP protona internamente.** NO pre-protonar la pose del ligando antes de alimentar a PLIP.

**UDX no es input del pipeline.** UDX nunca se dockea, puntua, ni analiza en hit_validation. Datos pre-computados de `reference_docking` se importan via `reference_context` en campaign_config como contexto informativo en tablas de salida. El pipeline funciona identicamente con o sin ellos.

**Binding site** se define por lista de residuos o coordenadas (no por un ligando mol2 de referencia). Alternativamente, usar `precomputed_grids_dir` para saltar 00d + 01b completamente.

**Rutas del servidor:**
- DOCK6: `/opt/dock6/`
- ChimeraX: `/usr/bin/chimerax-daily`
- Working dir: `/home/bprieto/hit_validation/`