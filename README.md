# Computational analysis of chromatin architecture | Tachm et al. (2026)

**Talianidis Laboratory — research code and analysis workflows**

This repository contains analysis scripts, configuration examples, and computational modeling utilities associated with the chromatin architecture studies of **Tachm et al. (2026)**. The workflows cover Hi-C and promoter capture Hi-C (PCHi-C) processing, chromatin domain and loop analysis, transcription factor perturbation comparisons, and data-constrained three-dimensional (3D) chromatin modeling.

> **Publication information:** Full article title, author list, journal, DOI, and data accession numbers will be added when confirmed. The repository name identifies the associated study; it does not by itself establish publication status.

## Overview

The code supports investigation of chromatin organization across mouse liver developmental stages (**E14, E18, P0, and P60**) and comparisons involving **HNF4A** and **CEBPA** perturbations. Available scripts illustrate methods for processing contact matrices, evaluating topologically associating domains (TADs), visualizing contact enrichment near protein-binding sites, and constructing coarse-grained chromosome models.

**Reference genome:** Mouse mm10 (where applicable to the provided analysis workflows). Verify the assembly and chromosome naming conventions of your own inputs before running commands.

## Repository contents

| File | Purpose |
| --- | --- |
| `HiCupconfig.txt` | Example configuration for HiCUP processing. |
| `align_build_matrix_hic_explorer_example.sh` | Alignment and Hi-C matrix construction commands. |
| `hicexplorer.sh` | HiCExplorer commands for matrix merging, normalization, correction, TAD detection, comparisons, and visualization. |
| `agregate_plots_around_CHip_binding_sites.sh` | Aggregate Hi-C contact analysis around ChIP-binding regions. |
| `bam2chicago_V02.sh` | Processing commands for CHiCAGO-oriented PCHi-C analysis. |
| `chicago_Man_Settings1kb.txt` | Example CHiCAGO analysis settings. |
| `makeDesignFiles_py3.py` | Design-file preparation utility for capture Hi-C analysis. |
| `commands2_for_align_PCHiC.txt` | PCHi-C alignment workflow notes/commands. |
| `Rscript2TAD_cmp.R`, `Rscript_TAD_compare.r.txt`, `TAD Compare.txt` | TAD comparison scripts and notes. |
| `insulation_scores.ipynb` | Notebook for insulation-score analysis. |
| `track_insulation_Heatmap_Stages.iniversion*.txt` | Configuration variants for insulation and heatmap tracks. |
| `chr8_P60_biophysics_CONTACTS_LOOPS10K_ONLY.py` | Coarse-grained chromosome 8 modeling with Hi-C contacts and loop constraints. |
| `chr8_P60_biophysics_H3K27ac_NO_LAMIN.py` | Modeling variant incorporating H3K27ac without lamina constraints. |
| `chr8_P60_biophysics_H3K27ac_WITH_LAMIN_3Dpdb_create.py` | Modeling variant incorporating H3K27ac and lamina-related constraints. |
| `make_interactive_genomeflow_3d.py` | Interactive HTML comparison of chromosome models from PDB files. |
| `digest_genome.py`, `fixchrom.pl` | Supporting genome-processing utilities. |

The filenames above reflect the repository contents at the time this README was drafted. Some `.txt` files contain executable commands or scripts; inspect them before use.

## Main analysis workflow

### 1. Hi-C processing and matrix generation

The repository includes examples for read processing, contact matrix generation, and downstream HiCExplorer analysis. The `hicexplorer.sh` workflow includes commands such as:

```bash
hicMergeMatrixBins --matrix E14.h5 --numBins 2 --outFileName E14.20kb.h5

hicFindTADs \
  -m E14.norm.20kb.0_1.cor.h5 \
  --outPrefix E14.20kbv2 \
  --numberOfProcessors 8 \
  --correctForMultipleTesting fdr \
  --threshold 0.001
```

These are **illustrative commands**, not a standalone end-to-end pipeline: the input files must already exist and some upstream preparation steps are not bundled with the repository.

### 2. PCHi-C and CHiCAGO

Scripts and configuration files support preparation and analysis of promoter capture Hi-C data. Consult `bam2chicago_V02.sh`, `makeDesignFiles_py3.py`, and `chicago_Man_Settings1kb.txt` for the relevant input assumptions and parameters.

### 3. TADs, insulation, and aggregate contacts

R scripts, notebooks, and shell commands cover TAD comparisons, insulation profiles, and aggregate contact analysis near ChIP-seq binding sites. Input matrices and genomic interval files must use compatible reference assemblies and chromosome naming conventions.

### 4. Data-constrained 3D chromosome modeling

The chromosome 8 Python scripts explore coarse-grained polymer representations using combinations of Hi-C contacts, loop constraints, H3K27ac-associated information, and lamina-associated information. Model outputs can include PDB representations and geometric statistics such as radius of gyration, end-to-end distance, and convex-hull volume, depending on the script.

**Interpretation caveat:** These are computationally inferred, data-constrained chromatin conformations, not uniquely measured single-cell chromosome structures. Model-dependent geometric quantities should not be interpreted as direct physical measurements without calibration and validation.

### 5. Interactive comparison of 3D models

`make_interactive_genomeflow_3d.py` compares PDB models for HNF4A and P60 conditions, generates alignment/displacement visualizations, and reports structural metrics in an interactive HTML interface.

Example command adapted from the script's documented usage:

```bash
python make_interactive_genomeflow_3d.py \
  --hnf4a HNF4A_model.pdb \
  --p60 P60_model.pdb \
  --chrom chr8 \
  --start 104860000 \
  --end 106860000 \
  --out-prefix chr8_HNF4A_P60_interactive
```

Replace the example PDB names with your own model files. The resulting HTML visualization uses 3Dmol.js from an external URL and may require an internet connection.

## Installation and dependencies

Clone the repository:

```bash
git clone https://github.com/Orsalia/Talianidis_lab_Tachm_et_al_2026.git
cd Talianidis_lab_Tachm_et_al_2026
```

Create an optional Python environment for the modeling and visualization utilities:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install numpy pandas scipy
```

Other workflow components require tools installed separately, depending on the script:

- **Hi-C / PCHi-C:** HiCUP, HiCExplorer, CHiCAGO and their respective dependencies.
- **R analysis:** R and the packages called by the specific TAD comparison scripts.
- **Notebook:** Jupyter and the packages imported by `insulation_scores.ipynb`.
- **Supporting tools:** Perl, shell utilities, and alignment tools as required by the selected workflow.

**Note:** This repository does not currently provide a validated, fully pinned software environment. Inspect each script and record tool versions for reproducible execution.

## Input data and reproducibility

The repository primarily distributes **code and configuration examples**, not the complete sequencing datasets or precomputed contact matrices. Depending on the analysis, you may need FASTQ/BAM files, `.h5` or `.cool` contact matrices, genomic BED/bedGraph tracks, ChIP-seq peak files, and/or chromosome-model PDB files.

For reproducibility:

1. Use consistent genome builds and chromosome labels across all inputs.
2. Confirm contact matrix resolution and normalization before comparing conditions.
3. Record command-line parameters, software versions, and random seeds for simulations.
4. Preserve biological replicate identities and avoid treating modeled conformations as independent experimental replicates.
5. Refer to the associated manuscript and deposited data for the final validated analysis choices when those resources become available.

## Citation

If you use this repository, please cite the associated **Tachm et al. (2026)** publication once its bibliographic details are available, and reference this code repository:

> Talianidis Laboratory. *Talianidis_lab_Tachm_et_al_2026*. GitHub repository: https://github.com/Orsalia/Talianidis_lab_Tachm_et_al_2026

**Manuscript citation:** *To be completed with the verified title, authors, journal, year, and DOI.*

## Contact

For questions or issues related to the code, please open an issue in this GitHub repository.

## License

No license is specified here. Unless a license is added to the repository, do not assume the code is distributed under an open-source license.
