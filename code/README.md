Please replace the contents of this file with relevant instructions for your repository or remove this file entirely.

This directory would generally contain source code files that contain the core code to implement the method and various utility/auxiliary functions.

Scripts/code files that execute the overall workflow to carry out an analysis and generate results for the manuscript might be placed in the main directory.



## Julia code for Shapley MFM example

This folder contains a Julia project used to fit mixture-of-finite-mixture
(MFM) models to the Shapley galaxy velocities:

- `Project.toml`, `Manifest.toml`: Julia environment.
- `shapley_mfm.ji`: precompiled Julia image used for the original runs.
- `data/processed/julia_run/MFM/`: CSV files of posterior draws of the number
  of clusters and other derived quantities.

For reproducibility of the JASA paper, you only need the precomputed CSV
files in `data/processed/julia_run/MFM/`, which are used by
`code/examples/01_shapley_gmm.R`. Re-generating the MFM fits would require
the original Julia source script (not included here).



## STRUCTURE runs (brook trout)

### Software
- STRUCTURE (version: <fill in>)
- Binary location: ./structure  (or on PATH)

### Inputs
- Parameter files: code/structure/mainparams, code/structure/extraparams
- Genotype file: data/raw/brooktrout/brooktrout.txt
  (Note: mainparams must point to this file, or add -i if your STRUCTURE build requires it.)

### Outputs
- Written to: data/processed/structure_run/
- Naming: results_K{K}_rep{rep} (plus any STRUCTURE suffixes)

### Command used
```bash
mkdir -p data/processed/structure_run

for K in {1..10}; do
  for rep in {1..20}; do
    SEED=$((1000 * K + rep))
    ./structure \
      -K $K \
      -D $SEED \
      -m code/structure/mainparams \
      -e code/structure/extraparams \
      -o data/processed/structure_run/results_K${K}_rep${rep}
  done
done
