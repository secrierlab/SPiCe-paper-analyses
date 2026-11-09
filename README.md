# SPiCe-paper-analyses

This repository contains the analysis code accompanying the [SPiCe](https://github.com/secrierlab/SPiCe) manuscript. 

It has five components: 

(i) `merfish_analyses` implements the graph neural network (GNN) framework. For each sample of the MERFISH mouse cortex dataset, it builds a spatial k-nearest-neighbour cell graph (k = 20, Euclidean-distance edge weights). It then holds out one cell type at a time and trains a three-layer graph convolutional network (`PyTorch` Geometric) on one-hot cell-type features to predict which remaining cells sit next to the removed type. Evaluation is inductive, on held-out tissue sections, and performance is measured by ROC-AUC. Predictions are interpreted with Integrated Gradients (Captum) for node-level (cell-type) importance, tested against a permutation baseline with Mann–Whitney U tests and multiple-testing correction. `GNNExplainer` is used for edge-level (cell-type-pair) importance. Code is in `src/` (config, dataset, models, training, explanations), with `analysis.ipynb` as the entry point and `base_env_gnn.yml` as the conda environment.

(ii) `spatial_regression` contains Python and R notebooks and scripts (`gwr.ipynb`, `gwr_r.ipynb`, `YY_05_00_spatial_regression.R`, `ZZ_05_spatial_regression.R`) that model spatially varying relationships in Xenium breast cancer data. They regress an EMT hallmark score on tumour-microenvironment and CNV-derived predictors, such as distances to neighbouring cell populations and principal components. The notebooks fit spatial error models (spreg), geographically weighted regression (`mgwr`) and multiscale `GWR` (`GWmodel` in R), with VIF/correlation filtering and variable-specific bandwidth selection.

(iii) `figure_plotting` holds the notebook and conda environment that reproduce the paper’s GNN performance, node/edge importance and network figures from the saved outputs. The conda environment file `figures_paper_gnn.yml` specifies the dependencies. The other notebooks show analyses regarding the generalisability, contender frameworks, etc.

(iv) `supplementary_analyses` holds the notebooks behind the supplementary results, organised by topic:
- `simulations/`: simulated data and a representative simulation with new state labels.
- `simvi/`: simVI run and archetype analysis.
- `xenium_breast_5k/`: Xenium 5k breast cancer processing, EMT plasticity hallmark scoring and spatial regression.
- `xenium_lung_5k/`: Xenium 5k lung cancer processing, EMT plasticity hallmark scoring, centroid analysis and spatial regression.

(v) `xenium_breast_v1_analyses` contains the Xenium v1 breast cancer notebooks for EMT plasticity hallmark scoring and CNA input preparation.

The scripts and associated data can also been found at [Zenodo](https://doi.org/10.5281/zenodo.23884754).

## Citation

> Withnell E, Celik C, Secrier M. **SPiCe: integrative spatial modelling of cellular plasticity using graph neural networks and geostatistics.** Genome Biology (2026).
