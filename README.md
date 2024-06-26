# Anchor test of population continuity
<img src="anchor_icon.png" alt="anchor" width="200">
This repository contains a collection of notebooks and scripts used in a project investigating genetic continuity among temporally distributed genomes.

Included in Simulations/ are notebooks for running msprime simulations for various demographic models and to calculate anchor statistics and f-statistics on the resulting dataset.

Included in Anchor_test/ is a pipeline for running the Anchor test; a novel test for genetic continuity, with empirical datasets.

## Author
James McKenna (james.andrew.mckenna@hi.no) \
https://github.com/jammc313/Genetic-continuity

## Dependencies include:
* Python 3.12
* jupyterlab 4.1.0
* msprime 1.3.1
* numpy 1.26.4
* scikit-allel 1.3.8

## Getting started
* Create the environment: Run the following command to create a Conda environment using the downloaded YAML file:
```
conda env create -f anchor_env.yml
```
* Activate the new environment: 
```
conda activate anchor_env
```
