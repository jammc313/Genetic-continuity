# Anchor test of population continuity
<p align="center">
<img src="anchor_icon.png" alt="anchor" width="200">
</p>
This repository contains a collection of notebooks and scripts used in a project investigating genetic continuity among temporally distributed genomes.

Included in Simulations/ are scripts for running msprime simulations for various demographic models.

Simulation_set_1: demonstrates that polymorphic datasets simulated under a demographic model of population-specific drift, and a model of historical admixture, can produce very similar summary statistics (Patterson's f2 and PCA coordinates).

Simulation_set_2: simulate data under a model of population specific drift and a model of historical admixture. Plot FST matrices for resulting data and calculate the anchor statistic for each dataset. Includes an investigation of the sensitivity of the approach to limited data quality and quantity. These scripts also contain functions which take "perfect" simulated and introduce a genotyping error and limited coverages, and then down-sample the resulting dataset to subsets of anchor sites. 

Simulation_set_3: Scripts simulating under a model of historical admixture, and studying the power of the anchor approach to estimate the proportion of admixture from an unsampled "Ghost" population. Three parameters are allowed to vary: the admixture proportion, the divergence time of the admixing population, and the amount of drift separating the two individuals chosen as anchors.  

Included in Anchor_test/ is a pipeline for running the Anchor test on empirical data; a novel test for discriminating population continuity from historical admixture.

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
