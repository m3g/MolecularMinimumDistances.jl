[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://m3g.github.io/MolecularMinimumDistances.jl/stable)
[![Tests](https://img.shields.io/badge/build-passing-green)](https://github.com/m3g/MolecularMinimumDistances.jl/actions)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

# MolecularMinimumDistances

## Documentation

Go to: [https://m3g.github.io/MolecularMinimumDistances.jl](https://m3g.github.io/MolecularMinimumDistances.jl)

## Introduction

This package computes the minimum distance between *molecules*, which are represented as arrays of coordinates in two or three dimensions. 

To understand the utility and purpose of this package, consider the image below:

![nearest.png](/docs/src/assets/nearest.png)

Here, there is one *blue* molecule, with 6 atoms, and several *red* molecules, with 2 atoms each. The package has identified which are the molecules of the *red* set that have at leat one atom within a cutoff from the atoms of the *blue* molecule, and annotated the corresponding atoms and the distances.

## Features

- Fast and parallel [cell-list approach](https://github.com/m3g/CellListMap.jl), to compute minimum-distance for thousands, or millions of atoms. 
- General periodic boundary conditions supported. 
- Advanced mode for in-place calculations, for non-allocating iterative calls (for analysis of MD trajectories, for example).
- Modes for the calculation of minimum-distances in sets of molecules.

# Citation

If this package was useful, please cite the article describing the main
algorithms on which it is based:

L. Martínez, CellListMap.jl: Efficient and customizable cell list implementation for calculation of pairwise particle properties within a cutoff. Computer Physics Communications, 279, 108452, 2022. https://doi.org/10.1016/j.cpc.2022.108452
