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

## Most typical use: Understanding solvation

This package was designed as the backend for computing [minimum distance distribution functions](http://m3g.github.io/ComplexMixtures.jl), which are useful for understanding solute-solvent interactions when the molecules have complex shapes. 

The most typical scenario is that of a protein, or another macromolecule, in a box of solvent. For example, here we download a frame of a protein which was simulated in a mixture of water and TMAO: 

```julia
julia> using PDBTools

julia> system = MolecularMinimumDistances.download_example()
   Array{Atoms,1} with 62026 atoms with fields:
   index name resname chain   resnum  residue        x        y        z  beta occup model segname index_pdb
       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  1.00     1    PROT         1
       2  HT1     ALA     A        1        1  -10.048  -15.427   -5.569  0.00  0.00     1    PROT         2
       3  HT2     ALA     A        1        1   -9.488  -13.913   -5.295  0.00  0.00     1    PROT         3
                                                       ⋮ 
   62024  OH2    TIP3     C     9339    19638   13.485   -4.534  -34.438  0.00  1.00     1    WAT2     62024
   62025   H1    TIP3     C     9339    19638   13.218   -3.647  -34.453  0.00  1.00     1    WAT2     62025
   62026   H2    TIP3     C     9339    19638   12.618   -4.977  -34.303  0.00  1.00     1    WAT2     62026
```

Next, we extract the protein coordinates, and the TMAO coordinates:

```julia
julia> protein = coor(system,"protein")
1463-element Vector{SVector{3, Float64}}:
 [-9.229, -14.861, -5.481]
 [-10.048, -15.427, -5.569]
 [-9.488, -13.913, -5.295]
 ⋮
 [6.408, -12.034, -8.343]
 [6.017, -10.967, -9.713]

julia> tmao = coor(system,"resname TMAO")
2534-element Vector{SVector{3, Float64}}:
 [-23.532, -9.347, 19.545]
 [-23.567, -7.907, 19.381]
 [-22.498, -9.702, 20.497]
 ⋮
 [13.564, -16.517, 12.419]
 [12.4, -17.811, 12.052]
```

The system was simulated with periodic boundary conditions, with sides in this frame of `[83.115, 83.044, 83.063]`, and this information will be provided to the minimum-distance computation.

Finally, we find all the TMAO molecules having at least one atom closer than 12 Angstroms to the protein, using the current package (TMAO has 14 atoms):

```julia
julia> box = Box([83.115, 83.044, 83.063], 12.);

julia> list = minimum_distances(tmao, protein, 14, box)
181-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(false, 0, 0, Inf)
 MinimumDistance{Float64}(false, 0, 0, Inf)
 MinimumDistance{Float64}(false, 0, 0, Inf)
 ⋮
 MinimumDistance{Float64}(false, 0, 0, Inf)
 MinimumDistance{Float64}(true, 2526, 97, 9.652277658666891)

julia> count(x -> x.within_cutoff, list)
33
```

Thus, 33 TMAO molecules are within the cutoff distance from the protein, and the distances can be used to study the solvation of the protein.

## Performance

This package exists because this computation is fast. For example, let us choose the water molecules instead, and benchmark the time required to compute these set of distances:
```julia
julia> water = coor(system,"resname TIP3")
58014-element Vector{SVector{3, Float64}}:
 [-28.223, 19.92, -27.748]
 [-27.453, 20.358, -27.476]
 [-27.834, 19.111, -28.148]
 ⋮
 [13.218, -3.647, -34.453]
 [12.618, -4.977, -34.303]

julia> using BenchmarkTools

julia> @btime minimum_distances($water, $protein, 3, $box);
  4.726 ms (2748 allocations: 11.82 MiB)
```

To compare, a naive algorithm to compute the same thing takes:

```julia
julia> @btime MolecularMinimumDistances.naive_md($water, $protein, 3, $box);
  911.580 ms (2 allocations: 604.36 KiB)
```

And the computation can be made faster and in-place using the more advanced interface that allows preallocation of all necessary arrays:

```julia
julia> using CellListMap

julia> list = init_list(water, i -> mol_indices(i,3)); # 3 atoms per molecule

julia> cl = CellList(water,protein,box);

julia> list_threaded = [ copy(list) for i in 1:nbatches(cl) ];

julia> minimum_distances!(i -> mol_indices(i,3), list, box, cl; list_threaded = list_threaded)
19338-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(false, 0, 0, Inf)
 MinimumDistance{Float64}(false, 0, 0, Inf)
 MinimumDistance{Float64}(false, 0, 0, Inf)
 ⋮
 MinimumDistance{Float64}(true, 58011, 383, 10.24673074692606)
 MinimumDistance{Float64}(false, 0, 0, Inf)
```

Allocations occur only for the launching of multiple threads:

```julia
julia> @btime minimum_distances!(
           $(i -> mol_indices(i,3)), 
           $list, $box, $cl; 
           parallel = false
        );
  12.723 ms (0 allocations: 0 bytes)

julia> @btime minimum_distances!(
          $(i -> mol_indices(i,3)), 
          $list, $box, $cl;
          list_threaded = $list_threaded,
          parallel = true # default
        );
  3.473 ms (76 allocations: 9.75 KiB)
```

## Details of the illustration

The initial illustration here consists of a toy solute-solvent example, where the solute is a approximatelly hexagonal molecule, and the solvent is composed by 40 diatomic molecules. The toy system is built as follows:

```julia
using MolecularMinimumDistances, StaticArrays
T = SVector{2,Float64}
# x will contain the "solvent", composed by 40 diatomic molecules
x = T[]
cmin = T(-20,-20)
for i in 1:40
    v = cmin .+ 40*rand(T)
    push!(x, v)
    theta = 2pi*rand()
    push!(x, v .+ T(sin(theta),cos(theta)))
end
# y will contain the "solute", composed by an approximate hexagonal molecule
y = [ T(1,1), T(1,-1), T(0,-1.5), T(-1,-1), T(-1,1), T(0,1.5) ]
```

Next, we compute the minimum distances between each molecule of `x` (the solvent)
and the solute. In the input we need to specify the number of atoms of each molecule
in `x`, and the cutoff up to which we want the distances to be computed:

```julia
julia> list = minimum_distances(x,y,2,10.0)
40-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(true, 2, 3, 1.0764931248364737)
 MinimumDistance{Float64}(false, 0, 0, Inf)
 MinimumDistance{Float64}(false, 0, 0, Inf)
 ⋮
 MinimumDistance{Float64}(true, 74, 5, 7.899981412729262)
 MinimumDistance{Float64}(false, 0, 0, Inf)
```

The output is a list of `MinimumDistance` data structures, one for each molecule in `x`. The `true` indicates that a distance smaller than the cutoff was found, and for these the indexes of the atoms in `x` and `y` associated are reported, along with the distance between them.

In this example, from the 40 molecules of `x`, eleven had atoms closer than the cutoff to some
atom of `y`:
```julia
julia> count(x -> x.within_cutoff, list)
11
```

We have an auxiliary function to plot the result, in this case where the "atoms" are bi-dimensional:

```julia
using Plots
import MolecularMinimumDistances: plot_md!
p = plot(lims=(-20,20),framestyle=:box,grid=false,aspect_ratio=1)
plot_md!(p, x, 2, y, 6, list, y_cycle=true)
```
will produce the illustration plot above, in which the nearest point between the two sets is identified.

# Citation

If this package was useful, please cite the article describing the main
algorithms on which it is based:

L. Martínez, CellListMap.jl: Efficient and customizable cell list implementation for calculation of pairwise particle properties within a cutoff. Computer Physics Communications, 279, 108452, 2022. https://doi.org/10.1016/j.cpc.2022.108452



