# User guide

Here the usage of the functions that allocate the list of distances will be described. Different running modes are available depending on the expected output.

## Installation and loading

To install `MolecularMinimumDistances`, first download and install Julia (1.6 or greater) from [https://julialang.org/downloads/](https://julialang.org/downloads). Install and run it. Then, use:

```julia-repl
julia> import Pkg; Pkg.add("MolecularMinimumDistances")

julia> using MolecularMinimumDistances
```

## Example input files

The input of the following tests can be obtained with:

```julia-repl
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
The system consists of a protein (with 1463 atoms), solvated by 181 TMAO molecules (with 14 atoms each), 19338 water molecules, and some ions. 

This is a snapshot of a simulation which was performed with cubic periodic boundary conditions, with a box side of `84.48` Angstrom. We will use periodic boundary conditions in the examples. 

The coordinates of each of the types of molecules can be extracted from the `system` array of atoms with (using `PDBTools`):

```julia-repl
julia> xprot = coor(system,"protein")
1463-element Vector{StaticArrays.SVector{3, Float64}}:
 [-9.229, -14.861, -5.481]
 [-10.048, -15.427, -5.569]
 [-9.488, -13.913, -5.295]
 ⋮
 [6.408, -12.034, -8.343]
 [6.017, -10.967, -9.713]

julia> tmao = coor(system,"resname TMAO")
2534-element Vector{StaticArrays.SVector{3, Float64}}:
 [-23.532, -9.347, 19.545]
 [-23.567, -7.907, 19.381]
 [-22.498, -9.702, 20.497]
 ⋮
 [13.564, -16.517, 12.419]
 [12.4, -17.811, 12.052]

julia> water = coor(system,"water")
58014-element Vector{StaticArrays.SVector{3, Float64}}:
 [-28.223, 19.92, -27.748]
 [-27.453, 20.358, -27.476]
 [-27.834, 19.111, -28.148]
 ⋮
 [13.218, -3.647, -34.453]
 [12.618, -4.977, -34.303]
```

Using these vectors of coordinates, we will illustrate the use of the current package.

## Two sets of molecules

The simplest usage consists of finding for each molecule of one set the atoms of the other set which are closer to them. For example, if we want the atoms of the proteins which are closer to each TMAO molecule (14 atoms), within a cutoff of `12.0` Angstroms, we do:

### Without periodic boundary conditions

```julia-repl
julia> list = minimum_distances(tmao, protein, 14, 12.0)
181-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(false, -1, -1, Inf)
 MinimumDistance{Float64}(false, -1, -1, Inf)
 MinimumDistance{Float64}(false, -1, -1, Inf)
 ⋮
 MinimumDistance{Float64}(false, -1, -1, Inf)
 MinimumDistance{Float64}(true, 2526, 97, 9.652277658666891)
```

The list obtained contains, for each molecule of TMAO, a structure of type `MinimumDistance`, containing:
1. A boolean marker, which is `true` if some protein atom was found within the desired cutoff, `false` otherwise (field `x.within_cutoff`).
2. The index of the TMAO atom (field `x.i`).
3. The index of the protein atom (field `x.j`).
4. The distance between these atoms (field `x.d`).

For instance, the number of molecules of TMAO having a protein atom within the cutoff are, here:
```julia-repl
julia> count(x -> x.within_cutoff, list)
33
```

For each molecule of water, he have, similarly:

```julia-repl
julia> list = minimum_distances(water,protein,3,12.)
19338-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(false, -1, -1, Inf)
 MinimumDistance{Float64}(false, -1, -1, Inf)
 MinimumDistance{Float64}(false, -1, -1, Inf)
 ⋮
 MinimumDistance{Float64}(true, 58011, 383, 10.24673074692606)
 MinimumDistance{Float64}(false, -1, -1, Inf)

julia> count(x -> x.within_cutoff, list)
2251
```

### With periodic boundary conditions

The example simulation was performed with cubic periodic boundary conditions. Let us provide the box information now. We will exemplify with the calculation of the nearest atoms of the water molecules. The interface here is that define by the `Box` constructor of `CellListMap.jl`, described in detail [here](https://m3g.github.io/CellListMap.jl/stable/pbc/). General periodic boundary conditions are supported. 

The box here is cubic, and we need to provide to the `Box` constructor the sides and the cutoff:

```julia-repl
julia> box = Box([84.48, 84.48, 84.48], 12.)
Box{CellListMap.OrthorhombicCell, 3, Float64, Float64, 9}
  unit cell matrix = [ 84.48, 0.0, 0.0; 0.0, 84.48, 0.0; 0.0, 0.0, 84.48 ]
  cutoff = 12.0
  number of computing cells on each dimension = [9, 9, 9]
  computing cell sizes = [12.06857142857143, 12.06857142857143, 12.06857142857143] (lcell: 1)
  Total number of cells = 729
```

And the `minimum_distance` function is called with the `box` instead of the `cutoff`:

```julia-repl
julia> list = minimum_distances(water,protein,3,box)
19338-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(false, -1, -1, Inf)
 MinimumDistance{Float64}(false, -1, -1, Inf)
 MinimumDistance{Float64}(false, -1, -1, Inf)
 ⋮
 MinimumDistance{Float64}(true, 58011, 383, 10.24673074692606)
 MinimumDistance{Float64}(false, -1, -1, Inf)
```








