# Advanced usage

## System build and update

If the molecular minimum distances will be computed many times for similar systems, it is possible
to construct the system and update its properties. The use of the interface of `CellListMap.PeriodicSystems`
is required (requires `CellListMap` version `0.7.24` or greater). 

For example, let us build one system with a protein and water:

```julia-repl
julia> using MolecularMinimumDistances, PDBTools

julia> system = MolecularMinimumDistances.download_example();

julia> protein = coor(system, "protein");

julia> water = coor(system, "water");
```

We now build the `CrossPairs`  type of system, instead of calling the `minimum_distances` function directly:

```julia-repl
julia> sys = CrossPairs(
           xpositions=water, # solvent
           ypositions=protein, # solute
           xn_atoms_per_molecule=3,
           cutoff=12.0,
           unitcell=[84.48, 84.48, 84.48]
       )
CrossPairs system with:

Number of atoms of set x: 58014
Number of molecules in set x: 19338
Number of atoms of target structure y: 1463
Cutoff: 12.0
unitcell: [84.48, 0.0, 0.0, 0.0, 84.48, 0.0, 0.0, 0.0, 84.48]
```

Now `sys`  contains the necessary arrays for computing the list of minimum distances. We use now the
`minimum_distances!`  function (with the `!`), to update that list:
```julia-repl
julia> minimum_distances!(sys)
19338-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(false, 0, 0, Inf)
 MinimumDistance{Float64}(false, 0, 0, Inf)
 MinimumDistance{Float64}(false, 0, 0, Inf)
 ⋮
 MinimumDistance{Float64}(true, 58011, 383, 10.24673074692606)
 MinimumDistance{Float64}(false, 0, 0, Inf)
```

The system can be now updated: the positions, cutoff, or unitcell can be modified, with the 
following interfaces:

## Updating positions

To update the positions, modify the `sys.xpositions` (or `ypositions`)  array. We will
boldy demonstrate this by making the first atom of the `x` set to be close to the first
atom of the protein, and recomputing the distances:
```julia-repl
julia> using StaticArrays

julia> sys.xpositions[2] = sys.ypositions[1] + SVector(1.0,0.0,0.0);

julia> minimum_distances!(sys)
19338-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(true, 2, 4, 0.9202923448556931)
 MinimumDistance{Float64}(false, 0, 0, Inf)
 MinimumDistance{Float64}(false, 0, 0, Inf)
 ⋮
 MinimumDistance{Float64}(true, 58011, 383, 10.24673074692606)
 MinimumDistance{Float64}(false, 0, 0, Inf)
```

## Updating the cutoff, unitcell and parallel flag

The `cutoff`, `unitcell` and `parallel` data of the `sys` objects can be modified 
directly. For example:
```julia-repl
julia> sys
CrossPairs system with:

Number of atoms of set x: 58014
Number of molecules in set x: 19338
Number of atoms of target structure y: 1463
Cutoff: 15.0
unitcell: [100.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 100.0]

julia> sys.cutoff = 10.0
10.0

julia> sys.unitcell = [84.4, 84.4, 84.4]
3-element Vector{Float64}:
 84.4
 84.4
 84.4

julia> sys.parallel = false
false

julia> sys
CrossPairs system with:

Number of atoms of set x: 58014
Number of molecules in set x: 19338
Number of atoms of target structure y: 1463
Cutoff: 10.0
unitcell: [84.4, 0.0, 0.0, 0.0, 84.4, 0.0, 0.0, 0.0, 84.4]
```

!!! note
    It is not possible to update the `unitcell` from a Orthorhombic to a general Triclinic cell. If the system
    will be Triclinic at any moment, the `unitcell` must be initialized with the full matrix instead of a 
    vector of sides.

## Index of molecules

Additionally, the low level interface allows the definition of more general groups of particles, in the sense that "molecule" can have different number of atoms in the same set. Therefore, one needs to provide *a function* that returns the index of the molecule of each atom, given the index of the atom. 

Briefly, if a set of atoms belong to molecules of the same number of atoms, one can compute the index of each molecule using
```julia
mol_indices(i,n) = div((i - 1), n) + 1
```
where `i` is the atom index in the array of coordinates, and `n` is the number of atoms per molecule. This is the default assumed in the basic interface, and can be called with:
```julia-repl
julia> using StaticArrays

julia> x = rand(SVector{3,Float64},9); # 3 water molecules

julia> mol_indices(2,3) # second atom belongs to first molecule
1

julia> mol_indices(4,3) # fourth atom belongs to second molecule
2
```

Typically, as we will show, this function will be used for setting up molecule indexes.

However, more general indexing can be used. For instance, let us suppose that the 9 atoms of the `x` array of coordinates above belong to `2` molecules, with `4` and `5` atoms each. Then, we could define, for example:

```julia-repl
julia> my_mol_indices(i) = i <= 4 ? 1 : 2
my_mol_indices (generic function with 1 method)

julia> my_mol_indices(4)
1

julia> my_mol_indices(5)
2
```

Since the function can close-over an array of molecular indexes, the definition can be completely general, that is:

```julia-repl
julia> molecular_indexes = [ 1, 3, 3, 2, 2, 1, 3, 1, 2 ];

julia> my_mol_indices(i) = molecular_indexes[i]
my_mol_indices (generic function with 1 method)

julia> my_mol_indices(1)
1

julia> my_mol_indices(5)
2
```

In summary, this function that given the index of the atom returns the index of the corresponding molecule must be provided in the advanced interface, and typically will be just a closure around the number of atoms per molecule, using the already available `mol_indices` function. 

#### Example

Let us mix water and TMAO molecules in the same set, and use a general function to compute the indices of the molecules of each atom: 

```julia-repl
julia> system = MolecularMinimumDistances.download_example();

julia> protein = coor(system, "protein");

julia> tmao_and_water = select(system, "resname TMAO or resname TIP3")
   Array{Atoms,1} with 60548 atoms with fields:
   index name resname chain   resnum  residue        x        y        z  beta occup model segname index_pdb
    1479    N    TMAO     A        1      120  -23.532   -9.347   19.545  0.00  1.00     1    TMAO      1479
    1480   C1    TMAO     A        1      120  -23.567   -7.907   19.381  0.00  1.00     1    TMAO      1480
    1481   C2    TMAO     A        1      120  -22.498   -9.702   20.497  0.00  1.00     1    TMAO      1481
                                                       ⋮ 
   62024  OH2    TIP3     C     9339    19638   13.485   -4.534  -34.438  0.00  1.00     1    WAT2     62024
   62025   H1    TIP3     C     9339    19638   13.218   -3.647  -34.453  0.00  1.00     1    WAT2     62025
   62026   H2    TIP3     C     9339    19638   12.618   -4.977  -34.303  0.00  1.00     1    WAT2     62026

julia> findfirst(at -> at.resname == "TIP3", tmao_and_water)
2535
```

Thus, the `tmao_and_water` atom array has two different types of molecules, TMAO with 14 atoms, and water with 3 atoms. 
The first atom of a water molecule is atom `2535` of the array. We extract the coordinates of the atoms with:
```julia-repl
julia> coor(tmao_and_water)
60548-element Vector{SVector{3, Float64}}:
 [-23.532, -9.347, 19.545]
 [-23.567, -7.907, 19.381]
 [-22.498, -9.702, 20.497]
 ⋮
 [13.218, -3.647, -34.453]
 [12.618, -4.977, -34.303]
```

 And now we define a function that, given the index of the atom, returns the molecule to which it belongs:

```julia-repl
julia> function mol_indices(i) 
           if i < 2535 # TMAO (14 atoms per molecule) 
               div(i-1,14) + 1 
           else # water (3 atoms per molecule)
               mol_indices(2534) + div(i-2534-1,3) + 1
           end
       end
mol_indices (generic function with 1 method)
```

The function above computes the molecular indices for TMAO in the standard way, and computes the water 
molecular indices by first summing the molecule index of the last TMAO molecule, and subtracting from the
atomic index of water the last index of the last TMAO atom. We can test this: 

```julia-repl
julia> mol_indices(14) # last atom of first TMAO
1

julia> mol_indices(15) # first atom of second TMAO
2

julia> mol_indices(2534) # last atom of last TMAO
181

julia> mol_indices(2535) # first atom of first water
182

julia> mol_indices(2537) # last atom of first water
182

julia> mol_indices(2538) # first atom of second water
183
```

With this function, we can construct the system using it instead of the `xn_atoms_per_molecule` integer variable,
to obtain the solvation of the protein by both TMAO and water in a single run:

```julia-repl
julia> sys = CrossPairs(
           xpositions=solvent, # solvent
           ypositions=protein, # solute
           xmol_indices = mol_indices,
           cutoff=12.0,
           unitcell=[84.48, 84.48, 84.48]
       )
CrossPairs system with:

Number of atoms of set x: 60548
Number of molecules in set x: 19519
Number of atoms of target structure y: 1463
Cutoff: 12.0
unitcell: [84.48, 0.0, 0.0, 0.0, 84.48, 0.0, 0.0, 0.0, 84.48]
```

As we can see, the number of molecules is correct (the sum of the number of water and tmao
molecules). And the list of minimum distances will retrive the information of the closest
protein atom to all solvent molecules of the set:

```julia-repl
julia> minimum_distances!(sys)
19519-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(false, 0, 0, Inf)
 MinimumDistance{Float64}(false, 0, 0, Inf)
 MinimumDistance{Float64}(false, 0, 0, Inf)
 ⋮
 MinimumDistance{Float64}(true, 60545, 383, 10.24673074692606)
 MinimumDistance{Float64}(false, 0, 0, Inf)
```




