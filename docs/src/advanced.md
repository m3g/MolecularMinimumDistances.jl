# Advanced usage

## Array preallocation

The advanced usage of this package exposes the interface of `CellListMap`, such that it is possible to iterate through sets varying coordinates and perhaps box types without need for reallocating lists of neighbors every time.  

Basically, preallocation is needed for:

- The cell lists.
- The resulting lists of minimum-distances. 

And for parallel runs (probably the most common ones):

- Auxiliary arrays for storing threaded versions of the cell lists.
- Auxiliary arrays for storing threaded versions of the resulting lists of minimum-distances. 

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

However, more general indexing can be used. For instance, let us suppose that the 9 atoms of the `x` array of coordinates above belong to `2` molecules, with `4` and `5` indexes each. Then, we could define, for example:

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

## Example

The example below illustrates the usage of the interface described above: 

```julia
using MolecularMinimumDistances
using CellListMap # required
using PDBTools # example

function iterate_lists(nsteps, water, protein, box)
    cl = CellList(protein,water,box)
    # Function that given the index of the atom, returns the index of the molecule. Here, all molecules are similar, and we use the standard `mol_indices` function:
    water_index(i) = mol_indices(i,3)
    # Initialize output array; i -> mol_indices(i,3) returns the molecule index for each water molecule
    list = init_list(water, water_index)
    # Initalize auxiliary arrays for multi-threading: The number of batches of multi-threading is obtained from a property of the cell-lists (see the CellListMap documentation for additional information and options).
    aux_cl = CellListMap.AuxThreaded(cl)
    list_threaded = [ copy(list) for _ in 1:nbatches(cl) ]
    # Now we are ready to iterate over many possible different sets of coordinates. Here the coordinates or the box won't change, but they could:
    for i in 1:nsteps
        # box = Box(...) the box could be redefined here.
        # Update the cell lists: water and protein coordinates, and box, could change here. 
        cl = UpdateCellList!(water, protein, box, cl, aux_cl)
        # Update the `list` array:
        minimum_distances!( water_index, list, box, cl, list_threaded = list_threaded) 
        # Perform whatever futher analysis using the `list` of minimum distances.
    end
end
```

Running the above example shows that allocations are small for each iteration:
```julia-repl
julia> using BenchmarkTools

julia> system = MolecularMinimumDistances.download_example();

julia> water = coor(system, "water"); # water coordinates

julia> protein = coor(system, "protein"); # protein coordinates

julia> box = Box([84.48, 84.48, 84.48], 12.);

julia> @btime iterate_lists(10, $water, $protein, $box)
  1.113 s (21471 allocations: 36.42 MiB)

julia> @btime iterate_lists(100, $water, $protein, $box)
  13.568 s (37946 allocations: 38.36 MiB)
```

In fact, if we opted to run the calculation in serial, with:

```julia
function iterate_lists_serial(nsteps, water, protein, box)
    cl = CellList(protein,water,box, parallel=false)
    water_index(i) = mol_indices(i,3)
    list = init_list(water, water_index)
    for i in 1:nsteps
        cl = UpdateCellList!(water, protein, box, cl, parallel=false)
        minimum_distances!(water_index, list, box, cl, parallel=false) 
        # Perform whatever futher analysis using the `list` of minimum distances.
    end
end
```

We can see that the updating of the cell lists and the computation of the minimum-distance lists is completely allocation free, such that the loop is allocation free:

```julia-repl
julia> @btime iterate_lists_serial(10, $water, $protein, $box)
  7.989 s (1797 allocations: 14.42 MiB)

julia> @btime iterate_lists_serial(20, $water, $protein, $box)
  16.716 s (1797 allocations: 14.42 MiB)
```

Thus, using a parallelization scheme at a upper level can be also an alternative.

!!! note
    The `Box` constructor is type-unstable because of the different possible dimensions and unit cells. This instability does not propagate to inner functions and is, thus, benign, but to avoid completely allocations in a loop using a predefined box, a function barrier may be required. 


