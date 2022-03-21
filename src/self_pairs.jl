#
# This file containst the functions for single-sets, that is for those cases where
# the list of minimum-distances is between the molecules of a single component.
#

function update_list!(
    i, j, d2,
    mol_index::F,
    list::AbstractVector{<:MinimumDistance}
) where F <: Function
    imol = mol_index(i)
    jmol = mol_index(j)
    if imol != jmol 
        d = sqrt(d2)
        if d < list[imol].d
            list[imol] = MinimumDistance(i,j,d)
        end
        if d < list[jmol].d
            list[jmol] = MinimumDistance(j,i,d)
        end
    end
    return list
end

"""

```
function minimum_distances!(
    mol_index::F,
    x_list::AbstractVector{<:MinimumDistance},
    box, cl;
    parallel = true
) where F <: Function
```

Compute the list of minimum distances given the precomputed cell lists, auxiliary vectors, and 
lists of indexes of molecules of each atom. Should be preferred for multiple calls of this function,
with the outer update of the cell lists. The `mol_index` function must return, for each atom index,
the corresponding molecular index. 

# Examples

```julia-repl
julia> using CellListMap, StaticArrays

julia> x = [ rand(SVector{2,Float32}) for _ in 1:100 ];

julia> box = Box([1,1],0.2f0);

julia> cl = CellList(x,box);

julia> minimum_distances!(i -> mol_index(i,5), init_list(x,5), box, cl)
20-element Vector{MinimumDistance{Float32}}:
 MinimumDistance{Float32}(79, 0.02097016f0)
 MinimumDistance{Float32}(32, 0.018027343f0)
 ⋮
 MinimumDistance{Float32}(55, 0.013549047f0)
 MinimumDistance{Float32}(61, 0.015638554f0)
```


"""
function minimum_distances!(
    mol_index::F,
    x_list::AbstractVector{<:MinimumDistance},
    box, cl;
    parallel = true
) where F <: Function
    map_pairwise!(
        (x, y, i, j, d2, list) -> update_list!(i, j, d2, mol_index, list),
        x_list, box, cl; 
        parallel=parallel,
        reduce = reduce_list!,
    )
    return x_list
end

"""

```
function minimum_distances(
    x, n_atoms_per_molecule_x, box;
    parallel = true
)
```

Compute the list of minimum distances given the vector of atomic coordinates `x`, 
and the number of atoms of the molecules, the box (of type `CellListMap.Box`). This is the simplest
interface, but it allocates a new list. For iterative computations use `minimum_distances!`
instead. 

## Example

```
julia> x = [ rand(SVector{3,Float64}) for _ in 1:100 ];

julia> box =  Box([1,1,1],0.2);

julia> x_list = minimum_distances(x, 5, box);

julia> x_list
20-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(51, 0.048400669109669024)
 MinimumDistance{Float64}(23, 0.08296253778070845)
 ⋮
 MinimumDistance{Float64}(71, 0.04853375996234046)
 MinimumDistance{Float64}(56, 0.03089814908795506)

```

Each entry of the output lists contains, for each atom (in `x` for example), the index of the atom of 
some other molecule in `x` that is closer to it, and the distance between these atoms. If no atom is found within
the cutoff of a given atom `MinimumDistance(-1,+inf)` will be returned. 

A more general interface, if the molecules do not have the same number of atoms, allows providing a function 
that returns the index of each molecule instead of the number of atoms of each molecule:

```
function minimum_distances(
    x, mol_index::F, box;
    parallel = true
) where F <: Function
```

See the user guide for further information.

"""
function minimum_distances(
    x, mol_index::F, box::Box;
    parallel=true
) where F<:Function
    cl = CellList(x,box,parallel=parallel)
    x_list = init_list(x, mol_index)
    minimum_distances!(
        mol_index, 
        x_list, box, cl; parallel=parallel
    )
    return x_list
end

minimum_distances(x, n_atoms_per_molecule_x::Int, box::Box; parallel=true) =
    minimum_distances(x,  i -> mol_index(i,n_atoms_per_molecule_x), box; parallel=parallel)
