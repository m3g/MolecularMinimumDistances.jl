module MolecularMinimumDistances

using CellListMap
export MinimumDistance
export minimum_distances, minimum_distances!
export Box

struct MinimumDistance{T}
    i::Int
    d::T
end

"""

```
mol_index(i_atom,n_atoms_per_molecule) = (i_atom-1) ÷ n_atoms_per_molecule + 1
```

Internal structure or function - interface may change. 

Sets the index of the molecule of an atom in the simples situation, in which all 
molecules have the same number of atoms. This is the default setting, and the 
`mol_index` parameter of the `minimum_distance` functions must be defined manually
in other situations. 

"""
mol_index(i,n_atoms_per_molecule) = (i-1) ÷ n_atoms_per_molecule + 1

"""

```
update_list!(i,j,d2,...)
```

Internal structure or function - interface may change. 

Function that updates one entry of the list minimum distances.

# Extended help

## Call for single-set of molecules

```
function update_list!(
    i, j, d2,
    mol_index::F,
    list::AbstractVector{<:MinimumDistance}
)  where F <: Function
```

## Call for disjoint set of molecules

```
function update_list!(
    i, j, d2,
    mol_index_i::Fi,
    mol_index_j::Fj,
    lists::Tuple{T,T}
) where {Fi<:Function, Fj<:Function, T<:AbstractVector{<:MinimumDistance}}
```

"""
# For a single set of molecules (self-md count)
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
            list[imol] = MinimumDistance(j,d)
        end
        if d < list[jmol].d
            list[jmol] = MinimumDistance(i,d)
        end
    end
    return list
end

# For disjoint sets of molecules
function update_list!(
    i, j, d2,
    mol_index_i::Fi,
    mol_index_j::Fj,
    lists::Tuple{T,T}
) where {Fi<:Function, Fj<:Function, T<:AbstractVector{<:MinimumDistance}}
    x_list = lists[1]
    y_list = lists[2]
    d = sqrt(d2)
    imol = mol_index_i(i)
    if d < x_list[imol].d
        x_list[imol] = MinimumDistance(j,d)
    end
    jmol = mol_index_j(j)
    if d < y_list[jmol].d
        y_list[jmol] = MinimumDistance(i,d)
    end
    return lists
end

"""
```
keep_best_list!(list1,list2)
```

Update `list1` considering the list of minimum distances given in `list2`.

Internal function or structure - interface may change. 

"""
function keep_best_list!(list1,list2)
    for i in eachindex(list1)
        if list2[i].d < list1[i].d
            list1[i] = list2[i]
        end
    end
end

"""

```
reduce_list!(list, list_threaded)
```

Update the final `list` of minimum-distances given the threaded list `list_threaded`.

Internal function or structure - interface may change. 

"""
function reduce_list!(list, list_threaded)
    list .= list_threaded[1]
    for it in 2:length(list_threaded)
        keep_best_list!(list, list_threaded[it])
    end
    return list
end

"""

```
reduce_list_pair!(lists, lists_threaded)
```

Update the final tuple `lists`, of minimum-distances of the two molecule sets given the threaded lists `lists_threaded`.

Internal function or structure - interface may change. 

"""
function reduce_list_pair!(lists, lists_threaded)
    lists[1] .= lists_threaded[1][1]
    lists[2] .= lists_threaded[1][2]
    for it in 2:length(lists_threaded)
        keep_best_list!(lists[1], lists_threaded[it][1])
        keep_best_list!(lists[2], lists_threaded[it][2])
    end
    return lists
end

"""

```
init_list(x, mol_index::F) where F<:Function
```

Initializes an array of type `Vector{MinimumDistance}` with length equal to the number of 
molecules. `x` must be provided so that the type of variable of the coordinates can be 
propagated to the distances, and `mol_index` is the function that given an atomic index `i`
returns the index of the molecule.

```
init_list(::Type{T}, number_of_molecules::Int) 
```

Given the type of the coordinates of the vector of atomic coordinates and the number of molecules,
returns the initialized array of minimum distances. 

# Examples

Providing the type of variable and the number of molecules:

```julia-repl
julia> x = [ rand(SVector{2,Float32}) for _ in 1:100 ]; # 50 molecules of 2 atoms

julia> MolecularMinimumDistances.init_list(Float32,50)
50-element Vector{MinimumDistance{Float32}}:
 MinimumDistance{Float32}(-1, Inf32)
 ⋮
 MinimumDistance{Float32}(-1, Inf32)
```

Providing the vector of coordinates and a function that returns the index of the
molecule of each element:

```julia-repl
julia> MolecularMinimumDistances.init_list(x, i -> (i-1)÷2 + 1)
50-element Vector{MinimumDistance{Float32}}:
 MinimumDistance{Float32}(-1, Inf32)
 ⋮
 MinimumDistance{Float32}(-1, Inf32)
```

The above annonymous function `i -> (i-1)÷2 + 1` is equivalent to `i -> MolecularMinimumDistances.mol_index(i,2)`,
and can be generalized if the the number of atoms of each molecule is not the same.

"""
function init_list(
    x::AbstractVector{<:AbstractVector}, 
    mol_index::F
) where F<:Function
    T = eltype(eltype(x))
    number_of_molecules = 0
    previous_molecule = 0
    for i in eachindex(x)
        imol = mol_index(i)
        if imol != previous_molecule
            number_of_molecules += 1
            previous_molecule = imol
        end
    end
    return init_list(T,number_of_molecules)
end
init_list(::Type{T}, number_of_molecules::Int) where T = 
    fill(MinimumDistance(-1,typemax(T)),number_of_molecules)

#
# Main function for a single set of molecules 
#
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

julia> minimum_distances!(
           i -> MolecularMinimumDistances.mol_index(i,5),
           MolecularMinimumDistances.init_list(x,5),
           box, cl
       )
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


#
# Main functions for two disjoint sets of molecules
#
"""

```
function minimum_distances!(
    molecule_of_i::AbstractVector{<:Integer},
    molecule_of_j::AbstractVector{<:Integer},
    x_list::AbstractVector{<:MinimumDistance},
    y_list::AbstractVector{<:MinimumDistance},
    box, cl;
    parallel = true
)
```

Compute the list of minimum distances given the precomputed cell lists, auxiliary vectors, and 
lists of indexes of molecules of each atom. Should be preferred for multiple calls of this function,
with the outer update of the cell lists.

"""
function minimum_distances!(
    mol_index_i::Fi,
    mol_index_j::Fj,
    x_list::AbstractVector{<:MinimumDistance},
    y_list::AbstractVector{<:MinimumDistance},
    box, cl;
    parallel = true
) where {Fi<:Function, Fj<:Function}
    lists = (x_list, y_list)
    map_pairwise!(
        (x, y, i, j, d2, lists) -> update_list!(i, j, d2, mol_index_i, mol_index_j, lists),
        lists, box, cl; 
        parallel=parallel,
        reduce = reduce_list_pair!,
    )
    return x_list, y_list
end

"""

```
function minimum_distances(
    x, y, 
    n_atoms_per_molecule_x,
    n_atoms_per_molecule_y,
    box;
    parallel = true
)
```

Compute the list of minimum distances given the vector of atomic coordinates `x` and `y`, 
the number of atoms of the molecules, the box (of type `CellListMap.Box`). This is the simplest
interface, but it allocates new arrays. For iterative computations use `minimum_distances!`
instead. 

## Example

```
julia> x = [ rand(SVector{3,Float64}) for _ in 1:100 ];

julia> y = [ rand(SVector{3,Float64}) for _ in 1:90 ];

julia> box =  Box([1,1,1],0.2);

julia> x_list, y_list = minimum_distances(x, y, 5, 3, box);

julia> x_list
20-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(51, 0.048400669109669024)
 MinimumDistance{Float64}(23, 0.08296253778070845)
 ⋮
 MinimumDistance{Float64}(71, 0.04853375996234046)
 MinimumDistance{Float64}(56, 0.03089814908795506)

```

Each entry of the output lists contains, for each atom (in `x` for example), the index of the atom of the other
set (i. e. `y`) that is closer to it, and the distance between these atoms. If no atom is found within
the cutoff of a given atom `MinimumDistance(-1,+inf)` will be returned. 

"""
function minimum_distances(
    x, y, n_atoms_per_molecule_x::Int, n_atoms_per_molecule_y::Int, box::Box;
    parallel=true
)
    parallel ? nbatches = (0,0) : nbatches=(1,1)
    cl = CellList(x,y,box,nbatches=nbatches)
    x_list = init_list(x, i -> mol_index(i,n_atoms_per_molecule_x))
    y_list = init_list(y, i -> mol_index(i,n_atoms_per_molecule_y))
    minimum_distances!(
        i -> mol_index(i,n_atoms_per_molecule_x),
        j -> mol_index(j,n_atoms_per_molecule_y),
        x_list, y_list, box, cl; 
        parallel=parallel
    )
    return x_list, y_list
end

#
# Testing routines
#
include("./testing.jl")

end # MolecularMinimumDistances



