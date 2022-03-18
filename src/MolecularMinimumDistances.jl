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
update_list!(i,j,d2,...)
```

Internal structure or function - interface may change. 

Function that updates one entry of the list minimum distances.

# Extended help

## Call for single-set of molecules

```
function update_list!(
    i, j, d2,
    molecule_of_i::AbstractVector{<:Integer},
    list::AbstractVector{<:MinimumDistance}
) 
```

## Call for disjoint set of molecules

```
function update_list!(
    i, j, d2,
    molecule_of_i::AbstractVector{<:Integer},
    molecule_of_j::AbstractVector{<:Integer},
    lists::Tuple{T,T}
) where T<:AbstractVector{<:MinimumDistance}
```

"""
# For a single set of molecules (self-md count)
function update_list!(
    i, j, d2,
    molecule_of_i::AbstractVector{<:Integer},
    list::AbstractVector{<:MinimumDistance}
) 
    imol = molecule_of_i[i]
    jmol = molecule_of_i[j]
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
    molecule_of_i::AbstractVector{<:Integer},
    molecule_of_j::AbstractVector{<:Integer},
    lists::Tuple{T,T}
) where T<:AbstractVector{<:MinimumDistance}
    x_list = lists[1]
    y_list = lists[2]
    d = sqrt(d2)
    imol = molecule_of_i[i]
    if d < x_list[imol].d
        x_list[imol] = MinimumDistance(j,d)
    end
    jmol = molecule_of_j[j]
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
init_list(x,n_atoms_per_molecule_x)  
```

Initializes an array of type `Vector{MinimumDistance}` with length equal to the number of 
molecules of `x`. 

"""
init_list(x,n_atoms_per_molecule_x) = 
    fill(MinimumDistance(-1,typemax(eltype(x[begin]))),length(x)÷n_atoms_per_molecule_x)

function molecule_indices(x,n_atoms_per_molecule)  
    if length(x)%n_atoms_per_molecule != 0
        error("The number of atoms must be a multiple of the number of atoms per molecule.")
    end
    nmols = length(x) ÷ n_atoms_per_molecule
    return collect(i for i in 1:nmols for j in 1:n_atoms_per_molecule)
end

#
# Main function for a single set of molecules 
#
"""

```
function minimum_distances!(
    molecule_of_i::AbstractVector{<:Integer},
    x_list::AbstractVector{<:MinimumDistance},
    box, cl;
    parallel = true
)
```

Compute the list of minimum distances given the precomputed cell lists, auxiliary vectors, and 
lists of indexes of molecules of each atom. Should be preferred for multiple calls of this function,
with the outer update of the cell lists.

"""
function minimum_distances!(
    molecule_of_i::AbstractVector{<:Integer},
    x_list::AbstractVector{<:MinimumDistance},
    box, cl;
    parallel = true
)
    map_pairwise!(
        (x, y, i, j, d2, list) -> update_list!(i, j, d2, molecule_of_i, list),
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
interface, but it allocates new arrays. For iterative computations use `minimum_distances!`
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

"""
function minimum_distances(
    x, n_atoms_per_molecule_x::Int, box::Box;
    parallel=true
)
    cl = CellList(x,box,parallel=parallel)
    molecule_of_i = molecule_indices(x, n_atoms_per_molecule_x)
    x_list = init_list(x, n_atoms_per_molecule_x)
    minimum_distances!(molecule_of_i, x_list, box, cl; parallel=parallel)
    return x_list
end


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
    molecule_of_i::AbstractVector{<:Integer},
    molecule_of_j::AbstractVector{<:Integer},
    x_list::AbstractVector{<:MinimumDistance},
    y_list::AbstractVector{<:MinimumDistance},
    box, cl;
    parallel = true
)
    lists = (x_list, y_list)
    map_pairwise!(
        (x, y, i, j, d2, lists) -> update_list!(i, j, d2, molecule_of_i, molecule_of_j, lists),
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
    molecule_of_i = molecule_indices(x, n_atoms_per_molecule_x)
    molecule_of_j = molecule_indices(y, n_atoms_per_molecule_y)
    x_list = init_list(x, n_atoms_per_molecule_x)
    y_list = init_list(y, n_atoms_per_molecule_y)
    minimum_distances!(molecule_of_i, molecule_of_j, x_list, y_list, box, cl ; parallel=parallel)
    return x_list, y_list
end

#
# Testing routines
#
include("./testing.jl")

end # MolecularMinimumDistances



