#
# Functions that return *all* pairs of minimum distances, that is, *two lists*, one
# of the atoms of one set closer to the molecules of the other set, and vice-versa.
#

"""

```
function update_list!(
    i, j, d2,
    mol_index_i::Fi,
    mol_index_j::Fj,
    lists::Tuple{T,T}
) where {Fi<:Function, Fj<:Function, T<:AbstractVector{<:MinimumDistance}}
```

Internal function or structure - interface may change. 

# Extended help

Function to update the minimum distance in the case where two lists are being constructed.

"""
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
        x_list[imol] = MinimumDistance(i,j,d)
    end
    jmol = mol_index_j(j)
    if d < y_list[jmol].d
        y_list[jmol] = MinimumDistance(j,i,d)
    end
    return lists
end

"""

```
function minimum_distances!(
    mol_index_i::Fi,
    mol_index_j::Fj
    x_list::AbstractVector{<:MinimumDistance},
    y_list::AbstractVector{<:MinimumDistance},
    box, cl;
    parallel = true
) where {Fi<:Function, Fj<:Function}
```

Compute the list of minimum distances given the precomputed cell lists, auxiliary vectors, and 
the functions that associate each atom to the corresponding molecule. Should be preferred for multiple calls of this function,
with the outer update of the cell lists.

## Example

```
julia> x = [ rand(SVector{3,Float64}) for _ in 1:12 ];

julia> y = [ rand(SVector{2,Float64}) for _ in 1:800 ];

julia> using CellListMap

julia> box = Box([1,1],0.2);

julia> cl = CellList(x,y,box);

julia> x_list = init_list(x, i -> mol_index(i,4)) # 4 atoms per molecule
3-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(-1, Inf)
 MinimumDistance{Float64}(-1, Inf)
 MinimumDistance{Float64}(-1, Inf)

julia> y_list = init_list(y, i -> mol_index(i,10)); # 10 atoms per molecule

julia> minimum_distances!(i -> mol_index(i,4), i -> mol_index(i,10), x_list, y_list, box, cl);

julia> x_list
3-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(53, 0.005824258092196386)
 MinimumDistance{Float64}(439, 0.003167356897471225)
 MinimumDistance{Float64}(467, 0.01260590534704902)

julia> y_list
80-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(5, 0.030187364306879957)
 MinimumDistance{Float64}(9, 0.06305749578995341)
 ⋮
 MinimumDistance{Float64}(9, 0.022385093830840475)
 MinimumDistance{Float64}(1, 0.07757423990370381)
```

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
