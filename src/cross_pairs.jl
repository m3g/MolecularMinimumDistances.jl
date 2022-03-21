#
# Functions that return the atoms of the second set that are closer to
# each molecule of the first set (only one list is returned).
#

function update_list_cross!(
    i, j, d2,
    mol_index_i::Fi,
    list::Vector{<:MinimumDistance}
) where {Fi<:Function}
    d = sqrt(d2)
    imol = mol_index_i(i)
    if d < list[imol].d
        list[imol] = MinimumDistance(i,j,d)
    end
    return list
end

"""

```
function minimum_distances!(
    mol_index_i::Fi,
    x_list::AbstractVector{<:MinimumDistance},
    y::AbstractVector,
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
    x_list::AbstractVector{<:MinimumDistance},
    box, cl::CellListMap.CellListPair;
    parallel = true
) where {Fi<:Function}
    map_pairwise!(
        (x, y, i, j, d2, x_list) -> update_list_cross!(i, j, d2, mol_index_i, x_list),
        x_list, box, cl; 
        parallel=parallel,
        reduce = reduce_list!,
    )
    return x_list
end

"""

```
function minimum_distances(
    x, y, 
    n_atoms_per_molecule_x,
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
    x::AbstractVector, y::AbstractVector, n_atoms_per_molecule_x::Int, box::Box;
    parallel=true
)
    parallel ? nbatches = (0,0) : nbatches=(1,1)
    cl = CellList(x,y,box,nbatches=nbatches)
    x_list = init_list(x, i -> mol_index(i,n_atoms_per_molecule_x))
    minimum_distances!(
        i -> mol_index(i,n_atoms_per_molecule_x),
        x_list, box, cl; 
        parallel=parallel
    )
    return x_list
end
