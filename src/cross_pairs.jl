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
 MinimumDistance{Float64}(-1, -1, Inf)
 MinimumDistance{Float64}(-1, -1, Inf)
 MinimumDistance{Float64}(-1, -1, Inf)

julia> y_list = init_list(y, i -> mol_index(i,10)); # 10 atoms per molecule

julia> minimum_distances!(i -> mol_index(i,4), i -> mol_index(i,10), x_list, y_list, box, cl);

julia> x_list
3-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(1, 516, 0.007047201657257146)
 MinimumDistance{Float64}(7, 415, 0.007128399163692169)
 MinimumDistance{Float64}(9, 511, 0.0033074993141050577)

julia> y_list
80-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(10, 11, 0.015032910410462463)
 MinimumDistance{Float64}(16, 7, 0.025198239285666908)
 MinimumDistance{Float64}(27, 4, 0.02006220925909916)
 ⋮
 MinimumDistance{Float64}(787, 11, 0.034239040614140855)
 MinimumDistance{Float64}(797, 7, 0.026912200186074015)

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
 MinimumDistance{Float64}(4, 39, 0.11311002375318155)
 MinimumDistance{Float64}(9, 42, 0.0752872244199588)
 MinimumDistance{Float64}(12, 8, 0.06411501358230641)
 ⋮
 MinimumDistance{Float64}(92, 88, 0.0468461712318602)
 MinimumDistance{Float64}(96, 50, 0.09478144847155845)

```

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
