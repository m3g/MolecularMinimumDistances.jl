"""

$(TYPEDEF)

$(INTERNAL)

$(TYPEDFIELDS)

"""
struct SelfPairs{T,F<:Function}
    system::T
    mol_indices::F
end

"""

```
SelfPairs(;
    positions::AbstractVector{<:SVector{N,T}},
    cutoff::T,
    unitcell::AbstractVecOrMat,
    n_atoms_per_molecule::Int,
    parallel::Bool=true
) where T<:Real
```

Initializes a particle system for the calculation of minimum distances
within a single set of molecules. The shortest distance of of each molecule
to any other molecule of the same set is computes.

Instead of the number of atoms per molecule, the user can also provide a 
more general `mol_index` function, which, for each atomic index, returns the 
corresponding
molecular index (which is `mol_index(i) = (i-1)%n + 1` where `n` is the
number of atoms per molecule if all molecules have the same number of
atoms and are continously stored in the array of positions). 

# Examples

"""
function SelfPairs(;
    positions::AbstractVector{<:SVector{N,T}},
    cutoff::T,
    unitcell::AbstractVecOrMat,
    n_atoms_per_molecule::Union{Nothing,Int}=nothing,
    mol_index::Union{Nothing,Function}=nothing,
    parallel::Bool=true
) where {N,T<:Real}
    mol_index = _get_mol_index(mol_index, n_atoms_per_molecule)
    system = PeriodicSystem(;
        positions=positions,
        cutoff=cutoff,
        unitcell=unitcell,
        output=init_list(positions, mol_index),
        parallel=parallel
    )
    return SelfPairs(system, mol_index)
end

#
# This file containst the functions for single-sets, that is for those cases where
# the list of minimum-distances is between the molecules of a single component.
#
function update_list!(i, j, d2, list::List, system::SelfPairs)
    imol = system.mol_index(i)
    jmol = system.mol_index(j)
    if imol != jmol
        d = sqrt(d2)
        if d < list[imol].d
            list[imol] = MinimumDistance(true, i, j, d)
        end
        if d < list[jmol].d
            list[jmol] = MinimumDistance(true, j, i, d)
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
    parallel = true,
    list_threaded = nothing
) where F <: Function
```

Compute the list of minimum distances given the precomputed cell lists, auxiliary vectors, and 
lists of indexes of molecules of each atom. Should be preferred for multiple calls of this function,
with the outer update of the cell lists. The `mol_index` function must return, for each atom index,
the corresponding molecular index. To obtain faster an minimally-allocating calls, 
preallocate the `list_threaded` variable.

# Examples

```julia-repl
julia> using CellListMap, StaticArrays

julia> x = [ rand(SVector{2,Float32}) for _ in 1:100 ];

julia> box = Box([1,1],0.2f0);

julia> cl = CellList(x,box);

julia> list = init_list(x, i -> mol_index(i,5));

julia> list_threaded = [ copy(list) for _ in 1:nbatches(cl) ]

julia> minimum_distances!(i -> mol_index(i,5), list, box, cl)
20-element Vector{MinimumDistance{Float32}}:
 MinimumDistance{Float32}(true, 5, 56, 0.033878796f0)
 MinimumDistance{Float32}(true, 9, 69, 0.009123626f0)
 MinimumDistance{Float32}(true, 15, 100, 0.054782398f0)
 ⋮
 MinimumDistance{Float32}(true, 93, 7, 0.03434341f0)
 MinimumDistance{Float32}(true, 96, 29, 0.03825064f0)

```

"""
function minimum_distances(
    positions::AbstractVector{<:SVector{N,T}},
    cutoff::T,
    unitcell::AbstractVecOrMat,
    mol_index::Function,
    parallel=true
) where {N,T}
    system = SelfPair(;
        positions=positions,
        cutoff=cutoff,
        unitcell=unitcell,
        mol_index=mol_index,
        parallel=parallel
    )
    minimum_distances!(system)
    return x_list
end

function minimum_distances(
    positions::AbstractVector{<:SVector{N,T}},
    cutoff::T,
    unitcell::AbstractVecOrMat,
    n_atoms_per_molecule::Int,
    parallel=true
) where {N,T}
    system = ParticleSystem(;
        positions=positions,
        cutoff=cutoff,
        unitcell=unitcell,
        mol_index=i -> mol_index(i, n_atoms_per_molecule),
        parallel=parallel
    )
    minimum_distances!(system)
    return x_list
end

function minimumdistances(system::SelfPairs)
    list = map_pairwise!(
        (x, y, i, j, d2, list) -> update_list!(i, j, d2, system.mol_index, list),
        system
    )
    return list
end


