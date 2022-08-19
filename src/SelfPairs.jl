"""

$(TYPEDEF)

$(INTERNAL)

$(TYPEDFIELDS)

"""
struct SelfPairs{T,F<:Function}
    system::T
    mol_indices::F
end

import Base.show
function Base.show(io::IO, mime::MIME"text/plain", sys::SelfPairs)
    print(io,"""
    SelfPairs system with:

    Number of atoms: $(length(sys.system.positions))
    Cutoff: $(sys.system._box.cutoff)
    unitcell: [$(join(CellListMap._uround.(sys.system._box.unit_cell.matrix),", "))]
    Number of molecules: $(_number_of_molecules(sys.mol_indices, sys.system.positions))""")
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
more general `mol_indices` function, which, for each atomic index, returns the 
corresponding
molecular index (which is `mol_indices(i) = (i-1)%n + 1` where `n` is the
number of atoms per molecule if all molecules have the same number of
atoms and are continously stored in the array of positions). 

# Examples

```julia-repl
julia> using MolecularMinimumDistances, StaticArrays

julia> sys = SelfPairs(
           positions=rand(SVector{3,Float64},10^5),
           cutoff=0.1,
           unitcell=[1,1,1],
           n_atoms_per_molecule=5
       )
SelfPairs system with:

Number of atoms: 100000
Cutoff: 0.1
unitcell: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
Number of molecules: 20000

julia> minimum_distances!(sys)
20000-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(true, 4, 24930, 0.008039482961077074)
 MinimumDistance{Float64}(true, 6, 74055, 0.0049818659155905255)
 ⋮
 MinimumDistance{Float64}(true, 99999, 75403, 0.0025051670801269433)
```

"""
function SelfPairs(;
    positions::AbstractVector{<:SVector{N,T}},
    cutoff::T,
    unitcell::AbstractVecOrMat,
    n_atoms_per_molecule::Union{Nothing,Int}=nothing,
    mol_indices::Union{Nothing,Function}=nothing,
    parallel::Bool=true
) where {N,T<:Real}
    mol_indices = _get_mol_indices(mol_indices, n_atoms_per_molecule)
    system = PeriodicSystem(;
        positions=positions,
        cutoff=cutoff,
        unitcell=unitcell,
        output=init_list(positions, mol_indices),
        parallel=parallel
    )
    return SelfPairs(system, mol_indices)
end

#
# This file containst the functions for single-sets, that is for those cases where
# the list of minimum-distances is between the molecules of a single component.
#
function update_list!(i, j, d2, list::List, system::SelfPairs)
    imol = system.mol_indices(i)
    jmol = system.mol_indices(j)
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
    mol_indices::F,
    x_list::AbstractVector{<:MinimumDistance},
    box, cl;
    parallel = true,
    list_threaded = nothing
) where F <: Function
```

Compute the list of minimum distances given the precomputed cell lists, auxiliary vectors, and 
lists of indexes of molecules of each atom. Should be preferred for multiple calls of this function,
with the outer update of the cell lists. The `mol_indices` function must return, for each atom index,
the corresponding molecular index. To obtain faster an minimally-allocating calls, 
preallocate the `list_threaded` variable.

# Examples

```julia-repl
julia> using CellListMap, StaticArrays

julia> x = [ rand(SVector{2,Float32}) for _ in 1:100 ];

julia> box = Box([1,1],0.2f0);

julia> cl = CellList(x,box);

julia> list = init_list(x, i -> mol_indices(i,5));

julia> list_threaded = [ copy(list) for _ in 1:nbatches(cl) ]

julia> minimum_distances!(i -> mol_indices(i,5), list, box, cl)
20-element Vector{MinimumDistance{Float32}}:
 MinimumDistance{Float32}(true, 5, 56, 0.033878796f0)
 MinimumDistance{Float32}(true, 9, 69, 0.009123626f0)
 MinimumDistance{Float32}(true, 15, 100, 0.054782398f0)
 ⋮
 MinimumDistance{Float32}(true, 93, 7, 0.03434341f0)
 MinimumDistance{Float32}(true, 96, 29, 0.03825064f0)

```

"""
function minimum_distances(;
    positions::AbstractVector{<:SVector{N,T}},
    cutoff::T,
    unitcell::AbstractVecOrMat,
    mol_indices::Union{Nothing,Function}=nothing,
    n_atoms_per_molecule::Union{Nothing,Int}=nothing,
    parallel=true
) where {N,T}
    mol_indices = _get_mol_indices(mol_indices, n_atoms_per_molecule)
    system = SelfPair(;
        positions=positions,
        cutoff=cutoff,
        unitcell=unitcell,
        mol_indices=mol_indices,
        parallel=parallel
    )
    return minimum_distances!(system)
end

function minimum_distances!(sys::SelfPairs)
    list = map_pairwise!(
        (x, y, i, j, d2, list) -> update_list!(i, j, d2, list, sys),
        sys.system
    )
    return list
end


