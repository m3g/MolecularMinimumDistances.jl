"""

$(TYPEDEF)

$(INTERNAL)

$(TYPEDFIELDS)

"""
struct SelfPairs{T,F<:Function} <: SystemPairs
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
within a single set of molecules. The shortest distance of each molecule
to any other molecule of the same set is computed.

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
