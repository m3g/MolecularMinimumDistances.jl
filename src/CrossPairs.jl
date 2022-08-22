"""

$(TYPEDEF)

$(INTERNAL)

$(TYPEDFIELDS)

"""
struct CrossPairs{T,F<:Function} <: SystemPairs
    system::T
    mol_indices::F
end

import Base.show
function Base.show(io::IO, mime::MIME"text/plain", sys::CrossPairs)
    print(io,"""
    CrossPairs system with:

    Number of atoms of set: $(length(sys.system.xpositions))
    Number of atoms of target structure: $(length(sys.system.ypositions))
    Cutoff: $(sys.system._box.cutoff)
    unitcell: [$(join(CellListMap._uround.(sys.system._box.unit_cell.matrix),", "))]
    Number of molecules in set: $(_number_of_molecules(sys.mol_indices, sys.system.xpositions))""")
end

"""

```
CrossPairs(;
    xpositions::AbstractVector{<:SVector{N,T}},
    ypositions::AbstractVector{<:SVector{N,T}},
    cutoff::T,
    unitcell::AbstractVecOrMat,
    xn_atoms_per_molecule::Int,
    parallel::Bool=true
) where T<:Real
```

Initializes a particle system for the calculation of minimum distances
between one molecule and a set of other molecules. Returns a list 
minimum distances (`MinimumDistance` type), containing for each
molecule of the set the information about the closest distance to the
reference molecule.

Instead of the number of atoms per molecule, the user can also provide a 
more general `mol_indices` function, which, for each atomic index, returns the 
corresponding
molecular index (which is `mol_indices(i) = (i-1)%n + 1` where `n` is the
number of atoms per molecule if all molecules have the same number of
atoms and are continously stored in the array of positions). 

# Examples

```julia-repl
julia> using MolecularMinimumDistances, StaticArrays

julia> sys = CrossPairs(
           xpositions=rand(SVector{3,Float64},100), # "solvent" (set of molecules)
           ypositions=rand(SVector{3,Float64},10^5), # "solute" (target structure)
           cutoff=0.1,
           unitcell=[1,1,1],
           xn_atoms_per_molecule=5 # of the "solvent"
       )
CrossPairs system with:

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
function CrossPairs(;
    xpositions::AbstractVector{<:SVector{N,T}},
    ypositions::AbstractVector{<:SVector{N,T}},
    cutoff::T,
    unitcell::AbstractVecOrMat,
    xn_atoms_per_molecule::Union{Nothing,Int}=nothing,
    xmol_indices::Union{Nothing,Function}=nothing,
    parallel::Bool=true
) where {N,T<:Real}
    xmol_indices = _get_mol_indices(xmol_indices, xn_atoms_per_molecule; flag="x")
    system = PeriodicSystem(;
        xpositions=xpositions,
        ypositions=ypositions,
        cutoff=cutoff,
        unitcell=unitcell,
        output=init_list(xpositions, xmol_indices),
        parallel=parallel
    )
    return CrossPairs(system, xmol_indices)
end

#
# Functions that return the atoms of the second set that are closer to
# each molecule of the first set (only one list is returned).
#
function update_list!(i, j, d2, list, system::CrossPairs)
    d = sqrt(d2)
    imol = system.mol_indices(i)
    if d < list[imol].d
        list[imol] = MinimumDistance(true, i, j, d)
    end
    return list
end