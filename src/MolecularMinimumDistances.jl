module MolecularMinimumDistances

using DocStringExtensions
const INTERNAL = "Internal function or structure - interface may change."

using StaticArrays
import CellListMap
using .CellListMap.PeriodicSystems

export MinimumDistance
export SelfPairs, CrossPairs, AllPairs
export minimum_distances, minimum_distances!
export getlist

"""

```
MinimumDistance{T}
```

The lists of minimum-distances are stored in arrays of type `Vector{MinimumDistance{T}}`. The index
of this vector corresponds to the index of the molecule in the original array.

`MinimumDistance{T}` is a simple structure that contains four fields: a boolean marker indicating
if the distance is within the cutoff, the indexes `i` and `j` of the atoms of the 
molecules that are closer to each other, and the distance `d`, with type `T`, which is
the same as that of the coordinates of the input vectors of coordinates.

## Example

```julia-repl
julia> md = MinimumDistance{Float32}(true, 2, 5, 1.f0)
MinimumDistance{Float32}(true, 2, 5, 1.0f0)

julia> md.i
2

julia> md.j
5

julia> md.d
1.0f0
```

"""
struct MinimumDistance{T}
    within_cutoff::Bool
    i::Int
    j::Int
    d::T
end
import Base: zero, copy
zero(::Type{MinimumDistance{T}}) where {T} = MinimumDistance(false, 0, 0, typemax(T))
copy(md::MinimumDistance) = MinimumDistance(md.within_cutoff, md.i, md.j, md.d)

# Simplify signature of arrays of MinimumDistance and its tuples.
List{T} = Vector{<:MinimumDistance{T}}
ListTuple{T} = Tuple{Vector{<:MinimumDistance{T}}, Vector{<:MinimumDistance{T}}} 

"""

```
_mol_indices(i_atom,n_atoms_per_molecule) = (i_atom-1) ÷ n_atoms_per_molecule + 1
```

$(INTERNAL)

# Extended help

Sets the index of the molecule of an atom in the simples situation, in which all 
molecules have the same number of atoms. This is the default setting, and the 
`_mol_indices` parameter of the `minimum_distance` functions must be defined manually
in other situations. 

"""
_mol_indices(i, n_atoms_per_molecule) = (i - 1) ÷ n_atoms_per_molecule + 1
_number_of_molecules(mol_indices, positions) = length(unique(mol_indices(i) for i in eachindex(positions)))

function _get_mol_indices(mol_indices, n_atoms_per_molecule; flag::String="")
    if (isnothing(n_atoms_per_molecule) && isnothing(mol_indices)) ||
       (!isnothing(n_atoms_per_molecule) && !isnothing(mol_indices))
        throw(ArgumentError("Please specify *either* $(flag)n_atoms_per_molecule *or* $(flag)mol_indices"))
    end
    if isnothing(mol_indices)
        mol_indices = i -> _mol_indices(i, n_atoms_per_molecule)
    end
    return mol_indices
end

"""

```
init_list(x, mol_indices::F) where F<:Function
```

$(INTERNAL)

# Extended help

Initializes an array of type `Vector{MinimumDistance}` with length equal to the number of 
molecules. `x` must be provided so that the type of variable of the coordinates can be 
propagated to the distances, and `mol_indices` is the function that given an atomic index `i`
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

julia> init_list(Float32,50)
50-element Vector{MinimumDistance{Float32}}:
 MinimumDistance{Float32}(false, -1, -1, Inf32)
 MinimumDistance{Float32}(false, -1, -1, Inf32)
 MinimumDistance{Float32}(false, -1, -1, Inf32)
 ⋮
 MinimumDistance{Float32}(false, -1, -1, Inf32)
 MinimumDistance{Float32}(false, -1, -1, Inf32)

```

Providing the vector of coordinates and a function that returns the index of the
molecule of each element:

```julia-repl
julia> init_list(x, i -> (i-1)÷2 + 1)
50-element Vector{MinimumDistance{Float32}}:
 MinimumDistance{Float32}(false, -1, -1, Inf32)
 MinimumDistance{Float32}(false, -1, -1, Inf32)
 MinimumDistance{Float32}(false, -1, -1, Inf32)
 ⋮
 MinimumDistance{Float32}(false, -1, -1, Inf32)
 MinimumDistance{Float32}(false, -1, -1, Inf32)

```

The above annonymous function `i -> (i-1)÷2 + 1` is equivalent to `i -> mol_indices(i,2)`,
and can be generalized if the the number of atoms of each molecule is not the same.

"""
function init_list(x::AbstractVector{<:AbstractVector}, mol_indices::F) where {F<:Function}
    T = eltype(eltype(x))
    number_of_molecules = 0
    previous_molecule = 0
    for i in eachindex(x)
        imol = mol_indices(i)
        if imol != previous_molecule
            number_of_molecules += 1
            previous_molecule = imol
        end
    end
    return init_list(T, number_of_molecules)
end
init_list(::Type{T}, n::Int) where {T} = fill(zero(MinimumDistance{T}), n)

#
# Functions required for threaded computations
#
import .PeriodicSystems: copy_output, reset_output!, reducer
copy_output(list::List) = copy(list) 
reset_output!(list::List) = list .= Ref(zero(eltype(list)))
reducer(md1::MinimumDistance{T}, md2::MinimumDistance{T}) where T = md1.d < md2.d ? md1 : md2

copy_output(list::ListTuple) = (copy_output(list[1]), copy_output(list[2]))
function reset_output!(list::ListTuple)
    reset_output!(list[1])
    reset_output!(list[2])
    return list
end
function reducer(l1::ListTuple, l2::ListTuple)
    for i in eachindex(l1[1], l2[1])
        l1[1][i] = reducer(l1[1][i], l2[1][i])
    end
    for i in eachindex(l1[2], l2[2])
        l1[2][i] = reducer(l1[2][i], l2[2][i])
    end
    return l1
end


# be careful because ListTuple is not mutable
function reduce_output!(list::ListTuple, list_threaded::Vector{<:ListTuple})
    list = reset_output!(list)
    for i in eachindex(list_threaded)
        list[1] = reduce_output!(list[1], list_threaded[1][i])
        list[2] = reduce_output!(list[2], list_threaded[2][i])
    end
    return list
end

"""

```
minimum_distances!(system)
```

General function that computes the minimum distances for an initialized system,
of `SelfPairs`, `CrossPairs`, or `AllPairs` types. Used as an advanced alternative
from preallocated system inputs.

# Examples

```julia-repl
```

"""
function minimum_distances!(sys)
    map_pairwise!(
        (x, y, i, j, d2, list) -> update_list!(i, j, d2, list, sys),
        sys.system
    )
    return getlist(sys)
end

"""

```
function minimum_distances(
   positions=rand(SVector{3,Float64},10^5),
   cutoff=0.1,
   unitcell=[1,1,1],
   n_atoms_per_molecule=5
)
```
# Examples

```julia-repl
```

"""
function minimum_distances(;
    positions::Union{Nothing,AbstractVector{<:SVector{N,T}}}=nothing,
    xpositions::Union{Nothing,AbstractVector{<:SVector{N,T}}}=nothing,
    ypositions::Union{Nothing,AbstractVector{<:SVector{N,T}}}=nothing,
    cutoff::T,
    unitcell::AbstractVecOrMat,
    mol_indices::Union{Nothing,Function}=nothing,
    xmol_indices::Union{Nothing,Function}=nothing,
    ymol_indices::Union{Nothing,Function}=nothing,
    n_atoms_per_molecule::Union{Nothing,Int}=nothing,
    xn_atoms_per_molecule::Union{Nothing,Int}=nothing,
    yn_atoms_per_molecule::Union{Nothing,Int}=nothing,
    parallel=true
) where {N,T}
    # SelfPairs
    if !isnothing(positions)
        mol_indices = _get_mol_indices(mol_indices, n_atoms_per_molecule)
        system = SelfPairs(;
            positions=positions,
            cutoff=cutoff,
            unitcell=unitcell,
            mol_indices=mol_indices,
            parallel=parallel
        )
        return minimum_distances!(system)
    end
    # CrossPairs
    if isnothing(ymol_indices) && isnothing(yn_atoms_per_molecule)
        xmol_indices = _get_mol_indices(xmol_indices, xn_atoms_per_molecule; flag="x")
        system = CrossPairs(;
            xpositions=xpositions,
            ypositions=ypositions,
            cutoff=cutoff,
            unitcell=unitcell,
            xmol_indices=xmol_indices,
            parallel=parallel
        )
        return minimum_distances!(system)
    end
    # AllPairs
    if !isnothing(xpositions) && (!isnothing(yn_atoms_per_molecule) || !isnothing(ymol_indices))
        xmol_indices = _get_mol_indices(xmol_indices, xn_atoms_per_molecule; flag="x")
        ymol_indices = _get_mol_indices(ymol_indices, yn_atoms_per_molecule; flag="y")
        system = AllPairs(;
            xpositions=xpositions,
            ypositions=ypositions,
            cutoff=cutoff,
            unitcell=unitcell,
            xmol_indices=xmol_indices,
            ymol_indices=ymol_indices,
            parallel=parallel
        )
        return minimum_distances!(system)
    end
    throw(ArgumentError("Incorrect set of keyword input parameters."))
end

abstract type SystemPairs end
getlist(sys::SystemPairs) = sys.system.output

#
# Functions for when the lists of minimum-distances is that of a single
# set of molecules (between the molecules of that set)
#
include("./SelfPairs.jl")

#
# Functions for when one wants the list of the atoms of the second set
# that are closer to each molecule of the first set (only one list is 
# returned)
#
include("./CrossPairs.jl")

# 
# Functions for when all pairs of minimum distances are desired,
# between two disjoint sets of molecules
#
include("./AllPairs.jl")

#
# Testing routines
#
include("./testing.jl")

end # MolecularMinimumDistances
