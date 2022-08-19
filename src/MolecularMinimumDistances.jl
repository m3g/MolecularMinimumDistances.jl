module MolecularMinimumDistances

using DocStringExtensions
const INTERNAL = "Internal function or structure - interface may change."

using StaticArrays
import CellListMap
using .CellListMap.PeriodicSystems

export MinimumDistance
export SelfPairs, CrossPairs, AllPairs
export minimum_distances, minimum_distances!

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
import Base: zero
zero(::Type{MinimumDistance{T}}) where {T} = MinimumDistance(false, 0, 0, typemax(T))

# Simplify signature of arrays of MinimumDistance and its tuples.
const List{T} = Vector{<:MinimumDistance{T}}
const ListTuple{T} = Tuple{<:List{T},<:List{T}}

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

function _get_mol_indices(mol_indices, n_atoms_per_molecule)
    if (isnothing(n_atoms_per_molecule) && isnothing(mol_indices)) ||
       (!isnothing(n_atoms_per_molecule) && !isnothing(mol_indices))
        throw(ArgumentError("Please specify *either* n_atoms_per_molecule *or* mol_indices"))
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
import .PeriodicSystems: copy_output, reset_output!, reducer, reduce_output!
copy_output(md::MinimumDistance) = md
reset_output!(md::MinimumDistance{T}) where {T} = zero(MinimumDistance{T})
reducer(md1::T, md2::T) where {T<:MinimumDistance} = md1.d < md2.d ? md1 : md2

copy_output(list::ListTuple) = (copy_output(list[1]), copy_output(list[2]))
function reset_output!(list::ListTuple)
    reset_output!(list[1])
    reset_output!(list[2])
    return list_tuple
end
function reducer(list::ListTuple, list_threaded::Vector{ListTuple})
    reset_output!(list)
    for i in eachindex(list_threaded)
        list = reduce_output!(list, list_threaded[i])
    end
    return list
end

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
#include("./cross_pairs.jl")

# 
# Functions for when all pairs of minimum distances are desired,
# between two disjoint sets of molecules
#
#include("./all_pairs.jl")

#
# Testing routines
#
include("./testing.jl")

end # MolecularMinimumDistances
