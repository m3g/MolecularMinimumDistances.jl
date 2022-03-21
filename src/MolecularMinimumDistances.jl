module MolecularMinimumDistances

using CellListMap
export MinimumDistance
export minimum_distances, minimum_distances!
export mol_index, init_list
export Box

"""

```
MinimumDistance{T}
```

The lists of minimum-distances are stored in arrays of type `Vector{MinimumDistance{T}}`. The index
of this vector corresponds to the index of the molecule in the original array.

`MinimumDistance{T}` is a simple structure that contains two fields: the indexes `i` and `j` of the atoms of the 
molecules that are closer to each other, and the distance `d`, with type `T`, which is
the same as that of the coordinates of the input vectors of coordinates.

## Example

```julia-repl
julia> md = MinimumDistance{Float32}(2, 5, 1.f0)
MinimumDistance{Float32}(2, 5, 1.0f0)

julia> md.i
2

julia> md.j
5

julia> md.d
1.0f0
```

"""
struct MinimumDistance{T}
    i::Int
    j::Int
    d::T
end

"""

```
mol_index(i_atom,n_atoms_per_molecule) = (i_atom-1) ÷ n_atoms_per_molecule + 1
```

Internal structure or function - interface may change. 

Sets the index of the molecule of an atom in the simples situation, in which all 
molecules have the same number of atoms. This is the default setting, and the 
`mol_index` parameter of the `minimum_distance` functions must be defined manually
in other situations. 

"""
mol_index(i,n_atoms_per_molecule) = (i-1) ÷ n_atoms_per_molecule + 1

"""

```
init_list(x, mol_index::F) where F<:Function
```

Initializes an array of type `Vector{MinimumDistance}` with length equal to the number of 
molecules. `x` must be provided so that the type of variable of the coordinates can be 
propagated to the distances, and `mol_index` is the function that given an atomic index `i`
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
 MinimumDistance{Float32}(-1, -1, Inf32)
 MinimumDistance{Float32}(-1, -1, Inf32)
 MinimumDistance{Float32}(-1, -1, Inf32)
 ⋮
 MinimumDistance{Float32}(-1, -1, Inf32)
 MinimumDistance{Float32}(-1, -1, Inf32)

```

Providing the vector of coordinates and a function that returns the index of the
molecule of each element:

```julia-repl
julia> init_list(x, i -> (i-1)÷2 + 1)
50-element Vector{MinimumDistance{Float32}}:
 MinimumDistance{Float32}(-1, -1, Inf32)
 MinimumDistance{Float32}(-1, -1, Inf32)
 MinimumDistance{Float32}(-1, -1, Inf32)
 ⋮
 MinimumDistance{Float32}(-1, -1, Inf32)
 MinimumDistance{Float32}(-1, -1, Inf32)

```

The above annonymous function `i -> (i-1)÷2 + 1` is equivalent to `i -> mol_index(i,2)`,
and can be generalized if the the number of atoms of each molecule is not the same.

"""
function init_list(
    x::AbstractVector{<:AbstractVector}, 
    mol_index::F
) where F<:Function
    T = eltype(eltype(x))
    number_of_molecules = 0
    previous_molecule = 0
    for i in eachindex(x)
        imol = mol_index(i)
        if imol != previous_molecule
            number_of_molecules += 1
            previous_molecule = imol
        end
    end
    return init_list(T,number_of_molecules)
end
init_list(::Type{T}, number_of_molecules::Int) where T = 
    fill(MinimumDistance(-1, -1,typemax(T)),number_of_molecules)

# Reduction functions for lists of minimum-distances
include("./reduction_functions.jl")

#
# Functions for when the lists of minimum-distances is that of a single
# set of molecules (between the molecules of that set)
#
include("./self_pairs.jl")

#
# Functions for when one wants the list of the atoms of the second set
# that are closer to each molecule of the first set (only one list is 
# returned)
#
include("./cross_pairs.jl")

# 
# Functions for when all pairs of minimum distances are desired,
# between two disjoint sets of molecules
#
include("./all_pairs.jl")

#
# Testing routines
#
include("./testing.jl")

end # MolecularMinimumDistances



