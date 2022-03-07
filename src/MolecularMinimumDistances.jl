module MolecularMinimumDistances

using CellListMap

export MinimumDistance
export minimum_distances, minimum_distances!

struct MinimumDistance{T}
    i::Int
    d::T
end

struct CellListParameters{Box,CellList,AuxThreaded}
    box::Box
    cl::CellList
    aux::AuxThreaded
end

function update_list!(
    i, j, d2,
    x_molecule_of_atom::AbstractVector{<:Integer},
    y_molecule_of_atom::AbstractVector{<:Integer},
    lists::NamedTuple{(:x_list, :y_list), <:Tuple}
)
    (;x_list, y_list) = lists
    d = sqrt(d2)
    imol = x_molecule_of_atom[i]
    if d < list[imol].d
        x_list[imol] = MinimumDistance(j,d)
    end
    jmol = y_molecule_of_atom[j]
    if d < list[jmol].d
        y_list[imol] = MinimumDistance(i,d)
    end
    return lists
end

function minimum_distances!(
    x, y,
    x_molecule_of_atom::AbstractVector{<:Integer},
    y_molecule_of_atom::AbstractVector{<:Integer},
    list_x::AbstractVector{<:MinimumDistance},
    list_y::AbstractVector{<:MinimumDistance},
    cell_list_parameters::CellListMapParameters,
)
    lists = (x_list = x_list, y_list = y_list)
    (; box, cl, aux = cell_list_parameters)
    map_pairwise!(
        (x, y, i, j,d2, lists) -> update_list!(i, j, d2, x_molecule_of_atom, y_molecule_of_atom, lists),
        lists, box, cl, aux
    )
    return list
end





end # module