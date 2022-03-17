module MolecularMinimumDistances

using Parameters
using CellListMap

export MinimumDistance
export minimum_distances, minimum_distances!

struct MinimumDistance{T}
    i::Int
    d::T
end

mutable struct CellListMapData{Box,CellList,AuxThreaded}
    box::Box
    cl::CellList
    aux::AuxThreaded
end

function update_list!(
    i, j, d2,
    molecule_of_i::AbstractVector{<:Integer},
    molecule_of_j::AbstractVector{<:Integer},
    lists
)
    @unpack x_list, y_list = lists
    d = sqrt(d2)
    imol = molecule_of_i[i]
    if d < list[imol].d
        x_list[imol] = MinimumDistance(j,d)
    end
    jmol = molecule_of_j[j]
    if d < list[jmol].d
        y_list[imol] = MinimumDistance(i,d)
    end
    return lists
end

function minimum_distances!(
    x, y,
    molecule_of_i::AbstractVector{<:Integer},
    molecule_of_j::AbstractVector{<:Integer},
    x_list::AbstractVector{<:MinimumDistance},
    y_list::AbstractVector{<:MinimumDistance},
    cell_list_data::CellListMapData;
    parallel = true
)
    @unpack box, cl, aux = cell_list_parameters
    lists = (x_list = x_list, y_list = y_list)
    # Update the cell lists with the new coordinates and box
    cl = UpdateCellList!(x, y, box, cl, aux)
    map_pairwise!(
        (x, y, i, j,d2, lists) -> update_list!(i, j, d2, molecule_of_i, molecule_of_j, lists),
        lists, box, cl, parallel=parallel
    )
    # Update the data container, for possible reuse
    cell_list_data.cl = cl
    call_list_data.aux = aux

    # Return the list of minimum distances
    return lists
end

function minimum_distances(x, n_atoms_per_molecule_x, y, n_atoms_per_molecule_y, box)
    cl = CellList(x,y,box)
    aux = CellListMap.AuxThreaded(cl)
    cl_data = CellListMapData(box,cl,aux)
    molecule_of_i = molecule_indices(x, n_atoms_per_molecule_x)
    molecule_of_j = molecule_indices(y, n_atoms_per_molecule_y)
    x_list = Vector{MinimumDistance}(undef, )
    y_lis




end

function molecule_indices(x,n_atoms_per_molecule)  
    nmols = div(length(x),n_atoms_per_molecule) 
    return collect(i for i in 1:nmols for j in 1:n_atoms_per_molecule)
end

end # module



