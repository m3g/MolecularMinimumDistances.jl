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
    if d < x_list[imol].d
        x_list[imol] = MinimumDistance(j,d)
    end
    jmol = molecule_of_j[j]
    if d < y_list[jmol].d
        y_list[jmol] = MinimumDistance(i,d)
    end
    return lists
end

function keep_best_list!(list1,list2)
    for i in eachindex(list1)
        if list2[i].d < list1[i].d
            list1[i] = list2[i]
        end
    end
end

function reduce_list!(list, list_threaded)
    list .= list_threaded[1]
    for it in 2:length(list_threaded)
        keep_best_list!(list, list_threaded[it])
    end
    return list
end

function reduce_list_pair!(lists, lists_threaded)
    lists[1] .= lists_threaded[1][1]
    lists[2] .= lists_threaded[1][2]
    for it in 2:length(lists_threaded)
        keep_best_list!(lists[1], lists_threaded[it][1])
        keep_best_list!(lists[2], lists_threaded[it][2])
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
    @unpack box, cl, aux = cell_list_data
    lists = (x_list = x_list, y_list = y_list)
    # Update the cell lists with the new coordinates and box
    #cl = UpdateCellList!(x, y, box, cl, aux)
    map_pairwise!(
        (x, y, i, j, d2, lists) -> update_list!(i, j, d2, molecule_of_i, molecule_of_j, lists),
        lists, box, cl; 
        parallel=parallel,
        reduce = reduce_list_pair!,
    )
    # Update the data container, for possible reuse
    cell_list_data.cl = cl
    cell_list_data.aux = aux

    # Return updated lists 
    return x_list, y_list
end

function minimum_distances(
    x, y, n_atoms_per_molecule_x, n_atoms_per_molecule_y, box;
    parallel=true
)
    parallel ? nbatches = (0,0) : nbatches=(1,1)
    cl = CellList(x,y,box,nbatches=nbatches)
    aux = CellListMap.AuxThreaded(cl)
    cell_list_data = CellListMapData(box,cl,aux)
    molecule_of_i = molecule_indices(x, n_atoms_per_molecule_x)
    molecule_of_j = molecule_indices(y, n_atoms_per_molecule_y)
    x_list = fill(MinimumDistance(-1,typemax(eltype(x[begin]))),length(x)÷n_atoms_per_molecule_x)
    y_list = fill(MinimumDistance(-1,typemax(eltype(y[begin]))),length(y)÷n_atoms_per_molecule_y)
    minimum_distances!(x, y, molecule_of_i, molecule_of_j, x_list, y_list, cell_list_data; parallel=parallel)
    return x_list, y_list
end

function molecule_indices(x,n_atoms_per_molecule)  
    nmols = div(length(x),n_atoms_per_molecule) 
    return collect(i for i in 1:nmols for j in 1:n_atoms_per_molecule)
end

#
# Testing routines
#
include("./testing.jl")

end # MolecularMinimumDistances



