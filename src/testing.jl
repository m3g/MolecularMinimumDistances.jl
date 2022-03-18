#
# Testing routines
#
function naive_md(x, y, n_atoms_per_molecule_x,n_atoms_per_molecule_y, box::Box)
    nmols_x = length(x) ÷ n_atoms_per_molecule_x
    nmols_y = length(y) ÷ n_atoms_per_molecule_y
    x_list = fill(MinimumDistance(-1,typemax(eltype(x[begin]))),nmols_x)
    y_list = fill(MinimumDistance(-1,typemax(eltype(y[begin]))),nmols_y)
    molecule_of_i = MolecularMinimumDistances.molecule_indices(x,n_atoms_per_molecule_x)
    molecule_of_j = MolecularMinimumDistances.molecule_indices(y,n_atoms_per_molecule_y)
    for (i,vx) in pairs(x)
        for (j,vy) in pairs(y)
            vy_wrapped = CellListMap.wrap_relative_to(vy, vx, box)
            d = sqrt(sum(abs2,vy_wrapped - vx))
            if d > box.cutoff 
                continue
            end
            imol = molecule_of_i[i]
            if d < x_list[imol].d 
                x_list[imol] = MinimumDistance(j,d)
            end
            jmol = molecule_of_j[j]
            if d < y_list[jmol].d
                y_list[jmol] = MinimumDistance(i,d)
            end
        end
    end
    return x_list, y_list
end

import Base.isapprox
function isapprox(
    list1::AbstractVector{<:MinimumDistance},
    list2::AbstractVector{<:MinimumDistance}
)
    length(list1) != length(list2) && return false
    for i in eachindex(list1)
        list1[i].i != list2[i].i && return false
        !(list1[i].d ≈ list2[i].d) && return false
    end
    return true
end
function isapprox(
    lists1::Tuple{T,T},
    lists2::Tuple{T,T}
) where T<:AbstractVector{<:MinimumDistance}
    !(lists1[1] ≈ lists2[1]) && return false
    !(lists1[2] ≈ lists2[2]) && return false
    return true
end

function plot(x,y)
    Main.scatter(Tuple.(x))
    Main.scatter(Tuple.(y))
    


end
