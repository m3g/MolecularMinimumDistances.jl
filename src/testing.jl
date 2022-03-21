#
# Testing routines: computing the lists with naive algorithms
#

# For a single set of molecules
function naive_md(x, n_atoms_per_molecule_x::Int, box::Box)
    x_list = init_list(x, i -> mol_index(i, n_atoms_per_molecule_x))
    for i in 1:length(x)-1
        vx = x[i]
        for j in i+1:length(x)
            vy = x[j]
            imol = mol_index(i, n_atoms_per_molecule_x)
            jmol = mol_index(j, n_atoms_per_molecule_x)
            if imol == jmol
                continue
            end
            vy_wrapped = CellListMap.wrap_relative_to(vy, vx, box)
            d = sqrt(sum(abs2,vy_wrapped - vx))
            if d > box.cutoff 
                continue
            end
            if d < x_list[imol].d 
                x_list[imol] = MinimumDistance(i,j,d)
            end
            if d < x_list[jmol].d
                x_list[jmol] = MinimumDistance(j,i,d)
            end
        end
    end
    return x_list
end

# For disjoint sets, returning only one list
function naive_md(x, y, n_atoms_per_molecule_x::Int, box::Box)
    x_list = init_list(x, i -> mol_index(i,n_atoms_per_molecule_x))
    for (i,vx) in pairs(x)
        for (j,vy) in pairs(y)
            vy_wrapped = CellListMap.wrap_relative_to(vy, vx, box)
            d = sqrt(sum(abs2,vy_wrapped - vx))
            if d > box.cutoff 
                continue
            end
            imol = mol_index(i,n_atoms_per_molecule_x)
            if d < x_list[imol].d 
                x_list[imol] = MinimumDistance(i,j,d)
            end
        end
    end
    return x_list
end

# For disjoint sets, returning both lists
function naive_md(x, y, n_atoms_per_molecule_x::Int, n_atoms_per_molecule_y::Int, box::Box)
    x_list = init_list(x, i -> mol_index(i,n_atoms_per_molecule_x))
    y_list = init_list(y, i -> mol_index(i,n_atoms_per_molecule_y))
    for (i,vx) in pairs(x)
        for (j,vy) in pairs(y)
            vy_wrapped = CellListMap.wrap_relative_to(vy, vx, box)
            d = sqrt(sum(abs2,vy_wrapped - vx))
            if d > box.cutoff 
                continue
            end
            imol = mol_index(i,n_atoms_per_molecule_x)
            if d < x_list[imol].d 
                x_list[imol] = MinimumDistance(i,j,d)
            end
            jmol = mol_index(j,n_atoms_per_molecule_y)
            if d < y_list[jmol].d
                y_list[jmol] = MinimumDistance(j,i,d)
            end
        end
    end
    return x_list, y_list
end

#
# Comparison functions
#
import Base.isapprox
function isapprox(
    list1::AbstractVector{<:MinimumDistance},
    list2::AbstractVector{<:MinimumDistance}
)
    length(list1) != length(list2) && return false
    for i in eachindex(list1)
        list1[i].i != list2[i].i && return false
        list1[i].j != list2[i].j && return false
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

#
# Plotting example functions
#
function plot(type=1)

    SVector = CellListMap.SVector

    if type == 1

        x = rand(SVector{2,Float64}, 30)
        for i in 1:2:length(x)
            x[i] = -10 .+ 20 * x[i]
            ϕ = rand(0:1e-5:2π)
            x[i+1] = x[i] + SVector{2,Float64}(sin(ϕ),cos(ϕ))
        end

        box = Box(limits(x),10.)
        x_list = minimum_distances(x, 2, box)

        plt = Main.plot()

        for pair in x_list
            if pair.d < 10. 
                xline = [ x[pair.i][1], x[pair.j][1] ]
                yline = [ x[pair.i][2], x[pair.j][2] ]
                Main.plot!(plt, xline, yline, label=:none, color=:black, alpha=0.25)
            end
        end

        for i in 1:2:length(x)
            xline = [ x[i][1], x[i+1][1] ]
            yline = [ x[i][2], x[i+1][2] ]
            Main.plot!(plt, xline, yline, label=:none, color=:black, linewidth=2)
        end

        Main.scatter!(plt,Tuple.(x),markercolor=:blue, markersize=5.0, label=:none)

        Main.plot!(plt, 
            framestyle=:box, 
            grid=:false, 
            aspect_ratio=1,
            lims=[-11,11]
        )

    end

    if type == 2

        x = [
            SVector{2,Float64}( -1, -1),
            SVector{2,Float64}( 0, 1 ),
            SVector{2,Float64}( 1, -1 ),
        ]

        y = rand(SVector{2,Float64}, 30)
        for i in 1:2:length(y)
            y[i] = -10 .+ 20 * y[i]
            ϕ = rand(0:1e-5:2π)
            y[i+1] = y[i] + SVector{2,Float64}(sin(ϕ),cos(ϕ))
        end

        box = Box(limits(x,y),10.)
        y_list = minimum_distances(y, x, 2, box)

        plt = Main.plot()

        for pair in y_list
            if pair.d < 10. 
                xline = [ y[pair.i][1], x[pair.j][1] ]
                yline = [ y[pair.i][2], x[pair.j][2] ]
                Main.plot!(plt, xline, yline, label=:none, color=:black, alpha=0.5, linestyle=:dash)
            end
        end

        for i in 1:length(x)-1
            xline = [ x[i][1], x[i+1][1] ]
            yline = [ x[i][2], x[i+1][2] ]
            Main.plot!(plt, xline, yline, label=:none, color=:black, linewidth=2)
        end

        for i in 1:2:length(y)
            xline = [ y[i][1], y[i+1][1] ]
            yline = [ y[i][2], y[i+1][2] ]
            Main.plot!(plt, xline, yline, label=:none, color=:black, linewidth=2)
        end

        Main.scatter!(plt,Tuple.(x),markercolor=:red, markersize=5.0, label=:none)
        Main.scatter!(plt,Tuple.(y),markercolor=:blue, markersize=5.0, label=:none)

        Main.plot!(plt, 
            framestyle=:box, 
            grid=:false, 
            aspect_ratio=1,
            lims=[-10,10]
        )

    end
    
    
    return plt
end
