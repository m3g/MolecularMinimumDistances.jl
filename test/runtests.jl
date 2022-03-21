using Test
using StaticArrays
using CellListMap

using MolecularMinimumDistances
import MolecularMinimumDistances: init_list, mol_index
import MolecularMinimumDistances: naive_md

@testset "MolecularMinimumDistances.jl" begin

    x = [ rand(SVector{3,Float64}) for _ in 1:12 ]
    #
    # Internal functions
    #
    @test length(init_list(x, i -> mol_index(i,4))) == 3
    @test length(init_list(Float64, 3)) == 3
    inds = [ 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2 ]
    @test length(init_list(x, i -> inds[i])) == 2

    # 
    # Single set of molecules
    #
    x = [ rand(SVector{3,Float64}) for _ in 1:135 ]

    for box in [ Box([1,1,1],0.2),
                 Box([1.0 0.2 0.0
                      0.2 1.0 0.2
                      0.0 0.0 1.0], 0.2) ]
   
        x_list_naive = naive_md(x,3,box)

        x_list = minimum_distances(x,3,box,parallel=false)
        @test x_list ≈ x_list_naive

        x_list = minimum_distances(x,3,box,parallel=true)
        @test x_list ≈ x_list_naive

        cl = CellList(x,box,parallel=false)
        x_list = init_list(x, i -> mol_index(i,3)) 
        minimum_distances!(i -> mol_index(i,3), x_list,box,cl,parallel=false)
        @test x_list ≈ x_list_naive

        cl = CellList(x,box,parallel=true)
        minimum_distances!(i -> mol_index(i,3),x_list,box,cl,parallel=true)
        @test x_list ≈ x_list_naive

    end

    #
    # Disjoint sets of molecules
    #
    x = [ rand(SVector{3,Float64}) for _ in 1:100 ]
    y = [ rand(SVector{3,Float64}) for _ in 1:90 ]

    for box in [ Box([1,1,1],0.2),
                 Box([1.0 0.2 0.0
                      0.2 1.0 0.2
                      0.0 0.0 1.0], 0.2) ]
   
        x_list_naive, y_list_naive = naive_md(x,y,5,3,box)

        x_list, y_list = minimum_distances(x,y,5,3,box,parallel=false)
        @test x_list ≈ x_list_naive
        @test y_list ≈ y_list_naive

        x_list, y_list = minimum_distances(x,y,5,3,box,parallel=true)
        @test x_list ≈ x_list_naive
        @test y_list ≈ y_list_naive

        cl = CellList(x,y,box,parallel=false)
        x_list = init_list(x, i -> mol_index(i,5)) 
        y_list = init_list(y, i -> mol_index(i,3)) 
        minimum_distances!(
            i -> mol_index(i,5),
            i -> mol_index(i,3),
            x_list,y_list,box,cl;
            parallel=false
        )
        @test x_list ≈ x_list_naive
        @test y_list ≈ y_list_naive

        cl = CellList(x,y,box,parallel=true)
        minimum_distances!(
            i -> mol_index(i,5),
            i -> mol_index(i,3),
            x_list,y_list,box,cl;
            parallel=true
        )
        @test x_list ≈ x_list_naive
        @test y_list ≈ y_list_naive

    end

end
