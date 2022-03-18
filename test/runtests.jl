using Test
using StaticArrays
using CellListMap
using MolecularMinimumDistances

@testset "MolecularMinimumDistances.jl" begin

    # 
    # Single set of molecules
    #
    x = [ rand(SVector{3,Float64}) for _ in 1:135 ]

    for box in [ Box([1,1,1],0.2),
                 Box([1.0 0.2 0.0
                      0.2 1.0 0.2
                      0.0 0.0 1.0], 0.2) ]
   
        x_list_naive = MolecularMinimumDistances.naive_md(x,3,box)

        x_list = minimum_distances(x,3,box,parallel=false)
        @test x_list ≈ x_list_naive

        x_list = minimum_distances(x,3,box,parallel=true)
        @test x_list ≈ x_list_naive

        cl = CellList(x,box,parallel=false)
        molecule_of_i = MolecularMinimumDistances.molecule_indices(x,3)
        x_list = MolecularMinimumDistances.init_list(x,3) 
        minimum_distances!(molecule_of_i,x_list,box,cl,parallel=false)
        @test x_list ≈ x_list_naive

        cl = CellList(x,box,parallel=true)
        minimum_distances!(molecule_of_i,x_list,box,cl,parallel=true)
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
   
        x_list_naive, y_list_naive = MolecularMinimumDistances.naive_md(x,y,5,3,box)

        x_list, y_list = minimum_distances(x,y,5,3,box,parallel=false)
        @test x_list ≈ x_list_naive
        @test y_list ≈ y_list_naive

        x_list, y_list = minimum_distances(x,y,5,3,box,parallel=true)
        @test x_list ≈ x_list_naive
        @test y_list ≈ y_list_naive

        cl = CellList(x,y,box,parallel=false)
        molecule_of_i = MolecularMinimumDistances.molecule_indices(x,5)
        molecule_of_j = MolecularMinimumDistances.molecule_indices(y,3)
        x_list = MolecularMinimumDistances.init_list(x,5) 
        y_list = MolecularMinimumDistances.init_list(y,3) 
        minimum_distances!(molecule_of_i,molecule_of_j,x_list,y_list,box,cl,parallel=false)
        @test x_list ≈ x_list_naive
        @test y_list ≈ y_list_naive

        cl = CellList(x,y,box,parallel=true)
        minimum_distances!(molecule_of_i,molecule_of_j,x_list,y_list,box,cl,parallel=true)
        @test x_list ≈ x_list_naive
        @test y_list ≈ y_list_naive

    end

end
