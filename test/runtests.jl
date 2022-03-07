using MolecularMinimumDistances
using Test

using StaticArrays
using CellListMap

@testset "MolecularMinimumDistances.jl" begin

    x = [ rand(SVector{3,Float64}) for i in 1:30 ]
    y = [ rand(SVector{3,Float64}) for i in 1:50 ]
    x_molecules = [ j for j in 1:10 for i in 1:3 ]
    y_molecules = [ j for j in 1:10 for i in 1:5 ]

    box = Box([1,1,1],0.1)
    cl = CellList(x,y,box)
    aux = AuxThreaded(cl)
    cell_list_parameters = MolecularMinimumDistances.CellListParameters(bo,cl,aux)







end
