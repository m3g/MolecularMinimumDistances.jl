using MolecularMinimumDistances
using Test

import PDBTools
using StaticArrays
using CellListMap


@testset "MolecularMinimumDistances.jl" begin

    x = [ rand(SVector{2,Float64}) for _ in 1:10 ]
    y = [ rand(SVector{2,Float64}) for _ in 1:8 ]

    box = Box([1,1,1],0.1)
    cl = CellList(x,y,box)
    aux = AuxThreaded(cl)
    cell_list_parameters = MolecularMinimumDistances.CellListParameters(box,cl,aux)







end
