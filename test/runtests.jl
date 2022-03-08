using MolecularMinimumDistances
using Test

import PDBTools
using StaticArrays
using CellListMap

@testset "MolecularMinimumDistances.jl" begin


    atoms = PDBTools.readPDB("./simple_test.pdb")
    x = PDBTools.select(atoms,"resname HOH") 
    y = PDBTools.select(atoms,"not resname HOH")


    box = Box([1,1,1],0.1)
    cl = CellList(x,y,box)
    aux = AuxThreaded(cl)
    cell_list_parameters = MolecularMinimumDistances.CellListParameters(box,cl,aux)







end
