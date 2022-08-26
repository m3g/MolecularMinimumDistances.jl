var documenterSearchIndex = {"docs":
[{"location":"basic/#Basic-user-guide","page":"Basic use","title":"Basic user guide","text":"","category":"section"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"Here the usage of the functions that allocate the list of distances will be described. Different running modes are available depending on the expected output.","category":"page"},{"location":"basic/#Installation-and-loading","page":"Basic use","title":"Installation and loading","text":"","category":"section"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"To install MolecularMinimumDistances, first download and install Julia (1.6 or greater) from https://julialang.org/downloads/. Install and run it. Then, use:","category":"page"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"julia> import Pkg; Pkg.add(\"MolecularMinimumDistances\")\n\njulia> using MolecularMinimumDistances","category":"page"},{"location":"basic/#Example-input-files","page":"Basic use","title":"Example input files","text":"","category":"section"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"The examples here use a molecular system, but the package actually only considers the coordinates of the atoms and the number of atoms of each molecule. Thus, more general distance problems can be tackled.","category":"page"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"The input atomic positions used in the following examples can be obtained with:","category":"page"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"julia> using PDBTools\n\njulia> system = MolecularMinimumDistances.download_example() \n   Array{Atoms,1} with 62026 atoms with fields:\n   index name resname chain   resnum  residue        x        y        z  beta occup model segname index_pdb\n       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  1.00     1    PROT         1\n       2  HT1     ALA     A        1        1  -10.048  -15.427   -5.569  0.00  0.00     1    PROT         2\n       3  HT2     ALA     A        1        1   -9.488  -13.913   -5.295  0.00  0.00     1    PROT         3\n                                                       ⋮ \n   62024  OH2    TIP3     C     9339    19638   13.485   -4.534  -34.438  0.00  1.00     1    WAT2     62024\n   62025   H1    TIP3     C     9339    19638   13.218   -3.647  -34.453  0.00  1.00     1    WAT2     62025\n   62026   H2    TIP3     C     9339    19638   12.618   -4.977  -34.303  0.00  1.00     1    WAT2     62026\n","category":"page"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"The system consists of a protein (with 1463 atoms), solvated by 181 TMAO molecules (with 14 atoms each), 19338 water molecules, and some ions. ","category":"page"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"These coordinates belong to a snapshot of a simulation which was performed with cubic periodic boundary conditions, with a box side of 84.48 Angstrom. ","category":"page"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"The coordinates of each of the types of molecules can be extracted from the system array of atoms with (using PDBTools - v0.13 or greater):","category":"page"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"julia> protein = coor(system,\"protein\")\n1463-element Vector{StaticArrays.SVector{3, Float64}}:\n [-9.229, -14.861, -5.481]\n [-10.048, -15.427, -5.569]\n [-9.488, -13.913, -5.295]\n ⋮\n [6.408, -12.034, -8.343]\n [6.017, -10.967, -9.713]\n\njulia> tmao = coor(system,\"resname TMAO\")\n2534-element Vector{StaticArrays.SVector{3, Float64}}:\n [-23.532, -9.347, 19.545]\n [-23.567, -7.907, 19.381]\n [-22.498, -9.702, 20.497]\n ⋮\n [13.564, -16.517, 12.419]\n [12.4, -17.811, 12.052]\n\njulia> water = coor(system,\"water\")\n58014-element Vector{StaticArrays.SVector{3, Float64}}:\n [-28.223, 19.92, -27.748]\n [-27.453, 20.358, -27.476]\n [-27.834, 19.111, -28.148]\n ⋮\n [13.218, -3.647, -34.453]\n [12.618, -4.977, -34.303]","category":"page"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"Using these vectors of coordinates, we will illustrate the use of the current package.","category":"page"},{"location":"basic/#Shortest-distances-from-a-solute","page":"Basic use","title":"Shortest distances from a solute","text":"","category":"section"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"The simplest usage consists of finding for each molecule of one set the atoms of the other set which are closer to them. For example, here we want the atoms of the proteins which are closer to each TMAO molecule (14 atoms), within a cutoff of 12.0 Angstroms.","category":"page"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"The simulations was performed with periodic boundary conditions, in a cubic box of sides [84.48, 84.48, 84.48]. We compute the minimum distances with:","category":"page"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"julia> list = minimum_distances(\n           xpositions=tmao, # solvent\n           ypositions=protein, # solute\n           xn_atoms_per_molecule=14,\n           cutoff=12.0,\n           unitcell=[84.48, 84.48, 84.48]\n       )\n181-element Vector{MinimumDistance{Float64}}:\n MinimumDistance{Float64}(false, 0, 0, Inf)\n MinimumDistance{Float64}(false, 0, 0, Inf)\n MinimumDistance{Float64}(false, 0, 0, Inf)\n ⋮\n MinimumDistance{Float64}(false, 0, 0, Inf)\n MinimumDistance{Float64}(true, 2526, 97, 9.652277658666891)","category":"page"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"The list contains, for each molecule of TMAO, a MinimumDistance object, containing the following fields,  exemplified by printing the last entry of the list:","category":"page"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"julia> list[end]\nMinimumDistance{Float64}(true, 2526, 97, 9.652277658666891)\n\nDistance within cutoff, within_cutoff = true\nx atom of pair, i = 2526\ny atom of pair, j = 97\nDistance found, d = 9.652277658666891","category":"page"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"The fields within_cutoff, i, j, and d show if a distance was found within the cutoff, the indexes of the atoms involved in the contact, and their distance.","category":"page"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"note: Note\nIf the solute has more than one molecule, this will not be taken into  consideration in this mode. All molecules will be considered as part of the same structure (the number of atoms per molecule of the protein is not a parameter here).","category":"page"},{"location":"basic/#All-shortest-distances","page":"Basic use","title":"All shortest distances","text":"","category":"section"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"A similar call of the previous section can be used to compute, for each molecule of a set of molecules, which is the closest atom of every other molecule of another set. ","category":"page"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"In the example, we can compute for each TMAO molecule, which is the closest atom of water, and vice-versa. The difference from the previous call is that now wee need to provide the number of atoms of both TMAO and water:","category":"page"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"julia> water_list, tmao_list = minimum_distances(\n           xpositions=water,\n           ypositions=tmao,\n           xn_atoms_per_molecule=3,\n           yn_atoms_per_molecule=14,\n           unitcell=[84.48, 84.48, 84.48],\n           cutoff=12.0\n       );\n\njulia> water_list\n19338-element Vector{MinimumDistance{Float64}}:\n MinimumDistance{Float64}(true, 2, 1512, 4.779476331147592)\n MinimumDistance{Float64}(true, 6, 734, 2.9413928673334357)\n MinimumDistance{Float64}(true, 8, 859, 5.701548824661595)\n ⋮\n MinimumDistance{Float64}(true, 58010, 1728, 3.942870781549911)\n MinimumDistance{Float64}(true, 58014, 2058, 2.2003220218867936)\n\njulia> tmao_list\n181-element Vector{MinimumDistance{Float64}}:\n MinimumDistance{Float64}(true, 12, 22520, 2.1985345118965056)\n MinimumDistance{Float64}(true, 20, 33586, 2.1942841657360606)\n MinimumDistance{Float64}(true, 37, 26415, 2.1992319113726926)\n ⋮\n MinimumDistance{Float64}(true, 2512, 37323, 2.198738501959709)\n MinimumDistance{Float64}(true, 2527, 33664, 2.1985044916943015)","category":"page"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"Two lists were returned, the first containing, for each water molecule, MinimumDistance data associated to the closest TMAO molecule (meaning the atoms involved in the contact and their distance). Similarly, the second list contains, for each TMAO molecule, the MinimumDistance data associated to each TMAO molecule. ","category":"page"},{"location":"basic/#Shortest-distances-within-molecules","page":"Basic use","title":"Shortest distances within molecules","text":"","category":"section"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"There is an interface to compute the shortest distances of molecules within a set of molecules. That is, given one group of molecules, compute for each molecule which is the shortest distance among the other molecules of the same type. ","category":"page"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"A typical call would be:","category":"page"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"julia> water_list = minimum_distances(\n           positions=water,\n           n_atoms_per_molecule=3,\n           unitcell=[84.48, 84.48, 84.48],\n           cutoff=12.0\n       )\n19338-element Vector{MinimumDistance{Float64}}:\n MinimumDistance{Float64}(true, 2, 33977, 2.1997806708851724)\n MinimumDistance{Float64}(true, 4, 43684, 2.1994928961012814)\n MinimumDistance{Float64}(true, 9, 28030, 2.1997583958244142)\n ⋮\n MinimumDistance{Float64}(true, 58010, 22235, 2.1992096307537414)\n MinimumDistance{Float64}(true, 58012, 9318, 2.20003227249056)","category":"page"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"Which contains for each water molecule the atoms involved in the closest contact to any other water molecule, and the distances (within the cutoff). A pictorial representation of a result of this type is, for a simpler system:","category":"page"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"(Image: self pairs)","category":"page"},{"location":"basic/","page":"Basic use","title":"Basic use","text":"This can be used for the identification of connectivity networks, for example, or for some types of clustering.","category":"page"},{"location":"advanced/#Advanced-usage","page":"Advanced use","title":"Advanced usage","text":"","category":"section"},{"location":"advanced/#System-build-and-update","page":"Advanced use","title":"System build and update","text":"","category":"section"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"If the molecular minimum distances will be computed many times for similar systems, it is possible to construct the system and update its properties. The use of the interface of CellListMap.PeriodicSystems is required (requires CellListMap version 0.7.24 or greater). ","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"For example, let us build one system with a protein and water:","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"julia> using MolecularMinimumDistances, PDBTools\n\njulia> system = MolecularMinimumDistances.download_example();\n\njulia> protein = coor(system, \"protein\");\n\njulia> water = coor(system, \"water\");","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"We now build the CrossPairs  type of system, instead of calling the minimum_distances function directly:","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"julia> sys = CrossPairs(\n           xpositions=water, # solvent\n           ypositions=protein, # solute\n           xn_atoms_per_molecule=3,\n           cutoff=12.0,\n           unitcell=[84.48, 84.48, 84.48]\n       )\nCrossPairs system with:\n\nNumber of atoms of set x: 58014\nNumber of molecules in set x: 19338\nNumber of atoms of target structure y: 1463\nCutoff: 12.0\nunitcell: [84.48, 0.0, 0.0, 0.0, 84.48, 0.0, 0.0, 0.0, 84.48]","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"Now sys  contains the necessary arrays for computing the list of minimum distances. We use now the minimum_distances!  function (with the !), to update that list:","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"julia> minimum_distances!(sys)\n19338-element Vector{MinimumDistance{Float64}}:\n MinimumDistance{Float64}(false, 0, 0, Inf)\n MinimumDistance{Float64}(false, 0, 0, Inf)\n MinimumDistance{Float64}(false, 0, 0, Inf)\n ⋮\n MinimumDistance{Float64}(true, 58011, 383, 10.24673074692606)\n MinimumDistance{Float64}(false, 0, 0, Inf)","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"The system can be now updated: the positions, cutoff, or unitcell can be modified, with the  following interfaces:","category":"page"},{"location":"advanced/#Updating-positions","page":"Advanced use","title":"Updating positions","text":"","category":"section"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"To update the positions, modify the sys.xpositions (or ypositions)  array. We will boldy demonstrate this by making the first atom of the x set to be close to the first atom of the protein, and recomputing the distances:","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"julia> using StaticArrays\n\njulia> sys.xpositions[2] = sys.ypositions[1] + SVector(1.0,0.0,0.0);\n\njulia> minimum_distances!(sys)\n19338-element Vector{MinimumDistance{Float64}}:\n MinimumDistance{Float64}(true, 2, 4, 0.9202923448556931)\n MinimumDistance{Float64}(false, 0, 0, Inf)\n MinimumDistance{Float64}(false, 0, 0, Inf)\n ⋮\n MinimumDistance{Float64}(true, 58011, 383, 10.24673074692606)\n MinimumDistance{Float64}(false, 0, 0, Inf)","category":"page"},{"location":"advanced/#Updating-the-cutoff,-unitcell-and-parallel-flag","page":"Advanced use","title":"Updating the cutoff, unitcell and parallel flag","text":"","category":"section"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"The cutoff, unitcell and parallel data of the sys objects can be modified  directly. For example:","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"julia> sys\nCrossPairs system with:\n\nNumber of atoms of set x: 58014\nNumber of molecules in set x: 19338\nNumber of atoms of target structure y: 1463\nCutoff: 15.0\nunitcell: [100.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 100.0]\n\njulia> sys.cutoff = 10.0\n10.0\n\njulia> sys.unitcell = [84.4, 84.4, 84.4]\n3-element Vector{Float64}:\n 84.4\n 84.4\n 84.4\n\njulia> sys.parallel = false\nfalse\n\njulia> sys\nCrossPairs system with:\n\nNumber of atoms of set x: 58014\nNumber of molecules in set x: 19338\nNumber of atoms of target structure y: 1463\nCutoff: 10.0\nunitcell: [84.4, 0.0, 0.0, 0.0, 84.4, 0.0, 0.0, 0.0, 84.4]","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"note: Note\nIt is not possible to update the unitcell from a Orthorhombic to a general Triclinic cell. If the system will be Triclinic at any moment, the unitcell must be initialized with the full matrix instead of a  vector of sides.","category":"page"},{"location":"advanced/#Index-of-molecules","page":"Advanced use","title":"Index of molecules","text":"","category":"section"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"Additionally, the low level interface allows the definition of more general groups of particles, in the sense that \"molecule\" can have different number of atoms in the same set. Therefore, one needs to provide a function that returns the index of the molecule of each atom, given the index of the atom. ","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"Briefly, if a set of atoms belong to molecules of the same number of atoms, one can compute the index of each molecule using","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"mol_indices(i,n) = div((i - 1), n) + 1","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"where i is the atom index in the array of coordinates, and n is the number of atoms per molecule. This is the default assumed in the basic interface, and can be called with:","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"julia> using StaticArrays\n\njulia> x = rand(SVector{3,Float64},9); # 3 water molecules\n\njulia> mol_indices(2,3) # second atom belongs to first molecule\n1\n\njulia> mol_indices(4,3) # fourth atom belongs to second molecule\n2","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"Typically, as we will show, this function will be used for setting up molecule indexes.","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"However, more general indexing can be used. For instance, let us suppose that the 9 atoms of the x array of coordinates above belong to 2 molecules, with 4 and 5 atoms each. Then, we could define, for example:","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"julia> my_mol_indices(i) = i <= 4 ? 1 : 2\nmy_mol_indices (generic function with 1 method)\n\njulia> my_mol_indices(4)\n1\n\njulia> my_mol_indices(5)\n2","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"Since the function can close-over an array of molecular indexes, the definition can be completely general, that is:","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"julia> molecular_indexes = [ 1, 3, 3, 2, 2, 1, 3, 1, 2 ];\n\njulia> my_mol_indices(i) = molecular_indexes[i]\nmy_mol_indices (generic function with 1 method)\n\njulia> my_mol_indices(1)\n1\n\njulia> my_mol_indices(5)\n2","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"In summary, this function that given the index of the atom returns the index of the corresponding molecule must be provided in the advanced interface, and typically will be just a closure around the number of atoms per molecule, using the already available mol_indices function. ","category":"page"},{"location":"advanced/#Example","page":"Advanced use","title":"Example","text":"","category":"section"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"Let us mix water and TMAO molecules in the same set, and use a general function to compute the indices of the molecules of each atom: ","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"julia> system = MolecularMinimumDistances.download_example();\n\njulia> protein = coor(system, \"protein\");\n\njulia> tmao_and_water = select(system, \"resname TMAO or resname TIP3\")\n   Array{Atoms,1} with 60548 atoms with fields:\n   index name resname chain   resnum  residue        x        y        z  beta occup model segname index_pdb\n    1479    N    TMAO     A        1      120  -23.532   -9.347   19.545  0.00  1.00     1    TMAO      1479\n    1480   C1    TMAO     A        1      120  -23.567   -7.907   19.381  0.00  1.00     1    TMAO      1480\n    1481   C2    TMAO     A        1      120  -22.498   -9.702   20.497  0.00  1.00     1    TMAO      1481\n                                                       ⋮ \n   62024  OH2    TIP3     C     9339    19638   13.485   -4.534  -34.438  0.00  1.00     1    WAT2     62024\n   62025   H1    TIP3     C     9339    19638   13.218   -3.647  -34.453  0.00  1.00     1    WAT2     62025\n   62026   H2    TIP3     C     9339    19638   12.618   -4.977  -34.303  0.00  1.00     1    WAT2     62026\n\njulia> findfirst(at -> at.resname == \"TIP3\", tmao_and_water)\n2535","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"Thus, the tmao_and_water atom array has two different types of molecules, TMAO with 14 atoms, and water with 3 atoms.  The first atom of a water molecule is atom 2535 of the array. We extract the coordinates of the atoms with:","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"julia> coor(tmao_and_water)\n60548-element Vector{SVector{3, Float64}}:\n [-23.532, -9.347, 19.545]\n [-23.567, -7.907, 19.381]\n [-22.498, -9.702, 20.497]\n ⋮\n [13.218, -3.647, -34.453]\n [12.618, -4.977, -34.303]","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"And now we define a function that, given the index of the atom, returns the molecule to which it belongs:","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"julia> function mol_indices(i) \n           if i < 2535 # TMAO (14 atoms per molecule) \n               div(i-1,14) + 1 \n           else # water (3 atoms per molecule)\n               mol_indices(2534) + div(i-2534-1,3) + 1\n           end\n       end\nmol_indices (generic function with 1 method)","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"The function above computes the molecular indices for TMAO in the standard way, and computes the water  molecular indices by first summing the molecule index of the last TMAO molecule, and subtracting from the atomic index of water the last index of the last TMAO atom. We can test this: ","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"julia> mol_indices(14) # last atom of first TMAO\n1\n\njulia> mol_indices(15) # first atom of second TMAO\n2\n\njulia> mol_indices(2534) # last atom of last TMAO\n181\n\njulia> mol_indices(2535) # first atom of first water\n182\n\njulia> mol_indices(2537) # last atom of first water\n182\n\njulia> mol_indices(2538) # first atom of second water\n183","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"With this function, we can construct the system using it instead of the xn_atoms_per_molecule integer variable, to obtain the solvation of the protein by both TMAO and water in a single run:","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"julia> sys = CrossPairs(\n           xpositions=solvent, # solvent\n           ypositions=protein, # solute\n           xmol_indices = mol_indices,\n           cutoff=12.0,\n           unitcell=[84.48, 84.48, 84.48]\n       )\nCrossPairs system with:\n\nNumber of atoms of set x: 60548\nNumber of molecules in set x: 19519\nNumber of atoms of target structure y: 1463\nCutoff: 12.0\nunitcell: [84.48, 0.0, 0.0, 0.0, 84.48, 0.0, 0.0, 0.0, 84.48]","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"As we can see, the number of molecules is correct (the sum of the number of water and tmao molecules). And the list of minimum distances will retrive the information of the closest protein atom to all solvent molecules of the set:","category":"page"},{"location":"advanced/","page":"Advanced use","title":"Advanced use","text":"julia> minimum_distances!(sys)\n19519-element Vector{MinimumDistance{Float64}}:\n MinimumDistance{Float64}(false, 0, 0, Inf)\n MinimumDistance{Float64}(false, 0, 0, Inf)\n MinimumDistance{Float64}(false, 0, 0, Inf)\n ⋮\n MinimumDistance{Float64}(true, 60545, 383, 10.24673074692606)\n MinimumDistance{Float64}(false, 0, 0, Inf)","category":"page"},{"location":"reference/#Citation","page":"Reference","title":"Citation","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"If this package was useful, please cite the article describing the main algorithms on which it is based:","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"L. Martínez, CellListMap.jl: Efficient and customizable cell list implementation for calculation of pairwise particle properties within a cutoff. Computer Physics Communications 279, 108452 (2022). ","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"DOI: 10.1016/j.cpc.2022.108452","category":"page"},{"location":"#MolecularMinimumDistances","page":"Home","title":"MolecularMinimumDistances","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package computes the minimum distance between molecules, which are represented as arrays of coordinates in two or three dimensions. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"To understand the utility and purpose of this package, consider the image below:","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: nearest.png)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Here, there is one blue molecule, with 6 atoms, and several red molecules, with 2 atoms each. The package has identified which are the molecules of the red set that have at leat one atom within a cutoff from the atoms of the blue molecule, and annotated the corresponding atoms and the distances.","category":"page"},{"location":"#Features","page":"Home","title":"Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Fast cell-list approach, to compute minimum-distance for thousands, or millions of atoms. \nGeneral periodic boundary conditions supported. \nAdvanced mode for in-place calculations, for non-allocating iterative calls (for analysis of MD trajectories, for example).\nModes for the calculation of minimum-distances in sets of molecules.","category":"page"},{"location":"#Most-typical-use:-Understanding-solvation","page":"Home","title":"Most typical use: Understanding solvation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package was designed as the backend for computing minimum distance distribution functions, which are useful for understanding solute-solvent interactions when the molecules have complex shapes. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"The most typical scenario is that of a protein, or another macromolecule, in a box of solvent. For example, here we download a frame of a protein which was simulated in a mixture of water and TMAO: ","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using PDBTools\n\njulia> system = MolecularMinimumDistances.download_example()\n   Array{Atoms,1} with 62026 atoms with fields:\n   index name resname chain   resnum  residue        x        y        z  beta occup model segname index_pdb\n       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  1.00     1    PROT         1\n       2  HT1     ALA     A        1        1  -10.048  -15.427   -5.569  0.00  0.00     1    PROT         2\n       3  HT2     ALA     A        1        1   -9.488  -13.913   -5.295  0.00  0.00     1    PROT         3\n                                                       ⋮ \n   62024  OH2    TIP3     C     9339    19638   13.485   -4.534  -34.438  0.00  1.00     1    WAT2     62024\n   62025   H1    TIP3     C     9339    19638   13.218   -3.647  -34.453  0.00  1.00     1    WAT2     62025\n   62026   H2    TIP3     C     9339    19638   12.618   -4.977  -34.303  0.00  1.00     1    WAT2     62026","category":"page"},{"location":"","page":"Home","title":"Home","text":"Next, we extract the protein coordinates, and the TMAO coordinates:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> protein = coor(system,\"protein\")\n1463-element Vector{SVector{3, Float64}}:\n [-9.229, -14.861, -5.481]\n [-10.048, -15.427, -5.569]\n [-9.488, -13.913, -5.295]\n ⋮\n [6.408, -12.034, -8.343]\n [6.017, -10.967, -9.713]\n\njulia> tmao = coor(system,\"resname TMAO\")\n2534-element Vector{SVector{3, Float64}}:\n [-23.532, -9.347, 19.545]\n [-23.567, -7.907, 19.381]\n [-22.498, -9.702, 20.497]\n ⋮\n [13.564, -16.517, 12.419]\n [12.4, -17.811, 12.052]","category":"page"},{"location":"","page":"Home","title":"Home","text":"The system was simulated with periodic boundary conditions, with sides in this frame of [83.115, 83.044, 83.063], and this information will be provided to the minimum-distance computation.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Finally, we find all the TMAO molecules having at least one atom closer than 12 Angstroms to the protein, using the current package (TMAO has 14 atoms):","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> list = minimum_distances(\n           xpositions=tmao, # solvent\n           ypositions=protein, # solute\n           xn_atoms_per_molecule=14,\n           cutoff=12.0,\n           unitcell=[83.115, 83.044, 83.063]\n       )\n181-element Vector{MinimumDistance{Float64}}:\n MinimumDistance{Float64}(false, 0, 0, Inf)\n MinimumDistance{Float64}(false, 0, 0, Inf)\n ⋮\n MinimumDistance{Float64}(true, 2526, 97, 9.652277658666891)\n\njulia> count(x -> x.within_cutoff, list)\n33","category":"page"},{"location":"","page":"Home","title":"Home","text":"Thus, 33 TMAO molecules are within the cutoff distance from the protein, and the distances can be used to study the solvation of the protein.","category":"page"},{"location":"#Performance","page":"Home","title":"Performance","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package exists because this computation is fast. For example, let us choose the water molecules instead, and benchmark the time required to compute this set of distances:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> water = coor(system,\"resname TIP3\")\n58014-element Vector{SVector{3, Float64}}:\n [-28.223, 19.92, -27.748]\n [-27.453, 20.358, -27.476]\n [-27.834, 19.111, -28.148]\n ⋮\n [13.218, -3.647, -34.453]\n [12.618, -4.977, -34.303]\n\njulia> using BenchmarkTools\n\njulia> @btime minimum_distances(\n           xpositions=$water, # solvent\n           ypositions=$protein, # solute\n           xn_atoms_per_molecule=3,\n           cutoff=12.0,\n           unitcell=[83.115, 83.044, 83.063]\n       );\n  6.288 ms (3856 allocations: 13.03 MiB)","category":"page"},{"location":"","page":"Home","title":"Home","text":"To compare, a naive algorithm to compute the same thing takes roughly 400x more for this system size:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> @btime MolecularMinimumDistances.naive_md($water, $protein, 3, [83.115, 83.044, 83.063], 12.0);\n  2.488 s (97 allocations: 609.16 KiB)","category":"page"},{"location":"","page":"Home","title":"Home","text":"And the computation can be made faster and in-place using the more advanced interface that allows preallocation of main necessary arrays:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> sys = CrossPairs(\n           xpositions=water, # solvent\n           ypositions=protein, # solute\n           xn_atoms_per_molecule=3,\n           cutoff=12.0,\n           unitcell=[83.115, 83.044, 83.063]\n       )\nCrossPairs system with:\n\nNumber of atoms of set: 58014\nNumber of atoms of target structure: 1463\nCutoff: 12.0\nunitcell: [83.12, 0.0, 0.0, 0.0, 83.04, 0.0, 0.0, 0.0, 83.06]\nNumber of molecules in set: 4144\n\njulia> @btime minimum_distances!($sys);\n  2.969 ms (196 allocations: 22.80 KiB)","category":"page"},{"location":"","page":"Home","title":"Home","text":"The remaining allocations occur only for the launching of multiple threads:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> sys = CrossPairs(\n           xpositions=water, # solvent\n           ypositions=protein, # solute\n           xn_atoms_per_molecule=14,\n           cutoff=12.0,\n           unitcell=[83.115, 83.044, 83.063],\n           parallel=false # default is true\n       );\n\njulia> @btime minimum_distances!($sys);\n  15.249 ms (0 allocations: 0 bytes)","category":"page"},{"location":"#Details-of-the-illustration","page":"Home","title":"Details of the illustration","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The initial illustration here consists of a toy solute-solvent example, where the solute is an approximately hexagonal molecule, and the solvent is composed by 40 diatomic molecules. The toy system is built as follows:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using MolecularMinimumDistances, StaticArrays\n# x will contain the \"solvent\", composed by 40 diatomic molecules\nT = SVector{2,Float64}\nx = T[]\ncmin = T(-20,-20)\nfor i in 1:40\n    v = cmin .+ 40*rand(T)\n    push!(x, v)\n    theta = 2pi*rand()\n    push!(x, v .+ T(sin(theta),cos(theta)))\nend\n# y will contain the \"solute\", composed by an approximate hexagonal molecule\ny = [ T(1,1), T(1,-1), T(0,-1.5), T(-1,-1), T(-1,1), T(0,1.5) ]","category":"page"},{"location":"","page":"Home","title":"Home","text":"Next, we compute the minimum distances between each molecule of x (the solvent) and the solute. In the input we need to specify the number of atoms of each molecule in x, and the cutoff up to which we want the distances to be computed:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> list = minimum_distances(\n           xpositions=x,\n           ypositions=y,\n           xn_atoms_per_molecule=2,\n           unitcell=[40.0, 40.0],\n           cutoff=10.0\n       )\n40-element Vector{MinimumDistance{Float64}}:\n MinimumDistance{Float64}(true, 2, 3, 1.0764931248364737)\n MinimumDistance{Float64}(false, 0, 0, Inf)\n MinimumDistance{Float64}(false, 0, 0, Inf)\n ⋮\n MinimumDistance{Float64}(true, 74, 5, 7.899981412729262)\n MinimumDistance{Float64}(false, 0, 0, Inf)","category":"page"},{"location":"","page":"Home","title":"Home","text":"The output is a list of MinimumDistance data structures, one for each molecule in x. The true indicates that a distance smaller than the cutoff was found, and for these the indexes of the atoms in x and y associated are reported, along with the distance between them.","category":"page"},{"location":"","page":"Home","title":"Home","text":"In this example, from the 40 molecules of x, eleven had atoms closer than the cutoff to some atom of y:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> count(x -> x.within_cutoff, list)\n11","category":"page"},{"location":"","page":"Home","title":"Home","text":"We have an auxiliary function to plot the result, in this case where the \"atoms\" are bi-dimensional:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Plots\nimport MolecularMinimumDistances: plot_md!\np = plot(lims=(-20,20),framestyle=:box,grid=false,aspect_ratio=1)\nplot_md!(p, x, 2, y, 6, list, y_cycle=true)","category":"page"},{"location":"","page":"Home","title":"Home","text":"will produce the illustration plot above, in which the nearest point between the two sets is identified.","category":"page"}]
}
