using Test
using TestItemRunner
using StaticArrays
using MolecularMinimumDistances
import MolecularMinimumDistances: init_list, _mol_indices

@testitem "Initalization functions" begin
    using StaticArrays
    using MolecularMinimumDistances
    import MolecularMinimumDistances: init_list, _mol_indices

    x = [rand(SVector{3,Float64}) for _ in 1:12]
    #
    # Initialization functions
    #
    @test length(init_list(x, i -> _mol_indices(i, 4))) == 3
    @test length(init_list(Float64, 3)) == 3
    inds = [1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2]
    @test length(init_list(x, i -> inds[i])) == 2

end

@testitem "Single set" begin
    using StaticArrays
    using MolecularMinimumDistances
    import MolecularMinimumDistances: init_list, _mol_indices

    # 
    # Single set of molecules
    #
    x = [rand(SVector{3,Float64}) for _ in 1:135]

    for parallel in (true, false),
        unitcell in ([1, 1, 1], [1.0 0.2 0.0; 0.2 1.0 0.2; 0.0 0.0 1.0])

        cutoff = 0.2

        x_list_naive = MolecularMinimumDistances.naive_md(x, 3, unitcell, cutoff)

        x_list = minimum_distances(xpositions=x, xn_atoms_per_molecule=3, cutoff=cutoff, unitcell=unitcell, parallel=parallel)
        @test x_list ≈ x_list_naive

        x_list = minimum_distances(xpositions=x, cutoff=cutoff, mol_indices=i -> _mol_indices(i, 3), unitcell=unitcell, parallel=parallel)
        @test x_list ≈ x_list_naive

        sys = SelfPairs(xpositions=x, cutoff=cutoff, mol_indices=i -> _mol_indices(i, 3), unitcell=unitcell, parallel=parallel)
        minimum_distances!(sys)
        @test sys.minimum_distances ≈ x_list_naive

    end

end

@testitem "Disjoint: one list" begin
    using StaticArrays
    using MolecularMinimumDistances
    import MolecularMinimumDistances: init_list, _mol_indices

    #
    # Disjoint sets: return only one list
    #
    x = [rand(SVector{3,Float64}) for _ in 1:100]
    y = [rand(SVector{3,Float64}) for _ in 1:90]

    for parallel in (true, false),
        unitcell in ([1, 1, 1], [1.0 0.2 0.0; 0.2 1.0 0.2; 0.0 0.0 1.0])

        cutoff = 0.2
        x_list_naive = MolecularMinimumDistances.naive_md(x, y, 5, unitcell, cutoff)

        x_list = minimum_distances(xpositions=x, ypositions=y, cutoff=cutoff, unitcell=unitcell, xn_atoms_per_molecule=5, parallel=parallel)
        @test x_list ≈ x_list_naive

        x_list = minimum_distances(xpositions=x, ypositions=y, cutoff=cutoff, unitcell=unitcell, xmol_indices=i -> _mol_indices(i, 5), parallel=parallel)
        @test x_list ≈ x_list_naive

        sys = CrossPairs(xpositions=x, ypositions=y, cutoff=cutoff, unitcell=unitcell, xn_atoms_per_molecule=5, parallel=parallel)
        minimum_distances!(sys)
        @test sys.minimum_distances ≈ x_list_naive

    end

end

@testitem "Disjoint: two lists" begin
    using StaticArrays
    using MolecularMinimumDistances
    import MolecularMinimumDistances: init_list, _mol_indices

    #
    # Disjoint sets of molecules
    #
    x = [rand(SVector{3,Float64}) for _ in 1:100]
    y = [rand(SVector{3,Float64}) for _ in 1:90]

    for parallel in (true, false),
        unitcell in ([1, 1, 1], [1.0 0.2 0.0; 0.2 1.0 0.2; 0.0 0.0 1.0])

        cutoff = 0.2

        x_list_naive, y_list_naive = MolecularMinimumDistances.naive_md(x, y, 5, 3, unitcell, cutoff)

        x_list, y_list = minimum_distances(
            xpositions=x,
            ypositions=y,
            cutoff=cutoff,
            unitcell=unitcell,
            xn_atoms_per_molecule=5,
            yn_atoms_per_molecule=3,
            parallel=parallel
        )
        @test x_list ≈ x_list_naive
        @test y_list ≈ y_list_naive

        sys = AllPairs(
            xpositions=x,
            ypositions=y,
            cutoff=cutoff,
            unitcell=unitcell,
            xn_atoms_per_molecule=5,
            yn_atoms_per_molecule=3,
            parallel=parallel
        )
        x_list, y_list = minimum_distances!(sys)
        @test x_list ≈ x_list_naive
        @test y_list ≈ y_list_naive

        sys = AllPairs(
            xpositions=x,
            ypositions=y,
            cutoff=cutoff,
            unitcell=unitcell,
            xmol_indices=i -> _mol_indices(i, 5),
            ymol_indices=i -> _mol_indices(i, 3),
            parallel=parallel
        )
        x_list, y_list = minimum_distances!(sys)
        @test x_list ≈ x_list_naive
        @test y_list ≈ y_list_naive

    end

end # @testitem

@run_package_tests