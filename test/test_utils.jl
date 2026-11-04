# Reference values computed with numpy 1.x and scipy for the routines Aeolis reproduces.
@testset "numpy and scipy conventions" begin
    @test Aeolis.pad_odd_reflect([1.0, 2.0, 4.0], 5) ==
          [-6.0, -5.0, -4.0, -2.0, 0.0, 1.0, 2.0, 4.0, 6.0, 7.0, 8.0, 10.0, 12.0]
    @test Aeolis.fftfreq(5, 0.5) ≈ [0.0, 0.4, 0.8, -0.8, -0.4]
    @test Aeolis.fftfreq(6, 2.0) ≈ [0.0, 1 / 12, 1 / 6, -0.25, -1 / 6, -1 / 12]
    @test Aeolis.np_arange(-1.2, 3.1, 0.7) ≈ [-1.2, -0.5, 0.2, 0.9, 1.6, 2.3, 3.0]
    @test length(Aeolis.np_arange(-1.2, 3.1, 0.7)) == 7
    A = Aeolis.diags(([1.0, 2.0], [3.0, 4.0, 5.0], [6.0, 7.0, 8.0, 9.0],
                      [10.0, 11.0, 12.0]),
                     (-2, -1, 0, 1), 4)
    @test Matrix(A) ==
          [6.0 10.0 0.0 0.0; 3.0 7.0 11.0 0.0; 1.0 4.0 8.0 12.0; 0.0 2.0 5.0 9.0]
    @test_throws ArgumentError Aeolis.diags(([1.0],), (1,), 4)
    # RegularGridInterpolator((y, x), z[y, x]) at (y, x) points; here z is (x, y)
    z = permutedims([1.0 2.0 3.0; 4.0 5.0 6.0])
    v = Aeolis.regular_grid_interp([0.0, 1.0, 2.0], [0.0, 2.0], z,
                                   [0.5, 1.5, 2.0, 1.0, 0.0], [0.5, 1.0, 2.0, 3.0, 0.0],
                                   -9.0)
    @test v ≈ [2.25, 4.0, 6.0, -9.0, 1.0]
    @test Aeolis.gaussian_filter1d([0.0, 0.0, 1.0, 0.0, 0.0, 3.0], 1.0, 1) ≈
          [0.05842298904073567, 0.2425067681550597, 0.4126405460900358, 0.41724041277880775,
           0.9418788466526198, 1.9273104372827414]
    g = Aeolis.gaussian_filter(permutedims([1.0 2.0 3.0; 4.0 5.0 6.0]), 0.7)
    @test permutedims(g) ≈ [1.9085786602254202 2.674088116655677 3.4395975730859334;
                            3.5604024269140675 4.325911883344324 5.0914213397745804]
    @test Aeolis.gaussian_filter([1.0, 2.0], 0.0) == [1.0, 2.0]
    @test Aeolis.pyrange(10, -3, -1) == 8:9
    @test Aeolis.pyrange(10, 2, 5) == 3:5
    @test isempty(Aeolis.pyrange(10, -1, -2))
    # rotation: a positive angle turns points clockwise
    x, y = Aeolis.rotate(1.0, 0.0, 90.0)
    @test x ≈ 0.0 atol = 1e-15
    @test y ≈ -1.0
    X, Y = Aeolis.rotate([1.0 2.0], [0.0 0.0], 90.0; origin=(1.0, 0.0))
    @test X ≈ [1.0 1.0] && Y ≈ [0.0 -1.0]
    @test Aeolis.normalize_along([1.0 3.0; 0.0 0.0], 2) == [0.25 0.75; 0.0 0.0]
    a = [-1e-12, -1.0, 1.0]
    Aeolis.prevent_tiny_negatives!(a, 1e-10)
    @test a == [0.0, -1.0, 1.0]
    m = ones(2, 2)
    Aeolis.apply_mask!(m, [2.0+1.0im 1.0; 0.0 1.0])
    @test m == [3.0 1.0; 0.0 1.0]
end

@testset "time series" begin
    c = Aeolis.CircularSeries([0.0, 100.0], [0.0, 1.0])
    @test c(50.0) ≈ 0.5
    @test c(100.0) ≈ 1.0
    # AeoLiS's interp_circular has period tmax - tmin + 1
    @test c(150.0) ≈ c(150.0 - 101.0)
    w = Aeolis.CircularStepSeries([0.0, 10.0, 20.0], [1.0, 2.0, 3.0])
    @test w(5.0) == 1.0
    @test w(10.0) == 2.0
    @test w(25.0) == w(5.0)
    @test Aeolis.clamped_interp(-1.0, [0.0, 1.0], [2.0, 4.0]) == 2.0
    @test Aeolis.clamped_interp(0.25, [0.0, 1.0], [2.0, 4.0]) == 2.5
    @test Aeolis.clamped_interp(5.0, [0.0, 1.0], [2.0, 4.0]) == 4.0
end
