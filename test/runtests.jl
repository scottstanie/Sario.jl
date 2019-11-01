try
    using Dates
    import Sario
catch
    println("import error")
    include("../src/Sario.jl")
end
using Test

@testset "Sario.jl" begin
    # Write your own tests here.

    demfile = joinpath(@__DIR__, "elevation.dem")

    dem = Sario.load(demfile)
    rows, cols = (5, 6)

    @test dem == reshape(collect(Int16, 1:30), rows, cols)
    #  1   6  11  16  21  26
    #  2   7  12  17  22  27
    #  3   8  13  18  23  28
    #  4   9  14  19  24  29
    #  5  10  15  20  25  30

    @test Sario.load(demfile, (3, 4)) == dem[3, 4] == 18

    @test_throws DomainError Sario.load(demfile, (30, 40))

    @test Sario.load(demfile, (2:4, 2:4)) == dem[2:4, 2:4]
    @test Sario.load(demfile, (2:4, 2:4)) == Sario.load(demfile)[2:4, 2:4]

    expected_ints = [(Date(2015,1,1), Date(2016,1,1)),
                     (Date(2016,1,1), Date(2017,1,1))]

    @test Sario.parse_intlist_strings([("20150101", "20160101"),
                                       ("20160101", "20170101")]) == expected_ints
    @test Sario.parse_intlist_strings(["20150101_20160101.int", "20160101_20170101.int"]) == expected_ints


    # Take looks test
    a = [.1  0.01  2 ; 3  4  1 + 1im]

    @test all(Sario.take_looks(a, 2, 1) .≈ [1.55  2.005  1.5 + 0.5im])
    @test all(Sario.take_looks(a, 1, 2) .≈ [0.055; 3.5])

    out = zeros(eltype(a), 2, 1)
    @test all(Sario.take_looks!(out, a, 1, 2, 2, 1) .≈ [0.055; 3.5])

    b = reshape(collect(1:12), 3, 4)
    @test all(Sario.take_looks(b, 1, 2) .≈ [0.055; 3.5])
end
