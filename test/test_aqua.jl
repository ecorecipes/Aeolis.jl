using Aqua

@testset "Aqua" begin
    Aqua.test_all(Aeolis; ambiguities=false, piracies=false)
end
