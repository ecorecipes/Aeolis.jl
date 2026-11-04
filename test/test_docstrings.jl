# Every exported name this package owns has a docstring. Names re-exported from another
# package are that package's to document.
function owned_undocumented(M; allow=Symbol[])
    return filter(n -> Base.which(M, n) === M && n ∉ allow, Base.Docs.undocumented_names(M))
end

@testset "docstrings on every exported name" begin
    @test isempty(owned_undocumented(Aeolis))
end
