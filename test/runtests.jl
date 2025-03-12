using Hakai_Project
using Test
using Aqua

@testset "Hakai_Project.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(Hakai_Project; ambiguities = false,)
    end
    # Write your tests here.
end
