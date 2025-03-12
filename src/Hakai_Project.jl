# should you ask why the last line of the docstring looks like that:
# it will show the package path when help on the package is invoked like     help?> Hakai_Project
# but will interpolate to an empty string on CI server, preventing appearing the path in the documentation built there

"""
    Package Hakai_Project v$(pkgversion(Hakai_Project))

This package refactor code from https://github.com/yozoyugen/HAKAI-fem for user to use it.

$(isnothing(get(ENV, "CI", nothing)) ? ("\n" * "Package local path: " * pathof(Hakai_Project)) : "") 
"""
module Hakai_Project
export main

using LinearAlgebra
using StaticArrays
using Printf
using Plots
using Random
using Quadmath
using Base.Threads
using FLoops
using CUDA
# Write your package code here.

include("ReadInp/readinp.jl")
include("Hakai/hakai.jl")

"""
    drawGraph(output)

Plot the output (e.g. a time history) using Plots.
"""
function drawGraph(output)
    plot(output)
    gui()
end


"""
    main()

Read the input file (from ARGS), run the simulation and produce output.
"""
function main(file)
    res = run_simulation(file)
end
end
