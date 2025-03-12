mutable struct IntegDataType
    index::Int
    integ_stress::Array{Float64,2}
    integ_strain::Array{Float64,2}
    integ_plastic_strain::Array{Float64,2}
    integ_eq_plastic_strain::Array{Float64,1}
    integ_triax_stress::Array{Float64,1}
end

mutable struct NodeDataType
    node_stress::Array{Float64,2}
    node_strain::Array{Float64,2}
    node_plastic_strain::Array{Float64,2}
    node_eq_plastic_strain::Array{Float64,1}
    node_mises_stress::Array{Float64,1}
    node_triax_stress::Array{Float64,1}
end

mutable struct ContactValue
    A::Array{Float64,2}
    b::Array{Float64,1}
    n::Array{Float64,1}
    x::Array{Float64,1}
    p::Array{Float64,1}
    q0::Array{Float64,1}
    q1::Array{Float64,1}
    q2::Array{Float64,1}
    v1::Array{Float64,1}
    v2::Array{Float64,1}
    ve::Array{Float64,1}
    vs::Array{Float64,1}
    f::Array{Float64,1}
    fric::Array{Float64,1}
    v::Array{Float64,1}
end

mutable struct ContactTriangle
    c_nodes_i::Array{Int,1}
    c_nodes_j::Array{Int,1}
    c_triangles_eleid::Array{Int,1}
    c_triangles::Array{Int,2}
    young::Float64
end