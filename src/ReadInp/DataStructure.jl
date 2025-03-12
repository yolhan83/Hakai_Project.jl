
mutable struct NsetType
    name::String
    instance_name::String
    instance_id::Int
    part_name::String
    part_id::Int
    nodes::Array{Int,1}
end

mutable struct ELsetType
    name::String
    instance_name::String
    instance_id::Int
    part_name::String
    part_id::Int
    elements::Array{Int,1}
end

mutable struct SurfaceType
    name::String
    elset_name::Array{String,1}
    instance_id::Int
    elements::Array{Int,1}
end

mutable struct PartType
    name::String
    nNode::Int
    coordmat::Array{Float64,2}
    nElement::Int
    elementmat::Array{Int,2}
    NSET::Array{NsetType,1}
    material_name::String
    material_id::Int
end

mutable struct InstanceType
    name::String
    part_name::String
    part_id::Int
    material_id::Int
    translate::Array{String,1}
    node_offset::Int
    nNode::Int
    element_offset::Int
    nElement::Int
    elements::Array{Int,1}
    surfaces::Array{Int,2}
    sorted_surfaces::Array{Int,2}
    surfaces_eleid::Array{Int,1}
    c_triangles::Array{Int,2}
    c_triangles_eleid::Array{Int,1}
    c_nodes::Array{Int,1}
end

mutable struct AmplitudeType
    name::String
    time::Array{Float64,1}
    value::Array{Float64,1}
end

mutable struct MaterialType
    name::String
    density::Float64
    young::Float64
    poisson::Float64
    plastic::Array{Float64,2}
    Hd::Array{Float64,1}
    fracture_flag::Int
    failure_stress::Float64
    ductile::Array{Float64,2}
    G::Float64
    Dmat::Array{Float64,2}
end

mutable struct BCType
    Nset_name::String
    dof::Array{Array{Int,1},1}
    value::Vector{Float64} 
    amp_name::String
    amplitude::AmplitudeType
end

mutable struct ICType
    Nset_name::String
    type::String
    dof::Array{Array{Int,1},1}
    value::Vector{Float64}
end

mutable struct CPType
    name::String
    surface_name_1::String
    surface_name_2::String
    instance_id_1::Int
    instance_id_2::Int
    elements_1::Array{Int,1}
    elements_2::Array{Int,1}
    c_triangles_1::Array{Int,2}
    c_triangles_2::Array{Int,2}
    c_triangles_eleid_1::Array{Int,1}
    c_triangles_eleid_2::Array{Int,1}
    c_nodes_1::Array{Int,1}
    c_nodes_2::Array{Int,1}
end

mutable struct ModelType
    PART::Array{PartType,1}
    INSTANCE::Array{InstanceType,1}
    NSET::Array{NsetType,1}
    ELSET::Array{ELsetType,1}
    SURFACE::Array{SurfaceType,1}
    AMPLITUDE::Array{AmplitudeType,1}
    MATERIAL::Array{MaterialType,1}
    BC::Array{BCType,1}
    IC::Array{ICType,1}
    CP::Array{CPType,1}
    nNode::Int
    coordmat::Array{Float64,2}
    nElement::Int
    elementmat::Array{Int,2}
    element_material::Vector{Int}
    element_instance::Vector{Int}
    d_time::Float64
    end_time::Float64
    mass_scaling::Float64
    contact_flag::Int
end
