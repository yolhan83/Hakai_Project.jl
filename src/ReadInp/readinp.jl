# =============================================================================
# Type Definitions
# =============================================================================

"""
    NsetType

Represents a node set with its name, associated instance and part information, and the node indices.
"""
mutable struct NsetType
    name::String
    instance_name::String
    instance_id::Int
    part_name::String
    part_id::Int
    nodes::Vector{Int}
end

"""
    ELsetType

Represents an element set with its name, associated instance and part information, and the element indices.
"""
mutable struct ELsetType
    name::String
    instance_name::String
    instance_id::Int
    part_name::String
    part_id::Int
    elements::Vector{Int}
end

"""
    SurfaceType

Represents a surface defined by a name, a list of element set names, an instance id, and the surface element indices.
"""
mutable struct SurfaceType
    name::String
    elset_name::Vector{String}
    instance_id::Int
    elements::Vector{Int}
end

"""
    PartType

Represents a part with nodes, elements, node sets, and material information.

Fields:
- `name`: Part name.
- `nNode`: Number of nodes.
- `coordmat`: Node coordinate matrix (stored column‐major).
- `nElement`: Number of elements.
- `elementmat`: Element connectivity matrix (stored column‐major).
- `NSET`: Vector of node sets defined in the part.
- `material_name`: Name of the material assigned.
- `material_id`: Material id (assigned later).
"""
mutable struct PartType
    name::String
    nNode::Int
    coordmat::Array{Float64,2}
    nElement::Int
    elementmat::Array{Int,2}
    NSET::Vector{NsetType}
    material_name::String
    material_id::Int
end

"""
    InstanceType

Represents an instance of a part with its own translation and offset information, and connectivity indices.

Fields:
- `name`: Instance name.
- `part_name`: Name of the part to which the instance refers.
- `part_id`: Index of the part in the PART array.
- `material_id`: Material id assigned from the part.
- `translate`: List of translation/rotation directives.
- `node_offset`: Global offset for nodes.
- `nNode`: Number of nodes in the instance.
- `element_offset`: Global offset for elements.
- `nElement`: Number of elements in the instance.
- `elements`: Element indices (local numbering).
- `surfaces`, `sorted_surfaces`, `surfaces_eleid`, `c_triangles`, `c_triangles_eleid`, `c_nodes`: Additional contact/surface data.
"""
mutable struct InstanceType
    name::String
    part_name::String
    part_id::Int
    material_id::Int
    translate::Vector{String}
    node_offset::Int
    nNode::Int
    element_offset::Int
    nElement::Int
    elements::Vector{Int}
    surfaces::Array{Int,2}
    sorted_surfaces::Array{Int,2}
    surfaces_eleid::Vector{Int}
    c_triangles::Array{Int,2}
    c_triangles_eleid::Vector{Int}
    c_nodes::Vector{Int}
end

"""
    AmplitudeType

Represents an amplitude curve with its name, time values, and corresponding amplitude values.
"""
mutable struct AmplitudeType
    name::String
    time::Vector{Float64}
    value::Vector{Float64}
end

"""
    MaterialType

Represents material properties including density, elasticity, plasticity (if any), ductile failure, etc.

Fields:
- `plastic` and `ductile` matrices store data; note that large heterogeneous arrays (e.g. Vector{Any}) are avoided.
"""
mutable struct MaterialType
    name::String
    density::Float64
    young::Float64
    poisson::Float64
    plastic::Array{Float64,2}
    Hd::Vector{Float64}
    fracture_flag::Int
    failure_stress::Float64
    ductile::Array{Float64,2}
    G::Float64
    Dmat::Array{Float64,2}
end

"""
    BCType

Represents a boundary condition with the node set name, degrees of freedom (dof) constraints,
values, amplitude name and an associated amplitude curve.
"""
mutable struct BCType
    Nset_name::String
    dof::Vector{Vector{Int}}
    value::Vector{Float64}
    amp_name::String
    amplitude::AmplitudeType
end

"""
    ICType

Represents an initial condition with the node set name, type, degrees of freedom, and prescribed values.
"""
mutable struct ICType
    Nset_name::String
    type::String
    dof::Vector{Vector{Int}}
    value::Vector{Float64}
end

"""
    CPType

Represents a contact pair (interaction) between two surfaces with related element, triangle, and node data.
"""
mutable struct CPType
    name::String
    surface_name_1::String
    surface_name_2::String
    instance_id_1::Int
    instance_id_2::Int
    elements_1::Vector{Int}
    elements_2::Vector{Int}
    c_triangles_1::Array{Int,2}
    c_triangles_2::Array{Int,2}
    c_triangles_eleid_1::Vector{Int}
    c_triangles_eleid_2::Vector{Int}
    c_nodes_1::Vector{Int}
    c_nodes_2::Vector{Int}
end

"""
    ModelType

Represents the global model, containing all parts, instances, node/element sets, surfaces, amplitudes,
materials, boundary conditions, initial conditions, contact pairs, and overall connectivity and dynamic parameters.
"""
mutable struct ModelType
    PART::Vector{PartType}
    INSTANCE::Vector{InstanceType}
    NSET::Vector{NsetType}
    ELSET::Vector{ELsetType}
    SURFACE::Vector{SurfaceType}
    AMPLITUDE::Vector{AmplitudeType}
    MATERIAL::Vector{MaterialType}
    BC::Vector{BCType}
    IC::Vector{ICType}
    CP::Vector{CPType}
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

# =============================================================================
# Helper Functions
# =============================================================================

"""
    read_file_lines(fname::String) -> Vector{String}

Read lines from the file named `fname` and return them as a vector of strings.
"""
function read_file_lines(fname::String)
    open(fname, "r") do f
        return readlines(f)
    end
end

"""
    find_indices(lines::Vector{String}, pattern::String) -> Vector{Int}

Return the indices of lines that contain the specified `pattern`.
"""
function find_indices(lines::Vector{String}, pattern::String)
    indices = Int[]
    for (i, line) in enumerate(lines)
        if occursin(pattern, line)
            push!(indices, i)
        end
    end
    return indices
end

# -----------------------------------------------------------------------------
# Parsing Functions for Each Section
# -----------------------------------------------------------------------------

"""
    parse_parts(lines::Vector{String}) -> Vector{PartType}

Parse the part sections from `lines` and return a vector of `PartType`.
"""
function parse_parts(lines::Vector{String})
    n = length(lines)
    part_indices = find_indices(lines, "*Part, name=")
    n_part = length(part_indices)
    parts = PartType[]
    for k in 1:n_part
        # Create a default part object
        part = PartType("", 0, zeros(Float64, 1,1), 0, zeros(Int,1,1), NsetType[], "", 0)
        push!(parts, part)
        # Parse part name
        s = replace(lines[part_indices[k]], " " => "")
        ss = split(s, ",", keepempty=false)
        parts[k].name = ss[2][findfirst("name=", ss[2]).stop+1:end]
        # --- Node parsing ---
        node_index = 0
        for i in part_indices[k]:n
            if occursin("*Node", lines[i])
                node_index = i
                break
            end
        end
        n_node = 0
        for i in node_index+1:n
            if occursin("*", lines[i])
                break
            end
            n_node += 1
        end
        parts[k].nNode = n_node
        coordmat = zeros(Float64, n_node, 3)
        for i in 1:n_node
            if occursin("*", lines[node_index+i])
                break
            end
            s_node = replace(lines[node_index+i], " " => "")
            ss_node = split(s_node, ",", keepempty=false)
            coordmat[i,1] = parse(Float64, ss_node[2])
            coordmat[i,2] = parse(Float64, ss_node[3])
            coordmat[i,3] = parse(Float64, ss_node[4])
        end
        parts[k].coordmat = coordmat'  # store as column major
        # --- Element parsing ---
        elem_index = 0
        for i in part_indices[k]:n
            if occursin("*Element", lines[i])
                elem_index = i
                break
            end
        end
        n_elem = 0
        for i in elem_index+1:n
            if occursin("*", lines[i])
                break
            end
            n_elem += 1
        end
        parts[k].nElement = n_elem
        elementmat = zeros(Int, n_elem, 8)
        for i in 1:n_elem
            if occursin("*", lines[elem_index+i])
                break
            end
            s_elem = replace(lines[elem_index+i], " " => "")
            ss_elem = split(s_elem, ",", keepempty=false)
            for j in 1:8
                elementmat[i,j] = parse(Int, ss_elem[1+j])
            end
        end
        parts[k].elementmat = elementmat'
        # --- Nset parsing (within part) ---
        nset_idxs = Int[]
        for i in part_indices[k]:n
            if occursin("*Nset", lines[i]) && occursin("generate", lines[i])
                push!(nset_idxs, i)
            end
            if occursin("*End Part", lines[i])
                break
            end
        end
        for i in eachindex(nset_idxs)
            ns = NsetType("", "", 0, "", 0, zeros(Int,1))
            push!(parts[k].NSET, ns)
            s_nset = replace(lines[nset_idxs[i]], " " => "")
            ss_nset = split(s_nset, ",", keepempty=false)
            parts[k].NSET[end].name = ss_nset[2][findfirst("nset=", ss_nset[2]).stop+1:end]
            s_nodes = replace(lines[nset_idxs[i]+1], " " => "")
            ss_nodes = split(s_nodes, ",", keepempty=false)
            a = [j for j in parse(Int, ss_nodes[1]):parse(Int, ss_nodes[3]):parse(Int, ss_nodes[2])]
            parts[k].NSET[end].nodes = a
        end
        # --- Material assignment parsing ---
        for i in part_indices[k]:n
            if occursin("*Solid Section", lines[i])
                s_sec = replace(lines[i], " " => "")
                ss_sec = split(s_sec, ",", keepempty=false)
                for token in ss_sec
                    if occursin("material=", token)
                        parts[k].material_name = token[findfirst("material=", token).stop+1:end]
                        break
                    end
                end
                break
            end
        end
    end
    return parts
end

"""
    parse_instances(lines::Vector{String}, parts::Vector{PartType}) -> Vector{InstanceType}

Parse instance sections from `lines` and associate each with its part.
"""
function parse_instances(lines::Vector{String}, parts::Vector{PartType})
    n = length(lines)
    instance_idxs = find_indices(lines, "*Instance")
    instance_num = length(instance_idxs)
    instances = InstanceType[]
    for k in 1:instance_num
        inst = InstanceType("", "", 0, 0, String[], 0, 0, 0, 0, zeros(Int,0), zeros(Int,0,0), zeros(Int,0,0), zeros(Int,0), zeros(Int,0,0), zeros(Int,0), zeros(Int,0))
        push!(instances, inst)
        s = replace(lines[instance_idxs[k]], " " => "")
        ss = split(s, ",", keepempty=false)
        instances[end].name = ss[2][findfirst("name=", ss[2]).stop+1:end]
        instances[end].part_name = ss[3][findfirst("part=", ss[3]).stop+1:end]
        for (i, part) in enumerate(parts)
            if part.name == instances[end].part_name
                instances[end].part_id = i
                break
            end
        end
        translate = String[]
        for i in instance_idxs[k]+1:n
            if occursin("*End Instance", lines[i])
                break
            end
            push!(translate, replace(lines[i], " " => ""))
        end
        instances[end].translate = translate
    end
    return instances
end

"""
    parse_global_nsets(lines::Vector{String}, instances::Vector{InstanceType}) -> Vector{NsetType}

Parse global node set sections (with instance information) from `lines`.
"""
function parse_global_nsets(lines::Vector{String}, instances::Vector{InstanceType})
    n = length(lines)
    nset_idxs = Int[]
    for i in 1:n
        if occursin("*Nset", lines[i]) && occursin("instance=", lines[i])
            push!(nset_idxs, i)
        end
    end
    global_nsets = NsetType[]
    for k in eachindex(nset_idxs)
        ns = NsetType("", "", 0, "", 0, zeros(Int,1))
        push!(global_nsets, ns)
        s = replace(lines[nset_idxs[k]], " " => "")
        ss = split(s, ",", keepempty=false)
        global_nsets[end].name = ss[2][findfirst("nset=", ss[2]).stop+1:end]
        global_nsets[end].instance_name = ss[3][findfirst("instance=", ss[3]).stop+1:end]
        for (i, inst) in enumerate(instances)
            if inst.name == global_nsets[end].instance_name
                global_nsets[end].part_name = inst.part_name
                global_nsets[end].part_id = inst.part_id
                global_nsets[end].instance_id = i
                break
            end
        end
        a = Int[]
        index = nset_idxs[k]
        ss_parts = split(s, ",", keepempty=false)
        if length(ss_parts) == 4 && ss_parts[4] == "generate"
            s_line = replace(lines[index+1], " " => "")
            ss_line = split(s_line, ",", keepempty=false)
            a = [j for j in parse(Int, ss_line[1]):parse(Int, ss_line[3]):parse(Int, ss_line[2])]
        else
            for i in index+1:n
                if occursin("*", lines[i])
                    break
                end
                s_line = replace(lines[i], " " => "")
                ss_line = split(s_line, ",", keepempty=false)
                for token in ss_line
                    if !isempty(token)
                        push!(a, parse(Int, token))
                    end
                end
            end
        end
        global_nsets[end].nodes = a
    end
    return global_nsets
end

"""
    parse_elsets(lines::Vector{String}, instances::Vector{InstanceType}) -> Vector{ELsetType}

Parse element set sections from `lines` and return a vector of `ELsetType`.
"""
function parse_elsets(lines::Vector{String}, instances::Vector{InstanceType})
    n = length(lines)
    elset_idxs = Int[]
    for i in 1:n
        if occursin("*Elset", lines[i]) && occursin("instance=", lines[i])
            push!(elset_idxs, i)
        end
    end
    elsets = ELsetType[]
    for k in eachindex(elset_idxs)
        elset = ELsetType("", "", 0, "", 0, zeros(Int,1))
        push!(elsets, elset)
        s = replace(lines[elset_idxs[k]], " " => "")
        ss = split(s, ",", keepempty=false)
        elsets[end].name = ss[2][findfirst("elset=", ss[2]).stop+1:end]
        if length(ss) >= 3 && occursin("instance=", ss[3])
            elsets[end].instance_name = ss[3][findfirst("instance=", ss[3]).stop+1:end]
        elseif length(ss) >= 4 && occursin("instance=", ss[4])
            elsets[end].instance_name = ss[4][findfirst("instance=", ss[4]).stop+1:end]
        end
        for (i, inst) in enumerate(instances)
            if inst.name == elsets[end].instance_name
                elsets[end].part_name = inst.part_name
                elsets[end].part_id = inst.part_id
                elsets[end].instance_id = i
                break
            end
        end
        a = Int[]
        index = elset_idxs[k]
        if length(ss) == 4 && ss[4] == "generate"
            s_line = replace(lines[index+1], " " => "")
            ss_line = split(s_line, ",", keepempty=false)
            a = [j for j in parse(Int, ss_line[1]):parse(Int, ss_line[3]):parse(Int, ss_line[2])]
        elseif length(ss) == 5 && ss[3] == "internal" && ss[5] == "generate"
            s_line = replace(lines[index+1], " " => "")
            ss_line = split(s_line, ",", keepempty=false)
            a = [j for j in parse(Int, ss_line[1]):parse(Int, ss_line[3]):parse(Int, ss_line[2])]
        elseif length(ss) == 4 && ss[3] == "internal"
            for i in index+1:n
                if occursin("*", lines[i])
                    break
                end
                s_line = replace(lines[i], " " => "")
                ss_line = split(s_line, ",", keepempty=false)
                for token in ss_line
                    if !isempty(token)
                        push!(a, parse(Int, token))
                    end
                end
            end
        end
        elsets[end].elements = a
    end
    return elsets
end

"""
    parse_surfaces(lines::Vector{String}, elsets::Vector{ELsetType}) -> Vector{SurfaceType}

Parse surface sections from `lines` and return a vector of `SurfaceType`.
"""
function parse_surfaces(lines::Vector{String}, elsets::Vector{ELsetType})
    n = length(lines)
    surface_idxs = find_indices(lines, "*Surface,")
    surfaces = SurfaceType[]
    for k in eachindex(surface_idxs)
        surf = SurfaceType("", String[], 0, zeros(Int,1))
        push!(surfaces, surf)
        s = replace(lines[surface_idxs[k]], " " => "")
        ss = split(s, ",", keepempty=false)
        surfaces[end].name = ss[3][findfirst("name=", ss[3]).stop+1:end]
        index = surface_idxs[k]
        a = Int[]
        c = 1
        for i in index+1:n
            if occursin("*", lines[i])
                break
            end
            s_line = replace(lines[i], " " => "")
            ss_line = split(s_line, ",", keepempty=false)
            push!(surfaces[end].elset_name, ss_line[1])
            for el in elsets
                if surfaces[end].elset_name[c] == el.name
                    surfaces[end].instance_id = el.instance_id
                    push!(a, el.elements)
                end
            end
            c += 1
        end
        surfaces[end].elements = unique(sort(a))
    end
    return surfaces
end

"""
    parse_amplitudes(lines::Vector{String}) -> Vector{AmplitudeType}

Parse amplitude sections from `lines` and return a vector of `AmplitudeType`.
"""
function parse_amplitudes(lines::Vector{String})
    n = length(lines)
    amplitude_idxs = find_indices(lines, "*Amplitude")
    amplitudes = AmplitudeType[]
    for k in eachindex(amplitude_idxs)
        amp = AmplitudeType("", zeros(Float64,1), zeros(Float64,1))
        push!(amplitudes, amp)
        s = replace(lines[amplitude_idxs[k]], " " => "")
        ss = split(s, ",", keepempty=false)
        amplitudes[end].name = ss[2][findfirst("name=", ss[2]).stop+1:end]
        index = amplitude_idxs[k]
        for i in index+1:n
            if occursin("*", lines[i])
                break
            end
            s_line = replace(lines[i], " " => "")
            ss_line = split(s_line, ",", keepempty=false)
            time_vals = Float64[]
            value_vals = Float64[]
            for j in 1:Int(length(ss_line)/2)
                push!(time_vals, parse(Float64, ss_line[2*j-1]))
                push!(value_vals, parse(Float64, ss_line[2*j]))
            end
            amplitudes[end].time = time_vals
            amplitudes[end].value = value_vals
        end
    end
    return amplitudes
end

"""
    parse_materials(lines::Vector{String}, parts::Vector{PartType}) -> Vector{MaterialType}

Parse material sections from `lines` and return a vector of `MaterialType`. Also assign material IDs to parts.
"""
function parse_materials(lines::Vector{String}, parts::Vector{PartType})
    n = length(lines)
    material_idxs = find_indices(lines, "*Material")
    materials = MaterialType[]
    for k in eachindex(material_idxs)
        mat = MaterialType("", 0.0, 0.0, 0.0, zeros(Float64, 0,0), zeros(Float64, 0), 0, 0.0, zeros(Float64, 0,0), 0.0, zeros(1,1))
        push!(materials, mat)
        s = replace(lines[material_idxs[k]], " " => "")
        ss = split(s, ",", keepempty=false)
        materials[end].name = ss[2][findfirst("name=", ss[2]).stop+1:end]
        index = material_idxs[k]
        plastic_index = 1
        ductile_index = 1
        for i in index+1:n
            if occursin("*Material", lines[i]) || occursin("**", lines[i])
                break
            end
            if occursin("*Density", lines[i])
                s_line = replace(lines[i+1], " " => "")
                ss_line = split(s_line, ",", keepempty=false)
                materials[end].density = parse(Float64, ss_line[1])
            end
            if occursin("*Elastic", lines[i])
                s_line = replace(lines[i+1], " " => "")
                ss_line = split(s_line, ",", keepempty=false)
                materials[end].young = parse(Float64, ss_line[1])
                materials[end].poisson = parse(Float64, ss_line[2])
            end
            if occursin("*Plastic", lines[i])
                plastic_index = i
            end
            if occursin("*Damage Initiation", lines[i]) && occursin("criterion=DUCTILE", lines[i])
                ductile_index = i
                materials[end].fracture_flag = 1
            end
            if occursin("*Tensile Failure", lines[i])
                s_line = replace(lines[i+1], " " => "")
                ss_line = split(s_line, ",", keepempty=false)
                materials[end].failure_stress = parse(Float64, ss_line[1])
                materials[end].fracture_flag = 1
            end
        end
        plastic = Any[]
        if plastic_index > material_idxs[k]
            for i in plastic_index+1:n
                if occursin("*", lines[i])
                    break
                end
                s_line = replace(lines[i], " " => "")
                ss_line = split(s_line, ",", keepempty=false)
                push!(plastic, [parse(Float64, ss_line[1])  parse(Float64, ss_line[2])])
            end
            m_ = zeros(Float64, length(plastic), 2)
            for i in eachindex(plastic)
                m_[i,1] = plastic[i][1]
                m_[i,2] = plastic[i][2]
            end
            materials[end].plastic = m_
        end
        Hd = Float64[]
        for i in 1:size(materials[end].plastic, 1)-1
            v = (materials[end].plastic[i+1,1] - materials[end].plastic[i,1]) /
                (materials[end].plastic[i+1,2] - materials[end].plastic[i,2])
            push!(Hd, v)
        end
        materials[end].Hd = Hd
        ductile = Any[]
        if ductile_index > material_idxs[k]
            for i in ductile_index+1:n
                if occursin("*", lines[i])
                    break
                end
                s_line = replace(lines[i], " " => "")
                ss_line = split(s_line, ",", keepempty=false)
                push!(ductile, [parse(Float64, ss_line[1])  parse(Float64, ss_line[2])  parse(Float64, ss_line[3])])
            end
            m_d = zeros(Float64, length(ductile), 3)
            for i in eachindex(ductile)
                m_d[i,1] = ductile[i][1]
                m_d[i,2] = ductile[i][2]
                m_d[i,3] = ductile[i][3]
            end
            materials[end].ductile = m_d
        end
    end
    # Assign material IDs to parts
    for (i, part) in enumerate(parts)
        for (j, material) in enumerate(materials)
            if part.material_name == material.name
                part.material_id = j
                break
            end
        end
    end
    return materials
end

"""
    parse_global_model_data(lines::Vector{String}, parts::Vector{PartType}, instances::Vector{InstanceType})
        -> Tuple{Int, Array{Float64,2}, Int, Array{Int,2}, Vector{Int}, Vector{Int}}

Parse the global model data including node coordinates, element connectivity, and assign element material
and instance information. Returns a tuple:
(nNode, coordmat, nElement, elementmat, element_material, element_instance)
"""
function parse_global_model_data(lines::Vector{String}, parts::Vector{PartType}, instances::Vector{InstanceType})
    n = length(lines)
    n_node_total = 0
    n_elem_total = 0
    global_coordmat = Float64[]
    global_elementmat = Int[]
    element_material = Int[]
    element_instance = Int[]
    for i in eachindex(instances)
        part_id = instances[i].part_id
        coordmat_i = parts[part_id].coordmat
        instances[i].node_offset = n_node_total
        instances[i].element_offset = n_elem_total
        instances[i].nNode = parts[part_id].nNode
        instances[i].nElement = parts[part_id].nElement
        instances[i].elements = collect(1:instances[i].nElement)
        T = Matrix(I, 3, 3)
        for j in length(instances[i].translate):-1:1
            s_trans = instances[i].translate[j]
            ss_trans = split(s_trans, ",", keepempty=false)
            if length(ss_trans) == 3
                offset_ = [parse(Float64, ss_trans[1]),
                           parse(Float64, ss_trans[2]),
                           parse(Float64, ss_trans[3])]'  # row vector
                coordmat_i = coordmat_i + offset_ * ones(1, size(coordmat_i, 2))
            elseif length(ss_trans) == 7
                nv = [parse(Float64, ss_trans[4]) - parse(Float64, ss_trans[1]),
                      parse(Float64, ss_trans[5]) - parse(Float64, ss_trans[2]),
                      parse(Float64, ss_trans[6]) - parse(Float64, ss_trans[3])]
                nv = nv / norm(nv)
                n1, n2, n3 = nv
                d = parse(Float64, ss_trans[7]) / 180.0 * pi
                T = [n1*n1*(1-cos(d))+cos(d)   n1*n2*(1-cos(d))-n3*sin(d)   n1*n3*(1-cos(d))+n2*sin(d);
                     n1*n2*(1-cos(d))+n3*sin(d)   n2*n2*(1-cos(d))+cos(d)   n2*n3*(1-cos(d))-n1*sin(d);
                     n1*n3*(1-cos(d))-n2*sin(d)   n2*n3*(1-cos(d))+n1*sin(d)   n3*n3*(1-cos(d))+cos(d)]
                coordmat_i = T * coordmat_i
            end
        end
        elementmat_i = parts[part_id].elementmat .+ n_node_total
        if i == 1
            global_coordmat = coordmat_i
            global_elementmat = elementmat_i
        else
            global_coordmat = [global_coordmat  coordmat_i]
            global_elementmat = [global_elementmat  elementmat_i]
        end
        n_node_total += parts[part_id].nNode
        n_elem_total += parts[part_id].nElement
        mat_id = fill(parts[part_id].material_id, parts[part_id].nElement)
        append!(element_material, mat_id)
        inst_ids = fill(i, parts[part_id].nElement)
        push!(element_instance, inst_ids)
    end
    return (n_node_total, global_coordmat, n_elem_total, global_elementmat, element_material, element_instance)
end

"""
    parse_dynamic_step(lines::Vector{String}) -> Tuple{Float64, Float64}

Parse the dynamic step parameters (d_time and end_time) from `lines`.
"""
function parse_dynamic_step(lines::Vector{String})
    d_time = 0.0
    end_time = 0.0
    n = length(lines)
    for i in 1:n
        if occursin("*Dynamic, Explicit", lines[i])
            s_line = replace(lines[i+1], " " => "")
            ss_line = split(s_line, ",", keepempty=false)
            d_time = parse(Float64, ss_line[1])
            end_time = parse(Float64, ss_line[2])
            break
        end
    end
    return (d_time, end_time)
end

"""
    parse_mass_scaling(lines::Vector{String}) -> Float64

Parse the mass scaling factor from `lines`.
"""
function parse_mass_scaling(lines::Vector{String})
    mass_scaling = 1.0
    n = length(lines)
    for i in 1:n
        if occursin("*Fixed Mass Scaling", lines[i])
            s_line = replace(lines[i], " " => "")
            ss_line = split(s_line, ",", keepempty=false)
            ss_line = ss_line[2]
            factor_str = ss_line[findfirst("factor=", ss_line).stop+1:end]
            mass_scaling = parse(Float64, factor_str)
            break
        end
    end
    return mass_scaling
end

"""
    parse_bc(lines::Vector{String}, global_nsets::Vector{NsetType},
             amplitudes::Vector{AmplitudeType}, instances::Vector{InstanceType},
             parts::Vector{PartType}) -> Vector{BCType}

Parse boundary condition sections from `lines` and return a vector of `BCType`.
"""
function parse_bc(lines::Vector{String}, global_nsets::Vector{NsetType},
                  amplitudes::Vector{AmplitudeType}, instances::Vector{InstanceType},
                  parts::Vector{PartType})
    n = length(lines)
    bc_idxs = find_indices(lines, "*Boundary")
    bc_conditions = BCType[]
    for k in eachindex(bc_idxs)
        bc = BCType("", Vector{Vector{Int}}(), Float64[], "", AmplitudeType("", zeros(Float64,1), zeros(Float64,1)))
        push!(bc_conditions, bc)
        index = bc_idxs[k]
        s = replace(lines[index], " " => "")
        ss = split(s, ",", keepempty=false)
        if length(ss) == 2 && occursin("amplitude=", ss[2])
            token = ss[2]
            amp_name = token[findfirst("amplitude=", token).stop+1:end]
            bc_conditions[end].amp_name = amp_name
            for amp in amplitudes
                if amp.name == amp_name
                    bc_conditions[end].amplitude = amp
                    break
                end
            end
        end
        for i in index+1:n
            if occursin("*Boundary", lines[i]) || occursin("**", lines[i])
                break
            end
            s_line = replace(lines[i], " " => "")
            ss_line = split(s_line, ",", keepempty=false)
            bc_conditions[end].Nset_name = ss_line[1]
            nodes = Int[]
            if occursin(".", bc_conditions[end].Nset_name)
                sss = split(bc_conditions[end].Nset_name, ".", keepempty=false)
                instance_name = sss[1]
                nset_name = sss[2]
                instance_id = 0
                part_id = 0
                for j in eachindex(instances)
                    if instances[j].name == instance_name
                        instance_id = j
                        part_id = instances[j].part_id
                        break
                    end
                end
                for j in eachindex(parts[part_id].NSET)
                    if parts[part_id].NSET[j].name == nset_name
                        nodes_j = parts[part_id].NSET[j].nodes .+ instances[instance_id].node_offset
                        push!(nodes, nodes_j)
                        break
                    end
                end
            else
                for ns in global_nsets
                    if bc_conditions[end].Nset_name == ns.name
                        nodes_j = ns.nodes .+ instances[ns.instance_id].node_offset
                        push!(nodes, nodes_j)
                    end
                end
            end
            dof = Int[]
            if length(ss_line) == 2 && occursin("ENCASTRE", ss_line[2])
                push!(dof, nodes*3 .- 2)
                push!(dof, nodes*3 .- 1)
                push!(dof, nodes*3)
                push!(bc_conditions[end].dof, dof)
                bc_conditions[end].value = [0.0]
            elseif length(ss_line) == 3
                dir = parse(Int, ss_line[3])
                if dir <= 3
                    dof = nodes*3 .- (3 - dir)
                    push!(bc_conditions[end].dof, dof)
                    push!(bc_conditions[end].value, 0.0)
                end
            elseif length(ss_line) == 4
                dir = parse(Int, ss_line[3])
                value = parse(Float64, ss_line[4])
                if dir <= 3
                    dof = nodes*3 .- (3 - dir)
                    push!(bc_conditions[end].dof, dof)
                    push!(bc_conditions[end].value, value)
                end
            end
        end
    end
    return bc_conditions
end

"""
    parse_ic(lines::Vector{String}, global_nsets::Vector{NsetType},
             instances::Vector{InstanceType}, parts::Vector{PartType}) -> Vector{ICType}

Parse initial condition sections from `lines` and return a vector of `ICType`.
"""
function parse_ic(lines::Vector{String}, global_nsets::Vector{NsetType},
                  instances::Vector{InstanceType}, parts::Vector{PartType})
    n = length(lines)
    ic_idxs = find_indices(lines, "*Initial Conditions")
    ics = ICType[]
    for k in eachindex(ic_idxs)
        ic = ICType("", "", Vector{Vector{Int}}(), Float64[])
        push!(ics, ic)
        index = ic_idxs[k]
        s = replace(lines[index], " " => "")
        ss = split(s, ",", keepempty=false)
        ics[end].type = ss[2][findfirst("type=", ss[2]).stop+1:end]
        for i in index+1:n
            if occursin("*Initial Conditions", lines[i]) || occursin("**", lines[i])
                break
            end
            s_line = replace(lines[i], " " => "")
            ss_line = split(s_line, ",", keepempty=false)
            ics[end].Nset_name = ss_line[1]
            nodes = Int[]
            if occursin(".", ics[end].Nset_name)
                sss = split(ics[end].Nset_name, ".", keepempty=false)
                instance_name = sss[1]
                nset_name = sss[2]
                instance_id = 0
                part_id = 0
                for j in eachindex(instances)
                    if instances[j].name == instance_name
                        instance_id = j
                        part_id = instances[j].part_id
                        break
                    end
                end
                for j in eachindex(parts[part_id].NSET)
                    if parts[part_id].NSET[j].name == nset_name
                        nodes = parts[part_id].NSET[j].nodes .+ instances[instance_id].node_offset
                        break
                    end
                end
            else
                for ns in global_nsets
                    if ics[end].Nset_name == ns.name
                        nodes = ns.nodes .+ instances[ns.instance_id].node_offset
                        break
                    end
                end
            end
            dir = parse(Int, ss_line[2])
            dof = nodes*3 .- (3 - dir)
            push!(ics[end].dof, dof)
            v = parse(Float64, ss_line[3])
            push!(ics[end].value, v)
        end
    end
    return ics
end

"""
    parse_contact_and_cp(lines::Vector{String}, surfaces::Vector{SurfaceType})
        -> Tuple{Int, Vector{CPType}}

Parse contact flags and contact pair sections from `lines`. Returns a tuple containing the contact flag
and a vector of `CPType`.
"""
function parse_contact_and_cp(lines::Vector{String}, surfaces::Vector{SurfaceType})
    n = length(lines)
    contact_flag = 0
    for i in 1:n
        if occursin("*Contact", lines[i])
            contact_flag = 1
            break
        end
    end
    for i in 1:n
        if occursin("*Contact Inclusions", lines[i]) && occursin("HAKAIoption=self-contact", lines[i])
            contact_flag = 2
            break
        end
    end
    cp_idxs = Int[]
    for i in 1:n
        if occursin("*Contact Pair,", lines[i])
            push!(cp_idxs, i)
        end
    end
    cp_conditions = CPType[]
    for k in eachindex(cp_idxs)
        cp = CPType("", "", "", 0, 0, Int[], Int[],
                    zeros(Int,0,0), zeros(Int,0,0), zeros(Int,0), zeros(Int,0),
                    zeros(Int,0), zeros(Int,0))
        push!(cp_conditions, cp)
        index = cp_idxs[k]
        s = replace(lines[index], " " => "")
        ss = split(s, ",", keepempty=false)
        cp_conditions[end].name = ss[4][findfirst("cpset=", ss[4]).stop+1:end]
        s_next = replace(lines[index+1], " " => "")
        ss_next = split(s_next, ",", keepempty=false)
        cp_conditions[end].surface_name_1 = ss_next[1]
        cp_conditions[end].surface_name_2 = ss_next[2]
        for srf in surfaces
            if cp_conditions[end].surface_name_1 == srf.name
                cp_conditions[end].instance_id_1 = srf.instance_id
                cp_conditions[end].elements_1 = srf.elements
            end
            if cp_conditions[end].surface_name_2 == srf.name
                cp_conditions[end].instance_id_2 = srf.instance_id
                cp_conditions[end].elements_2 = srf.elements
            end
        end
    end
    return (contact_flag, cp_conditions)
end

# =============================================================================
# Main Function
# =============================================================================

"""
    read_inp_file(fname::String) -> ModelType

Read an Abaqus inp file from `fname` and return a `ModelType` structure containing
all parsed parts, instances, sets, surfaces, amplitudes, materials, BC, IC, contact pairs,
and global model data.
"""
function read_inp_file(fname::String)
    println("readInpFile: ", fname)
    lines = read_file_lines(fname)
    # Parse each section of the inp file
    parts         = parse_parts(lines)
    instances     = parse_instances(lines, parts)
    global_nsets  = parse_global_nsets(lines, instances)
    elsets        = parse_elsets(lines, instances)
    surfaces      = parse_surfaces(lines, elsets)
    amplitudes    = parse_amplitudes(lines)
    materials     = parse_materials(lines, parts)
    n_node, coordmat, n_elem, elementmat, element_material, element_instance =
                        parse_global_model_data(lines, parts, instances)
    d_time, end_time = parse_dynamic_step(lines)
    mass_scaling  = parse_mass_scaling(lines)
    bc            = parse_bc(lines, global_nsets, amplitudes, instances, parts)
    ic            = parse_ic(lines, global_nsets, instances, parts)
    contact_flag, cp = parse_contact_and_cp(lines, surfaces)
    # Assemble and return the global model
    model = ModelType(parts, instances, global_nsets, elsets, surfaces, amplitudes,
                      materials, bc, ic, cp, n_node, coordmat, n_elem, elementmat,
                      element_material, element_instance, d_time, end_time, mass_scaling,
                      contact_flag)
    return model
end