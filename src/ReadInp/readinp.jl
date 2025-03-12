#= 
HAKAI - A 3-dimensional finite element program
Copyright (c) 2024 Yozo Yugen

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see http://www.gnu.org/licenses/.
=#

using LinearAlgebra

#########################################################
#               Data structures (unchanged)             #
#########################################################

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

#########################################################
#          Helper Functions for Reading/Parsing         #
#########################################################

@doc raw"""
    readFile(fname::String) -> Vector{String}

Reads the entire input file into a vector of lines.
"""
function readFile(fname::String)
    f = open(fname, "r")
    lines = readlines(f)
    close(f)
    return lines
end

@doc raw"""
    findLineIndices(lines::Vector{String}, pattern::String) -> Vector{Int}

Returns all line indices in `lines` that contain the given `pattern`.
"""
function findLineIndices(lines::Vector{String}, pattern::String)
    idx = Int[]
    for (i, line) in enumerate(lines)
        if occursin(pattern, line)
            push!(idx, i)
        end
    end
    return idx
end

@doc raw"""
    parseParts(lines::Vector{String}) -> Vector{PartType}

Parses all *Part sections (with nodes/elements) into an array of `PartType`.
"""
function parseParts(lines::Vector{String})
    n = length(lines)
    # Find lines with "*Part, name="
    part_index = Int[]
    for i in 1:n
        if occursin("*Part, name=", lines[i])
            push!(part_index, i)
        end
    end
    nPart = length(part_index)

    PART = PartType[]
    for k in 1:nPart
        p = PartType("", 0, zeros(1,1), 0, zeros(Int,1,1), NsetType[], "", 0)
        push!(PART, p)

        # Extract part name
        s = replace(lines[part_index[k]], " " => "")
        ss = split(s, ",", keepempty=false)
        ss = ss[2]
        name_ = ss[ findfirst("name=", ss).stop+1 : end ]
        PART[k].name = name_

        #----- Find and parse *Node block -----
        index_node = 0
        for i in part_index[k]:n
            if occursin("*Node", lines[i])
                index_node = i
                break
            end
        end

        # Count how many nodes
        nNode = 0
        for i in (index_node+1):n
            if i > n || occursin("*", lines[i])
                break
            end
            nNode += 1
        end
        PART[k].nNode = nNode
        coordmat = zeros(nNode, 3)
        for i in 1:nNode
            if occursin("*", lines[index_node + i])
                break
            end
            s = replace(lines[index_node + i], " " => "")
            ss = split(s, ",", keepempty=false)
            coordmat[i,1] = parse(Float64, ss[2])
            coordmat[i,2] = parse(Float64, ss[3])
            coordmat[i,3] = parse(Float64, ss[4])
        end
        PART[k].coordmat = coordmat'  # make it 3 x nNode

        #----- Find and parse *Element block -----
        index_elem = 0
        for i in part_index[k]:n
            if occursin("*Element", lines[i])
                index_elem = i
                break
            end
        end

        nElement = 0
        for i in (index_elem+1):n
            if i > n || occursin("*", lines[i])
                break
            end
            nElement += 1
        end
        PART[k].nElement = nElement
        elementmat = zeros(Int, nElement, 8)

        for i in 1:nElement
            if occursin("*", lines[index_elem + i])
                break
            end
            s = replace(lines[index_elem + i], " " => "")
            ss = split(s, ",", keepempty=false)
            for j in 1:8
                elementmat[i, j] = parse(Int, ss[1 + j])
            end
        end
        PART[k].elementmat = elementmat'

        #----- Find and parse any *Nset generate block inside Part -----
        nset_index_local = Int[]
        for i in part_index[k]:n
            if occursin("*Nset", lines[i]) && occursin("generate", lines[i])
                push!(nset_index_local, i)
            end
            if occursin("*End Part", lines[i])
                break
            end
        end

        for idx_ns in nset_index_local
            ns = NsetType("", "", 0, "", 0, zeros(Int,1))
            push!(PART[k].NSET, ns)
            s = replace(lines[idx_ns], " " => "")
            ss = split(s, ",", keepempty=false)
            ss = ss[2]
            name_ = ss[ findfirst("nset=", ss).stop+1 : end ]
            ns.name = name_
            s = replace(lines[idx_ns + 1], " " => "")
            ss_gen = split(s, ",", keepempty=false)
            a = [j for j in parse(Int, ss_gen[1]):parse(Int, ss_gen[3]):parse(Int, ss_gen[2])]
            ns.nodes = a
        end

        #----- Extract material name from *Solid Section -----
        for i in part_index[k]:n
            if occursin("*Solid Section", lines[i])
                s = replace(lines[i], " " => "")
                ss = split(s, ",", keepempty=false)
                for sss in ss
                    if occursin("material=", sss)
                        name_ = sss[ findfirst("material=", sss).stop+1 : end ]
                        PART[k].material_name = name_
                        break
                    end
                end
                break
            end
        end
    end

    return PART
end

@doc raw"""
    parseInstances(lines, PART) -> Vector{InstanceType}

Finds all *Instance blocks and populates an array of InstanceType.
"""
function parseInstances(lines::Vector{String}, PART::Vector{PartType})
    n = length(lines)
    instance_index = findLineIndices(lines, "*Instance")
    instance_num = length(instance_index)

    INSTANCE = InstanceType[]
    for k in 1:instance_num
        inst = InstanceType("", "", 0, 0, String[], 0, 0, 0, 0, Int[], zeros(Int,0,0), zeros(Int,0,0), Int[], zeros(Int,0,0), Int[], Int[])
        push!(INSTANCE, inst)

        # parse instance name, part name
        s = replace(lines[instance_index[k]], " " => "")
        ss = split(s, ",", keepempty=false)
        # name=
        name_ = ss[2]
        name_ = name_[ findfirst("name=", name_).stop+1 : end ]
        inst.name = name_
        # part=
        part_ = ss[3]
        part_ = part_[ findfirst("part=", part_).stop+1 : end ]
        inst.part_name = part_

        # find matching part index
        for i in 1:length(PART)
            if PART[i].name == inst.part_name
                inst.part_id = i
                break
            end
        end

        # parse translate lines until *End Instance
        translate = String[]
        for line_i in (instance_index[k]+1):n
            if occursin("*End Instance", lines[line_i])
                break
            end
            push!(translate, replace(lines[line_i], " " => ""))
        end
        inst.translate = translate
    end
    return INSTANCE
end

@doc raw"""
    parseGlobalNsets(lines, INSTANCE) -> Vector{NsetType}

Parses any *Nset with "instance=" at the global level.
"""
function parseGlobalNsets(lines::Vector{String}, INSTANCE::Vector{InstanceType})
    nset_index = Int[]
    n = length(lines)
    for i in 1:n
        if occursin("*Nset", lines[i]) && occursin("instance=", lines[i])
            push!(nset_index, i)
        end
    end

    nset_num = length(nset_index)
    NSET = NsetType[]
    for k in 1:nset_num
        ns = NsetType("", "", 0, "", 0, Int[])
        push!(NSET, ns)

        s = replace(lines[nset_index[k]], " " => "")
        ss = split(s, ",", keepempty=false)
        # nset name
        name_ = ss[2]
        name_ = name_[ findfirst("nset=", name_).stop+1 : end ]
        ns.name = name_
        # instance name
        inst_ = ss[3]
        inst_ = inst_[ findfirst("instance=", inst_).stop+1 : end ]
        ns.instance_name = inst_

        # find instance ID
        for i in 1:length(INSTANCE)
            if ns.instance_name == INSTANCE[i].name
                ns.part_name = INSTANCE[i].part_name
                ns.part_id = INSTANCE[i].part_id
                ns.instance_id = i
            end
        end

        # parse nodes
        a = Int[]
        idx = nset_index[k]
        # check generate
        sgen = split(s, ",", keepempty=false)
        if length(sgen) == 4 && sgen[4] == "generate"
            line_ = replace(lines[idx+1], " " => "")
            ss2 = split(line_, ",", keepempty=false)
            a = [j for j in parse(Int, ss2[1]):parse(Int, ss2[3]):parse(Int, ss2[2])]
        else
            # read subsequent lines until next '*'
            for i in (idx+1):n
                if i>n || occursin("*", lines[i])
                    break
                end
                l = replace(lines[i], " " => "")
                tokens = split(l, ",", keepempty=false)
                for t in tokens
                    if !isempty(t)
                        push!(a, parse(Int, t))
                    end
                end
            end
        end
        ns.nodes = a
    end
    return NSET
end

@doc raw"""
    parseGlobalElsets(lines, INSTANCE) -> Vector{ELsetType}

Parses any *Elset with "instance=" at the global level.
"""
function parseGlobalElsets(lines::Vector{String}, INSTANCE::Vector{InstanceType})
    n = length(lines)
    elset_index = Int[]
    for i in 1:n
        if occursin("*Elset", lines[i]) && occursin("instance=", lines[i])
            push!(elset_index, i)
        end
    end

    ELSET = ELsetType[]
    for idx in elset_index
        elset_obj = ELsetType("", "", 0, "", 0, Int[])
        push!(ELSET, elset_obj)

        s = replace(lines[idx], " " => "")
        ss = split(s, ",", keepempty=false)
        # name
        name_ = ss[2]
        name_ = name_[ findfirst("elset=", name_).stop+1 : end ]
        elset_obj.name = name_

        # instance name can appear in either ss[3] or ss[4]
        # Check
        if length(ss) >= 3 && occursin("instance=", ss[3])
            name_ = ss[3]
            name_ = name_[ findfirst("instance=", name_).stop+1 : end ]
            elset_obj.instance_name = name_
        elseif length(ss) >= 4 && occursin("instance=", ss[4])
            name_ = ss[4]
            name_ = name_[ findfirst("instance=", name_).stop+1 : end ]
            elset_obj.instance_name = name_
        end

        # find instance ID
        for i in 1:length(INSTANCE)
            if elset_obj.instance_name == INSTANCE[i].name
                elset_obj.part_name = INSTANCE[i].part_name
                elset_obj.part_id   = INSTANCE[i].part_id
                elset_obj.instance_id = i
            end
        end

        # gather elements
        a = Int[]
        ssplit = split(s, ",", keepempty=false)
        idx_plus1 = idx + 1

        is_generate = (length(ssplit) >= 4 && ssplit[end] == "generate") ||
                      (length(ssplit) >= 5 && ssplit[3] == "internal" && ssplit[end] == "generate")

        if is_generate
            gen_line = replace(lines[idx_plus1], " " => "")
            ss2 = split(gen_line, ",", keepempty=false)
            a = [j for j in parse(Int, ss2[1]):parse(Int, ss2[3]):parse(Int, ss2[2])]
        else
            # read subsequent lines until next '*'
            for i in idx_plus1: length(lines)
                if i>length(lines) || occursin("*", lines[i])
                    break
                end
                line_ = replace(lines[i], " " => "")
                tokens = split(line_, ",", keepempty=false)
                for t in tokens
                    if !isempty(t)
                        push!(a, parse(Int, t))
                    end
                end
            end
        end
        elset_obj.elements = a
    end
    return ELSET
end

@doc raw"""
    parseSurfaces(lines, ELSET) -> Vector{SurfaceType}

Parses *Surface blocks and populates SurfaceType.
"""
function parseSurfaces(lines::Vector{String}, ELSET::Vector{ELsetType})
    n = length(lines)
    surface_index = findLineIndices(lines, "*Surface,")
    SURFACE = SurfaceType[]

    for idx in surface_index
        surf = SurfaceType("", String[], 0, Int[])
        push!(SURFACE, surf)

        s = replace(lines[idx], " " => "")
        ss = split(s, ",", keepempty=false)
        # name
        name_ = ss[3]
        name_ = name_[ findfirst("name=", name_).stop+1 : end ]
        surf.name = name_

        # read subsequent lines until next '*'
        a = Int[]
        c = 1
        line_idx = idx + 1
        while line_idx <= n && !occursin("*", lines[line_idx])
            l = replace(lines[line_idx], " " => "")
            tokens = split(l, ",", keepempty=false)
            # tokens[1] is elset name
            push!(surf.elset_name, tokens[1])

            # find matching elset
            for j in 1:length(ELSET)
                if surf.elset_name[c] == ELSET[j].name
                    surf.instance_id = ELSET[j].instance_id
                    append!(a, ELSET[j].elements)
                end
            end
            c += 1
            line_idx += 1
        end
        surf.elements = unique(sort(a))
    end
    return SURFACE
end

@doc raw"""
    buildGlobalModel!(PART, INSTANCE) -> (nNode, coordmat, nElement, elementmat)

Given the PART and INSTANCE arrays, iterates over each instance to:
  - apply translations/rotations
  - assemble global coord and element arrays
Returns the global node count, global coordinate matrix, global element count, and global element matrix.
"""
function buildGlobalModel!(PART::Vector{PartType}, INSTANCE::Vector{InstanceType})
    nNode = 0
    nElement = 0
    coordmat_global = Float64[]
    elementmat_global = Int[]

    for (i, inst) in enumerate(INSTANCE)
        part_id = inst.part_id
        coordmat_i = copy(PART[part_id].coordmat)
        inst.node_offset = nNode
        inst.element_offset = nElement
        inst.nNode    = PART[part_id].nNode
        inst.nElement = PART[part_id].nElement
        inst.elements = collect(1:inst.nElement)

        # Compose transformations from end to front
        T = Matrix(I, 3, 3)
        for tr in reverse(inst.translate)
            ss = split(tr, ",", keepempty=false)
            # translation
            if length(ss) == 3
                offset_ = [parse(Float64, ss[1]),
                           parse(Float64, ss[2]),
                           parse(Float64, ss[3])]
                coordmat_i .= coordmat_i .+ offset_ * ones(1, size(coordmat_i, 2))

            # rotation
            elseif length(ss) == 7
                # axis from (x1,y1,z1) to (x2,y2,z2), plus angle
                x1, y1, z1 = parse.(Float64, ss[1:3])
                x2, y2, z2 = parse.(Float64, ss[4:6])
                angle_deg  = parse(Float64, ss[7])
                nv = [x2 - x1, y2 - y1, z2 - z1]
                nv ./= norm(nv)
                n1, n2, n3 = nv
                d = angle_deg / 180.0 * π

                R = [
                    n1*n1*(1-cos(d))+cos(d)   n1*n2*(1-cos(d))-n3*sin(d)  n1*n3*(1-cos(d))+n2*sin(d);
                    n1*n2*(1-cos(d))+n3*sin(d)  n2*n2*(1-cos(d))+cos(d)   n2*n3*(1-cos(d))-n1*sin(d);
                    n1*n3*(1-cos(d))-n2*sin(d)  n2*n3*(1-cos(d))+n1*sin(d)  n3*n3*(1-cos(d))+cos(d)
                ]
                coordmat_i .= R * coordmat_i
            end
        end

        # element offset
        elementmat_i = PART[part_id].elementmat .+ nNode

        # append to global
        if i == 1
            coordmat_global = copy(vec(coordmat_i))   # store as vector, then reshape
            elementmat_global = copy(vec(elementmat_i))
        else
            coordmat_global = vcat(coordmat_global, vec(coordmat_i))
            elementmat_global = vcat(elementmat_global, vec(elementmat_i))
        end

        nNode    += PART[part_id].nNode
        nElement += PART[part_id].nElement
    end

    # Reshape to final form
    # Global coords are stored in a big vector. We want 3 x nNode
    coordmat_final = reshape(coordmat_global, 3, nNode)
    elementmat_final = reshape(elementmat_global, size(PART[1].elementmat,1), nElement)

    return (nNode, coordmat_final, nElement, elementmat_final)
end

@doc raw"""
    parseAmplitudes(lines) -> Vector{AmplitudeType}

Parses *Amplitude blocks.
"""
function parseAmplitudes(lines::Vector{String})
    amplitude_index = findLineIndices(lines, "*Amplitude")
    AMPLITUDE = AmplitudeType[]

    for idx in amplitude_index
        amp = AmplitudeType("", Float64[], Float64[])
        push!(AMPLITUDE, amp)
        s = replace(lines[idx], " " => "")
        ss = split(s, ",", keepempty=false)
        # name
        name_ = ss[2]
        name_ = name_[ findfirst("name=", name_).stop+1 : end ]
        amp.name = name_

        # parse data lines until next '*'
        i = idx + 1
        time_, value_ = Float64[], Float64[]
        while i <= length(lines) && !occursin("*", lines[i])
            line_ = replace(lines[i], " " => "")
            tokens = split(line_, ",", keepempty=false)
            # tokens come in pairs: time, value
            for j in 1:2:length(tokens)
                push!(time_, parse(Float64, tokens[j]))
                push!(value_, parse(Float64, tokens[j+1]))
            end
            i += 1
        end
        amp.time = time_
        amp.value = value_
    end

    return AMPLITUDE
end

@doc raw"""
    parseMaterials(lines) -> Vector{MaterialType}

Parses *Material blocks.
"""
function parseMaterials(lines::Vector{String})
    material_index = findLineIndices(lines, "*Material")
    MATERIAL = MaterialType[]

    for idx in material_index
        mat = MaterialType("", 0.0, 0.0, 0.0, zeros(Float64,0,0), Float64[], 0, 0.0, zeros(Float64,0,0), 0.0, zeros(1,1))
        push!(MATERIAL, mat)

        s = replace(lines[idx], " " => "")
        ss = split(s, ",", keepempty=false)
        # name
        name_ = ss[2]
        name_ = name_[ findfirst("name=", name_).stop+1 : end ]
        mat.name = name_

        # keep track of lines for plastic/ductile blocks
        plastic_index   = 0
        ductile_index   = 0

        i = idx + 1
        while i <= length(lines)
            if i>length(lines) || occursin("*Material", lines[i]) || occursin("**", lines[i])
                break
            end
            if occursin("*Density", lines[i])
                s2 = replace(lines[i+1], " " => "")
                ss2 = split(s2, ",", keepempty=false)
                mat.density = parse(Float64, ss2[1])
            elseif occursin("*Elastic", lines[i])
                s2 = replace(lines[i+1], " " => "")
                ss2 = split(s2, ",", keepempty=false)
                mat.young = parse(Float64, ss2[1])
                mat.poisson = parse(Float64, ss2[2])
            elseif occursin("*Plastic", lines[i])
                plastic_index = i
            elseif occursin("*Damage Initiation", lines[i]) && occursin("criterion=DUCTILE", lines[i])
                ductile_index = i
                mat.fracture_flag = 1
            elseif occursin("*Tensile Failure", lines[i])
                s2 = replace(lines[i+1], " " => "")
                ss2 = split(s2, ",", keepempty=false)
                mat.failure_stress = parse(Float64, ss2[1])
                mat.fracture_flag = 1
            end
            i += 1
        end

        # parse plastic
        if plastic_index > 0
            plastic_data = Float64[]
            j = plastic_index + 1
            while j <= length(lines)
                if occursin("*", lines[j])
                    break
                end
                s3 = replace(lines[j], " " => "")
                ss3 = split(s3, ",", keepempty=false)
                push!(plastic_data, parse(Float64, ss3[1]))
                push!(plastic_data, parse(Float64, ss3[2]))
                j += 1
            end
            # reshape into Nx2
            n_pairs = length(plastic_data) ÷ 2
            mat.plastic = reshape(plastic_data, 2, n_pairs)'  # shape is n_pairs x 2

            # compute hardening slopes
            Hd_ = Float64[]
            for i in 1:(size(mat.plastic,1)-1)
                num   = mat.plastic[i+1,1] - mat.plastic[i,1]
                denom = mat.plastic[i+1,2] - mat.plastic[i,2]
                push!(Hd_, num/denom)
            end
            mat.Hd = Hd_
        end

        # parse ductile
        if ductile_index > 0
            ductile_data = Float64[]
            j = ductile_index + 1
            while j <= length(lines)
                if occursin("*", lines[j])
                    break
                end
                s3 = replace(lines[j], " " => "")
                ss3 = split(s3, ",", keepempty=false)
                push!(ductile_data, parse(Float64, ss3[1]))
                push!(ductile_data, parse(Float64, ss3[2]))
                push!(ductile_data, parse(Float64, ss3[3]))
                j += 1
            end
            # reshape into Nx3
            n_triples = length(ductile_data) ÷ 3
            mat.ductile = reshape(ductile_data, 3, n_triples)'
        end
    end

    return MATERIAL
end

@doc raw"""
    assignElementMaterial!(PART, INSTANCE, MATERIAL) -> (element_material, element_instance)

Given PART, INSTANCE, and MATERIAL arrays, figure out the material assignments per element
and also record which instance each element belongs to.
"""
function assignElementMaterial!(PART::Vector{PartType},
                                INSTANCE::Vector{InstanceType},
                                MATERIAL::Vector{MaterialType})
    element_material = Int[]
    element_instance = Int[]
    for i in 1:length(INSTANCE)
        part_id = INSTANCE[i].part_id
        # match part's material with the one in the global MATERIAL list
        for j in 1:length(MATERIAL)
            if PART[part_id].material_name == MATERIAL[j].name
                PART[part_id].material_id = j
                INSTANCE[i].material_id   = j
                break
            end
        end
        mat_id = fill(PART[part_id].material_id, PART[part_id].nElement)
        append!(element_material, mat_id)
        inst_id = fill(i, PART[part_id].nElement)
        append!(element_instance, inst_id)
    end
    return (element_material, element_instance)
end

@doc raw"""
    parseStep(lines) -> (d_time, end_time)

Parses *Dynamic, Explicit parameters.
"""
function parseStep(lines::Vector{String})
    d_time, end_time = 0.0, 0.0
    for i in 1:length(lines)
        if occursin("*Dynamic, Explicit", lines[i])
            s = replace(lines[i+1], " " => "")
            ss = split(s, ",", keepempty=false)
            d_time   = parse(Float64, ss[1])
            end_time = parse(Float64, ss[2])
            break
        end
    end
    return (d_time, end_time)
end

@doc raw"""
    parseMassScaling(lines) -> Float64

Parses *Fixed Mass Scaling to obtain mass scaling factor.
"""
function parseMassScaling(lines::Vector{String})
    mass_scaling = 1.0
    for i in 1:length(lines)
        if occursin("*Fixed Mass Scaling", lines[i])
            s = replace(lines[i], " " => "")
            ss = split(s, ",", keepempty=false)
            # factor=
            fac_ = ss[2]
            fac_ = fac_[ findfirst("factor=", fac_).stop+1 : end ]
            mass_scaling = parse(Float64, fac_)
            break
        end
    end
    return mass_scaling
end

@doc raw"""
    parseBC(lines, AMPLITUDE, INSTANCE, PART, NSET) -> Vector{BCType}

Parses *Boundary blocks, linking to amplitude if specified.
"""
function parseBC(lines::Vector{String},
                 AMPLITUDE::Vector{AmplitudeType},
                 INSTANCE::Vector{InstanceType},
                 PART::Vector{PartType},
                 NSET::Vector{NsetType})

    bc_index = findLineIndices(lines, "*Boundary")
    BC = BCType[]

    for idx in bc_index
        bc = BCType("", [], Float64[], "", AmplitudeType("", Float64[], Float64[]))
        push!(BC, bc)

        # check if amplitude= is in the same line
        line_ = replace(lines[idx], " " => "")
        ss = split(line_, ",", keepempty=false)
        if length(ss) == 2 && occursin("amplitude=", ss[2])
            amp_s = ss[2]
            amp_name = amp_s[ findfirst("amplitude=", amp_s).stop+1 : end ]
            bc.amp_name = amp_name
            # link to amplitude
            for a in AMPLITUDE
                if a.name == amp_name
                    bc.amplitude = a
                    break
                end
            end
        end

        # read subsequent lines until next "*Boundary" or "**"
        line_i = idx + 1
        while line_i <= length(lines)
            if line_i>length(lines) || occursin("*Boundary", lines[line_i]) || occursin("**", lines[line_i])
                break
            end
            line2 = replace(lines[line_i], " " => "")
            tokens = split(line2, ",", keepempty=false)
            bc.Nset_name = tokens[1]
            # gather nodes
            nodes = Int[]

            # check for local part nset or global nset
            if occursin(".", bc.Nset_name)
                # parse "instanceName.nsetName"
                s2 = split(bc.Nset_name, '.', keepempty=false)
                instance_name = s2[1]
                nset_name     = s2[2]
                instance_id   = 0
                part_id       = 0
                for j in 1:length(INSTANCE)
                    if INSTANCE[j].name == instance_name
                        instance_id = j
                        part_id     = INSTANCE[j].part_id
                        break
                    end
                end
                # look up the nset inside that part
                for ns_ in PART[part_id].NSET
                    if ns_.name == nset_name
                        nodes_j = ns_.nodes .+ INSTANCE[instance_id].node_offset
                        append!(nodes, nodes_j)
                        break
                    end
                end
            else
                # global nset
                for ns_ in NSET
                    if bc.Nset_name == ns_.name
                        nodes_j = ns_.nodes .+ INSTANCE[ns_.instance_id].node_offset
                        append!(nodes, nodes_j)
                    end
                end
            end

            # interpret the DOFs and values
            if length(tokens) == 2 && occursin("ENCASTRE", tokens[1]) || occursin("ENCASTRE", tokens[2])
                # ENCASTRE => fix all dofs in [1..3]
                dof_ = Int[]
                append!(dof_, nodes*3 .- 2)
                append!(dof_, nodes*3 .- 1)
                append!(dof_, nodes*3)
                push!(bc.dof, sort(dof_))
                push!(bc.value, 0.0)
            elseif length(tokens) == 3
                # (nsetName, iDir, iComp)
                dir = parse(Int, tokens[2])
                comp = parse(Int, tokens[3])
                if comp <= 3
                    dof_ = nodes*3 .- (3-comp)
                    push!(bc.dof, dof_)
                    push!(bc.value, 0.0)
                end
            elseif length(tokens) == 4
                # (nsetName, iDir, iComp, value)
                dir = parse(Int, tokens[2])
                comp = parse(Int, tokens[3])
                val_ = parse(Float64, tokens[4])
                if comp <= 3
                    dof_ = nodes*3 .- (3-comp)
                    push!(bc.dof, dof_)
                    push!(bc.value, val_)
                end
            end

            line_i += 1
        end
    end

    return BC
end

@doc raw"""
    parseIC(lines, INSTANCE, PART, NSET) -> Vector{ICType}

Parses *Initial Conditions blocks.
"""
function parseIC(lines::Vector{String},
                 INSTANCE::Vector{InstanceType},
                 PART::Vector{PartType},
                 NSET::Vector{NsetType})

    ic_index = findLineIndices(lines, "*Initial Conditions")
    IC = ICType[]

    for idx in ic_index
        ic = ICType("", "", [], Float64[])
        push!(IC, ic)

        # parse "type="
        s = replace(lines[idx], " " => "")
        ss = split(s, ",", keepempty=false)
        if length(ss) >= 2 && occursin("type=", ss[2])
            type_ = ss[2]
            type_ = type_[ findfirst("type=", type_).stop+1 : end ]
            ic.type = type_
        end

        # read subsequent lines until next *Initial Conditions or "**"
        line_i = idx + 1
        while line_i <= length(lines)
            if line_i>length(lines) || occursin("*Initial Conditions", lines[line_i]) || occursin("**", lines[line_i])
                break
            end
            line2 = replace(lines[line_i], " " => "")
            tokens = split(line2, ",", keepempty=false)
            ic.Nset_name = tokens[1]
            # gather nodes
            nodes = Int[]

            if occursin(".", ic.Nset_name)
                # parse "instanceName.nsetName"
                s2 = split(ic.Nset_name, '.', keepempty=false)
                instance_name = s2[1]
                nset_name     = s2[2]
                instance_id   = 0
                part_id       = 0
                for j in 1:length(INSTANCE)
                    if INSTANCE[j].name == instance_name
                        instance_id = j
                        part_id     = INSTANCE[j].part_id
                        break
                    end
                end
                for ns_ in PART[part_id].NSET
                    if ns_.name == nset_name
                        nodes_j = ns_.nodes .+ INSTANCE[instance_id].node_offset
                        nodes   = nodes_j
                        break
                    end
                end
            else
                # global
                for ns_ in NSET
                    if ic.Nset_name == ns_.name
                        nodes_j = ns_.nodes .+ INSTANCE[ns_.instance_id].node_offset
                        nodes = nodes_j
                        break
                    end
                end
            end

            # parse the dof and value
            dir = parse(Int, tokens[2])
            val_ = parse(Float64, tokens[3])
            dof_ = nodes*3 .- (3-dir)
            push!(ic.dof, dof_)
            push!(ic.value, val_)

            line_i += 1
        end
    end

    return IC
end

@doc raw"""
    parseContactFlag(lines) -> Int

Checks if *Contact or *Contact Inclusions (with self-contact) is present.
Returns:
  0 = no contact
  1 = general contact
  2 = self-contact
"""
function parseContactFlag(lines::Vector{String})
    contact_flag = 0
    for line_ in lines
        if occursin("*Contact", line_)
            contact_flag = 1
            break
        end
    end

    # check for self-contact
    for line_ in lines
        if occursin("*Contact Inclusions", line_) && occursin("HAKAIoption=self-contact", line_)
            contact_flag = 2
            break
        end
    end

    return contact_flag
end

@doc raw"""
    parseContactPairs(lines, SURFACE) -> Vector{CPType}

Parses *Contact Pair blocks into CPType.
"""
function parseContactPairs(lines::Vector{String}, SURFACE::Vector{SurfaceType})
    cp_index = findLineIndices(lines, "*Contact Pair,")
    CP = CPType[]

    for idx in cp_index
        cp = CPType("", "", "", 0, 0, Int[], Int[],
                    zeros(Int,0,0), zeros(Int,0,0), Int[], Int[], Int[], Int[])
        push!(CP, cp)

        s = replace(lines[idx], " " => "")
        ss = split(s, ",", keepempty=false)
        # cpset
        name_ = ss[4]
        name_ = name_[ findfirst("cpset=", name_).stop+1 : end ]
        cp.name = name_

        # next line => surfaces
        s2 = replace(lines[idx+1], " " => "")
        ss2 = split(s2, ",", keepempty=false)
        cp.surface_name_1 = ss2[1]
        cp.surface_name_2 = ss2[2]

        # find surfaces
        for surf in SURFACE
            if cp.surface_name_1 == surf.name
                cp.instance_id_1 = surf.instance_id
                cp.elements_1    = surf.elements
            end
            if cp.surface_name_2 == surf.name
                cp.instance_id_2 = surf.instance_id
                cp.elements_2    = surf.elements
            end
        end
    end
    return CP
end

#########################################################
#                  Master read function                 #
#########################################################

@doc raw"""
    readInpFile(fname::String) -> ModelType

Top-level function that reads all lines of the .inp file, then parses out
all relevant sections into a `ModelType`.
"""
function readInpFile(fname::String)
    println("readInpFile: ", fname)

    # 1) Read lines
    lines = readFile(fname)

    # 2) Parts
    PART = parseParts(lines)

    # 3) Instances
    INSTANCE = parseInstances(lines, PART)

    # 4) Global Nsets
    NSET = parseGlobalNsets(lines, INSTANCE)

    # 5) Global Elsets
    ELSET = parseGlobalElsets(lines, INSTANCE)

    # 6) Surfaces
    SURFACE = parseSurfaces(lines, ELSET)

    # 7) Build global model (node/element assembly)
    (nNode, coordmat, nElement, elementmat) = buildGlobalModel!(PART, INSTANCE)

    # 8) Amplitudes
    AMPLITUDE = parseAmplitudes(lines)

    # 9) Materials
    MATERIAL = parseMaterials(lines)

    # 10) Assign element-material
    (element_material, element_instance) = assignElementMaterial!(PART, INSTANCE, MATERIAL)

    # 11) Dynamic step
    (d_time, end_time) = parseStep(lines)

    # 12) Mass scaling
    mass_scaling = parseMassScaling(lines)

    # 13) Boundary conditions
    BC = parseBC(lines, AMPLITUDE, INSTANCE, PART, NSET)

    # 14) Initial conditions
    IC = parseIC(lines, INSTANCE, PART, NSET)

    # 15) Contact
    contact_flag = parseContactFlag(lines)

    # 16) Contact pairs
    CP = parseContactPairs(lines, SURFACE)

    # Create model
    MODEL = ModelType(
        PART,
        INSTANCE,
        NSET,
        ELSET,
        SURFACE,
        AMPLITUDE,
        MATERIAL,
        BC,
        IC,
        CP,
        nNode,
        coordmat,
        nElement,
        elementmat,
        element_material,
        element_instance,
        d_time,
        end_time,
        mass_scaling,
        contact_flag
    )

    return MODEL
end
