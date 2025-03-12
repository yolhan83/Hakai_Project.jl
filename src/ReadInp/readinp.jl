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

include("DataStructure.jl")
include("Parsing.jl")


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
