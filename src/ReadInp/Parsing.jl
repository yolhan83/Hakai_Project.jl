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