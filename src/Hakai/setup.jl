
@doc raw"""
    prepareMaterialProperties!(MODEL)

Loops over MODEL.MATERIAL to compute and store properties like G, Dmat,
and detect if there's any fracture. Returns a flag indicating whether
fracture is enabled in the model.
"""
function prepareMaterialProperties!(MODEL)
    fracture_flag = 0
    for i in 1:length(MODEL.MATERIAL)
        mat = MODEL.MATERIAL[i]
        young   = mat.young
        poisson = mat.poisson

        # Shear modulus
        G = young / (2.0 * (1.0 + poisson))
        mat.G = G

        # Build Dmat (6x6)
        d1 = 1.0 - poisson
        d2 = poisson
        d3 = (1.0 - 2.0*poisson) / 2.0
        fac = young / ((1.0 + poisson)*(1.0 - 2.0*poisson))

        Dmat_i = fac * [
            d1 d2 d2  0.0 0.0 0.0
            d2 d1 d2  0.0 0.0 0.0
            d2 d2 d1  0.0 0.0 0.0
            0.0 0.0 0.0 d3  0.0 0.0
            0.0 0.0 0.0 0.0 d3  0.0
            0.0 0.0 0.0 0.0 0.0 d3
        ]
        mat.Dmat = Dmat_i

        # Check fracture
        if length(mat.failure_stress) > 0
            fracture_flag = 1
        end
        if size(mat.ductile,1) > 0
            fracture_flag = 1
        end
    end

    return fracture_flag
end


@doc raw"""
    computeElementVolumes(MODEL, Pusai_mat; integ_num=8)

Given the shape function matrix (Pusai_mat) and the MODEL's global coords,
compute the volume for each element (summing over integration points).
Returns an array of volumes, one per element.
"""
function computeElementVolumes(MODEL, Pusai_mat; integ_num=8)
    nElement   = MODEL.nElement
    coordmat   = MODEL.coordmat
    elementmat = MODEL.elementmat

    elementVolume = zeros(nElement)
    for e in 1:nElement
        # gather element node coords
        e_pos = coordmat[:, elementmat[:,e]]
        vol = 0.0
        for i in 1:integ_num
            J = Pusai_mat[i] * e_pos'
            vol += my3det(J)
        end
        elementVolume[e] = vol
    end
    return elementVolume
end


@doc raw"""
    buildMassDamping(MODEL, elementVolume; mass_scaling=1.0)

Constructs the diagonal mass and damping arrays diag_M, diag_C
based on element densities and volumes. Applies mass scaling.
"""
function buildMassDamping(MODEL, elementVolume; mass_scaling=1.0)
    nNode = MODEL.nNode
    fn = nNode * 3

    diag_M = zeros(fn)
    diag_C = zeros(fn)

    elementmat      = MODEL.elementmat
    elementMaterial = MODEL.element_material
    nElement        = MODEL.nElement

    @views for i in 1:nElement
        mat_id  = elementMaterial[i]
        density = MODEL.MATERIAL[mat_id].density
        vol     = elementVolume[i]
        node_mass = density * vol / 8.0

        # add to diagonal mass
        nodes_i = elementmat[:, i]
        diag_M[(nodes_i .-1).*3 .+ 1] .+= node_mass
        diag_M[(nodes_i .-1).*3 .+ 2] .+= node_mass
        diag_M[(nodes_i .-1).*3 .+ 3] .+= node_mass
    end

    # apply mass scaling
    diag_M .*= mass_scaling

    # for now, damping is just zero or some fraction:
    diag_C .= diag_M .* 0.0   # or set some fraction of M

    return diag_M, diag_C
end


@doc raw"""
    initDisplacementVelocity(MODEL, nNode, d_time, diag_M, mass_scaling)

Creates and returns arrays for disp, disp_new, disp_pre, velo, and position.
Also applies initial conditions from MODEL.IC if present.
"""
function initDisplacementVelocity(MODEL, nNode, d_time, diag_M, mass_scaling)
    fn       = nNode * 3
    coordmat = MODEL.coordmat

    disp     = zeros(fn)
    disp_new = zeros(fn)
    disp_pre = zeros(fn)
    velo     = zeros(fn)
    position = copy(coordmat)

    # apply initial velocity or displacement from MODEL.IC
    for ic in MODEL.IC
        for j in 1:length(ic.dof)
            dof_list = ic.dof[j]  # array of DOFs
            val      = ic.value[j]

            # set disp_pre and velo
            disp_pre[dof_list] .= -val * d_time
            velo[dof_list]     .=  val
        end
    end

    return disp, disp_new, disp_pre, velo, position
end


@doc raw"""
    setupContact!(MODEL)

Handles logic to build up contact pairs, surfaces, etc.
Returns a struct or NamedTuple with contact data (CT, instance_pair, cp_index, etc.).
"""
function setupContact!(MODEL)
    contact_flag = MODEL.contact_flag
    if contact_flag < 1
        return nothing  # no contact
    end

    # We replicate the logic from 'hakai' that sets up surfaces in each instance
    for i in 1:length(MODEL.INSTANCE)
        faces, faces_eleid, sorted_faces = get_element_face(MODEL, i)
        MODEL.INSTANCE[i].surfaces       = faces
        MODEL.INSTANCE[i].surfaces_eleid = faces_eleid
        MODEL.INSTANCE[i].sorted_surfaces= sorted_faces
    end

    # If CP is empty, do the "all exterior" logic, otherwise parse the existing CP
    all_exterior_flag = 0
    if length(MODEL.CP) == 0
        all_exterior_flag = 1
        # build a CP for all instance pairs
        newCP = CPType[]
        ni = length(MODEL.INSTANCE)
        c = 1
        if ni > 1
            for i in 1:ni
                js = i+1
                if contact_flag == 2
                    js = i
                end
                for j in js:ni
                    cp = CPType("", "", "", 0, 0, Int[], Int[],
                                zeros(Int,0,0), zeros(Int,0,0),
                                zeros(Int,0), zeros(Int,0), zeros(Int,0), zeros(Int,0))
                    cp.instance_id_1 = i
                    cp.instance_id_2 = j
                    cp.elements_1    = collect(1:MODEL.INSTANCE[i].nElement)
                    cp.elements_2    = collect(1:MODEL.INSTANCE[j].nElement)
                    push!(newCP, cp)
                    c+=1
                end
            end
        else
            # single instance => self-contact
            cp = CPType("", "", "", 0, 0, Int[], Int[],
                        zeros(Int,0,0), zeros(Int,0,0),
                        zeros(Int,0), zeros(Int,0), zeros(Int,0), zeros(Int,0))
            cp.instance_id_1 = 1
            cp.instance_id_2 = 1
            cp.elements_1 = collect(1:MODEL.INSTANCE[1].nElement)
            cp.elements_2 = collect(1:MODEL.INSTANCE[1].nElement)
            push!(newCP, cp)
        end
        MODEL.CP = newCP
    end

    # Now fill c_triangles_1, c_triangles_2, c_nodes_1, c_nodes_2 for each CP
    for i in 1:length(MODEL.CP)
        cp = MODEL.CP[i]
        iinst = cp.instance_id_1
        jinst = cp.instance_id_2

        array_i = collect(1:MODEL.INSTANCE[iinst].nElement)
        array_j = collect(1:MODEL.INSTANCE[jinst].nElement)

        # build surface triangles
        tri_i, tri_eleid_i, c_nodes_i =
            get_surface_triangle(MODEL.INSTANCE[iinst], array_i, cp.elements_1)
        tri_j, tri_eleid_j, c_nodes_j =
            get_surface_triangle(MODEL.INSTANCE[jinst], array_j, cp.elements_2)

        cp.c_triangles_1       = tri_i
        cp.c_triangles_eleid_1 = tri_eleid_i
        cp.c_nodes_1           = c_nodes_i

        cp.c_triangles_2       = tri_j
        cp.c_triangles_eleid_2 = tri_eleid_j
        cp.c_nodes_2           = c_nodes_j
    end

    # Build instance_pair + cp_index
    instance_pair = Vector{Vector{Int}}()
    cp_index      = Int[]
    for (idx, cp) in enumerate(MODEL.CP)
        if cp.instance_id_1 == cp.instance_id_2
            push!(instance_pair, [cp.instance_id_1, cp.instance_id_2])
            push!(cp_index, idx)
        else
            push!(instance_pair, [cp.instance_id_1, cp.instance_id_2])
            push!(instance_pair, [cp.instance_id_2, cp.instance_id_1])
            push!(cp_index, idx)
            push!(cp_index, idx)
        end
    end

    # Build the array of ContactTriangle objects (CT)
    CT = ContactTriangle[]
    for c in 1:length(cp_index)
        cc = cp_index[c]
        newCT = ContactTriangle([], [], [], zeros(Int,0,0), 0.0)
        push!(CT, newCT)

        iinst = instance_pair[c][1]
        jinst = instance_pair[c][2]

        # figure out whether c_nodes_1/c_nodes_2 correspond to iinst or jinst
        if MODEL.CP[cc].instance_id_1 == iinst
            newCT.c_nodes_i = MODEL.CP[cc].c_nodes_1 .+ MODEL.INSTANCE[iinst].node_offset
            newCT.c_nodes_j = MODEL.CP[cc].c_nodes_2 .+ MODEL.INSTANCE[jinst].node_offset
            newCT.c_triangles = MODEL.CP[cc].c_triangles_2 .+ MODEL.INSTANCE[jinst].node_offset
            newCT.c_triangles_eleid = MODEL.CP[cc].c_triangles_eleid_2 .+ MODEL.INSTANCE[jinst].element_offset
        else
            newCT.c_nodes_i = MODEL.CP[cc].c_nodes_2 .+ MODEL.INSTANCE[iinst].node_offset
            newCT.c_nodes_j = MODEL.CP[cc].c_nodes_1 .+ MODEL.INSTANCE[jinst].node_offset
            newCT.c_triangles = MODEL.CP[cc].c_triangles_1 .+ MODEL.INSTANCE[jinst].node_offset
            newCT.c_triangles_eleid = MODEL.CP[cc].c_triangles_eleid_1 .+ MODEL.INSTANCE[jinst].element_offset
        end

        # Use the jinst’s material for 'young'
        mat_id = MODEL.INSTANCE[jinst].material_id
        newCT.young = MODEL.MATERIAL[mat_id].young
    end

    return (
        contactFlag       = contact_flag,
        all_exterior_flag = all_exterior_flag,
        instance_pair     = instance_pair,
        cp_index          = cp_index,
        CT                = CT
    )
end