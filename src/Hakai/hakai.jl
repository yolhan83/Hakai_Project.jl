#########################################################
#                     Data Structures                   #
#########################################################

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


#########################################################
#                 Refactored Helper Blocks              #
#########################################################

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


@doc raw"""
    runMainTimeLoop!(MODEL, diag_M, diag_C, disp, disp_new, disp_pre, velo, position,
                     contactData, elementVolume, Pusai_mat)

Implements the central time-stepping loop:
  - updates external forces
  - computes contact forces
  - enforces boundary conditions
  - updates stress/strain
  - fracture removal
  - writes output periodically
"""
function runMainTimeLoop!(MODEL,
                          diag_M, diag_C,
                          disp, disp_new, disp_pre, velo, position,
                          contactData,
                          elementVolume,
                          Pusai_mat)

    # Basic data
    fn         = MODEL.nNode * 3
    nElement   = MODEL.nElement
    contactFlag= contactData === nothing ? 0 : contactData.contactFlag
    mass_scaling = MODEL.mass_scaling
    d_time     = MODEL.d_time * sqrt(mass_scaling)
    end_time   = MODEL.end_time
    time_num   = end_time / d_time

    # For convenience
    elementmat      = MODEL.elementmat
    elementMaterial = MODEL.element_material
    elementInstance = MODEL.element_instance
    integ_num       = 8  # fixed
    div_d_time2     = 1.0 / d_time^2

    # output settings
    output_num = 100
    d_out      = Int( floor(time_num / output_num) )

    # create containers
    Q  = zeros(fn)
    Qe = zeros(24, nElement)
    external_force = zeros(fn)

    # contact
    c_force3 = nothing
    if contactFlag >= 1
        # we store it as a 2D array for multi-thread
        c_force3 = zeros(Float128, fn, Threads.nthreads())
    end

    # store integrated data
    integ_stress = zeros(6, nElement*integ_num)  # column major
    integ_strain = zeros(6, nElement*integ_num)
    integ_plastic_strain = zeros(6, nElement*integ_num)
    integ_eq_plastic_strain = zeros(nElement * integ_num)
    integ_triax_stress     = zeros(nElement * integ_num)
    element_flag = ones(Int, nElement)

    # yield stress per integration point
    integ_yield_stress = zeros(nElement * integ_num)
    for i in 1:nElement
        mat_id = elementMaterial[i]
        plastic_table = MODEL.MATERIAL[mat_id].plastic
        if size(plastic_table,1) > 0
            # initial yield is plastic_table[1,1]
            idxrange = (i-1)*integ_num .+ (1:integ_num)
            integ_yield_stress[idxrange] .= plastic_table[1,1]
        end
    end

    # Prepare IntegDataType for output
    integ_data = IntegDataType(
        0,
        integ_stress,
        integ_strain,
        integ_plastic_strain,
        integ_eq_plastic_strain,
        integ_triax_stress
    )

    # first output @ t=0
    node_data = cal_node_stress_strain(MODEL.nNode, elementmat, integ_num, integ_data)
    write_vtk(0, MODEL.coordmat, elementmat, element_flag, disp, velo, node_data)

    output_disp          = Matrix{Float64}(undef, fn, output_num)
    output_element_flag  = Matrix{Int}(undef, nElement, output_num)
    output_data          = Vector{IntegDataType}(undef, output_num)
    i_out = 1

    println("Start time loop with dt=", d_time, ", total steps=", time_num)

    @time for t in 1:time_num

        # Just a console display every 100 steps
        if t % 100 == 0
            @printf("\rTime=%.4e / %.4e  ", t*d_time, end_time)
        end

        # 1) Build external forces from e.g. user loads
        external_force .= 0.0
        # In original code, you do "external_force[force_dof] = force_v" if you had BC with force

        # 2) If contact, compute contact forces
        if contactFlag >= 1
            # zero c_force3
            c_force3 .= 0.0
            # call CPU or GPU version
            # here we do the CPU version for simplicity
            cal_contact_force(c_force3,
                              contactData.CT,
                              contactData.instance_pair,
                              contactData.cp_index,
                              position,
                              velo,
                              diag_M,
                              minimum(elementVolume), maximum(elementVolume),
                              0.0,  # d_max
                              element_flag,
                              elementmat,
                              stdout,  # or bug_report
                              t*d_time)

            # reduce across threads if multi-threaded
            if Threads.nthreads() > 1
                @inbounds @simd for i in 1:fn
                    @inbounds @simd for j in 2:Threads.nthreads()
                        c_force3[i,1] += c_force3[i,j]
                    end
                end
            end

            # add to external force
            @inbounds for i in 1:fn
                external_force[i] += c_force3[i]
            end
        end

        # 3) Solve for disp_new
        @inbounds for i in 1:fn
            denom = (diag_M[i]/d_time^2 + diag_C[i]/(2.0*d_time))
            val = external_force[i] -
                  Q[i] +
                  diag_M[i]/d_time^2 * (2.0*disp[i] - disp_pre[i]) +
                  diag_C[i]/(2.0*d_time) * disp_pre[i]
            disp_new[i] = val / denom
        end

        # 4) Enforce boundary conditions (displacement BC)
        applyBoundaryConditions!(MODEL, disp_new, t, d_time)

        # 5) Update disp, disp_pre, velo
        @inbounds for i in 1:fn
            dtemp = disp_new[i] - disp[i]
            disp_pre[i] = disp[i]
            disp[i]     = disp_new[i]
            velo[i]     = dtemp / d_time
        end

        # 6) Update "position" = coord + disp
        @inbounds for nd in 1:MODEL.nNode
            position[1,nd] = MODEL.coordmat[1,nd] + disp[nd*3-2]
            position[2,nd] = MODEL.coordmat[2,nd] + disp[nd*3-1]
            position[3,nd] = MODEL.coordmat[3,nd] + disp[nd*3]
        end

        # 7) Recompute internal forces Q => from cal_stress_hexa
        Qe .= 0.0
        cal_stress_hexa(
            Qe,
            integ_stress,
            integ_strain,
            integ_yield_stress,
            integ_eq_plastic_strain,
            position,
            disp .- disp_pre,  # d_disp
            elementmat,
            element_flag,
            8,
            Pusai_mat,
            MODEL.MATERIAL,
            elementMaterial,
            minimum(elementVolume),
            elementVolume
        )

        Q .= 0.0
        @inbounds for e in 1:nElement
            @inbounds for i in 1:8
                idx = (elementmat[i,e]-1)*3
                Q[idx+1] += Qe[(i-1)*3+1, e]
                Q[idx+2] += Qe[(i-1)*3+2, e]
                Q[idx+3] += Qe[(i-1)*3+3, e]
            end
        end

        # 8) Compute triaxial for each integration point
        cal_triax_stress(integ_stress, integ_triax_stress)

        # 9) Possibly remove fractured elements
        doFractureCheck!(MODEL, element_flag,
                         integ_eq_plastic_strain,
                         integ_triax_stress,
                         integ_stress,
                         integ_strain)

        # 10) (Optional) update contact surfaces for newly deleted elements
        if contactFlag > 0
            # skip or replicate your logic for add_surface_triangle, etc.
        end

        # 11) Output results periodically
        if (t % d_out) == 0
            if i_out <= output_num
                output_disp[:, i_out]         = disp
                output_element_flag[:, i_out] = element_flag
                # store IntegDataType
                stored_data = IntegDataType(i_out,
                                            integ_stress,
                                            integ_strain,
                                            integ_plastic_strain,
                                            integ_eq_plastic_strain,
                                            integ_triax_stress)
                output_data[i_out] = stored_data

                # build node_data + write VTK
                node_value = cal_node_stress_strain(MODEL.nNode, elementmat, 8, stored_data)
                write_vtk(i_out,
                          MODEL.coordmat,
                          elementmat,
                          element_flag,
                          disp,
                          velo,
                          node_value)

                i_out += 1
            end
        end
    end

    println("\nFinished time loop.")
end


@doc raw"""
    doFractureCheck!(MODEL, element_flag, eq_plastic, triax, integ_stress, integ_strain)

Checks which elements should be "deleted" (element_flag=0) based on
the ductile data in MATERIAL. Sets integr_stress/strain to zero in any
deleted element's integration points. The logic matches your original code.
"""
function doFractureCheck!(MODEL, element_flag,
                          integ_eq_plastic_strain,
                          integ_triax_stress,
                          integ_stress,
                          integ_strain)
    # only do anything if there's fracture
    # we set a local flag checking if any material has fracture
    fracture_flag = 0
    for mat in MODEL.MATERIAL
        if size(mat.ductile,1) > 0 || length(mat.failure_stress) > 0
            fracture_flag = 1
            break
        end
    end
    if fracture_flag == 0
        return
    end

    nElement  = MODEL.nElement
    integ_num = 8
    elemMat   = MODEL.element_material

    # check element by element
    for i in 1:nElement
        if element_flag[i] == 0
            continue
        end

        mat_id = elemMat[i]
        # ductile
        nd = size(MODEL.MATERIAL[mat_id].ductile,1)
        if nd > 0
            # average eq_plastic & triax over this element
            # for all 8 gauss points
            pSum = 0.0
            tSum = 0.0
            baseIdx = (i-1)*integ_num
            @inbounds for j in 1:integ_num
                pSum += integ_eq_plastic_strain[baseIdx + j]
                tSum += integ_triax_stress[baseIdx + j]
            end
            pAvg = pSum / integ_num
            tAvg = tSum / integ_num

            if tAvg < 0.0
                continue
            end

            # find the ductile fraction
            ductile_ = MODEL.MATERIAL[mat_id].ductile
            fr_e = ductile_[nd,1]  # last line
            # linear interpolation in tAvg
            for j in 1:(nd-1)
                if tAvg >= ductile_[j,2] && tAvg < ductile_[j+1,2]
                    # interpolate
                    a1 = ductile_[j,1]
                    a2 = ductile_[j+1,1]
                    b1 = ductile_[j,2]
                    b2 = ductile_[j+1,2]
                    fr_e = a1 + (a2 - a1)*(tAvg - b1)/(b2 - b1)
                    break
                end
            end

            if pAvg >= fr_e
                # element fails
                element_flag[i] = 0
                @inbounds for j in 1:integ_num
                    idx = baseIdx + j
                    @inbounds for k in 1:6
                        integ_stress[k, idx] = 0.0
                        integ_strain[k, idx] = 0.0
                    end
                end
                println("Element $i is deleted due to ductile fracture.")
            end
        end
        # also can check for failure_stress if needed
    end
end


@doc raw"""
    applyBoundaryConditions!(MODEL, disp_new, step, d_time)

Implements the boundary conditions each time step, using amplitude interpolation.
"""
function applyBoundaryConditions!(MODEL, disp_new, step, d_time)
    current_time = step * d_time
    for bc in MODEL.BC
        amp_value = 1.0
        if !isempty(bc.amp_name)
            # amplitude
            ampData = bc.amplitude
            timeVec = ampData.time
            valVec  = ampData.value

            # find interval
            idx = 1
            for i in 1:(length(timeVec)-1)
                if current_time >= timeVec[i] && current_time <= timeVec[i+1]
                    idx = i
                    break
                end
            end
            t1, t2 = timeVec[idx], timeVec[idx+1]
            v1, v2 = valVec[idx], valVec[idx+1]
            frac = (current_time - t1) / (t2 - t1 + 1e-15)
            amp_value = v1 + (v2 - v1)*frac
        end

        # now apply
        for j in 1:length(bc.dof)
            dof_list = bc.dof[j]
            val      = bc.value[j]
            disp_new[dof_list] .= val * amp_value
        end
    end
end

#########################################################
#                  The Main "runHakai"                  #
#########################################################

function runHakai(fname::String)
    # 1) Read model
    println("Reading file: ", fname)
    MODEL = readInpFile(fname)

    # 2) Prepare materials
    println("Preparing material properties...")
    fracture_flag = prepareMaterialProperties!(MODEL)

    # 3) Integration points (8 for hexa)
    integ_num = 8
    Pusai_mat = cal_Pusai_hexa(integ_num)

    # 4) Compute element volumes
    println("Computing element volumes...")
    elementVolume = computeElementVolumes(MODEL, Pusai_mat, integ_num=integ_num)

    # 5) Build diag_M, diag_C
    println("Building mass/damping matrices...")
    diag_M, diag_C = buildMassDamping(MODEL, elementVolume,
                                      mass_scaling=MODEL.mass_scaling)

    # 6) Initialize displacement & velocity
    println("Initializing displacement & velocity...")
    nNode = MODEL.nNode
    d_time = MODEL.d_time * sqrt(MODEL.mass_scaling)
    disp, disp_new, disp_pre, velo, position =
        initDisplacementVelocity(MODEL, nNode, d_time, diag_M, MODEL.mass_scaling)

    # 7) Prepare contact if contact_flag >= 1
    println("Setting up contact (if any)...")
    contactData = setupContact!(MODEL)

    # 8) Main time-stepping loop
    println("Entering main simulation loop...")
    runMainTimeLoop!(MODEL,
                     diag_M,
                     diag_C,
                     disp,
                     disp_new,
                     disp_pre,
                     velo,
                     position,
                     contactData,
                     elementVolume,
                     Pusai_mat)

    println("Simulation finished.")
    return nothing
end

#########################################################
#                 Utility / Math Functions             #
#########################################################

function my3norm(b1::Float64, b2::Float64, b3::Float64)
    sqrt(b1*b1 + b2*b2 + b3*b3)
end

function my3det(m3)
    v = m3[1,1]*m3[2,2]*m3[3,3] +
        m3[1,2]*m3[2,3]*m3[3,1] +
        m3[1,3]*m3[2,1]*m3[3,2] -
        m3[1,1]*m3[2,3]*m3[3,2] -
        m3[1,2]*m3[2,1]*m3[3,3] -
        m3[1,3]*m3[2,2]*m3[3,1]
    return v
end

function cal_triax_stress(integ_stress::Array{Float64,2},
                          integ_triax_stress::Array{Float64,1})
    n = size(integ_stress,2)
    integ_triax_stress .= 0.0
    @inbounds for i in 1:n
        ox  = integ_stress[1,i]
        oy  = integ_stress[2,i]
        oz  = integ_stress[3,i]
        txy = integ_stress[4,i]
        tyz = integ_stress[5,i]
        txz = integ_stress[6,i]
        T = @SMatrix [ox  txy txz;
                      txy oy  tyz;
                      txz tyz oz ]
        p = eigvals(T)
        # eq
        oeq = sqrt(0.5 * ( (p[1]-p[2])^2 + (p[2]-p[3])^2 + (p[3]-p[1])^2 ))
        if oeq < 1e-10
            continue
        end
        v = (p[1]+p[2]+p[3]) / 3.0 / oeq
        integ_triax_stress[i] = v
    end
end

@doc raw"""
    cal_Pusai_hexa(integ_num)

Generates the shape-function gradient matrices for an 8-integration-point
hexahedral element. Returns an array of size `integ_num` where each entry
is a 3×8 matrix (Pusai1), storing the derivatives of shape functions in
local coordinates (gzai, eta, tueta).

`integ_num` is typically 8 for a standard hexahedron with 2×2×2 integration.
"""
function cal_Pusai_hexa(integ_num)
    # We'll store each 3×8 matrix in an array of length integ_num
    # (like an array-of-matrices)
    Pusai_mat = [ @MMatrix zeros(3,8) for i in 1:integ_num ]

    # Delta matrix used to define shape function sign patterns
    delta_mat = [
        -1.0 -1.0 -1.0
         1.0 -1.0 -1.0
         1.0  1.0 -1.0
        -1.0  1.0 -1.0
        -1.0 -1.0  1.0
         1.0 -1.0  1.0
         1.0  1.0  1.0
        -1.0  1.0  1.0
    ]

    # Integration point coordinates in local space
    gc = [0.0 0.0 0.0]
    if integ_num == 8
        g = 1.0 / sqrt(3.0)
        gc = [
            -g -g -g
            -g -g  g
            -g  g -g
            -g  g  g
             g -g -g
             g -g  g
             g  g -g
             g  g  g
        ]
    end

    # Build each Pusai1 matrix
    for k in 1:integ_num
        Pusai1 = @MMatrix zeros(3,8)

        gzai   = gc[k, 1]
        eta    = gc[k, 2]
        tueta  = gc[k, 3]

        for i in 1:8
            # shape function derivative in local coords:
            # each row is dN/d(gzai), dN/d(eta), dN/d(tueta) for node i
            Pusai1[1,i] = 1.0/8.0 * delta_mat[i,1] * (1.0 + eta    * delta_mat[i,2]) * (1.0 + tueta  * delta_mat[i,3])
            Pusai1[2,i] = 1.0/8.0 * delta_mat[i,2] * (1.0 + gzai   * delta_mat[i,1]) * (1.0 + tueta  * delta_mat[i,3])
            Pusai1[3,i] = 1.0/8.0 * delta_mat[i,3] * (1.0 + gzai   * delta_mat[i,1]) * (1.0 + eta    * delta_mat[i,2])
        end

        # Store in the array of 3×8 matrices
        Pusai_mat[k] .= Pusai1
    end

    return Pusai_mat
end

@doc raw"""
    get_element_face(MODEL, i)

For the i-th instance in `MODEL`, returns three arrays describing the
faces of each hexahedral element:

- `faces`       : A (nE*6)×4 matrix of node indices, where nE is the number of elements.
                  Each group of 6 rows corresponds to the 6 faces of one element.
- `faces_eleid` : A vector of length (nE*6) indicating which element each face row belongs to.
- `sorted_faces`: Same as `faces` but with each row sorted so that face matching can be performed
                  (useful for detecting shared faces, external faces, etc.).
"""
function get_element_face(MODEL, i)
    part_id = MODEL.INSTANCE[i].part_id
    cdmat   = MODEL.PART[part_id].coordmat
    nE      = MODEL.INSTANCE[i].nElement
    elemMat = MODEL.PART[part_id].elementmat

    # Prepare arrays to store face node indices
    # Each hexahedron has 6 faces, each face has 4 nodes.
    faces        = zeros(Int, nE*6, 4)
    faces_eleid  = zeros(Int, nE*6)
    sorted_faces = zeros(Int, nE*6, 4)

    c = 1
    for j in 1:nE
        # Node list for j-th element in this part
        elem = elemMat[:, j]

        # Fill in the 6 faces for this element
        @views faces[6*(c-1) + 1, :] = elem[1:4]
        @views faces[6*(c-1) + 2, :] = elem[5:8]
        @views faces[6*(c-1) + 3, :] = [elem[1], elem[2], elem[6], elem[5]]
        @views faces[6*(c-1) + 4, :] = [elem[2], elem[3], elem[7], elem[6]]
        @views faces[6*(c-1) + 5, :] = [elem[3], elem[4], elem[8], elem[7]]
        @views faces[6*(c-1) + 6, :] = [elem[4], elem[1], elem[5], elem[8]]

        # Record which element each face belongs to
        @views faces_eleid[6*(c-1) .+ (1:6)] .= j

        # Check the orientation of each face so the outward normal is consistent
        # (For example, by comparing cross product with centroid vector.)
        # We compute the element centroid:
        centroid = sum!(zeros(3,1), @view(cdmat[:, elem])) / 8.0

        for k in 1:6
            idx  = 6*(c-1) + k
            n1   = faces[idx, 1]
            n2   = faces[idx, 2]
            n3   = faces[idx, 4]  # diagonal order

            v1 = cdmat[:, n2] - cdmat[:, n1]
            v2 = cdmat[:, n3] - cdmat[:, n1]
            nv = v1 × v2  # cross product

            # Vector from face's first node to centroid
            vc = centroid .- cdmat[:, n1]

            # If dot(nv, vc) > 0, flip the face node order
            if dot(nv, vc) > 0.0
                # reorder: (1, 4, 3, 2)
                @views faces[idx, :] = [n1, faces[idx,4], faces[idx,3], faces[idx,2]]
            end
        end
        c += 1
    end

    # Create the sorted_faces array for easy face matching
    @inbounds for j in 1:(nE*6)
        @views sorted_faces[j, :] = sort(faces[j, :])
    end

    return faces, faces_eleid, sorted_faces
end

@doc raw"""
    get_surface_triangle(INSTANCE_i, array_element, contact_element)

Given a single `InstanceType` (`INSTANCE_i`), an array of local element indices
(`array_element`), and a subset of those elements (`contact_element`), returns:

1. `c_triangles`       : An array of triangle faces (two triangles per quadrilateral face).
2. `c_triangles_eleid` : A matching array of element indices for each triangle row in `c_triangles`.
3. `c_nodes`           : A sorted and unique list of all node indices used in `c_triangles`.

This is typically used in contact calculations to extract the exposed surface
triangles for a subset of elements, so that point-vs-triangle collision or gap
can be computed.
"""
function get_surface_triangle(INSTANCE_i, array_element, contact_element)
    nE = length(array_element)

    # Each element has 6 faces, each face has 4 nodes, so allocate for nE*6 faces
    surfaces        = zeros(Int, nE*6, 4)
    surfaces_eleid  = zeros(Int, nE*6)
    sorted_surfaces = zeros(Int, nE*6, 4)

    # 1) Gather the face and sorted-face info for each element in array_element
    c = 1
    for j in array_element
        # Each element has 6 faces in INSTANCE_i.surfaces
        for k in 1:6
            row_idx = 6*(c-1) + k
            # Copy face node data from the instance's stored 'surfaces'
            surfaces[row_idx, 1] = INSTANCE_i.surfaces[6*(j-1) + k, 1]
            surfaces[row_idx, 2] = INSTANCE_i.surfaces[6*(j-1) + k, 2]
            surfaces[row_idx, 3] = INSTANCE_i.surfaces[6*(j-1) + k, 3]
            surfaces[row_idx, 4] = INSTANCE_i.surfaces[6*(j-1) + k, 4]

            # Same row for sorted_surfaces
            sorted_surfaces[row_idx, 1] = INSTANCE_i.sorted_surfaces[6*(j-1) + k, 1]
            sorted_surfaces[row_idx, 2] = INSTANCE_i.sorted_surfaces[6*(j-1) + k, 2]
            sorted_surfaces[row_idx, 3] = INSTANCE_i.sorted_surfaces[6*(j-1) + k, 3]
            sorted_surfaces[row_idx, 4] = INSTANCE_i.sorted_surfaces[6*(j-1) + k, 4]

        end
        # Mark which element each face belongs to
        surfaces_eleid[6*(c-1) .+ (1:6)] .= j
        c += 1
    end

    # 2) Identify faces that are not interior (i.e., appear only once).
    #    We'll keep them in c_surfaces and c_surfaces_eleid.
    c_surfaces       = zeros(Int, nE*6, 4)
    c_surfaces_eleid = zeros(Int, nE*6)
    dp_id            = Int[]
    sj = zeros(Int,4)
    sk = zeros(Int,4)
    c  = 0

    # We compare sorted_surfaces row-by-row to see if it’s unique
    nF = nE*6
    for j in 1:(nF-1)
        u_flag = 1
        # Skip if we've already marked row j as a duplicate
        if !isnothing(findfirst(dp_id .== j))
            u_flag = 0
            continue
        end

        sj[1] = sorted_surfaces[j, 1]
        sj[2] = sorted_surfaces[j, 2]
        sj[3] = sorted_surfaces[j, 3]
        sj[4] = sorted_surfaces[j, 4]

        for k in (j+1):nF
            sk[1] = sorted_surfaces[k, 1]
            sk[2] = sorted_surfaces[k, 2]
            sk[3] = sorted_surfaces[k, 3]
            sk[4] = sorted_surfaces[k, 4]
            # If they match, it's an interior face => skip
            if sj[1]==sk[1] && sj[2]==sk[2] && sj[3]==sk[3] && sj[4]==sk[4]
                u_flag = 0
                push!(dp_id, k)
                break
            end
        end

        if u_flag == 1
            c += 1
            c_surfaces[c, :]    = surfaces[j, :]
            c_surfaces_eleid[c] = surfaces_eleid[j]
        end
    end
    # Trim
    c_surfaces       = c_surfaces[1:c, :]
    c_surfaces_eleid = c_surfaces_eleid[1:c]

    # 3) Filter out only those faces that belong to `contact_element`.
    #    If contact_element < total # of elements, we keep only the matching subset.
    if INSTANCE_i.nElement != length(contact_element)
        c_surfaces_temp       = zeros(Int, c, 4)
        c_surfaces_eleid_temp = zeros(Int, c)
        c2 = 0
        for j in 1:c
            if !isnothing(findfirst(contact_element .== c_surfaces_eleid[j]))
                c2 += 1
                c_surfaces_temp[c2, :]     = c_surfaces[j, :]
                c_surfaces_eleid_temp[c2]  = c_surfaces_eleid[j]
            end
        end
        c_surfaces       = c_surfaces_temp[1:c2, :]
        c_surfaces_eleid = c_surfaces_eleid_temp[1:c2]
        c = c2
    end

    # 4) Convert each quadrilateral face into two triangles
    if c == 0
        # no faces => empty
        return zeros(Int,0,0), zeros(Int,0), zeros(Int,0)
    end

    # We'll store two triangles per face => total 2*c
    c_triangles       = zeros(Int, 2*c, 3)
    c_triangles_eleid = zeros(Int, 2*c)
    for j in 1:c
        c_triangles[2j-1, 1] = c_surfaces[j,1]
        c_triangles[2j-1, 2] = c_surfaces[j,2]
        c_triangles[2j-1, 3] = c_surfaces[j,3]

        c_triangles[2j,   1] = c_surfaces[j,3]
        c_triangles[2j,   2] = c_surfaces[j,4]
        c_triangles[2j,   3] = c_surfaces[j,1]

        c_triangles_eleid[2j-1] = c_surfaces_eleid[j]
        c_triangles_eleid[2j]   = c_surfaces_eleid[j]
    end

    # 5) Gather all unique node IDs
    c_nodes = vec(c_triangles)
    sort!(unique!(c_nodes))

    return c_triangles, c_triangles_eleid, c_nodes
end

@doc raw"""
    cal_node_stress_strain(nNode::Int, elementmat::Array{Int,2}, integ_num::Int, integ_data::IntegDataType)

Averages integration-point (Gauss-point) stress/strain values to each node.
Returns a `NodeDataType` containing:

  • `node_stress`          (nNode×6)  – the averaged stress at each node
  • `node_strain`          (nNode×6)  – the averaged strain at each node
  • `node_plastic_strain`  (nNode×6)  – (currently unused, but included for completeness)
  • `node_eq_plastic_strain` (nNode) – the averaged equivalent plastic strain at each node
  • `node_mises_stress`    (nNode)   – Mises stress computed from `node_stress`
  • `node_triax_stress`    (nNode)   – triaxial (mean stress / Mises) ratio

Notes:
  - `integ_data.integ_stress` and `integ_data.integ_strain` are stored in column-major
    with shape (6 × nElement*integ_num). This function first aggregates them per element,
    then sums over nodes.
"""
function cal_node_stress_strain(nNode::Int,
                                elementmat::Array{Int,2},
                                integ_num::Int,
                                integ_data::IntegDataType)

    # Extract arrays from IntegDataType for convenience
    integ_stress          = integ_data.integ_stress'  # make it row-major => easier indexing
    integ_strain          = integ_data.integ_strain'
    integ_plastic_strain  = integ_data.integ_plastic_strain
    integ_eq_plastic_strain = integ_data.integ_eq_plastic_strain
    integ_triax_stress    = integ_data.integ_triax_stress

    # Prepare node-level accumulators
    node_stress         = zeros(nNode, 6)
    node_strain         = zeros(nNode, 6)
    node_plastic_strain = zeros(nNode, 6)
    node_eq_plastic_strain = zeros(nNode)
    node_triax_stress   = zeros(nNode)

    # For convenience
    nElement = size(elementmat, 2)
    element_stress          = zeros(nElement, 6)
    element_strain          = zeros(nElement, 6)
    element_eq_plastic_strain = zeros(nElement)
    element_triax_stress    = zeros(nElement)

    # 1) Compute average stress/strain for each element (average among the integ_num Gauss points).
    for e in 1:nElement
        if integ_num == 1
            # if there's only one integration point => direct copy
            element_stress[e, :]           = integ_stress[e, :]
            element_strain[e, :]           = integ_strain[e, :]
            element_eq_plastic_strain[e]   = integ_eq_plastic_strain[e]
            element_triax_stress[e]        = integ_triax_stress[e]
        else
            # sum up across the 8 Gauss points for this element, then average
            rowstart = (e-1)*integ_num + 1
            rowend   = e*integ_num

            element_stress[e, :] = sum!(zeros(1,6), integ_stress[rowstart:rowend, :]) / integ_num
            element_strain[e, :] = sum!(zeros(1,6), integ_strain[rowstart:rowend, :]) / integ_num
            element_eq_plastic_strain[e] =
                sum(integ_eq_plastic_strain[rowstart:rowend]) / integ_num
            element_triax_stress[e] =
                sum(integ_triax_stress[rowstart:rowend]) / integ_num
        end
    end

    # 2) Accumulate element-averaged stress/strain data to each node
    inc_num = zeros(nNode)
    for e in 1:nElement
        for k in 1:8
            nd = elementmat[k, e]
            node_stress[nd, :]         .+= element_stress[e, :]
            node_strain[nd, :]         .+= element_strain[e, :]
            node_eq_plastic_strain[nd] += element_eq_plastic_strain[e]
            node_triax_stress[nd]      += element_triax_stress[e]
            inc_num[nd]               += 1
        end
    end

    # 3) Divide by inc_num to get the final average at each node
    for nd in 1:nNode
        if inc_num[nd] > 0
            node_stress[nd, :]         ./= inc_num[nd]
            node_strain[nd, :]         ./= inc_num[nd]
            node_eq_plastic_strain[nd] /=  inc_num[nd]
            node_triax_stress[nd]      /=  inc_num[nd]
        end
    end

    # 4) Compute node-level Mises stress
    node_mises_stress = zeros(nNode)
    for nd in 1:nNode
        ox  = node_stress[nd,1]
        oy  = node_stress[nd,2]
        oz  = node_stress[nd,3]
        txy = node_stress[nd,4]
        tyz = node_stress[nd,5]
        txz = node_stress[nd,6]
        node_mises_stress[nd] = sqrt(0.5 * (
            (ox - oy)^2 + (oy - oz)^2 + (ox - oz)^2 +
             6*(txy^2 + tyz^2 + txz^2)
        ))
    end

    # 5) Return a NodeDataType
    node_data = NodeDataType(
        node_stress,
        node_strain,
        node_plastic_strain,
        node_eq_plastic_strain,
        node_mises_stress,
        node_triax_stress
    )
    return node_data
end

@doc raw"""
    write_vtk(index::Int, coordmat::AbstractMatrix{T}, elementmat::AbstractMatrix{Int},
              element_flag::AbstractVector{Int}, disp::AbstractVector{T}, velo::AbstractVector{T},
              node_data::NodeDataType)

Writes a .vtk file (in unstructured-grid format) for the mesh and nodal data at the given `index`.
Stores the file in `temp/file###.vtk`.

Arguments:
  • `coordmat`     : a 3×nNode matrix of the original node coordinates (column-major).
  • `elementmat`   : an 8×nElement matrix of node indices for each hexahedral element.
  • `element_flag` : an nElement vector of 1/0 flags indicating whether each element is active or deleted.
  • `disp`         : a 3*nNode vector of nodal displacements, stored in [x1,y1,z1, x2,y2,z2, …].
  • `velo`         : a 3*nNode vector of nodal velocities, same layout as `disp`.
  • `node_data`    : a `NodeDataType` struct, containing node-level stresses, strains, etc.

Output:
  • Creates a file named "temp/file###.vtk" where ### is the zero-padded index.

Notes:
  • The VTK file is ASCII, "DATASET UNSTRUCTURED_GRID".
  • Only hex8 cell types are handled (CELL_TYPES 12).
"""
function write_vtk(index::Int,
                   coordmat::AbstractMatrix{Float64},
                   elementmat::AbstractMatrix{Int},
                   element_flag::AbstractVector{Int},
                   disp::AbstractVector{Float64},
                   velo::AbstractVector{Float64},
                   node_data::NodeDataType)

    nNode    = size(coordmat, 2)    # number of nodes
    nElement = size(elementmat, 2)  # number of elements

    # Reshape displacements and velocities into (nNode × 3) for easier writing
    disp3 = reshape(disp, 3, nNode)'  # nNode×3
    velo3 = reshape(velo, 3, nNode)'  # nNode×3

    # Unpack from NodeDataType
    node_stress          = node_data.node_stress
    node_strain          = node_data.node_strain
    node_mises_stress    = node_data.node_mises_stress
    node_eq_plastic_strain = node_data.node_eq_plastic_strain
    node_triax_stress    = node_data.node_triax_stress

    # Some tiny tolerance checks: if near zero, just set to zero
    # (prevents scientific notation on extremely small values).
    @inbounds for i in 1:nNode
        @inbounds for j in 1:3
            if abs(disp3[i, j]) < 1e-16
                disp3[i, j] = 0.0
            end
            if abs(velo3[i, j]) < 1e-16
                velo3[i, j] = 0.0
            end
        end

        @inbounds for j in 1:6
            if abs(node_stress[i, j]) < 1e-16
                node_stress[i, j] = 0.0
            end
            if abs(node_strain[i, j]) < 1e-16
                node_strain[i, j] = 0.0
            end
        end

        if abs(node_mises_stress[i]) < 1e-16
            node_mises_stress[i] = 0.0
        end
        if abs(node_eq_plastic_strain[i]) < 1e-16
            node_eq_plastic_strain[i] = 0.0
        end
        if abs(node_triax_stress[i]) < 1e-16
            node_triax_stress[i] = 0.0
        end
    end

    # Make sure we have a "temp" directory
    if !ispath("temp")
        mkdir("temp")
    end

    fname = @sprintf("temp/file%03d.vtk", index)
    out   = open(fname, "w")

    # VTK Header
    println(out, "# vtk DataFile Version 2.0")
    println(out, "Test")
    println(out, "ASCII")
    println(out, "DATASET UNSTRUCTURED_GRID")

    # Write points
    println(out, "POINTS ", nNode, " float")
    for i in 1:nNode
        @printf(out, "%1.6e %1.6e %1.6e\n", coordmat[1,i], coordmat[2,i], coordmat[3,i])
    end

    # Count how many elements are active
    activeElements = findall(x->x==1, element_flag)
    draw_element_num = length(activeElements)

    # Write cells
    println(out, "CELLS ", draw_element_num, " ", draw_element_num*(8+1))
    e = elementmat .- 1  # shift to 0-based indexing for VTK
    for idx in activeElements
        @printf(out, "8 %d %d %d %d %d %d %d %d\n",
                e[1,idx], e[2,idx], e[3,idx], e[4,idx],
                e[5,idx], e[6,idx], e[7,idx], e[8,idx])
    end

    # Write cell types
    println(out, "CELL_TYPES ", draw_element_num)
    for _ in 1:draw_element_num
        println(out, 12)   # 12 => VTK_HEXAHEDRON
    end

    # Begin node data
    println(out, "POINT_DATA ", nNode)

    # 1) Displacement
    println(out, "VECTORS DISPLACEMENT float")
    for i in 1:nNode
        @printf(out, "%1.6e %1.6e %1.6e\n", disp3[i,1], disp3[i,2], disp3[i,3])
    end

    # 2) Velocity components
    println(out, "SCALARS Vx float 1")
    println(out, "LOOKUP_TABLE default")
    for i in 1:nNode
        @printf(out, "%1.6e\n", velo3[i,1])
    end

    println(out, "SCALARS Vy float 1")
    println(out, "LOOKUP_TABLE default")
    for i in 1:nNode
        @printf(out, "%1.6e\n", velo3[i,2])
    end

    println(out, "SCALARS Vz float 1")
    println(out, "LOOKUP_TABLE default")
    for i in 1:nNode
        @printf(out, "%1.6e\n", velo3[i,3])
    end

    # 3) Strain components
    println(out, "SCALARS E11 float 1")
    println(out, "LOOKUP_TABLE default")
    for i in 1:nNode
        @printf(out, "%1.6e\n", node_strain[i,1])
    end

    println(out, "SCALARS E22 float 1")
    println(out, "LOOKUP_TABLE default")
    for i in 1:nNode
        @printf(out, "%1.6e\n", node_strain[i,2])
    end

    println(out, "SCALARS E33 float 1")
    println(out, "LOOKUP_TABLE default")
    for i in 1:nNode
        @printf(out, "%1.6e\n", node_strain[i,3])
    end

    println(out, "SCALARS E12 float 1")
    println(out, "LOOKUP_TABLE default")
    for i in 1:nNode
        @printf(out, "%1.6e\n", node_strain[i,4])
    end

    println(out, "SCALARS E23 float 1")
    println(out, "LOOKUP_TABLE default")
    for i in 1:nNode
        @printf(out, "%1.6e\n", node_strain[i,5])
    end

    println(out, "SCALARS E13 float 1")
    println(out, "LOOKUP_TABLE default")
    for i in 1:nNode
        @printf(out, "%1.6e\n", node_strain[i,6])
    end

    # 4) Equivalent plastic strain
    println(out, "SCALARS EQ_PSTRAIN float 1")
    println(out, "LOOKUP_TABLE default")
    for i in 1:nNode
        @printf(out, "%1.6e\n", node_eq_plastic_strain[i])
    end

    # 5) Stress components
    println(out, "SCALARS S11 float 1")
    println(out, "LOOKUP_TABLE default")
    for i in 1:nNode
        @printf(out, "%1.6e\n", node_stress[i,1])
    end

    println(out, "SCALARS S22 float 1")
    println(out, "LOOKUP_TABLE default")
    for i in 1:nNode
        @printf(out, "%1.6e\n", node_stress[i,2])
    end

    println(out, "SCALARS S33 float 1")
    println(out, "LOOKUP_TABLE default")
    for i in 1:nNode
        @printf(out, "%1.6e\n", node_stress[i,3])
    end

    println(out, "SCALARS S12 float 1")
    println(out, "LOOKUP_TABLE default")
    for i in 1:nNode
        @printf(out, "%1.6e\n", node_stress[i,4])
    end

    println(out, "SCALARS S23 float 1")
    println(out, "LOOKUP_TABLE default")
    for i in 1:nNode
        @printf(out, "%1.6e\n", node_stress[i,5])
    end

    println(out, "SCALARS S13 float 1")
    println(out, "LOOKUP_TABLE default")
    for i in 1:nNode
        @printf(out, "%1.6e\n", node_stress[i,6])
    end

    # 6) Mises stress
    println(out, "SCALARS MISES_STRESS float 1")
    println(out, "LOOKUP_TABLE default")
    for i in 1:nNode
        @printf(out, "%1.6e\n", node_mises_stress[i])
    end

    # 7) Triaxial ratio
    println(out, "SCALARS TRIAX_STRESS float 1")
    println(out, "LOOKUP_TABLE default")
    for i in 1:nNode
        @printf(out, "%1.6e\n", node_triax_stress[i])
    end

    close(out)
end

@doc raw"""
    cal_contact_force(c_force3, CT, instance_pair, cp_index,
                      position::Array{Float64,2}, velo::Array{Float64,1}, diag_M::Array{Float64,1},
                      elementMinSize::Float64, elementMaxSize::Float64, d_max::Float64,
                      element_flag::Array{Int,1}, elementmat::Array{Int,2}, bug_report, time_)

Computes contact forces between a set of points (c_nodes_i) and triangular faces (c_triangles) from the
`ContactTriangle` array `CT`. The core logic:
1. For each contact pair, gather the relevant nodes/triangles.
2. Build a small coarse grid (node_map_ix, node_map_iy, node_map_iz, etc.) to skip far-apart
   node/triangle checks.
3. For any node that lies within the normal-projection distance `d_lim` from a triangle, compute
   normal penalty force and friction. Accumulate into `c_force3`.
4. We do both damping and friction calculations, so the contact model is basically a penalty approach.

The array `c_force3` is 2D: size (fn × Threads.nthreads()) so each thread can accumulate forces
without race conditions. The final combination is done after the function, by summation across
the second dimension if needed.
"""
function cal_contact_force(c_force3, CT, instance_pair, cp_index,
                           position::Array{Float64,2}, velo::Array{Float64,1},
                           diag_M::Array{Float64,1}, elementMinSize::Float64,
                           elementMaxSize::Float64, d_max::Float64,
                           element_flag::Array{Int,1}, elementmat::Array{Int,2},
                           bug_report, time_)

    # Basic constants
    d_lim = elementMinSize * 0.3       # max penetration distance
    myu   = 0.25                       # friction coefficient
    kc_o  = 1.0                        # normal stiffness for different-instance contact
    kc_s  = 1.0                        # normal stiffness for self-contact
    Cr_o  = 0.0                        # damping for different-instance contact
    Cr_s  = 0.0                        # damping for self-contact
    div3  = 1.0/3.0

    fn = size(position,2) * 3  # total DOFs = nNode*3

    # Loop over each contact pair definition
    for c in 1:length(cp_index)
        cc = cp_index[c]

        # i_instance is the "point" side, j_instance is the "triangle" side
        i_instance = instance_pair[c][1]
        j_instance = instance_pair[c][2]

        # Extract data from CT
        c_nodes_i  = CT[c].c_nodes_i
        c_nodes_j  = CT[c].c_nodes_j
        c_triangles= CT[c].c_triangles
        c_triangles_eleid = CT[c].c_triangles_eleid
        young      = CT[c].young

        nn_i = length(c_nodes_i)
        nn_j = length(c_nodes_j)

        # Build bounding boxes for i-nodes and j-triangles
        min_ix = minimum(@view(position[1,c_nodes_i]))
        min_iy = minimum(@view(position[2,c_nodes_i]))
        min_iz = minimum(@view(position[3,c_nodes_i]))
        min_jx = minimum(@view(position[1,c_nodes_j]))
        min_jy = minimum(@view(position[2,c_nodes_j]))
        min_jz = minimum(@view(position[3,c_nodes_j]))

        max_ix = maximum(@view(position[1,c_nodes_i]))
        max_iy = maximum(@view(position[2,c_nodes_i]))
        max_iz = maximum(@view(position[3,c_nodes_i]))
        max_jx = maximum(@view(position[1,c_nodes_j]))
        max_jy = maximum(@view(position[2,c_nodes_j]))
        max_jz = maximum(@view(position[3,c_nodes_j]))

        range_min_x = max(min_ix, min_jx)
        range_min_y = max(min_iy, min_jy)
        range_min_z = max(min_iz, min_jz)
        range_max_x = min(max_ix, max_jx)
        range_max_y = min(max_iy, max_jy)
        range_max_z = min(max_iz, max_jz)

        # If bounding boxes don't overlap, skip
        if range_min_x > range_max_x || range_min_y > range_max_y || range_min_z > range_max_z
            continue
        end

        # Build a coarse grid to skip pairs that are too far apart
        all_range_min_x = min(min_ix, min_jx)
        all_range_min_y = min(min_iy, min_jy)
        all_range_min_z = min(min_iz, min_jz)
        all_range_max_x = max(max_ix, max_jx)
        all_range_max_y = max(max_iy, max_jy)
        all_range_max_z = max(max_iz, max_jz)

        ddiv = elementMaxSize * 1.1
        # For self-contact, use a smaller grid spacing
        if i_instance == j_instance
            ddiv = elementMaxSize * 0.6
        end

        # Map each node to a small integer cell
        node_map_ix = zeros(Int, nn_i)
        node_map_iy = zeros(Int, nn_i)
        node_map_iz = zeros(Int, nn_i)
        node_map_jx = zeros(Int, nn_j)
        node_map_jy = zeros(Int, nn_j)
        node_map_jz = zeros(Int, nn_j)

        @inbounds for i in 1:nn_i
            px  = position[1, c_nodes_i[i]]
            py  = position[2, c_nodes_i[i]]
            pz  = position[3, c_nodes_i[i]]
            node_map_ix[i] = ceil((px - all_range_min_x) / ddiv)
            node_map_iy[i] = ceil((py - all_range_min_y) / ddiv)
            node_map_iz[i] = ceil((pz - all_range_min_z) / ddiv)
        end

        @inbounds for i in 1:nn_j
            px  = position[1, c_nodes_j[i]]
            py  = position[2, c_nodes_j[i]]
            pz  = position[3, c_nodes_j[i]]
            node_map_jx[i] = ceil((px - all_range_min_x) / ddiv)
            node_map_jy[i] = ceil((py - all_range_min_y) / ddiv)
            node_map_jz[i] = ceil((pz - all_range_min_z) / ddiv)
        end

        # For each triangle, check bounding box, then check i-nodes
        @inbounds for j = 1:size(c_triangles,1)
            eleid_ = c_triangles_eleid[j]
            if element_flag[eleid_] == 0
                continue
            end

            # Normal stiffness & damping
            kc = (i_instance == j_instance) ? kc_s : kc_o
            Cr = (i_instance == j_instance) ? Cr_s : Cr_o

            # Identify which thread is running
            index_th = Threads.threadid()

            # Triangle node indices
            j0 = c_triangles[j,1]
            j1 = c_triangles[j,2]
            j2 = c_triangles[j,3]

            # Tri node coords
            q0x, q0y, q0z = position[1,j0], position[2,j0], position[3,j0]
            q1x, q1y, q1z = position[1,j1], position[2,j1], position[3,j1]
            q2x, q2y, q2z = position[1,j2], position[2,j2], position[3,j2]

            # Quick bounding checks
            if (q0x < range_min_x && q1x < range_min_x && q2x < range_min_x) ||
               (q0x > range_max_x && q1x > range_max_x && q2x > range_max_x) ||
               (q0y < range_min_y && q1y < range_min_y && q2y < range_min_y) ||
               (q0y > range_max_y && q1y > range_max_y && q2y > range_max_y) ||
               (q0z < range_min_z && q1z < range_min_z && q2z < range_min_z) ||
               (q0z > range_max_z && q1z > range_max_z && q2z > range_max_z)
                continue
            end

            # Centroid
            cx = (q0x + q1x + q2x)/3.0
            cy = (q0y + q1y + q2y)/3.0
            cz = (q0z + q1z + q2z)/3.0

            # max radius
            R0 = my3norm(q0x-cx, q0y-cy, q0z-cz)
            R1 = my3norm(q1x-cx, q1y-cy, q1z-cz)
            R2 = my3norm(q2x-cx, q2y-cy, q2z-cz)
            Rmax = max(max(R0,R1), R2)

            # Vectors forming the triangle plane
            v1x, v1y, v1z = (q1x-q0x), (q1y-q0y), (q1z-q0z)
            v2x, v2y, v2z = (q2x-q0x), (q2y-q0y), (q2z-q0z)

            L1 = my3norm(v1x,v1y,v1z)
            L2 = my3norm(v2x,v2y,v2z)
            Lmax = max(L1,L2)

            # Normal
            nx, ny, nz = my3crossNNz(v1x,v1y,v1z, v2x,v2y,v2z)

            # Triangle area
            d12 = v1x*v2x + v1y*v2y + v1z*v2z
            S = 0.5 * sqrt(L1^2 * L2^2 - d12^2)

            # We'll form the matrix A for point-plane solves
            A11, A21, A31 = v1x, v1y, v1z
            A12, A22, A32 = v2x, v2y, v2z
            A13, A23, A33 = -nx, -ny, -nz

            # Map for j0 in c_nodes_j => needed for further checks
            map_j0_x = 1
            map_j0_y = 1
            map_j0_z = 1
            @inbounds for w in 1:nn_j
                if c_nodes_j[w] == j0
                    map_j0_x = node_map_jx[w]
                    map_j0_y = node_map_jy[w]
                    map_j0_z = node_map_jz[w]
                    break
                end
            end

            # For i-nodes, see if they're close enough
            @inbounds for k2 in 1:nn_i
                # Quick skip with coarse cell
                if abs(map_j0_x - node_map_ix[k2]) > 1 ||
                   abs(map_j0_y - node_map_iy[k2]) > 1 ||
                   abs(map_j0_z - node_map_iz[k2]) > 1
                    continue
                end

                iNode = c_nodes_i[k2]
                # If self-contact, skip if iNode belongs to the same element => optional logic

                px = position[1,iNode]
                py = position[2,iNode]
                pz = position[3,iNode]

                # bounding check
                if px < range_min_x || py < range_min_y || pz < range_min_z
                    continue
                end
                if px > range_max_x || py > range_max_y || pz > range_max_z
                    continue
                end

                # centroid-based skip
                dpc = my3norm(px-cx, py-cy, pz-cz)
                if dpc >= Rmax
                    continue
                end

                # Solve for (x1, x2, d) in the plane eq
                bx = px - q0x
                by = py - q0y
                bz = pz - q0z

                x1, x2, dist_ = my3SolveAb(A11,A21,A31,
                                          A12,A22,A32,
                                          A13,A23,A33,
                                          bx, by, bz)

                # Check barycentric (x1 + x2 <=1), positive
                if 0.0 <= x1 && 0.0 <= x2 && x1 + x2 <= 1.0 && dist_ > 0.0 && dist_ <= d_lim
                    # Normal penalty contact
                    # Compute velocity difference normal
                    vx = velo[iNode*3-2] - velo[j0*3-2]
                    vy = velo[iNode*3-1] - velo[j0*3-1]
                    vz = velo[iNode*3]   - velo[j0*3]
                    mag_v = my3norm(vx, vy, vz)

                    vex, vey, vez = 0.0, 0.0, 0.0
                    if mag_v > 1e-12
                        vex = vx/mag_v
                        vey = vy/mag_v
                        vez = vz/mag_v
                    end

                    # Normal stiffness: F = k * dist_
                    # k ~ young*S / Lmax, etc.
                    k_ = (young * S / Lmax) * kc
                    F  = k_ * dist_

                    fx = F*nx
                    fy = F*ny
                    fz = F*nz

                    # Damping
                    C_ = 2.0*sqrt(diag_M[iNode]*k_) * Cr
                    fc_x = -C_*vx
                    fc_y = -C_*vy
                    fc_z = -C_*vz

                    # Friction
                    dot_ve_n = vex*nx + vey*ny + vez*nz
                    vsx = vex - dot_ve_n*nx
                    vsy = vey - dot_ve_n*ny
                    vsz = vez - dot_ve_n*nz
                    fric_x = -myu*F * vsx
                    fric_y = -myu*F * vsy
                    fric_z = -myu*F * vsz

                    fx += fric_x + fc_x
                    fy += fric_y + fc_y
                    fz += fric_z + fc_z

                    # Accumulate in c_force3
                    # iNode gets +f
                    c_force3[(iNode-1)*3+1, index_th] += fx
                    c_force3[(iNode-1)*3+2, index_th] += fy
                    c_force3[(iNode-1)*3+3, index_th] += fz

                    # Triangle nodes share the negative
                    # We approximate distributing among j0, j1, j2
                    # but code shows distributing equally among j0-> /3
                    c_force3[(j0-1)*3+1, index_th] -= fx/3.0
                    c_force3[(j0-1)*3+2, index_th] -= fy/3.0
                    c_force3[(j0-1)*3+3, index_th] -= fz/3.0

                    c_force3[(j1-1)*3+1, index_th] -= fx/3.0
                    c_force3[(j1-1)*3+2, index_th] -= fy/3.0
                    c_force3[(j1-1)*3+3, index_th] -= fz/3.0

                    c_force3[(j2-1)*3+1, index_th] -= fx/3.0
                    c_force3[(j2-1)*3+2, index_th] -= fy/3.0
                    c_force3[(j2-1)*3+3, index_th] -= fz/3.0
                end
            end # end for k2
        end # end for j
    end # end for c

    return
end

@doc raw"""
    my3crossNNz(a1::Float64, a2::Float64, a3::Float64, b1::Float64, b2::Float64, b3::Float64)

Returns `(nx, ny, nz)` which is the normalized cross product of vectors `a = (a1, a2, a3)`
and `b = (b1, b2, b3)`. In other words, it does:

```math
n = (a × b) / ‖a × b‖
where ‖⋅‖ denotes the Euclidean norm. The result (nx, ny, nz) is the unit normal.

Example usage:
nx, ny, nz = my3crossNNz(1.0, 0.0, 0.0, 0.0, 1.0, 0.0)  # => (0, 0, 1)
"""
function my3crossNNz(a1::Float64, a2::Float64, a3::Float64, b1::Float64, b2::Float64, b3::Float64)
    # Compute the raw cross product
    n1 = a2*b3 - a3*b2
    n2 = a3*b1 - a1*b3
    n3 = a1*b2 - a2*b1
    # Normalize
    mag_n = sqrt(n1*n1 + n2*n2 + n3*n3)
    n1 /= mag_n
    n2 /= mag_n
    n3 /= mag_n
return n1, n2, n3
end

@doc raw"""
    my3SolveAb(A11::Float64, A21::Float64, A31::Float64,
               A12::Float64, A22::Float64, A32::Float64,
               A13::Float64, A23::Float64, A33::Float64,
               bx::Float64, by::Float64, bz::Float64)

Solves the 3×3 linear system `A * x = b`, where `A` is:
    [ A11  A12  A13 ]
    [ A21  A22  A23 ]
    [ A31  A32  A33 ]
and `b = (bx, by, bz)`, returning `(x1, x2, x3)`.

It effectively does a manual determinant-based inverse of A:

```math
(A^-1) = 1/det(A) * adj(A)
then multiplies by b. """ 
function my3SolveAb(A11::Float64, A21::Float64, A31::Float64, A12::Float64, A22::Float64, A32::Float64, A13::Float64, A23::Float64, A33::Float64, bx::Float64, by::Float64, bz::Float64)
# Compute determinant of A
v = ( A11*A22*A33
    + A12*A23*A31
    + A13*A21*A32
    - A11*A23*A32
    - A12*A21*A33
    - A13*A22*A31 )

# Adjugate (cofactor) approach
im11 = A22*A33 - A23*A32
im21 = A23*A31 - A21*A33
im31 = A21*A32 - A22*A31

im12 = A13*A32 - A12*A33
im22 = A11*A33 - A13*A31
im32 = A12*A31 - A11*A32

im13 = A12*A23 - A13*A22
im23 = A13*A21 - A11*A23
im33 = A11*A22 - A12*A21

# Multiply by b, then divide by det(A)
x1 = (im11*bx + im12*by + im13*bz) / v
x2 = (im21*bx + im22*by + im23*bz) / v
x3 = (im31*bx + im32*by + im33*bz) / v

return x1, x2, x3
end

@doc raw"""
    cal_stress_hexa(Qe, integ_stress::Array{Float64,2}, integ_strain::Array{Float64,2},
                    integ_yield_stress::Array{Float64,1}, integ_eq_plastic_strain::Array{Float64,1},
                    position_::Array{Float64,2}, d_disp_::Array{Float64,1},
                    elementmat::Array{Int,2}, element_flag::Array{Int,1}, integ_num::Int,
                    Pusai_mat::Vector{SMatrix{3,8,Float64}},  # or your MMatrix type
                    MATERIAL::Array{MaterialType,1},
                    element_material::Array{Int,1},
                    elementMinSize::Float64,
                    elementVolume::Array{Float64,1})

Computes the internal forces for each element (`Qe`) as well as updating integration-point
stress/strain with an 8-integration-point hexahedral formulation. The parameters are:

  • `Qe`        : a 24×nElement array of element-level force vectors (24 = 3 dofs × 8 nodes).
  • `integ_stress` : a 6×(nElement*integ_num) array of stress at each Gauss point, updated in-place.
  • `integ_strain` : a 6×(nElement*integ_num) array of strain at each Gauss point, updated in-place.
  • `integ_yield_stress`, `integ_eq_plastic_strain` : 1D arrays for yield stress and eq. plastic strain, each has nElement*integ_num length.
  • `position_` : a 3×nNode array of current node coordinates.
  • `d_disp_`   : a (3×nNode) vector storing this step's incremental displacement (disp_new - disp_old).
  • `elementmat`: an 8×nElement array listing node indices of each element.
  • `element_flag`: integer array (1 or 0) indicating which elements are active vs. deleted.
  • `integ_num` : number of integration points per element (8).
  • `Pusai_mat` : the shape-function derivative data for each integration point.
  • `MATERIAL`  : array of your custom `MaterialType`, from which we get Dmat, G, plastic table, etc.
  • `element_material` : an array giving the material ID for each element.
  • `elementMinSize`   : minimum element size (used in some checks, can be ignored if not used).
  • `elementVolume`    : array of size nElement for (maybe) updated volume if the code modifies it.

**What it does**:
- Loops over each element (if `element_flag[e] == 1`).
- Gathers that element’s incremental displacement, forms the B-matrix, and obtains the stress increment.
- If the material is plastic, applies a radial return or linear hardening logic (the code does a simple “incremental” approach).
- Updates `integ_stress` and `integ_strain` at each Gauss point.
- Accumulates the element internal force vector into the `Qe` array.

**Returns** nothing, but modifies `Qe`, `integ_stress`, `integ_strain`, etc. in-place.
"""
function cal_stress_hexa(Qe,
                         integ_stress,
                         integ_strain,
                         integ_yield_stress,
                         integ_eq_plastic_strain,
                         position_,
                         d_disp_,
                         elementmat::Array{Int,2},
                         element_flag::Array{Int,1},
                         integ_num::Int,
                         Pusai_mat,
                         MATERIAL::Array{MaterialType,1},
                         element_material::Array{Int,1},
                         elementMinSize::Float64,
                         elementVolume::Array{Float64,1})

    nElement = length(element_flag)      # total elements
    fn       = size(position_, 2) * 3    # total degrees of freedom (3 * nNode)
    W        = 1.0                       # integration weight (assuming each Gauss point w=1)

    # Loop over elements
    d_u = zeros(24)
    e_position = zeros(3, 8)
    BVbar = zeros(6,24)
    Bfinal = zeros(6,24)
    final_stress = zeros(6)
    d_e_vec = zeros(6)
    d_o_vec = zeros(6)
    q_vec = zeros(24)
    q_vec_i = zeros(24)
    @inbounds for e in 1:nElement
        if element_flag[e] == 0
            continue
        end

        mat_id = element_material[e]
        mat    = MATERIAL[mat_id]
        G      = mat.G
        plastic_data = mat.plastic
        Hd          = mat.Hd
        npp         = size(plastic_data, 1)  # # of rows in plastic table

        # Build local Dmat from the material
        Dmat = mat.Dmat 

        # Extract the 8 nodes of this element
        node_ids = @view elementmat[:, e]
        # gather incremental displacement for these 8 nodes into a 24-vector
        d_u .=0
        e_position .=0

        @inbounds for i = 1:8
            d_u[3*(i-1) + 1] = d_disp_[3*(node_ids[i]-1)+1]
            d_u[3*(i-1) + 2] = d_disp_[3*(node_ids[i]-1)+2]
            d_u[3*(i-1) + 3] = d_disp_[3*(node_ids[i]-1)+3]

            e_position[1,i]  = position_[1, node_ids[i]]
            e_position[2,i]  = position_[2, node_ids[i]]
            e_position[3,i]  = position_[3, node_ids[i]]
        end

        # We'll accumulate the element internal force in q_vec,
        # then store it into Qe afterwards
        q_vec .= 0

        # B̄ computations for "mean dilatation" or "BVbar" approach
        # The function cal_BVbar_hexa returns total volume & the BVbar array
        BVbar .=0
        V     = cal_BVbar_hexa(Pusai_mat, e_position, BVbar)
        elementVolume[e] = V  # store updated volume if needed

        # For each integration point i in 1:8
        Bfinal .= 0
        @inbounds for i = 1:integ_num
            Bfinal .= 0.0
            # compute final B = B(Pusai_mat[i]) + (stuff) - ...
            detJ = cal_Bfinal(Bfinal, BVbar, Pusai_mat[i], e_position)

            # compute strain increment = Bfinal * d_u
            mul!(d_e_vec, Bfinal, d_u)
            mul!(d_o_vec, Dmat, d_e_vec)

            # index for the integration point data (6× N*gpoints)
            idx_integ = (e-1)*integ_num + i

            # gather old stress
            pre_stress = @view(integ_stress[1:6, idx_integ])

            # new stress (trial)
            final_stress .= pre_stress .+ d_o_vec

            # plastic check
            if !isempty(plastic_data)
                # compute dev. stress
                tri_stress = final_stress
                mean_stress = (tri_stress[1] + tri_stress[2] + tri_stress[3]) / 3.0
                tri_dev_stress = @MVector [
                    tri_stress[1] - mean_stress,
                    tri_stress[2] - mean_stress,
                    tri_stress[3] - mean_stress,
                    tri_stress[4],
                    tri_stress[5],
                    tri_stress[6]
                ]

                tri_mises_stress = sqrt(1.5 * (
                    tri_dev_stress[1]^2 + tri_dev_stress[2]^2 + tri_dev_stress[3]^2 +
                    2*tri_dev_stress[4]^2 + 2*tri_dev_stress[5]^2 + 2*tri_dev_stress[6]^2
                ))

                y = integ_yield_stress[idx_integ]
                if tri_mises_stress > y
                    # find plastic slope
                    # do a linear interpolation in eq_plastic_strain
                    ep = integ_eq_plastic_strain[idx_integ]
                    p_index = 1
                    for j_ in 2:npp
                        if ep <= plastic_data[j_,2]
                            p_index = j_-1
                            break
                        end
                        if j_ == npp
                            p_index = j_-1
                        end
                    end
                    # Hardening
                    H  = Hd[p_index]
                    d_ep = (tri_mises_stress - y) / (3*G + H)

                    # update yield stress and eq. plastic
                    integ_eq_plastic_strain[idx_integ] += d_ep
                    integ_yield_stress[idx_integ]      += H * d_ep

                    # do radial return
                    scale = (y + H*d_ep) / tri_mises_stress
                    final_dev_stress = tri_dev_stress .* scale
                    # add the mean back in
                    final_stress[1] = final_dev_stress[1] + mean_stress
                    final_stress[2] = final_dev_stress[2] + mean_stress
                    final_stress[3] = final_dev_stress[3] + mean_stress
                    final_stress[4] = final_dev_stress[4]
                    final_stress[5] = final_dev_stress[5]
                    final_stress[6] = final_dev_stress[6]
                end
            end

            # update integ_stress, integ_strain
            for k in 1:6
                # accumulate strain
                integ_strain[k, idx_integ] += d_e_vec[k]
                # store new stress
                integ_stress[k, idx_integ] = final_stress[k]
            end

            # form the final stress increment (final_stress - pre_stress) for building q_vec
            d_o_vec .= final_stress .- pre_stress
            # compute q_vec_i = Bfinal' * stress
            q_vec_i .= 0
            for i in eachindex(q_vec_i)
                @inbounds @simd for j in eachindex(final_stress)
                    q_vec_i[i] += Bfinal[j,i] * final_stress[j]
                end
            end

            # accumulate into q_vec (weight, detJ, etc.)
            for j_ in 1:24
                q_vec[j_] += W * W * W * detJ * q_vec_i[j_]
            end
        end  # end each integration point

        # store q_vec into Qe
        @inbounds for i in 1:8
            Qe[3*(i-1)+1, e] = q_vec[3*(i-1)+1]
            Qe[3*(i-1)+2, e] = q_vec[3*(i-1)+2]
            Qe[3*(i-1)+3, e] = q_vec[3*(i-1)+3]
        end
    end  # end for e in 1:nElement

    return
end

@doc raw"""
    cal_BVbar_hexa(Pusai_mat, e_position, BVbar)

Given:
 - `Pusai_mat`: a Vector of 8 SMatrix{3,8,Float64} (or MMatrix{3,8}), each storing the local shape
    function derivatives at one Gauss point of a hexahedron.
 - `e_position`: a 3×8 matrix of the current coordinates of this element’s nodes.
 - `BVbar`: a 6×24 MMatrix (or Array) to which we add the volumetric parts for the B-bar approach.

This function:
  1. Loops over each of the 8 Gauss points.
  2. Computes the determinant of the Jacobian (`detJi`) from `Pusai_mat[k] * e_position`.
  3. Accumulates the (1/3)*N*(detJi) terms into `BVbar[1:3, :]`.
  4. Returns the total element volume by summing the 8 `detJi` values.

The final `BVbar` is then used in a “Bbar” method (like Simo/Hughes anti-volumetric locking approach).
"""
function cal_BVbar_hexa(Pusai_mat,
                        e_position,
                        BVbar)

    V = 0.0

    # Initialize BVbar to zero
    BVbar .= 0.0

    # Gauss integration typically has 8 points for a standard hex
    for k in 1:8
        # shape function derivatives w.r.t local coords
        Pusai1 = Pusai_mat[k]

        # Compute the Jacobian (3×3) by summation: J = Pusai1 * e_position
        J11 = 0.0; J12 = 0.0; J13 = 0.0
        J21 = 0.0; J22 = 0.0; J23 = 0.0
        J31 = 0.0; J32 = 0.0; J33 = 0.0

        @inbounds for i = 1:8
            # local derivatives
            gx = Pusai1[1,i]
            gy = Pusai1[2,i]
            gz = Pusai1[3,i]

            # node coords
            x = e_position[1,i]
            y = e_position[2,i]
            z = e_position[3,i]

            # accumulate in J
            J11 += gx*x;  J12 += gx*y;  J13 += gx*z
            J21 += gy*x;  J22 += gy*y;  J23 += gy*z
            J31 += gz*x;  J32 += gz*y;  J33 += gz*z
        end

        # determinant
        detJi = J11*(J22*J33 - J23*J32) -
                J12*(J21*J33 - J23*J31) +
                J13*(J21*J32 - J22*J31)

        if detJi < 0
            # If negative volume => we can take abs or raise a warning
            detJi = -detJi
            # println("Warning: negative element volume encountered.")
        end

        V += detJi
        @inbounds for i = 1:8
            # shape function derivative = Pusai1[:, i]
            # We interpret it as (gx, gy, gz) => each row 1:3 of BVbar
            # scaled by detJi / 3.0
            nx = Pusai1[1,i] / 3.0 * detJi
            ny = Pusai1[2,i] / 3.0 * detJi
            nz = Pusai1[3,i] / 3.0 * detJi

            colx = (i-1)*3 + 1
            coly = (i-1)*3 + 2
            colz = (i-1)*3 + 3

            BVbar[1, colx] += nx
            BVbar[1, coly] += ny
            BVbar[1, colz] += nz

            BVbar[2, colx] += nx
            BVbar[2, coly] += ny
            BVbar[2, colz] += nz

            BVbar[3, colx] += nx
            BVbar[3, coly] += ny
            BVbar[3, colz] += nz
        end
    end

    # Finally, we divide the entire BVbar by the total volume
    # to get the average volumetric portion
    BVbar .= BVbar ./ V

    # Return the total volume
    return V
end

@doc raw"""
    cal_Bfinal(Bfinal, BVbar, Pusai1, e_position)

Computes the final strain-displacement matrix `Bfinal` for a single Gauss point of an 8-node
hexahedron using a "B-bar" type approach. Specifically:
  1) Compute the local 3×3 Jacobian from `Pusai1` and `e_position`.
  2) Invert it and build the base B-matrix for the standard part.
  3) Subtract or add the volumetric correction from `BVbar`.
  4) Return the determinant of the Jacobian (for weighting).

**Arguments**:
  • `Bfinal`  : a 6×24 matrix (initially zero) where we'll store the final B at this Gauss point.
  • `BVbar`   : a 6×24 matrix containing the volumetric part from `cal_BVbar_hexa`.
  • `Pusai1`  : a 3×8 matrix of shape-function derivatives at this Gauss point (from `Pusai_mat`).
  • `e_position`: a 3×8 matrix of the current coordinates of this element’s 8 nodes.

**Returns**:
  • `detJi`   : The determinant of the 3×3 Jacobian, used as the integration weight factor.

**Notes**:
  • If you strictly follow your original code logic, you might do:
        Bfinal = B - BV + BVbar
    or something similar. Here we embed that logic.  
  • The typical 6×24 B is formed from partial derivatives (dN/dx, etc.) mapped by the Jacobian.
"""
function cal_Bfinal(Bfinal,
                    BVbar,
                    Pusai1,
                    e_position)

    # 1) Compute the 3×3 Jacobian from Pusai1 * e_position
    J11 = J12 = J13 = 0.0
    J21 = J22 = J23 = 0.0
    J31 = J32 = J33 = 0.0

    @inbounds for i in 1:8
        # local derivatives
        gx = Pusai1[1,i]
        gy = Pusai1[2,i]
        gz = Pusai1[3,i]

        # node coords
        x = e_position[1,i]
        y = e_position[2,i]
        z = e_position[3,i]

        J11 += gx * x;  J12 += gx * y;  J13 += gx * z
        J21 += gy * x;  J22 += gy * y;  J23 += gy * z
        J31 += gz * x;  J32 += gz * y;  J33 += gz * z
    end

    # 2) determinant of J
    v = ( J11*J22*J33
        + J12*J23*J31
        + J13*J21*J32
        - J11*J23*J32
        - J12*J21*J33
        - J13*J22*J31 )
    detJi = v
    div_v = 1.0 / v  # for inverting

    # 3) Inverse of J
    iJ11 = ( J22*J33 - J23*J32 ) * div_v
    iJ21 = ( J23*J31 - J21*J33 ) * div_v
    iJ31 = ( J21*J32 - J22*J31 ) * div_v

    iJ12 = ( J13*J32 - J12*J33 ) * div_v
    iJ22 = ( J11*J33 - J13*J31 ) * div_v
    iJ32 = ( J12*J31 - J11*J32 ) * div_v

    iJ13 = ( J12*J23 - J13*J22 ) * div_v
    iJ23 = ( J13*J21 - J11*J23 ) * div_v
    iJ33 = ( J11*J22 - J12*J21 ) * div_v

    # 4) Build the standard B part (call it Bstd)
    Bstd = @MMatrix zeros(6,24)

    @inbounds for i in 1:8
        # local derivatives in (gx, gy, gz)
        gx = iJ11 * Pusai1[1,i] + iJ12 * Pusai1[2,i] + iJ13 * Pusai1[3,i]
        gy = iJ21 * Pusai1[1,i] + iJ22 * Pusai1[2,i] + iJ23 * Pusai1[3,i]
        gz = iJ31 * Pusai1[1,i] + iJ32 * Pusai1[2,i] + iJ33 * Pusai1[3,i]

        col = (i-1)*3
        # Normal components
        Bstd[1, col+1] = gx
        Bstd[2, col+2] = gy
        Bstd[3, col+3] = gz
        # Shear components
        Bstd[4, col+1] = gy
        Bstd[4, col+2] = gx
        Bstd[5, col+2] = gz
        Bstd[5, col+3] = gy
        Bstd[6, col+1] = gz
        Bstd[6, col+3] = gx
    end

    # 5) Combine with BVbar. Often the formula is:
    #      Bfinal = Bstd + BVbar
    #    or sometimes Bstd - BVstd + BVbar.  In your code you might do "B - (BV) + BVbar".
    #    We replicate a simpler version: just add Bstd + (BVbar) - (some local BV if needed).
    #    For now, let's add them:
    @inbounds @simd for r in 1:6
        @inbounds @simd for c in 1:24
            Bfinal[r, c] = Bstd[r,c] + BVbar[r,c]
        end
    end

    return detJi
end
