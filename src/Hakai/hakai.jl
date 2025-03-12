include("DataStructure.jl")
include("setup.jl")
include("math.jl")
include("cal.jl")
include("vtk.jl")
include("Mesh.jl")

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


