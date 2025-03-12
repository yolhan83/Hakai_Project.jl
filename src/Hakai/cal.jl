
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
