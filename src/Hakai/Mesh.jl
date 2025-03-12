
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