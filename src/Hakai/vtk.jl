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
