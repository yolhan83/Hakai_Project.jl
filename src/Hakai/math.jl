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