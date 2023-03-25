using LinearAlgebra

print("hello")

x = [1.0 0.0 2.0]
# x1 = [1.0 0.0 2.0]
# x2 = [2.0 1.0 2.0]
# x3 = [5.0 2.0 2.0]
y = [0.0 0.0 2.0]


"""
 Calculate distance |r1-r2|. A small numerical cutoff (eps)
 may help avoiding masses going to infinity.
"""
function dist(r1::Array{Float64,2}, r2::Array{Float64,2}, eps::Float64)::Float64
    return sqrt(sum((r1 - r2) * (r1 - r2)') + eps)
end

"""
Calculate versor n_{ij} with appropriate direction (r1-r2) for force from
body 1 to 2.

Parameters
----------
r1, r2: np.ndarray
    Position vectors with shape (3).

Returns
-------
versor : np.ndarray
    versor n_{12} with shape (3)

"""
function nvec(r1::Array{Float64,2}, r2::Array{Float64,2}, eps::Float64)::Array{Float64,2}
    return (r1 - r2) / dist(r1, r2, eps)
end

# ------------------ Force terms  ------------------------


"""
Newtonian gravitational force.

Parameters
----------
x: np.ndarray
    Vectors for positions. Shape (3)

Returns
-------
f_pn0 : np.ndarray
    Newtonian contribution. Shape (3).
"""
function PN0(
    m2::Float64,
    m3::Float64,
    x::Vector{Matrix{Float64}},
    eps::Float64,
)::Array{Float64,2}
    x1, x2, x3 = x
    f12 = m2 * nvec(x2, x1, eps) / dist(x1, x2, eps)^2
    f13 = m3 * nvec(x3, x1, eps) / dist(x1, x3, eps)^2
    f_pn0 = f12 + f13
    return f_pn0
end


x1 = [1.0 0.0 2.0]
x2 = [2.0 1.0 3.0]
x3 = [5.0 2.0 1.3]
v1 = [1.0 0.0 2.0] * 2
v2 = [2.0 1.0 3.0] * 2
v3 = [5.0 2.0 1.3] * 2
xt = [x1, x2, x3]
vt = [v1, v2, v3]
print(x)
m2 = 1.0
m3 = 3.0
PN0(m2, m3, xt, 0.0)
a = rand(1:10, (3, 3, 4))
println("teste")


"""
First order relativistic correction (order c^{-2}) for newtonian gravitational force.

Parameters
----------
x, v: np.ndarray
    Vectors positions and velocities. Shape (3).

Returns
-------
f_pn1 : np.ndarray
    PN correction. Shape (3).

"""
global c = 1
function PN1(
    m1::Float64,
    m2::Float64,
    m3::Float64,
    x::Vector{Matrix{Float64}},
    v::Vector{Matrix{Float64}},
    eps::Float64
)::Array{Float64,2}
    x1, x2, x3 = x
    v1, v2, v3 = v
    f_pn1 =
        (
            (dot(nvec(x1, x2, eps), (4.0 * v1 - 3.0 * v2)) * m2 * (v1 - v2)) /
            dist(x1, x2, eps)^2 + 
            (dot(nvec(x1, x3, eps), 4.0 * v1 - 3.0 * v3) * m3 * (v1 - v3)) /
            dist(x1, x3, eps)^2 +
            (
                m2 *
                nvec(x1, x2, eps) .*
                (
                    4.0 * dot(v1, v2) .+ (3.0 * dot(v2, nvec(x1, x2, eps))^2) / 2.0 .-
                    v1 * v1' .- 2.0 * v2 * v2' .+
                    (5.0 * m1) / dist(x1, x2, eps) .+
                    (4.0 * m2) / dist(x1, x2, eps) .+
                    (4.0 * m3) / dist(x1, x3, eps) .-
                    (dot(nvec(x1, x2, eps), nvec(x2, x3, eps)) * m3 * dist(x1, x2, eps)) /
                    (2.0 * dist(x2, x3, eps)^2) .+ m3 / dist(x2, x3, eps)
                )
               ) / dist(x1, x2, eps)^2 -
            (7.0 * m2 * m3 * nvec(x2, x3, eps)) /
            (2.0 * dist(x1, x2, eps) * dist(x2, x3, eps)^2) +
            (
                m3 *
                nvec(x1, x3, eps) .*
                (
                    4.0 * dot(v1, v3) .+ (3.0 * dot(v3, nvec(x1, x3, eps))^2) / 2.0 .-
                    v1 * v1' .- 2.0 * v3 * v3' .+
                    (4.0 * m2) / dist(x1, x2, eps) .+
                    (5.0 * m1) / dist(x1, x3, eps) .+
                    (4.0 * m3) / dist(x1, x3, eps) .-
                    (dot(nvec(x1, x3, eps), nvec(x3, x2, eps)) * m2 * dist(x1, x3, eps)) /
                    (2.0 * dist(x3, x2, eps)^2) .+ m2 / dist(x3, x2, eps)
                )
            ) / dist(x1, x3, eps)^2 -
            (7.0 * m2 * m3 * nvec(x3, x2, eps)) /
            (2.0 * dist(x1, x3, eps) * dist(x3, x2, eps)^2)
           ) / c^2
    return f_pn1
end


# """
# Second order relativistic correction (order c^{-4}) for newtonian gravitational force.
#
# Parameters
# ----------
# x, v: np.ndarray
#     Vectors positions and velocities. Shape (3).
#
# Returns
# -------
# f_pn2 : np.ndarray
#     PN2 correction. Shape (3).
#
# """
# function PN2(
#     m1::Float64,
#     m2::Float64,
#     m3::Float64,
#     x::Vector{Matrix{Float64}},
#     v::Vector{Matrix{Float64}},
#     eps::Float64
#    )::Array{Float64,2}
#     x1, x2, x3 = x
#     v1, v2, v3 = v
#
#     term1 = (
#         m2
#         * nvec(x1, x2)
#         * (
#             -2 * dot(v1, v2) ** 2
#             - 6 * dot(v1, v2) * dot(nvec(x1, x2), v2) ** 2
#             - (15 * dot(nvec(x1, x2), v2) ** 4) / 8.0
#             + (3 * dot(nvec(x1, x2), v2) ** 2 * v1**2) / 2.0
#             + 4 * dot(v1, v2) * v2**2
#             + (9 * dot(nvec(x1, x2), v2) ** 2 * v2**2) / 2.0
#             - 2 * v2**4
#             - (57 * m1**2) / (4.0 * dist(x1, x2) ** 2)
#             - (69 * m1 * m2) / (2.0 * dist(x1, x2) ** 2)
#             - (9 * m2**2) / dist(x1, x2) ** 2
#             + (
#                 m1
#                 * (
#                     (-5 * dot(v1, v2)) / 2.0
#                     + (39 * dot(nvec(x1, x2), v1) ** 2) / 2.0
#                     - 39 * dot(nvec(x1, x2), v1) * dot(nvec(x1, x2), v2)
#                     + (17 * dot(nvec(x1, x2), v2) ** 2) / 2.0
#                     - (15 * v1**2) / 4.0
#                     + (5 * v2**2) / 4.0
#                 )
#             )
#             / dist(x1, x2)
#             + (
#                 m2
#                 * (
#                     -8 * dot(v1, v2)
#                     + 2 * dot(nvec(x1, x2), v1) ** 2
#                     - 4 * dot(nvec(x1, x2), v1) * dot(nvec(x1, x2), v2)
#                     - 6 * dot(nvec(x1, x2), v2) ** 2
#                     + 4 * v2**2
#                 )
#             )
#             / dist(x1, x2)
#         )
#     ) / (c**4 * dist(x1, x2) ** 2) + (
#         m3
#         * nvec(x1, x3)
#         * (
#             -2 * dot(v1, v3) ** 2
#             - 6 * dot(v1, v3) * dot(nvec(x1, x3), v3) ** 2
#             - (15 * dot(nvec(x1, x3), v3) ** 4) / 8.0
#             + (3 * dot(nvec(x1, x3), v3) ** 2 * v1**2) / 2.0
#             + 4 * dot(v1, v3) * v3**2
#             + (9 * dot(nvec(x1, x3), v3) ** 2 * v3**2) / 2.0
#             - 2 * v3**4
#             - (57 * m1**2) / (4.0 * dist(x1, x3) ** 2)
#             - (69 * m1 * m3) / (2.0 * dist(x1, x3) ** 2)
#             - (9 * m3**2) / dist(x1, x3) ** 2
#             + (
#                 m1
#                 * (
#                     (-5 * dot(v1, v3)) / 2.0
#                     + (39 * dot(nvec(x1, x3), v1) ** 2) / 2.0
#                     - 39 * dot(nvec(x1, x3), v1) * dot(nvec(x1, x3), v3)
#                     + (17 * dot(nvec(x1, x3), v3) ** 2) / 2.0
#                     - (15 * v1**2) / 4.0
#                     + (5 * v3**2) / 4.0
#                 )
#             )
#             / dist(x1, x3)
#             + (
#                 m3
#                 * (
#                     -8 * dot(v1, v3)
#                     + 2 * dot(nvec(x1, x3), v1) ** 2
#                     - 4 * dot(nvec(x1, x3), v1) * dot(nvec(x1, x3), v3)
#                     - 6 * dot(nvec(x1, x3), v3) ** 2
#                     + 4 * v3**2
#                 )
#             )
#             / dist(x1, x3)
#         )
#     ) / (
#         c**4 * dist(x1, x3) ** 2
#     )
#     term2 = (
#         m2
#         * (
#             -6 * dot(nvec(x1, x2), v1) * dot(nvec(x1, x2), v2) ** 2
#             + (9 * dot(nvec(x1, x2), v2) ** 3) / 2.0
#             - 4 * dot(v1, v2) * dot(nvec(x1, x2), (v1 - v2))
#             + dot(nvec(x1, x2), v2) * v1**2
#             + 4 * dot(nvec(x1, x2), v1) * v2**2
#             - 5 * dot(nvec(x1, x2), v2) * v2**2
#             + (
#                 (
#                     (-63 * dot(nvec(x1, x2), v1)) / 4.0
#                     + (55 * dot(nvec(x1, x2), v2)) / 4.0
#                 )
#                 * m1
#             )
#             / dist(x1, x2)
#             - (2 * (dot(nvec(x1, x2), v1) + dot(nvec(x1, x2), v2)) * m2)
#             / dist(x1, x2)
#         )
#         * (v1 - v2)
#     ) / (c**4 * dist(x1, x2) ** 2) + (
#         m3
#         * (
#             -6 * dot(nvec(x1, x3), v1) * dot(nvec(x1, x3), v3) ** 2
#             + (9 * dot(nvec(x1, x3), v3) ** 3) / 2.0
#             - 4 * dot(v1, v3) * dot(nvec(x1, x3), (v1 - v3))
#             + dot(nvec(x1, x3), v3) * v1**2
#             + 4 * dot(nvec(x1, x3), v1) * v3**2
#             - 5 * dot(nvec(x1, x3), v3) * v3**2
#             + (
#                 (
#                     (-63 * dot(nvec(x1, x3), v1)) / 4.0
#                     + (55 * dot(nvec(x1, x3), v3)) / 4.0
#                 )
#                 * m1
#             )
#             / dist(x1, x3)
#             - (2 * (dot(nvec(x1, x3), v1) + dot(nvec(x1, x3), v3)) * m3)
#             / dist(x1, x3)
#         )
#         * (v1 - v3)
#     ) / (
#         c**4 * dist(x1, x3) ** 2
#     )
#
#     f_pn2 = term1 + term2
#     return f_pn2
# end


# """
# Second order .5 relativistic correction (order c^{-5}) for newtonian gravitational force.
#
# Parameters
# ----------
# x, v: np.ndarray
#     Vectors positions and velocities. Shape (3).
#
# Returns
# -------
# f_pn2p5 : np.ndarray
#     PN2.5 correction. Shape (3).
#
# """
# function PN25(
#     m1::Float64,
#     m2::Float64,
#     m3::Float64,
#     x::Vector{Matrix{Float64}},
#     v::Vector{Matrix{Float64}},
#     eps::Float64
#    )::Array{Float64,2}
#     x1, x2, x3 = x
#     v1, v2, v3 = v
#
#     f_n2p5 = (
#         4
#         * m1
#         * m2
#         * (
#             (v1 - v2)
#             * ((2 * m1) / dist(x1, x2) - (8 * m2) / dist(x1, x2) - (v1 - v2) ** 2)
#             + dot(nvec(x1, x2), (v1 - v2))
#             * nvec(x1, x2)
#             * (
#                 (-6 * m1) / dist(x1, x2)
#                 + (52 * m2) / (3.0 * dist(x1, x2))
#                 + 3 * (v1 - v2) ** 2
#             )
#         )
#     ) / (5.0 * c**5 * dist(x1, x2) ** 3) + (
#         4
#         * m1
#         * m3
#         * (
#             (v1 - v3)
#             * ((2 * m1) / dist(x1, x3) - (8 * m3) / dist(x1, x3) - (v1 - v3) ** 2)
#             + dot(nvec(x1, x3), (v1 - v3))
#             * nvec(x1, x3)
#             * (
#                 (-6 * m1) / dist(x1, x3)
#                 + (52 * m3) / (3.0 * dist(x1, x3))
#                 + 3 * (v1 - v3) ** 2
#             )
#         )
#     ) / (
#         5.0 * c**5 * dist(x1, x3) ** 3
#     )
#     return f_n2p5
# end
#
#
# # ------------------ Solvers  ------------------------
