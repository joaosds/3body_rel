import numpy as np
from diffeqpy import de
import sys
import matplotlib.pyplot as plt
import os
import time

# diffeqpy.install()
c = 1
anispeed = 0.08  # 0.0002
t_f = 108
# TD:
# Testar intervalo de integração até convergência para outros métodos

m = np.array([1, 1, 1])
m1, m2, m3 = m


def dist(r1: np.ndarray, r2: np.ndarray) -> float:
    # eps = numerical cutoff of 1e-4 to avoid masses going to infinity
    eps = 0
    return np.sqrt(np.sum(r1 - r2) ** 2 + eps)


def nvec(r1: np.ndarray, r2: np.ndarray) -> np.ndarray:
    # Já inclui a direção
    return (r1 - r2) / dist(r1, r2)


def PN1(m: np.ndarray, x: np.ndarray, v: np.ndarray) -> np.ndarray:
    """
    First order relativistic correction (order c^{-2}) for newtonian gravitational force.
    """
    m1, m2, m3 = m
    x1, x2, x3 = x
    v1, v2, v3 = v
    cu = (
        (np.dot(nvec(x1, x2), (4 * v1 - 3 * v2)) * m2 * (v1 - v2)) / dist(x1, x2) ** 2
        + (np.dot(nvec(x1, x3), 4 * v1 - 3 * v3) * m3 * (v1 - v3)) / dist(x1, x3) ** 2
        + (
            m2
            * nvec(x1, x2)
            * (
                4 * np.dot(v1, v2)
                + (3 * np.dot(v2, nvec(x1, x2)) ** 2) / 2
                - v1**2
                - 2 * v2**2
                + (5 * m1) / dist(x1, x2)
                + (4 * m2) / dist(x1, x2)
                + (4 * m3) / dist(x1, x3)
                - (np.dot(nvec(x1, x2), nvec(x2, x3)) * m3 * dist(x1, x2))
                / (2.0 * dist(x2, x3) ** 2)
                + m3 / dist(x2, x3)
            )
        )
        / dist(x1, x2) ** 2
        - (7 * m2 * m3 * nvec(x2, x3)) / (2.0 * dist(x1, x2) * dist(x2, x3) ** 2)
        + (
            m3
            * nvec(x1, x3)
            * (
                4 * np.dot(v1, v3)
                + (3 * np.dot(v3, nvec(x1, x3)) ** 2) / 2
                - v1**2
                - 2 * v3**2
                + (4 * m2) / dist(x1, x2)
                + (5 * m1) / dist(x1, x3)
                + (4 * m3) / dist(x1, x3)
                - (np.dot(nvec(x1, x3), nvec(x3, x2)) * m2 * dist(x1, x3))
                / (2.0 * dist(x3, x2) ** 2)
                + m2 / dist(x3, x2)
            )
        )
        / dist(x1, x3) ** 2
        - (7 * m2 * m3 * nvec(x3, x2)) / (2.0 * dist(x1, x3) * dist(x3, x2) ** 2)
    ) / c**2
    return cu


def PN0(x: np.ndarray) -> np.ndarray:
    """
    Newtonian gravitational force.
    """
    # _, m2, m3 = m
    x1, x2, x3 = x
    f12 = m2 * nvec(x2, x1) / dist(x1, x2) ** 2
    f13 = m3 * nvec(x3, x1) / dist(x1, x3) ** 2

    return f12 + f13


def PN2(m: np.ndarray, x: np.ndarray, v: np.ndarray) -> np.ndarray:

    """
    Second order relativistic correction (order c^{-4}) for newtonian gravitational force.
    """
    m1, m2, m3 = m
    x1, x2, x3 = x
    v1, v2, v3 = v

    term1 = (
        m2
        * nvec(x1, x2)
        * (
            -2 * np.dot(v1, v2) ** 2
            - 6 * np.dot(v1, v2) * np.dot(nvec(x1, x2), v2) ** 2
            - (15 * np.dot(nvec(x1, x2), v2) ** 4) / 8.0
            + (3 * np.dot(nvec(x1, x2), v2) ** 2 * v1**2) / 2.0
            + 4 * np.dot(v1, v2) * v2**2
            + (9 * np.dot(nvec(x1, x2), v2) ** 2 * v2**2) / 2.0
            - 2 * v2**4
            - (57 * m1**2) / (4.0 * dist(x1, x2) ** 2)
            - (69 * m1 * m2) / (2.0 * dist(x1, x2) ** 2)
            - (9 * m2**2) / dist(x1, x2) ** 2
            + (
                m1
                * (
                    (-5 * np.dot(v1, v2)) / 2.0
                    + (39 * np.dot(nvec(x1, x2), v1) ** 2) / 2.0
                    - 39 * np.dot(nvec(x1, x2), v1) * np.dot(nvec(x1, x2), v2)
                    + (17 * np.dot(nvec(x1, x2), v2) ** 2) / 2.0
                    - (15 * v1**2) / 4.0
                    + (5 * v2**2) / 4.0
                )
            )
            / dist(x1, x2)
            + (
                m2
                * (
                    -8 * np.dot(v1, v2)
                    + 2 * np.dot(nvec(x1, x2), v1) ** 2
                    - 4 * np.dot(nvec(x1, x2), v1) * np.dot(nvec(x1, x2), v2)
                    - 6 * np.dot(nvec(x1, x2), v2) ** 2
                    + 4 * v2**2
                )
            )
            / dist(x1, x2)
        )
    ) / (c**4 * dist(x1, x2) ** 2) + (
        m3
        * nvec(x1, x3)
        * (
            -2 * np.dot(v1, v3) ** 2
            - 6 * np.dot(v1, v3) * np.dot(nvec(x1, x3), v3) ** 2
            - (15 * np.dot(nvec(x1, x3), v3) ** 4) / 8.0
            + (3 * np.dot(nvec(x1, x3), v3) ** 2 * v1**2) / 2.0
            + 4 * np.dot(v1, v3) * v3**2
            + (9 * np.dot(nvec(x1, x3), v3) ** 2 * v3**2) / 2.0
            - 2 * v3**4
            - (57 * m1**2) / (4.0 * dist(x1, x3) ** 2)
            - (69 * m1 * m3) / (2.0 * dist(x1, x3) ** 2)
            - (9 * m3**2) / dist(x1, x3) ** 2
            + (
                m1
                * (
                    (-5 * np.dot(v1, v3)) / 2.0
                    + (39 * np.dot(nvec(x1, x3), v1) ** 2) / 2.0
                    - 39 * np.dot(nvec(x1, x3), v1) * np.dot(nvec(x1, x3), v3)
                    + (17 * np.dot(nvec(x1, x3), v3) ** 2) / 2.0
                    - (15 * v1**2) / 4.0
                    + (5 * v3**2) / 4.0
                )
            )
            / dist(x1, x3)
            + (
                m3
                * (
                    -8 * np.dot(v1, v3)
                    + 2 * np.dot(nvec(x1, x3), v1) ** 2
                    - 4 * np.dot(nvec(x1, x3), v1) * np.dot(nvec(x1, x3), v3)
                    - 6 * np.dot(nvec(x1, x3), v3) ** 2
                    + 4 * v3**2
                )
            )
            / dist(x1, x3)
        )
    ) / (
        c**4 * dist(x1, x3) ** 2
    )
    term2 = (
        m2
        * (
            -6 * np.dot(nvec(x1, x2), v1) * np.dot(nvec(x1, x2), v2) ** 2
            + (9 * np.dot(nvec(x1, x2), v2) ** 3) / 2.0
            - 4 * np.dot(v1, v2) * np.dot(nvec(x1, x2), (v1 - v2))
            + np.dot(nvec(x1, x2), v2) * v1**2
            + 4 * np.dot(nvec(x1, x2), v1) * v2**2
            - 5 * np.dot(nvec(x1, x2), v2) * v2**2
            + (
                (
                    (-63 * np.dot(nvec(x1, x2), v1)) / 4.0
                    + (55 * np.dot(nvec(x1, x2), v2)) / 4.0
                )
                * m1
            )
            / dist(x1, x2)
            - (2 * (np.dot(nvec(x1, x2), v1) + np.dot(nvec(x1, x2), v2)) * m2)
            / dist(x1, x2)
        )
        * (v1 - v2)
    ) / (c**4 * dist(x1, x2) ** 2) + (
        m3
        * (
            -6 * np.dot(nvec(x1, x3), v1) * np.dot(nvec(x1, x3), v3) ** 2
            + (9 * np.dot(nvec(x1, x3), v3) ** 3) / 2.0
            - 4 * np.dot(v1, v3) * np.dot(nvec(x1, x3), (v1 - v3))
            + np.dot(nvec(x1, x3), v3) * v1**2
            + 4 * np.dot(nvec(x1, x3), v1) * v3**2
            - 5 * np.dot(nvec(x1, x3), v3) * v3**2
            + (
                (
                    (-63 * np.dot(nvec(x1, x3), v1)) / 4.0
                    + (55 * np.dot(nvec(x1, x3), v3)) / 4.0
                )
                * m1
            )
            / dist(x1, x3)
            - (2 * (np.dot(nvec(x1, x3), v1) + np.dot(nvec(x1, x3), v3)) * m3)
            / dist(x1, x3)
        )
        * (v1 - v3)
    ) / (
        c**4 * dist(x1, x3) ** 2
    )

    return term1 + term2


def PN25(m: np.ndarray, x: np.ndarray, v: np.ndarray) -> np.ndarray:
    """
    2.5 order relativistic correction (order c^{-5}) for newtonian gravitational force.
    """
    m1, m2, m3 = m
    x1, x2, x3 = x
    v1, v2, v3 = v

    term = (
        4
        * m1
        * m2
        * (
            (v1 - v2)
            * ((2 * m1) / dist(x1, x2) - (8 * m2) / dist(x1, x2) - (v1 - v2) ** 2)
            + np.dot(nvec(x1, x2), (v1 - v2))
            * nvec(x1, x2)
            * (
                (-6 * m1) / dist(x1, x2)
                + (52 * m2) / (3.0 * dist(x1, x2))
                + 3 * (v1 - v2) ** 2
            )
        )
    ) / (5.0 * c**5 * dist(x1, x2) ** 3) + (
        4
        * m1
        * m3
        * (
            (v1 - v3)
            * ((2 * m1) / dist(x1, x3) - (8 * m3) / dist(x1, x3) - (v1 - v3) ** 2)
            + np.dot(nvec(x1, x3), (v1 - v3))
            * nvec(x1, x3)
            * (
                (-6 * m1) / dist(x1, x3)
                + (52 * m3) / (3.0 * dist(x1, x3))
                + 3 * (v1 - v3) ** 2
            )
        )
    ) / (
        5.0 * c**5 * dist(x1, x3) ** 3
    )
    return term


## Figure 8 Initial Conditions ##
def init_cond(case: str) -> np.ndarray:

    if case == "fig8":
        Pos_0 = np.array(
            [[-0.97000436, 0.24308753, 0], [0, 0, 0], [0.97000436, -0.24308753, 0]]
        )

        Vel_0 = np.array(
            [
                [0.4662036850, 0.4323657300, 0],
                [-0.93240737, -0.86473146, 0],
                [0.4662036850, 0.4323657300, 0],
            ]
        )

    return np.concatenate((Pos_0, Vel_0))


r1, r2, r3, v1, v2, v3 = init_cond("fig8")


def ThreeBodyEquations2(u, p, t):
    """
    Three-body classical trajectories.
    u - > array coordinates, velocity
    """

    # dv1_dt = PN0(u[:3])
    # dv2_dt = PN0(u[[1, 2, 0]])
    # dv3_dt = PN0(u[[2, 0, 1]])
    # dr1_dt = u[3]
    # dr2_dt = u[4]
    # dr3_dt = u[5]

    du_dt = np.zeros((6, 3))
    # r1, r2, r3 = u[:3]
    # v1, v2, v3 = u[3:6]

    dv1_dt = PN0(u[:3])
    dv2_dt = PN0(u[[1, 2, 0]])
    dv3_dt = PN0(u[[2, 0, 1]])
    dr1_dt = u[3]
    dr2_dt = u[4]
    dr3_dt = u[5]

    r12_derivs = np.concatenate((dr1_dt, dr2_dt))
    r_derivs = np.concatenate((r12_derivs, dr3_dt))
    v12_derivs = np.concatenate((dv1_dt, dv2_dt))
    v_derivs = np.concatenate((v12_derivs, dv3_dt))
    derivs = np.concatenate((r_derivs, v_derivs))
    return derivs


def PENIS0_EQUATION(u, p, t):
    """
    Three-body classical trajectories.
    u - > array coordinates, velocity
    """

    # print(u)
    # dv1_dt = PN0(u[:3])
    # dv2_dt = PN0(u[[1, 2, 0]])
    # dv3_dt = PN0(u[[2, 0, 1]])
    # dr1_dt = u[3]
    # dr2_dt = u[4]
    # dr3_dt = u[5]

    du_dt = np.zeros((6, 3))
    # r1, r2, r3 = u[:3]
    # v1, v2, v3 = u[3:6]

    du_dt[0] = u[3]
    du_dt[1] = u[4]
    du_dt[2] = u[5]
    du_dt[3] = PN0(u[:3])
    du_dt[4] = PN0(u[[1, 2, 0]])
    du_dt[5] = PN0(u[[2, 0, 1]])
    # print(du_dt)

    return du_dt


p = 1
t = 1
PENIS0_EQUATION(init_cond("fig8"), p, t)

tspan = (0.0, t_f)
prob = de.ODEProblem(PENIS0_EQUATION, init_cond("fig8"), tspan)
sol = de.solve(prob, de.RK4(), saveat=anispeed)

solt = np.array(sol.u)
r1 = solt[:, 0, :].T
r2 = solt[:, 1, :].T
r3 = solt[:, 2, :].T
v1 = solt[:, 3, :].T
v2 = solt[:, 4, :].T
v3 = solt[:, 5, :].T

plt.plot(*r1[:2], c='black')
# plt.plot(*r2[:2], c='red')
# plt.plot(*r3[:2], c='green')
plt.show()

# print(np.array(sol.u).shape)
# print(r1.shape, v1.shape)
#
# start = time.time()
# sys.path.append("/home/jass/phd/projects/3body/")
# import TBODYplot as TBODYplot
#
# # time interval and the step size
# t = np.arange(0, t_f, anispeed)

# render = 60
# cu_fade = TBODYplot.plot2Dfade_solution(
#     render,
#     m1,
#     m2,
#     m3,
#     *r1,
#     t,
#     *r2,
#     *r3,
#     *r3,
#     *r3,
#     *r3,
#     *r3,
#     *r2,
#     *r1,
#     2,
# )
#
# cu_fade.save("/home/jass/phd/projects/3body/figure8_25PN.mp4", writer="ffmpeg", fps=30)
#
#
# end = time.time()
# elapsed = end - start
# print("Animation time =", elapsed)
# # sol = de.solve(prob, de.RK4, saveat=0.1)
# # three_body_sol = solve_ivp(
# #     dydx,
# #     [0, t_f],
# #     u0,
# #     t_eval=t,
# #     method="Radau",
#     max_step=1e-2,
#     args=(m1, m2, m3),
#     atol=1e-8,
#     rtol=1e-8,
# )
