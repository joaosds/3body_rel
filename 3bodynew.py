# TD:
# 1. Testar intervalo de integração até convergência para outros métodos, e se necessário
# adcionar uma integração simplética com termos de velocidade.
# 2. Add all initial conditions and constants in a separate file.
# 3. Think the most efficient way to document the code. I don't think the parameters/ returns
# is necessary

import numpy as np
from diffeqpy import de
import sys

# import matplotlib.pyplot as plt
# import os
import time

# For first time using diffeqpy
# diffeqpy.install()

# Global variables
c = 1  # speed of light
anispeed = 0.08  # animation step
t_f = 108  # final time for evolution of the system
tspan = (0.0, t_f)
m = np.array([1, 1, 1])  # masses of three bodies, global constants
m1, m2, m3 = m


# ------------------ Auxiliary Functions ------------------------


def init_cond(case: str) -> np.ndarray:
    """
    Predefined initial conditions (IC) for a three-body configuration. Options for case: fig8, triquette...
    # add other initial configs.

    Parameters
    ----------
    Pos_0, Vel_0: np.ndarray
        Position/Velocity vectors for 3 bodies shape (3, 3).

    Returns
    -------
    ic : np.ndarray
        Vector containing the IC with shape (6, 3).
    """
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

    ic = np.concatenate((Pos_0, Vel_0))
    return ic


def dist(r1: np.ndarray, r2: np.ndarray, eps=0) -> float:
    """
    Calculate distance |r1-r2|. A small numerical cutoff (eps)
    may help avoiding masses going to infinity.
    """
    return np.sqrt(np.sum((r1 - r2) ** 2) + eps)


def nvec(r1: np.ndarray, r2: np.ndarray) -> np.ndarray:
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
        Versos n_{12} with shape (3)

    """
    versor = (r1 - r2) / dist(r1, r2)
    return versor


# ------------------ Force terms  ------------------------


def PN0(x: np.ndarray) -> np.ndarray:
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
    x1, x2, x3 = x
    f12 = m2 * nvec(x2, x1) / dist(x1, x2) ** 2
    f13 = m3 * nvec(x3, x1) / dist(x1, x3) ** 2
    f_pn0 = f12 + f13

    return f_pn0


def PN1(x: np.ndarray, v: np.ndarray) -> np.ndarray:
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
    x1, x2, x3 = x
    v1, v2, v3 = v
    f_pn1 = (
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
    return f_pn1


def PN2(x: np.ndarray, v: np.ndarray) -> np.ndarray:
    """
    Second order relativistic correction (order c^{-4}) for newtonian gravitational force.

    Parameters
    ----------
    x, v: np.ndarray
        Vectors positions and velocities. Shape (3).

    Returns
    -------
    f_pn2 : np.ndarray
        PN2 correction. Shape (3).

    """
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

    f_pn2 = term1 + term2
    return f_pn2


def PN25(x: np.ndarray, v: np.ndarray) -> np.ndarray:
    """
    Second order .5 relativistic correction (order c^{-5}) for newtonian gravitational force.

    Parameters
    ----------
    x, v: np.ndarray
        Vectors positions and velocities. Shape (3).

    Returns
    -------
    f_pn2p5 : np.ndarray
        PN2.5 correction. Shape (3).

    """
    x1, x2, x3 = x
    v1, v2, v3 = v

    f_n2p5 = (
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
    return f_n2p5


# ------------------ Solvers  ------------------------


def de_newton(u, p, t):
    """
    Three-body classical trajectories. Newtonian solution of differential
    equations (DE).

    Parameters
    ----------
    u: np.ndarray
        Vector containing the right part of the DEs.
        Here, dv/dt is indexed in (u[1:3, :]) and dx/dt in (u[3:6, :]),
        such that u has a shape of (6,3). The 2nd index refers to each body.

    Returns
    -------
    du_dt: np.ndarray
        Array with all DEs.
    """

    du_dt = np.zeros((6, 3))

    du_dt[0] = u[3]
    du_dt[1] = u[4]
    du_dt[2] = u[5]

    du_dt[3] = PN0(u[:3])
    du_dt[4] = PN0(u[[1, 2, 0]])
    du_dt[5] = PN0(u[[2, 0, 1]])

    return du_dt

runtime = time.time()

# In de.odeProblem we call the julia function with first argument as the DE,
# second with, and number of solutions (tspan)
# TD create a small functiopn for these definitions (in order to use to other corrections)

prob = de.ODEProblem(de_newton, init_cond("fig8"), tspan)
sol = de.solve(prob, de.Vern9(), saveat=anispeed, abstol=1e-9, reltol=1e-9)
# In sol, the second argument is the method
print(time.time() - runtime)

# TD: think on how to save this in a smarter way
sol = np.array(sol.u)
r1 = sol[:, 0, :].T
r2 = sol[:, 1, :].T
r3 = sol[:, 2, :].T
v1 = sol[:, 3, :].T
v2 = sol[:, 4, :].T
v3 = sol[:, 5, :].T

# ------------------ Animation  ------------------------

# TD: make path general to current folder
start = time.time()
sys.path.append("/home/jass/phd/projects/3body/")
import TBODYplot as TBODYplot

# time interval and the step size
t = np.arange(0, t_f, anispeed)
render = 60

ani = TBODYplot.plot2Dfade_solution(
    render,
    m1,
    m2,
    m3,
    *r1,
    t,
    *r2,
    *r3,
    *r3,
    *r3,
    *r3,
    *r3,
    *r2,
    *r1,
    2,
)

# TD: Add option to already save as gif
ani.save("/home/jass/phd/projects/3body/figure8_25PN.mp4", writer="ffmpeg", fps=30)


end = time.time()
elapsed = end - start
