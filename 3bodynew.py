import numpy as np

c = 1


def dist(r1: np.ndarray, r2: np.ndarray) -> float:
    # eps = numerical cutoff of 1e-4 to avoid masses going to infinity
    eps = 1e-2
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


def PN0(m: np.ndarray, x: np.ndarray) -> np.ndarray:
    """
    Newtonian gravitational force.
    """
    _, m2, m3 = m
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
