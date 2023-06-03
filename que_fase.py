#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 11:06:56 2023

@author: feradofogo
"""

# TD:
# 1. Testar intervalo de integração até convergência para outros métodos, e se necessário
# adcionar uma integração simplética com termos de velocidade.
# 2. Add all initial conditions and constants in a separate file.
# 3. Think the most efficient way to document the code. I don't think the parameters/ returns
# is necessary

import numpy as np
#from diffeqpy import de
import sys

import matplotlib.pyplot as plt
# import os
import time

np.random.seed(6699)



# For first time using diffeqpy
# diffeqpy.install()

# Global variables
c = 1  # speed of light
G = 1
anispeed = 0.008  # animation step
t_f = 85 # final time for evolution of the system
tspan = (0.0, t_f)

eps = 1e-4

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
    #global m 
    if case == "fig8":
        m = np.array([1, 1, 1])  # masses of three bodies, global constants
        m1, m2, m3 = m
        Masses = np.array([[m1],
              [m2],
              [m3]])
        P_0 = np.array(
            [[-0.97000436, 0.24308753, 0],
             [0, 0, 0],
             [0.97000436, -0.24308753, 0]]
        )

        V_0 = np.array(
            [
                [0.4662036850, 0.4323657300, 0],
                [-0.93240737, -0.86473146, 0],
                [0.4662036850, 0.4323657300, 0],
            ]
        )

    #ic = np.concatenate((P_0, V_0))
    
    

    if case == "planaru1":
        m = np.array([1, 1, 1])  # masses of three bodies, global constants
        m1, m2, m3 = m
        Masses = np.array([[m1],
                [m2],
                [m3]])
        
        P_0 = np.array([[-1, 0, 0.0], [1, 0, 0.0], [0.00001, 0, 0.0]])
        V_0 = np.array(
            [
                [0.3420307307, 0.1809369236, 0.0],
                [0.3420307307, 0.1809369236, 0.0],
                [-2 * 0.3420307307 * m1 / m3, -2 * 0.1809369236 * m2 / m3, 0.0],
            ]
        )
        
    if case == "planar_m3":
        m1 = 1 
        m2 = m1 
        m3 = 0.5
        
        m = np.array([m1, m2, m3])
        Masses = np.array([[m1],
               [m2],
               [m3]])
        
        r1=[-1,0,0.0]
        r2=[1,0,0.0]
        r3=[0.00001,0,0.0]
        
        v1=[0.3420307307, 0.1809369236, 0.0]
        v2=v1
        v3=[-2*v1[0]*m1/m3, -2*v1[1]*m2/m3, 0.0]
        
        P_0 = np.array([r1, r2, r3])
        V_0 = np.array([v1, v2, v3])    
    
    if case == "planar_m1_bowtie":
        m1 = 1
        m2 = m1
        m3 = 4
        
        m = np.array([m1, m2, m3])
        Masses = np.array([[m1],
               [m2],
               [m3]])
        #print("Masses dentro", Masses)
        r1 = [-1, 0, 0.0]
        r2 = [1, 0, 0.0]
        r3 = [0.00001, 0, 0.0]
        
        v1 = [0.991198122, 0.711947212, 0.0]
        v2 = v1
        v3 = [-2 * v1[0] * m1 / m3, -2 * v1[1] * m2 / m3, 0.0]
        
        P_0 = np.array([r1, r2, r3])
        V_0 = np.array([v1, v2, v3]) 
    
    if case == "planar_m2_butterfly":
        m1 = 1
        m2 = m1
        m3 = 2
        
        m = np.array([m1, m2, m3])
        Masses = np.array([[m1],
               [m2],
               [m3]])
        
        r1 = [-1, 0, 0.0]
        r2 = [1, 0, 0.0]
        r3 = [0.00001, 0, 0.0]
        
        v1 = [0.664910758, 0.832416786, 0.0]
        v2 = v1
        v3 = [-2 * v1[0] * m1 / m3, -2 * v1[1] * m2 / m3, 0.0]
        
        P_0 = np.array([r1, r2, r3])
        V_0 = np.array([v1, v2, v3])
    
    if case == "circle":
        m = np.array([1, 1, 1])  # masses of three bodies, global constants
        m1, m2, m3 = m
        Masses = np.array([[m1],
               [m2],
               [m3]])
        
        L = 5
        #Define initial position vectors
        r1=[-L/2,-((np.sqrt(3)*L/2) - L/np.sqrt(3)),0.0]
        r2=[0.0,(np.sqrt(3)*L/2)-((np.sqrt(3)*L/2) - L/np.sqrt(3)),0.0]
        r3=[L/2,-((np.sqrt(3)*L/2) -L/np.sqrt(3)),0.0]
        v1=[np.sqrt(G/L)*1/2, -np.sqrt(G/L)*np.sqrt(3)/2, 0.0]
        v2=[-np.sqrt(G/L), 0.0,0.0]
        v3=[np.sqrt(G/L)*1/2, np.sqrt(G/L)*np.sqrt(3)/2, 0.0]
        
        P_0 = np.array([r1, r2, r3])
        V_0 = np.array([v1, v2, v3])
    
    if case == "triquette":
        m1 = 2.5
        m2 = m1
        m3 = m1
        m = np.array([m1, m2, m3])
        L = 5
        # Define initial position vectors
        r1 = [-L / 2, -((np.sqrt(3) * L / 2) - L / np.sqrt(3)), 0.0]
        r2 = [0.0, (np.sqrt(3) * L / 2) - ((np.sqrt(3) * L / 2) - L / np.sqrt(3)), 0.0]
        r3 = [L / 2, -((np.sqrt(3) * L / 2) - L / np.sqrt(3)), 0.0]
    
        v1 = [np.sqrt(G / L) * 1 / 2, -np.sqrt(G / L) * np.sqrt(3) / 2, 0.0]
        v2 = [-np.sqrt(G / L), 0.0, 0.0]
        v3 = [np.sqrt(G / L) * 1 / 2, np.sqrt(G / L) * np.sqrt(3) / 2, 0.0]
    
        P_0 = np.array([r1, r2, r3])
        V_0 = np.array([v1, v2, v3])
    
    if case == "satanico":
        m1 = 1
        m2 = m1
        m3 = 0.5
        m = np.array([m1,m2,m3])
        Masses = np.array([[m1],
                [m2],
                [m3]])
        
        P_0 = np.array([[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        v1 = 0.2009656237
        v2 = 0.2431076328
        V_0 = np.array(
            [
                [v1, v2 , 0.0],
                [v1, v2 , 0.0],
                [-2 * v1 * m1 / m3, -2 * v2 * m2 / m3, 0.0],
            ]
        )
    
    if case == "random":
        m = np.array([1, 1, 1])  # masses of three bodies, global constants
        m1, m2, m3 = m
        Masses = np.array([[m1],
               [m2],
               [m3]])
        P_0 = 1.3*np.random.randn(3,3)

        V_0 = 0.5*np.random.randn(3,3)

    ic = np.concatenate((P_0, V_0))
    return ic, m

icond, m = init_cond("random")






#m = np.array([1, 1, 1])  # masses of three bodies, global constants
m1, m2, m3 = m
Masses = np.array([[m1],
        [m2],
        [m3]])
print("Masses global", Masses)
# ------------------ Auxiliary Functions ------------------------




def dist(r1: np.ndarray, r2: np.ndarray, eps=eps) -> float:
    """
    Calculate distance |r1-r2|. A small numerical cutoff (eps)
    may help avoiding masses going to infinity.
    """
    return np.sqrt(np.sum((r1 - r2) ** 2) + eps)

def dist2(r1: np.ndarray, r2: np.ndarray, eps= eps ) -> float:
    """
    Calculate distance |r1-r2|. A small numerical cutoff (eps)
    may help avoiding masses going to infinity.
    """
    return np.sqrt(np.sum((r1 - r2) ** 2, axis = 0) + eps)



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



t = np.arange(0, t_f, anispeed)

def Sapao(DE,u0,ti):
    Pos_0 = u0[:3]
    Vel_0 = u0[3:6]

    Pos = Pos_0
    Vel = Vel_0
    
    Accel = DE(u0,0,0)[3:6]
    
    ti = 0
    
    Vel -= np.mean(Masses*Vel,0)/np.mean(Masses)

    sol = []

    dt = np.diff(t)[0]
    t_end = t[-1]
    Nt = int(t_end/dt)
    
    for i in range(Nt):
        # half kick
        Vel += Accel*dt/(2.0) 
        
        # drift
        Pos += Vel*dt
        
        # Update accelerations
        #Accel = AccelCalc(Pos, Masses)
        Accel = DE(np.vstack((Pos,Vel)),0,i*dt)[3:6]
        
        # half kick
        Vel += Accel*dt/(2.0)

        # update time
        ti += dt
        
        sol.append(np.vstack((Pos,Vel)))
        

    sol = np.array(sol)
    
    return sol

start = time.time()
solt = Sapao(de_newton, icond, t).T
end = time.time()

print("Sapao time = ", (end - start))

r1 = solt[:, 0, :]
r2 = solt[:, 1, :]
r3 = solt[:, 2, :]
v1 = solt[:, 3, :]
v2 = solt[:, 4, :]
v3 = solt[:, 5, :]

# K1 = 0.5*m1*    

def Energy(pos,vel):
    
    
    """
    Calculates the total newtonian energy of the system E = T + V
    
    """
    # Sum of the kinetic energy of all particles
    T = 0.5*( np.sum(m.reshape((3,1))*(( vel**2).sum(axis = 0)), axis = 0)  )
    
    # Sum of the potential energy of all particles
    x1 = pos[0]
    x2 = pos[1]
    x3 = pos[2]
    V = -G*m1*m2/dist2(x1, x2) - G*m2*m3/dist2(x2, x3) - G*m3*m1/dist2(x3, x1)
    
    return T,V

def Energy_cu(pos,vel,p):
    
    
    """
    Calculates the total newtonian energy of the system E = T + V
    
    """
    # Sum of the kinetic energy of all particles
    T = 0.5*m[p]*( vel[:,p,:]**2).sum(axis = 0)
    
    # Sum of the potential energy of all particles
    V = -G*m[p]*( m[(p-1)]/dist2(pos[:,p,:], pos[:,(p-1),:]) 
                 + m[(p+1) % 3]/dist2(pos[:,p,:], pos[:,(p+1) % 3,:]) )
    
    return T,V
    
print("T = ", Energy(solt[:, 0:3 ,:],solt[:, 3:6 ,:])[0] )
print("T shape = ", Energy(solt[:, 0:3 ,:],solt[:, 3:6 ,:])[0].shape)
print("Initial T = ", Energy(solt[:, 0:3 ,:],solt[:, 3:6 ,:])[0][0])

print("V = ", Energy(solt[:, 0:3 ,:],solt[:, 3:6 ,:])[1] )
print("V shape = ", Energy(solt[:, 0:3 ,:],solt[:, 3:6 ,:])[1].shape)
print("Initial V = ", Energy(solt[:, 0:3 ,:],solt[:, 3:6 ,:])[1][0])

print("E = ", Energy(solt[:, 0:3 ,:],solt[:, 3:6 ,:])[0] + 
      Energy(solt[:, 0:3 ,:],solt[:, 3:6 ,:])[1])
print("E shape = ", (Energy(solt[:, 0:3 ,:],solt[:, 3:6 ,:])[0] 
                     + Energy(solt[:, 0:3 ,:],solt[:, 3:6 ,:])[1]).shape)
print("Initial E = ", Energy(solt[:, 0:3 ,:],solt[:, 3:6 ,:])[0][0]
      + Energy(solt[:, 0:3 ,:],solt[:, 3:6 ,:])[1][0]
      )

T = Energy(solt[:, 0:3 ,:],solt[:, 3:6 ,:])[0]
V = Energy(solt[:, 0:3 ,:],solt[:, 3:6 ,:])[1]

Total = T + V

plt.figure(figsize = (10,8))
plt.plot(t[:-1], T/abs(Total[0]),
         label = 'Kinetic Energy')

plt.plot(t[:-1], V/abs(Total[0]),
         label = 'Potentital Energy')


plt.plot(t[:-1], Total[:]/abs(Total[0]),
         label = 'Total Energy')

plt.plot(t[:-1], -np.ones_like(t[:-1]), linestyle = '--', color = 'black'
         , label = 'Initial Energy')
plt.legend()
plt.show()

plt.figure(figsize = (10,10))
plt.plot(r1[0,:], r1[1,:],label = '1')
plt.scatter(r1[0,0], r1[1,0],label = '1')
plt.plot(r2[0,:], r2[1,:],label = '2')
plt.scatter(r2[0,0], r2[1,0],label = '2')

plt.plot(r3[0,:], r3[1,:],label = '3')
plt.scatter(r3[0,0], r3[1,0],label = '3')
plt.legend()
plt.show()

plt.figure(figsize = (10,8))
plt.plot(t[:-1], r1[0,:])
plt.show()


plt.figure(figsize = (10,8))
plt.plot(t[:-1], Energy_cu(solt[:, 0:3 ,:], solt[:, 3:6 ,:], 0)[0], label = '1'     )
plt.plot(t[:-1], Energy_cu(solt[:, 0:3 ,:], solt[:, 3:6 ,:], 1)[0], label = '2'     )
plt.plot(t[:-1], Energy_cu(solt[:, 0:3 ,:], solt[:, 3:6 ,:], 2)[0], label = '3'     )
plt.title("Kinetic Energy")
plt.legend()
plt.show()


plt.figure(figsize = (10,8))
plt.plot(t[:-1], Energy_cu(solt[:, 0:3 ,:], solt[:, 3:6 ,:], 0)[1], label = '1'     )
plt.plot(t[:-1], Energy_cu(solt[:, 0:3 ,:], solt[:, 3:6 ,:], 1)[1], label = '2'     )
plt.plot(t[:-1], Energy_cu(solt[:, 0:3 ,:], solt[:, 3:6 ,:], 2)[1], label = '3'     )
plt.title("Potential Energy")
plt.legend()
plt.show()
# ------------------ Animation  ------------------------

# TD: make path general to current folder
#start = time.time()
#sys.path.append("/home/feradofogo/3body_rel")
#import TBODYplot as TBODYplot

#render = 60

#ani = TBODYplot.plot2Dfade_solution_perspectives_energy(
#    render,
#    m1,
#    m2,
#    m3,
#    *r1,
#    t,
#    *r2,
#    *r3,
#    r1[1],r1[2],r1[0],
#    r2[1],r2[2],r2[0],
#    r3[1],r3[2],r3[0],
#    r1[2],r1[0],r1[1],
#    r2[2],r2[0],r2[1],
#    r3[2],r3[0],r3[1],
#    10,T,V,Total
#)

# TD: Add option to already save as gif
#ani.save("/home/feradofogo/3body_rel_phase_space/figure8_Frog2.mp4", writer="ffmpeg", fps=30)


#end = time.time()
#print("Animation time: ", (end - start))
