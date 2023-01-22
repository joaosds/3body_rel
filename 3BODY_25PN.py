#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 01:42:08 2022

@author: feradofogo
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 17:36:39 2022

@author: feradofogo
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 19:34:05 2022

@author: feradofogo
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
import sys
import matplotlib
#matplotlib.use('Qt5Agg')
#from PyQt5 import QtCore, QtWidgets
#from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import scipy as sp
import scipy.integrate
from scipy.integrate import solve_ivp
#from PyQt5.QtCore import QFile, QTextStream
import time


G,c=1,1#6.67408e-11 #N-m2/kg2

m1 = 1
m2 = m1
m3 = m1

fig = plt.figure(figsize=(12,6), dpi=300)    
fig.set_facecolor('#021226')

## Figure 8 Initial Conditions ##
Pos_0 = np.array([ [-0.97000436, 0.24308753,0],
                  [0,0,0],
                  [0.97000436, -0.24308753,0]   ])



Vel_0 = np.array([ [ 0.4662036850, 0.4323657300,0],[-0.93240737, -0.86473146,0],
                                      [ 0.4662036850, 0.4323657300,0]])



r1,r2,r3 = Pos_0
v1,v2,v3 = Vel_0


## Circle Initial Conditions ##
# L = 10
# #Define initial position vectors
# r1=[-L/2,-((np.sqrt(3)*L/2) - L/np.sqrt(3)),0.0]
# r2=[0.0,(np.sqrt(3)*L/2)-((np.sqrt(3)*L/2) - L/np.sqrt(3)),0.0]
# r3=[L/2,-((np.sqrt(3)*L/2) -L/np.sqrt(3)),0.0]

# #Convert pos vectors to arrays
# r1=np.array(r1,dtype="float64")
# r2=np.array(r2,dtype="float64")
# r3=np.array(r3,dtype="float64")

# #Define initial velocities
# v1=[np.sqrt(G/L)*1/2, -np.sqrt(G/L)*np.sqrt(3)/2, 0.0]
# v2=[-np.sqrt(G/L), 0.0,0.0]
# v3=[np.sqrt(G/L)*1/2, np.sqrt(G/L)*np.sqrt(3)/2, 0.0]
# #Convert velocity vectors to arrays
# v1=np.array(v1,dtype="float64")
# v2=np.array(v2,dtype="float64")
# v3=np.array(v3,dtype="float64")

# Pos_0 = np.array([r1,r2,r3])
# Vel_0 = np.array([v1,v2,v3])

# m1 = m2 = m3 = 1

# #Define initial position vectors
# r1=[-1,0,0.0]
# r2=[1,0,0.0]
# r3=[0.00001,0,0.0]
# Pos_0=np.array([r1,r2,r3])
# #Convert pos vectors to arrays
# r1=np.array(r1,dtype="float64")
# r2=np.array(r2,dtype="float64")
# r3=np.array(r3,dtype="float64")
# #Define initial velocities
# v1=[0.3420307307, 0.1809369236, 0.0]
# v2=v1
# v3=[-2*v1[0]*m1/m3, -2*v1[1]*m2/m3, 0.0]
# Vel_0=np.array([v1,v2,v3])
# v1=np.array(v1,dtype="float64")
# v2=np.array(v2,dtype="float64")
# v3=np.array(v3,dtype="float64")

# Pos_0 = np.array([r1,r2,r3])
# Vel_0 = np.array([v1,v2,v3])

    
def plot_solution(fig, x_1, y_1, z_1, t, x_2, y_2, z_2, x_3, y_3, z_3, 
         x_4, y_4, z_4, x_5, y_5, z_5, x_6, y_6, z_6, x_7, y_7, z_7, x_8, y_8, z_8,x_9, y_9, z_9):
    
    ax1 = fig.add_subplot(1, 3, 1, projection='3d')
    #ax1.figure.clear()
    ax1.grid(False)
    line, = ax1.plot([], [], [], color='#A325D9', linewidth=0.5, label='$m_1 = {}$'.format(m1))
    line2, = ax1.plot([], [], [], color='#04AD7B', linewidth=0.5, label='$m_2 = {}$'.format(m2))
    line3, = ax1.plot([], [], [], color='#D93425', linewidth=0.5, label='$m_3 = {}$'.format(m3))
    
    point, = ax1.plot([], [], [], marker='o', color='#A325D9', markersize=4)
    point2, = ax1.plot([], [], [], marker='o', color='#04AD7B', markersize=4)
    point3, = ax1.plot([], [], [], marker='o', color='#D93425', markersize=4)
    
    ax1.set_facecolor('#021226')
    ax1.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.2))
    ax1.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.2))
    ax1.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.2))
    # Take out the axis
    #ax1.set_axis_off()
    
    # set point-of-view: specified by (altitude degrees, azimuth degree)
    ax1.view_init(elev=23, azim=-131)
    
    # set limits 
    ax1.set_xlim(-10, 10)
    ax1.set_ylim(-10, 10)
    ax1.set_zlim(-10, 10)
    
    # remove tick labels 
    ax1.set_yticklabels([])
    ax1.set_xticklabels([])
    ax1.set_zticklabels([])
        
    # legend
    ax1.legend(loc=(0.2, -0.08), fancybox=True, facecolor='white', edgecolor='black', frameon=True)
    
    # set title
    ax1.set_title('Newtonian', fontsize=14, color = 'white')
        
    # plot x_2 post newtonian

    ax2 = fig.add_subplot(1, 3, 2, projection='3d')
    #ax1.figure.clear()
    ax2.grid(False)

    line4, = ax2.plot([], [], color = '#A325D9', linewidth=0.8)
    line5, = ax2.plot([], [], color = '#04AD7B', linewidth=0.8)
    line6, = ax2.plot([], [], color = '#D93425', linewidth=0.8)
    
    point4, = ax2.plot([], [], marker='o', color = '#A325D9', markersize=4)
    point5, = ax2.plot([], [], marker='o', color = '#04AD7B', markersize=4)
    point6, = ax2.plot([], [], marker='o', color = '#D93425', markersize=4)

    ax2.set_facecolor('#021226')
    ax2.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.2))
    ax2.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.2))
    ax2.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.2))
    # Take out the axis
    #ax1.set_axis_off()
    
    # set point-of-view: specified by (altitude degrees, azimuth degree)
    ax2.view_init(elev=23, azim=-131)
    
    # set limits 
    ax2.set_xlim(-10, 10)
    ax2.set_ylim(-10, 10)
    ax2.set_zlim(-10, 10)
    
    # # set axis labels
    # ax1.set_xlabel('X', fontsize=10)
    # ax1.set_ylabel('Y', fontsize=10)
    # ax1.set_zlabel('Z', fontsize=10)
    
    # remove tick labels 
    ax2.set_yticklabels([])
    ax2.set_xticklabels([])
    ax2.set_zticklabels([])
        
    # legend
    #ax2.text(0.2, -0.08, s= "Feel free to rotate the plot!")
    
    # set title
    ax2.set_title('1.0 Post-Newtonian ', fontsize=14, color = 'white')
    
    ax3 = fig.add_subplot(1, 3, 3, projection='3d')
    #ax1.figure.clear()
    ax3.grid(False)

    line7, = ax3.plot([], [], color = '#A325D9', linewidth=0.8)
    line8, = ax3.plot([], [], color = '#04AD7B', linewidth=0.8)
    line9, = ax3.plot([], [], color = '#D93425', linewidth=0.8)
    
    point7, = ax3.plot([], [], marker='o', color = '#A325D9', markersize=4)
    point8, = ax3.plot([], [], marker='o', color = '#04AD7B', markersize=4)
    point9, = ax3.plot([], [], marker='o', color = '#D93425', markersize=4)

    ax3.set_facecolor('#021226')
    ax3.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.2))
    ax3.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.2))
    ax3.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.2))
    # Take out the axis
    #ax1.set_axis_off()
    
    # set point-of-view: specified by (altitude degrees, azimuth degree)
    ax3.view_init(elev=23, azim=-131)
    
    # set limits 
    ax3.set_xlim(-10, 10)
    ax3.set_ylim(-10, 10)
    ax3.set_zlim(-10, 10)
    
    # # set axis labels
    # ax1.set_xlabel('X', fontsize=10)
    # ax1.set_ylabel('Y', fontsize=10)
    # ax1.set_zlabel('Z', fontsize=10)
    
    # remove tick labels 
    ax3.set_yticklabels([])
    ax3.set_xticklabels([])
    ax3.set_zticklabels([])
        
    # legend
    #ax2.text(0.2, -0.08, s= "Feel free to rotate the plot!")
    
    # set title
    ax3.set_title('2.5 Post-Newtonian ', fontsize=14, color = 'white')


    def update(i):
        i = 1*i
        if i >= x_1.size:
            i = x_1.size-1
        xq_1 = x_1[0:i]
        yq_1 = y_1[0:i]
        zq_1 = z_1[0:i]

        xq_2 = x_2[0:i]
        yq_2 = y_2[0:i]
        zq_2 = z_2[0:i]
        
        xq_3 = x_3[0:i]
        yq_3 = y_3[0:i]
        zq_3 = z_3[0:i]

        xq_4 = x_4[0:i]
        yq_4 = y_4[0:i]
        zq_4 = z_4[0:i]

        xq_5 = x_5[0:i]
        yq_5 = y_5[0:i]
        zq_5 = z_5[0:i]
        
        xq_6 = x_6[0:i]
        yq_6 = y_6[0:i]
        zq_6 = z_6[0:i]
        
        xq_7 = x_7[0:i]
        yq_7 = y_7[0:i]
        zq_7 = z_7[0:i]

        xq_8 = x_8[0:i]
        yq_8 = y_8[0:i]
        zq_8 = z_8[0:i]
        
        xq_9 = x_9[0:i]
        yq_9 = y_9[0:i]
        zq_9 = z_9[0:i]
        
        tq = t[0:i]
        
        line.set_data(xq_1, yq_1)
        line.set_3d_properties(zq_1)
        point.set_data(x_1[i], y_1[i])
        point.set_3d_properties(z_1[i])
        
        line2.set_data(xq_2, yq_2)
        line2.set_3d_properties(zq_2)
        point2.set_data(x_2[i], y_2[i])
        point2.set_3d_properties(z_2[i])
        
        line3.set_data(xq_3, yq_3)
        line3.set_3d_properties(zq_3)
        point3.set_data(x_3[i], y_3[i])
        point3.set_3d_properties(z_3[i])

        line4.set_data(xq_4, yq_4)
        line4.set_3d_properties(zq_4)
        point4.set_data(x_4[i], y_4[i])
        point4.set_3d_properties(z_4[i])
        
        line5.set_data(xq_5, yq_5)
        line5.set_3d_properties(zq_5)
        point5.set_data(x_5[i], y_5[i])
        point5.set_3d_properties(z_5[i])
        
        line6.set_data(xq_6, yq_6)
        line6.set_3d_properties(zq_6)
        point6.set_data(x_6[i], y_6[i])
        point6.set_3d_properties(z_6[i])
        
        line7.set_data(xq_7, yq_7)
        line7.set_3d_properties(zq_7)
        point7.set_data(x_7[i], y_7[i])
        point7.set_3d_properties(z_7[i])
        
        line8.set_data(xq_8, yq_8)
        line8.set_3d_properties(zq_8)
        point8.set_data(x_8[i], y_8[i])
        point8.set_3d_properties(z_8[i])
        
        line9.set_data(xq_9, yq_9)
        line9.set_3d_properties(zq_9)
        point9.set_data(x_9[i], y_9[i])
        point9.set_3d_properties(z_9[i])
        
        return line, line2, line3, line4, line5, line6, line7, line8, line9, point, point2, point3, point4, point5, point6, point7, point8, point9
    
    # instantiate the animator
    global ani
    ani = FuncAnimation(fig, update, frames=np.size(x_1), interval=0, blit=True)
    return ani








lw = 2
box = 2

def plot2D_solution(fig, x_1, y_1, z_1, t, x_2, y_2, z_2, x_3, y_3, z_3, 
         x_4, y_4, z_4, x_5, y_5, z_5, x_6, y_6, z_6, x_7, y_7, z_7, x_8, y_8, z_8,x_9, y_9, z_9):
    
    ax1 = fig.add_subplot(1, 3, 1)
    #ax1.figure.clear()
    ax1.grid(False)
    ax1.axis('off')
    line, = ax1.plot([], [], color='#A325D9', linewidth=lw, label='$m_1 = {}$'.format(m1))
    line2, = ax1.plot([], [], color='#04AD7B', linewidth=lw, label='$m_2 = {}$'.format(m2))
    line3, = ax1.plot([], [], color='#D93425', linewidth=lw, label='$m_3 = {}$'.format(m3))
    
    point, = ax1.plot([], [], marker='o', linewidth=lw, color='#A325D9', markersize=4)
    point2, = ax1.plot([], [], marker='o', linewidth=lw, color='#04AD7B', markersize=4)
    point3, = ax1.plot([], [], marker='o', linewidth=lw, color='#D93425', markersize=4)
    
    ax1.set_facecolor('#021226')
    #ax1.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.2))
    #ax1.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.2))
    #ax1.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.2))
    # Take out the axis
    #ax1.set_axis_off()
    
    # set point-of-view: specified by (altitude degrees, azimuth degree)
    #ax1.view_init(elev=23, azim=-131)
    
    # set limits 
    ax1.set_xlim(-box, box)
    ax1.set_ylim(-box, box)
    #ax1.set_zlim(-10, 10)
    
    # remove tick labels 
    ax1.set_yticklabels([])
    ax1.set_xticklabels([])
    #ax1.set_zticklabels([])
        
    # legend
    ax1.legend(loc=(0.2, -0.08), fancybox=True, facecolor='white', edgecolor='black', frameon=True)
    
    # set title
    ax1.set_title('Newtonian', fontsize=14, color = 'white')
        
    # plot x_2 post newtonian

    ax2 = fig.add_subplot(1, 3, 2)
    #ax1.figure.clear()
    ax2.grid(False)
    ax2.axis('off')

    line4, = ax2.plot([], [], color = '#A325D9', linewidth=lw)
    line5, = ax2.plot([], [], color = '#04AD7B', linewidth=lw)
    line6, = ax2.plot([], [], color = '#D93425', linewidth=lw)
    
    point4, = ax2.plot([], [], marker='o', color = '#A325D9',linewidth=lw, markersize=4)
    point5, = ax2.plot([], [], marker='o', color = '#04AD7B',linewidth=lw, markersize=4)
    point6, = ax2.plot([], [], marker='o', color = '#D93425',linewidth=lw, markersize=4)

    ax2.set_facecolor('#021226')
    #ax2.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.2))
    #ax2.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.2))
    #ax2.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.2))
    # Take out the axis
    #ax1.set_axis_off()
    
    # set point-of-view: specified by (altitude degrees, azimuth degree)
    #ax2.view_init(elev=23, azim=-131)
    
    # set limits 
    ax2.set_xlim(-box, box)
    ax2.set_ylim(-box, box)
    #ax2.set_zlim(-10, 10)
    
    # # set axis labels
    # ax1.set_xlabel('X', fontsize=10)
    # ax1.set_ylabel('Y', fontsize=10)
    # ax1.set_zlabel('Z', fontsize=10)
    
    # remove tick labels 
    ax2.set_yticklabels([])
    ax2.set_xticklabels([])
    #ax2.set_zticklabels([])
        
    # legend
    #ax2.text(0.2, -0.08, s= "Feel free to rotate the plot!")
    
    # set title
    ax2.set_title('1.0 Post-Newtonian ', fontsize=14, color = 'white')
    
    ax3 = fig.add_subplot(1, 3, 3)
    #ax1.figure.clear()
    ax3.grid(False)
    ax3.axis('off')

    line7, = ax3.plot([], [], color = '#A325D9', linewidth=lw)
    line8, = ax3.plot([], [], color = '#04AD7B', linewidth=lw)
    line9, = ax3.plot([], [], color = '#D93425', linewidth=lw)
    
    point7, = ax3.plot([], [], marker='o', color = '#A325D9',linewidth=lw, markersize=4)
    point8, = ax3.plot([], [], marker='o', color = '#04AD7B',linewidth=lw, markersize=4)
    point9, = ax3.plot([], [], marker='o', color = '#D93425',linewidth=lw, markersize=4)

    ax3.set_facecolor('#021226')
    #ax3.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.2))
    #ax3.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.2))
    #ax3.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.2))
    # Take out the axis
    #ax1.set_axis_off()
    
    # set point-of-view: specified by (altitude degrees, azimuth degree)
    #ax3.view_init(elev=23, azim=-131)
    
    # set limits 
    ax3.set_xlim(-box, box)
    ax3.set_ylim(-box, box)
    #ax3.set_zlim(-10, 10)
    
    # # set axis labels
    # ax1.set_xlabel('X', fontsize=10)
    # ax1.set_ylabel('Y', fontsize=10)
    # ax1.set_zlabel('Z', fontsize=10)
    
    # remove tick labels 
    ax3.set_yticklabels([])
    ax3.set_xticklabels([])
    #ax3.set_zticklabels([])
        
    # legend
    #ax2.text(0.2, -0.08, s= "Feel free to rotate the plot!")
    
    # set title
    ax3.set_title('2.5 Post-Newtonian ', fontsize=14, color = 'white')


    def update(i):
        i = 1*i
        if i >= x_1.size:
            i = x_1.size-1
        xq_1 = x_1[0:i]
        yq_1 = y_1[0:i]
        zq_1 = z_1[0:i]

        xq_2 = x_2[0:i]
        yq_2 = y_2[0:i]
        zq_2 = z_2[0:i]
        
        xq_3 = x_3[0:i]
        yq_3 = y_3[0:i]
        zq_3 = z_3[0:i]

        xq_4 = x_4[0:i]
        yq_4 = y_4[0:i]
        zq_4 = z_4[0:i]

        xq_5 = x_5[0:i]
        yq_5 = y_5[0:i]
        zq_5 = z_5[0:i]
        
        xq_6 = x_6[0:i]
        yq_6 = y_6[0:i]
        zq_6 = z_6[0:i]
        
        xq_7 = x_7[0:i]
        yq_7 = y_7[0:i]
        zq_7 = z_7[0:i]

        xq_8 = x_8[0:i]
        yq_8 = y_8[0:i]
        zq_8 = z_8[0:i]
        
        xq_9 = x_9[0:i]
        yq_9 = y_9[0:i]
        zq_9 = z_9[0:i]
        
        tq = t[0:i]
        
        line.set_data(xq_1, yq_1)
        #line.set_3d_properties(zq_1)
        point.set_data(x_1[i], y_1[i])
        #point.set_3d_properties(z_1[i])
        
        line2.set_data(xq_2, yq_2)
        #line2.set_3d_properties(zq_2)
        point2.set_data(x_2[i], y_2[i])
        #point2.set_3d_properties(z_2[i])
        
        line3.set_data(xq_3, yq_3)
        #line3.set_3d_properties(zq_3)
        point3.set_data(x_3[i], y_3[i])
        #point3.set_3d_properties(z_3[i])

        line4.set_data(xq_4, yq_4)
        #line4.set_3d_properties(zq_4)
        point4.set_data(x_4[i], y_4[i])
        #point4.set_3d_properties(z_4[i])
        
        line5.set_data(xq_5, yq_5)
        #line5.set_3d_properties(zq_5)
        point5.set_data(x_5[i], y_5[i])
        #point5.set_3d_properties(z_5[i])
        
        line6.set_data(xq_6, yq_6)
        #line6.set_3d_properties(zq_6)
        point6.set_data(x_6[i], y_6[i])
        #point6.set_3d_properties(z_6[i])
        
        line7.set_data(xq_7, yq_7)
        #line7.set_3d_properties(zq_7)
        point7.set_data(x_7[i], y_7[i])
        #point7.set_3d_properties(z_7[i])
        
        line8.set_data(xq_8, yq_8)
        #line8.set_3d_properties(zq_8)
        point8.set_data(x_8[i], y_8[i])
        #point8.set_3d_properties(z_8[i])
        
        line9.set_data(xq_9, yq_9)
        #line9.set_3d_properties(zq_9)
        point9.set_data(x_9[i], y_9[i])
        #point9.set_3d_properties(z_9[i])
        
        return line, line2, line3, line4, line5, line6, line7, line8, line9, point, point2, point3, point4, point5, point6, point7, point8, point9
    
    # instantiate the animator
    global ani
    ani = FuncAnimation(fig, update, frames=np.size(x_1), interval=0, blit=True)
    return ani
























init_params=np.array([r1,r2,r3,v1,v2,v3]) #Initial parameters
init_params=init_params.flatten() #Flatten to make 1D array

eps12 = 1e-8
eps13 = 1e-8
eps23 = 1e-8
anispeed = 0.08 #0.0002
#box_size = 30

t_f = 108

# time interval and the step size
t = np.arange(0, t_f, anispeed)

# vectors for the solutions
x_1 = np.zeros((len(t)))
y_1 = np.zeros((len(t)))
z_1 = np.zeros((len(t)))

x_2 = np.zeros((len(t)))
y_2 = np.zeros((len(t)))
z_2 = np.zeros((len(t)))


#fig = plt.figure(figsize=(10, 10), dpi=100)    
#fig.set_facecolor('#021226')
# ------------------------ Post Newtonian Solution -------------

def dist(x1, x2):
    
    d = np.sqrt(np.sum((x1- x2)**2)+1e-2) # numerical cutoff of 1e-4 to avoid masses going to infinity
    
    return d

def CenterDot(a,b):
    return np.dot(a,b)

def nvec(a,b):
    return (a-b)/dist(a,b)
    
def forcePN25(m1, m2, m3, x1, x2, x3, v1, v2, v3):
    
    term = (4*m1*m2*((v1 - v2)*((2*m1)/dist(x1, x2) - (8*m2)/dist(x1, x2) - (v1 - v2)**2) + CenterDot(nvec(x1,x2),(v1 - v2))*nvec(x1,x2)*((-6*m1)/dist(x1, x2) + (52*m2)/(3.*dist(x1, x2)) + 3*(v1 - v2)**2)))/(5.*c**5*dist(x1, x2)**3) + (4*m1*m3*((v1 - v3)*((2*m1)/dist(x1, x3) - (8*m3)/dist(x1, x3) - (v1 - v3)**2) + CenterDot(nvec(x1,x3),(v1 - v3))*nvec(x1,x3)*((-6*m1)/dist(x1, x3) + (52*m3)/(3.*dist(x1, x3)) + 3*(v1 - v3)**2)))/(5.*c**5*dist(x1, x3)**3)
    return term


def forcePN2_1(m1, m2, m3, x1, x2, x3, v1, v2, v3):
    
    term = (m2*nvec(x1,x2)*(-2*CenterDot(v1,v2)**2 - 6*CenterDot(v1,v2)*CenterDot(nvec(x1,x2),v2)**2 - (15*CenterDot(nvec(x1,x2),v2)**4)/8. + (3*CenterDot(nvec(x1,x2),v2)**2*v1**2)/2. + 4*CenterDot(v1,v2)*v2**2 + (9*CenterDot(nvec(x1,x2),v2)**2*v2**2)/2. - 2*v2**4 - (57*m1**2)/(4.*dist(x1,x2)**2) - (69*m1*m2)/(2.*dist(x1,x2)**2) - (9*m2**2)/dist(x1,x2)**2 + (m1*((-5*CenterDot(v1,v2))/2. + (39*CenterDot(nvec(x1,x2),v1)**2)/2. - 39*CenterDot(nvec(x1,x2),v1)*CenterDot(nvec(x1,x2),v2) + (17*CenterDot(nvec(x1,x2),v2)**2)/2. - (15*v1**2)/4. + (5*v2**2)/4.))/dist(x1,x2) + (m2*(-8*CenterDot(v1,v2) + 2*CenterDot(nvec(x1,x2),v1)**2 - 4*CenterDot(nvec(x1,x2),v1)*CenterDot(nvec(x1,x2),v2) - 6*CenterDot(nvec(x1,x2),v2)**2 + 4*v2**2))/dist(x1,x2)))/(c**4*dist(x1,x2)**2) + (m3*nvec(x1,x3)*(-2*CenterDot(v1,v3)**2 - 6*CenterDot(v1,v3)*CenterDot(nvec(x1,x3),v3)**2 - (15*CenterDot(nvec(x1,x3),v3)**4)/8. + (3*CenterDot(nvec(x1,x3),v3)**2*v1**2)/2. + 4*CenterDot(v1,v3)*v3**2 + (9*CenterDot(nvec(x1,x3),v3)**2*v3**2)/2. - 2*v3**4 - (57*m1**2)/(4.*dist(x1,x3)**2) - (69*m1*m3)/(2.*dist(x1,x3)**2) - (9*m3**2)/dist(x1,x3)**2 + (m1*((-5*CenterDot(v1,v3))/2. + (39*CenterDot(nvec(x1,x3),v1)**2)/2. - 39*CenterDot(nvec(x1,x3),v1)*CenterDot(nvec(x1,x3),v3) + (17*CenterDot(nvec(x1,x3),v3)**2)/2. - (15*v1**2)/4. + (5*v3**2)/4.))/dist(x1,x3) + (m3*(-8*CenterDot(v1,v3) + 2*CenterDot(nvec(x1,x3),v1)**2 - 4*CenterDot(nvec(x1,x3),v1)*CenterDot(nvec(x1,x3),v3) - 6*CenterDot(nvec(x1,x3),v3)**2 + 4*v3**2))/dist(x1,x3)))/(c**4*dist(x1,x3)**2)
    return term

def forcePN2_2(m1, m2, m3, x1, x2, x3, v1, v2, v3):
    
    term = (m2*(-6*CenterDot(nvec(x1,x2),v1)*CenterDot(nvec(x1,x2),v2)**2 + (9*CenterDot(nvec(x1,x2),v2)**3)/2. - 4*CenterDot(v1,v2)*CenterDot(nvec(x1,x2),(v1 - v2)) + CenterDot(nvec(x1,x2),v2)*v1**2 + 4*CenterDot(nvec(x1,x2),v1)*v2**2 - 5*CenterDot(nvec(x1,x2),v2)*v2**2 + (((-63*CenterDot(nvec(x1,x2),v1))/4. + (55*CenterDot(nvec(x1,x2),v2))/4.)*m1)/dist(x1,x2) - (2*(CenterDot(nvec(x1,x2),v1) + CenterDot(nvec(x1,x2),v2))*m2)/dist(x1,x2))*(v1 - v2))/(c**4*dist(x1,x2)**2) + (m3*(-6*CenterDot(nvec(x1,x3),v1)*CenterDot(nvec(x1,x3),v3)**2 + (9*CenterDot(nvec(x1,x3),v3)**3)/2. - 4*CenterDot(v1,v3)*CenterDot(nvec(x1,x3),(v1 - v3)) + CenterDot(nvec(x1,x3),v3)*v1**2 + 4*CenterDot(nvec(x1,x3),v1)*v3**2 - 5*CenterDot(nvec(x1,x3),v3)*v3**2 + (((-63*CenterDot(nvec(x1,x3),v1))/4. + (55*CenterDot(nvec(x1,x3),v3))/4.)*m1)/dist(x1,x3) - (2*(CenterDot(nvec(x1,x3),v1) + CenterDot(nvec(x1,x3),v3))*m3)/dist(x1,x3))*(v1 - v3))/(c**4*dist(x1,x3)**2)    
    return term

def forcePN2(m1, m2, m3, x1, x2, x3, v1, v2, v3):
    return forcePN2_1(m1, m2, m3, x1, x2, x3, v1, v2, v3) + forcePN2_2(m1, m2, m3, x1, x2, x3, v1, v2, v3)

def force(m2, x1, x2):
    
    f0 = m2*(x2 - x1)/dist(x1, x2)**3
    
    return f0

def force2(m1, m2, m3, x1, x2, x3):
    
    dis = m2/dist(x1, x2)
    f2 = force(m1, x2,x1) + force(m3, x2,x3)
    
    return (7/2)*dis*f2

def force3(m2, x1, x2, v1, v2):
    
    f = m2/(dist(x1,x2)**2)*(np.dot((x1 - x2)/dist(x1,x2),(4*v1 - 3*v2)) * (v1-v2))
    
    return f

def force4(m1, m2, m3, x1, x2, x3, v1, v2):
    
    f = force(m2, x1, x2)*(np.dot(v1, v1) + 2*np.dot(v2,v2) - 4*np.dot(v1,v2) - 
                        3/2 * np.dot((x1-x2)/dist(x1,x2), v2)**2 + 
                        1/2 * np.dot(x2 - x1, force(m1, x2, x1) + force(m3, x2, x3)))
    return f

def force5(m2, x1, x2, x3):
    
    r12 = 1/dist(x1, x2)
    r23 = 1/dist(x2, x3)
    r13 = 1/dist(x1, x3)
    f1 = force(m2, x1,x2)*(r12 + r13) 
    f2 = force(m2, x1,x2)*(r12 + r23)
    
    return (-4*f1 -f2)

def dydx(t, u, m1, m2, m3): ## 2.5 PN
    
    """
    x1, y1, z1, vx1, vy1, vz1 = u[:6]
    x2, y2, z2, vx2, vy2, vz2 = u[6:12]
    x3, y3, z3, vx3, vy3, vz3 = u[12:]
    """
    
    x1 = u[:3]
    x2 = u[6:9]
    x3 = u[12:15]
    
    v1 = u[3:6]
    v2 = u[9:12]
    v3 = u[15:]
    
    f01 = force(m2, x1, x2) + force(m3, x1, x3)
    f11 = force2(m1, m2, m3, x1, x2, x3) + force2(m1, m3, m2, x1 ,x3, x2)
    f21 = force3(m2, x1, x2, v1, v2) + force3(m3, x1, x3, v1, v3)
    f31 = force4(m1, m2, m3, x1, x2, x3, v1, v2) + force4(m1, m3, m2, x1, x3, x2, v1, v3)
    f41 = force5(m2, x1,x2,x3) + force5(m3, x1,x3,x2)
    f2PN_1 = forcePN2(m1, m2, m3, x1, x2, x3, v1, v2, v3) 
    f25PN_1 = forcePN25(m1, m2, m3, x1, x2, x3, v1, v2, v3)
    
    f02 = force(m1, x2,x1) + force(m3, x2, x3)
    f12 = force2(m2, m1, m3, x2, x1, x3) + force2(m2, m3, m1, x2, x3, x1)
    f22 = force3(m1, x2, x1, v2, v1) + force3(m3, x2, x3, v2, v3)
    f32 = force4(m2, m1, m3, x2, x1, x3, v2, v1) + force4(m2, m3, m1, x2, x3, x1, v2, v3)
    f42 = force5(m1, x2, x1, x3) + force5(m3, x2, x3, x1)
    f2PN_2 = forcePN2(m2, m3, m1, x2, x3, x1, v2, v3, v1) 
    f25PN_2 = forcePN25(m2, m3, m1, x2, x3, x1, v2, v3, v1)

    
    f03 = force(m1, x3, x1) + force(m2, x3, x2)
    f13 = force2(m3, m1, m2, x3, x1, x2) + force2(m3, m2, m1, x3, x2, x1)
    f23 = force3(m1, x3, x1, v3, v1) + force3(m2, x3, x2, v3, v2)
    f33 = force4(m3, m2, m1, x3, x2, x1, v3, v2) + force4(m3, m1, m2, x3, x1, x2, v3, v1)
    f43 = force5(m1, x3,x1,x2) + force5(m2, x3,x2,x1)
    f2PN_3 = forcePN2(m3, m1, m2, x3, x1, x2, v3, v1, v2) 
    f25PN_3 = forcePN25(m3, m1, m2, x3, x1, x2, v3, v1, v2)
    
    tol2 = 1e-2 #1e-3
    tol3 = 1e-2
    tol4 = 1e-2
    tol4 = 1e-2
    tol2PN = 1e-4
    tol25PN = 5e-5
    
    return np.array([*v1,*(f01 + tol2*f11 + tol3*f21 + tol4*f31 + tol4*f41 + tol2PN*f2PN_1 + tol25PN*f25PN_1  ),  ## 1,2 + 1,3
                      *v2,*(f02 + tol2*f12 + tol3*f22 + tol4*f32 + tol4*f42 + tol2PN*f2PN_2 + tol25PN*f25PN_2 ),  ## 2,1 + 2,3
                      *v3,*(f03 + tol2*f13 + tol3*f23 + tol4*f33 + tol4*f43 + tol2PN*f2PN_3 + tol25PN*f25PN_3 )])

G = M = 1

def dwdx(t, u, m1, m2, m3): ## 1PN
    
    """
    x1, y1, z1, vx1, vy1, vz1 = u[:6]
    x2, y2, z2, vx2, vy2, vz2 = u[6:12]
    x3, y3, z3, vx3, vy3, vz3 = u[12:]
    """
    
    x1 = u[:3]
    x2 = u[6:9]
    x3 = u[12:15]
    
    v1 = u[3:6]
    v2 = u[9:12]
    v3 = u[15:]
    
    f01 = force(m2, x1, x2) + force(m3, x1, x3)
    f11 = force2(m1, m2, m3, x1, x2, x3) + force2(m1, m3, m2, x1 ,x3, x2)
    f21 = force3(m2, x1, x2, v1, v2) + force3(m3, x1, x3, v1, v3)
    f31 = force4(m1, m2, m3, x1, x2, x3, v1, v2) + force4(m1, m3, m2, x1, x3, x2, v1, v3)
    f41 = force5(m2, x1,x2,x3) + force5(m3, x1,x3,x2)
    f2PN_1 = forcePN2(m1, m2, m3, x1, x2, x3, v1, v2, v3) 
    f25PN_1 = forcePN25(m1, m2, m3, x1, x2, x3, v1, v2, v3)
    
    f02 = force(m1, x2,x1) + force(m3, x2, x3)
    f12 = force2(m2, m1, m3, x2, x1, x3) + force2(m2, m3, m1, x2, x3, x1)
    f22 = force3(m1, x2, x1, v2, v1) + force3(m3, x2, x3, v2, v3)
    f32 = force4(m2, m1, m3, x2, x1, x3, v2, v1) + force4(m2, m3, m1, x2, x3, x1, v2, v3)
    f42 = force5(m1, x2, x1, x3) + force5(m3, x2, x3, x1)
    f2PN_2 = forcePN2(m2, m3, m1, x2, x3, x1, v2, v3, v1) 
    f25PN_2 = forcePN25(m2, m3, m1, x2, x3, x1, v2, v3, v1)

    
    f03 = force(m1, x3, x1) + force(m2, x3, x2)
    f13 = force2(m3, m1, m2, x3, x1, x2) + force2(m3, m2, m1, x3, x2, x1)
    f23 = force3(m1, x3, x1, v3, v1) + force3(m2, x3, x2, v3, v2)
    f33 = force4(m3, m2, m1, x3, x2, x1, v3, v2) + force4(m3, m1, m2, x3, x1, x2, v3, v1)
    f43 = force5(m1, x3,x1,x2) + force5(m2, x3,x2,x1)
    f2PN_3 = forcePN2(m3, m1, m2, x3, x1, x2, v3, v1, v2) 
    f25PN_3 = forcePN25(m3, m1, m2, x3, x1, x2, v3, v1, v2)
    
    tol2 = 1e-2
    tol3 = 1e-2
    tol4 = 1e-2
    tol4 = 1e-2
    tol2PN = 0
    tol25PN = 0
    
    return np.array([*v1,*(f01 + tol2*f11 + tol3*f21 + tol4*f31 + tol4*f41 + tol2PN*f2PN_1 + tol25PN*f25PN_1  ),  ## 1,2 + 1,3
                      *v2,*(f02 + tol2*f12 + tol3*f22 + tol4*f32 + tol4*f42 + tol2PN*f2PN_2 + tol25PN*f25PN_2 ),  ## 2,1 + 2,3
                      *v3,*(f03 + tol2*f13 + tol3*f23 + tol4*f33 + tol4*f43 + tol2PN*f2PN_3 + tol25PN*f25PN_3 )])


def ThreeBodyEquations(t,w,G,m1,m2,m3):
    r1=w[:3]
    r2=w[3:6]
    r3=w[6:9]
    v1=w[9:12]
    v2=w[12:15]
    v3=w[15:18]
    
    r12=sp.linalg.norm(r2-r1)+eps12
    r13=sp.linalg.norm(r3-r1)+eps13
    r23=sp.linalg.norm(r3-r2)+eps23
    
    dv1bydt=m2*(r2-r1)/r12**3+m3*(r3-r1)/r13**3
    dv2bydt=m1*(r1-r2)/r12**3+m3*(r3-r2)/r23**3
    dv3bydt=m1*(r1-r3)/r13**3+m2*(r2-r3)/r23**3
    dr1bydt=v1
    dr2bydt=v2
    dr3bydt=v3
    r12_derivs=np.concatenate((dr1bydt,dr2bydt))
    r_derivs=np.concatenate((r12_derivs,dr3bydt))
    v12_derivs=np.concatenate((dv1bydt,dv2bydt))
    v_derivs=np.concatenate((v12_derivs,dv3bydt))
    derivs=np.concatenate((r_derivs,v_derivs))
    return derivs
    
start = time.time()
three_body_sol=solve_ivp(ThreeBodyEquations, [0, t_f], init_params, t_eval = t, method = 'Radau', max_step = 1e-2, args=(G,m1,m2,m3))

end = time.time()

elapsed = end - start
print("Integration time: 0 PN =", elapsed)

r1_sol=three_body_sol.y[:3]
r2_sol=three_body_sol.y[3:6]
r3_sol=three_body_sol.y[6:9]

x_1 = r1_sol[0]
y_1 = r1_sol[1]
z_1 = r1_sol[2]

x_2 = r2_sol[0]
y_2 = r2_sol[1]
z_2 = r2_sol[2]

x_3 = r3_sol[0]
y_3 = r3_sol[1]
z_3 = r3_sol[2]

x_4 = np.zeros((len(t)))
y_4 = np.zeros((len(t)))
z_4 = np.zeros((len(t)))

x_5 = np.zeros((len(t)))
y_5 = np.zeros((len(t)))
z_5 = np.zeros((len(t)))

x_6 = np.zeros((len(t)))
y_6 = np.zeros((len(t)))
z_6 = np.zeros((len(t)))

u0 = np.zeros(18)

u0[:3] = Pos_0[0,:]
u0[3:6] = Vel_0[0,:]
u0[6:9] = Pos_0[1,:]
u0[9:12] = Vel_0[1,:]
u0[12:15] = Pos_0[2,:]
u0[15:] = Vel_0[2, :]

start = time.time()
three_body_sol=solve_ivp(dwdx, [0, t_f], u0, t_eval = t, method = 'Radau', max_step = 1e-2, args=(m1,m2,m3))
end = time.time()

elapsed = end - start
print("Integration time: 1.0 PN = ", elapsed)

r4_sol=three_body_sol.y[:3]
r5_sol=three_body_sol.y[6:9]
r6_sol=three_body_sol.y[12:15]

x_4 = r4_sol[0]
y_4 = r4_sol[1]
z_4 = r4_sol[2]

x_5 = r5_sol[0]
y_5 = r5_sol[1]
z_5 = r5_sol[2]

x_6 = r6_sol[0]
y_6 = r6_sol[1]
z_6 = r6_sol[2]

x_7 = np.zeros((len(t)))
y_7 = np.zeros((len(t)))
z_7 = np.zeros((len(t)))

x_8 = np.zeros((len(t)))
y_8 = np.zeros((len(t)))
z_8 = np.zeros((len(t)))

x_9 = np.zeros((len(t)))
y_9 = np.zeros((len(t)))
z_9 = np.zeros((len(t)))
##############################################
u0 = np.zeros(18)

u0[:3] = Pos_0[0,:]
u0[3:6] = Vel_0[0,:]
u0[6:9] = Pos_0[1,:]
u0[9:12] = Vel_0[1,:]
u0[12:15] = Pos_0[2,:]
u0[15:] = Vel_0[2, :]

start = time.time()
three_body_sol=solve_ivp(dydx, [0, t_f], u0, t_eval = t, method = 'Radau', max_step = 1e-2, args=(m1,m2,m3),atol=1e-8, rtol=1e-8)
end = time.time()

elapsed = end - start
print("Integration time: 2.5 PN =", elapsed)



r7_sol=three_body_sol.y[:3]
r8_sol=three_body_sol.y[6:9]
r9_sol=three_body_sol.y[12:15]

x_7 = r7_sol[0]
y_7 = r7_sol[1]
z_7 = r7_sol[2]

x_8 = r8_sol[0]
y_8 = r8_sol[1]
z_8 = r8_sol[2]

x_9 = r9_sol[0]
y_9 = r9_sol[1]
z_9 = r9_sol[2]

##### Testing with saving file with 2.5 PN simulation data #####
import time
start = time.time()
np.savez('SAVEtest_triquete.dat', 
         x_7=x_7, y_7=y_7, z_7=z_7,
         x_8=x_8, y_8=y_8, z_8=z_8,
         x_9=x_9, y_9=y_9, z_9=z_9)
end = time.time()
elapsed = end - start
print("Saving time 2.5 PN = ",elapsed)

#################################################################

start = time.time()
# cu = plot2D_solution(fig, x_1, y_1, z_1, t, x_2, y_2, z_2, x_3, y_3, z_3, 
#  x_4, y_4, z_4, x_5, y_5, z_5, x_6, y_6, z_6,x_7, y_7, z_7, x_8, y_8, z_8, x_9, y_9, z_9)

# cu.save('cu.mp4', writer = 'ffmpeg', fps = 30)

# pica = plot2D_solution(fig, x_1, y_1, z_1, t, x_2, y_2, z_2, x_3, y_3, z_3, 
#   x_4, y_4, z_4, x_5, y_5, z_5, x_6, y_6, z_6,x_7, y_7, z_7, x_8, y_8, z_8, x_9, y_9, z_9)

# pica.save('pica.mp4', writer = 'ffmpeg', fps = 30)


import sys, os
sys.path.append('/home/feradofogo/Relativistic 3b GW/')

import TBODYplot as TBODYplot

render = 60
cu_fade = TBODYplot.plot2Dfade_solution(render,m1,m2,m3,x_1, y_1, z_1, t, x_2, y_2, z_2, x_3, y_3, z_3, x_4, y_4, z_4, x_5, y_5, z_5, x_6, y_6, z_6, x_7, y_7, z_7, x_8, y_8, z_8,x_9, y_9, z_9, 2)

cu_fade.save('/home/feradofogo/Relativistic 3b GW/figure8_25PN.mp4', writer = 'ffmpeg', fps = 30)


end = time.time()
elapsed = end - start
print("Animation time =", elapsed)