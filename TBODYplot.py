#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 12:02:19 2022

@author: feradofogo
"""

# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# Created on Fri Nov 18 17:38:28 2022

# @author: feradofogo
# """

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
import matplotlib.colors as colors
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D


lw = 2

    
def plot_solution(fig,m1,m2,m3,x_1, y_1, z_1, t, x_2, y_2, z_2, x_3, y_3, z_3, 
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


def plot2D_solution(fig,m1,m2,m3,x_1, y_1, z_1, t, x_2, y_2, z_2, x_3, y_3, z_3, 
         x_4, y_4, z_4, x_5, y_5, z_5, x_6, y_6, z_6, x_7, y_7, z_7, x_8, y_8, z_8,x_9, y_9, z_9):
    
    ax1 = fig.add_subplot(1, 3, 1)
    #ax1.figure.clear()
    ax1.grid(False)
    ax1.axis('off')
    line, = ax1.plot([], [], color='#A325D9', linewidth=lw, label='$m_1 = {}$'.format(m1))
    #line, = ax1.plot([], [], cmap = 'inferno', linewidth=lw, label='$m_1 = {}$'.format(m1))
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


def plot2Dfade_solution(render,m1,m2,m3,x_1, y_1, z_1, t, x_2, y_2, z_2, x_3, y_3, z_3, 
         x_4, y_4, z_4, x_5, y_5, z_5, x_6, y_6, z_6, x_7, y_7, z_7, x_8, y_8, z_8,x_9, y_9, z_9,box):
        
    # n_points_to_render = 500
    # # opacity of the segments
    n_points_to_render = render
    alphas = np.linspace(0, 1, n_points_to_render)
    
    # create "solid" color maps with a varying opacity
    red = colors.to_rgb("#A325D9") + (0.0,)
    redfade = colors.LinearSegmentedColormap.from_list('my', [red, '#A325D9'])
    green = colors.to_rgb('#04AD7B') + (0.0,)
    greenfade = colors.LinearSegmentedColormap.from_list('my', [green, '#04AD7B'])
    blue = colors.to_rgb('#D93425') + (0.0,)
    bluefade = colors.LinearSegmentedColormap.from_list('my', [blue, '#D93425'])
    
    
    def get_segments(i):
        # LineCollection requires segments
        _x1 = x_1[i:i+n_points_to_render]
        _y1 = y_1[i:i+n_points_to_render]
        _x2 = x_2[i:i+n_points_to_render]
        _y2 = y_2[i:i+n_points_to_render]
        _x3 = x_3[i:i+n_points_to_render]
        _y3 = y_3[i:i+n_points_to_render]
        points1 = np.vstack((_x1, _y1)).T.reshape(-1, 1, 2)
        segments1 = np.hstack((points1[:-1], points1[1:]))
        points2 = np.vstack((_x2, _y2)).T.reshape(-1, 1, 2)
        segments2 = np.hstack((points2[:-1], points2[1:]))
        points3 = np.vstack((_x3, _y3)).T.reshape(-1, 1, 2)
        segments3 = np.hstack((points3[:-1], points3[1:]))
        
        _x4 = x_4[i:i+n_points_to_render]
        _y4 = y_4[i:i+n_points_to_render]
        _x5 = x_5[i:i+n_points_to_render]
        _y5 = y_5[i:i+n_points_to_render]
        _x6 = x_6[i:i+n_points_to_render]
        _y6 = y_6[i:i+n_points_to_render]
        points4 = np.vstack((_x4, _y4)).T.reshape(-1, 1, 2)
        segments4 = np.hstack((points4[:-1], points4[1:]))
        points5 = np.vstack((_x5, _y5)).T.reshape(-1, 1, 2)
        segments5 = np.hstack((points5[:-1], points5[1:]))
        points6 = np.vstack((_x6, _y6)).T.reshape(-1, 1, 2)
        segments6 = np.hstack((points6[:-1], points6[1:]))
        
        _x7 = x_7[i:i+n_points_to_render]
        _y7 = y_7[i:i+n_points_to_render]
        _x8 = x_8[i:i+n_points_to_render]
        _y8 = y_8[i:i+n_points_to_render]
        _x9 = x_9[i:i+n_points_to_render]
        _y9 = y_9[i:i+n_points_to_render]
        points7 = np.vstack((_x7, _y7)).T.reshape(-1, 1, 2)
        segments7 = np.hstack((points7[:-1], points7[1:]))
        points8 = np.vstack((_x8, _y8)).T.reshape(-1, 1, 2)
        segments8 = np.hstack((points8[:-1], points8[1:]))
        points9 = np.vstack((_x9, _y9)).T.reshape(-1, 1, 2)
        segments9 = np.hstack((points9[:-1], points9[1:]))
        
        return segments1, segments2, segments3, segments4, segments5, segments6, segments7, segments8, segments9

        
    fig, ax = plt.subplots(nrows = 1, ncols = 3, figsize = (18,6), dpi = 300)
    fig.set_facecolor('#021226')
    ax[0].set_ylim(-box, box)
    ax[0].set_xlim(-box, box)
    ax[0].grid(False)
    ax[0].axis('off')
    
    ax[1].set_ylim(-box, box)
    ax[1].set_xlim(-box, box)
    ax[1].grid(False)
    ax[1].axis('off')
    
    ax[2].set_ylim(-box, box)
    ax[2].set_xlim(-box, box)
    ax[2].grid(False)
    ax[2].axis('off')
    
    ax[0].set_title('Newtonian', fontsize=14, color = 'white')
    ax[1].set_title('1.0 Post-Newtonian', fontsize=14, color = 'white')
    ax[2].set_title('2.5 Post-Newtonian', fontsize=14, color = 'white')
    
    # plt.style.use("ggplot")
    point, = ax[0].plot([], [], marker='o', color = '#A325D9',linewidth=lw, markersize=4)
    point2, = ax[0].plot([], [], marker='o', color = '#04AD7B',linewidth=lw, markersize=4)
    point3, = ax[0].plot([], [], marker='o', color = '#D93425',linewidth=lw, markersize=4)
    
    point4, = ax[1].plot([], [], marker='o', color = '#A325D9',linewidth=lw, markersize=4)
    point5, = ax[1].plot([], [], marker='o', color = '#04AD7B',linewidth=lw, markersize=4)
    point6, = ax[1].plot([], [], marker='o', color = '#D93425',linewidth=lw, markersize=4)
    
    point7, = ax[2].plot([], [], marker='o', color = '#A325D9',linewidth=lw, markersize=4)
    point8, = ax[2].plot([], [], marker='o', color = '#04AD7B',linewidth=lw, markersize=4)
    point9, = ax[2].plot([], [], marker='o', color = '#D93425',linewidth=lw, markersize=4)
    
    # # create and add two LineCollections
    segments1, segments2, segments3, segments4, segments5, segments6, segments7, segments8, segments9 = get_segments(0)
    lc1 = LineCollection(segments1, array=alphas, cmap=redfade, lw=2)
    lc2 = LineCollection(segments2, array=alphas, cmap=greenfade, lw=2)
    lc3 = LineCollection(segments3, array=alphas, cmap=bluefade, lw=2)
    lc4 = LineCollection(segments4, array=alphas, cmap=redfade, lw=2)
    lc5 = LineCollection(segments5, array=alphas, cmap=greenfade, lw=2)
    lc6 = LineCollection(segments6, array=alphas, cmap=bluefade, lw=2)
    lc7 = LineCollection(segments7, array=alphas, cmap=redfade, lw=2)
    lc8 = LineCollection(segments8, array=alphas, cmap=greenfade, lw=2)
    lc9 = LineCollection(segments9, array=alphas, cmap=bluefade, lw=2)
    line1 = ax[0].add_collection(lc1)
    line2 = ax[0].add_collection(lc2)
    line3 = ax[0].add_collection(lc3)
    line4 = ax[1].add_collection(lc4)
    line5 = ax[1].add_collection(lc5)
    line6 = ax[1].add_collection(lc6)
    line7 = ax[2].add_collection(lc7)
    line8 = ax[2].add_collection(lc8)
    line9 = ax[2].add_collection(lc9)
    
    point.set_data(x_1[0+n_points_to_render - 1], y_1[0+n_points_to_render - 1])
    point2.set_data(x_2[0+n_points_to_render - 1], y_2[0+n_points_to_render - 1])
    point3.set_data(x_3[0+n_points_to_render - 1], y_3[0+n_points_to_render -1])
    
    point4.set_data(x_4[0+n_points_to_render - 1], y_4[0+n_points_to_render - 1])
    point5.set_data(x_5[0+n_points_to_render - 1], y_5[0+n_points_to_render - 1])
    point6.set_data(x_6[0+n_points_to_render - 1], y_6[0+n_points_to_render -1])
    
    point7.set_data(x_7[0+n_points_to_render - 1], y_7[0+n_points_to_render - 1])
    point8.set_data(x_8[0+n_points_to_render - 1], y_8[0+n_points_to_render - 1])
    point9.set_data(x_9[0+n_points_to_render - 1], y_9[0+n_points_to_render -1])
    
    def animate(i):
        segments1, segments2, segments3, segments4, segments5, segments6, segments7, segments8, segments9 = get_segments(i)
        line1.set_segments(segments1)
        line2.set_segments(segments2)
        line3.set_segments(segments3)
        line4.set_segments(segments4)
        line5.set_segments(segments5)
        line6.set_segments(segments6)
        line7.set_segments(segments7)
        line8.set_segments(segments8)
        line9.set_segments(segments9)
        
        point.set_data(x_1[i+n_points_to_render - 1], y_1[i+n_points_to_render -1])
        point2.set_data(x_2[i+n_points_to_render - 1], y_2[i+n_points_to_render -1])
        point3.set_data(x_3[i+n_points_to_render -1], y_3[i+n_points_to_render - 1])
        
        point4.set_data(x_4[i+n_points_to_render - 1], y_4[i+n_points_to_render -1])
        point5.set_data(x_5[i+n_points_to_render - 1], y_5[i+n_points_to_render -1])
        point6.set_data(x_6[i+n_points_to_render -1], y_6[i+n_points_to_render - 1])
        
        point7.set_data(x_7[i+n_points_to_render - 1], y_7[i+n_points_to_render -1])
        point8.set_data(x_8[i+n_points_to_render - 1], y_8[i+n_points_to_render -1])
        point9.set_data(x_9[i+n_points_to_render -1], y_9[i+n_points_to_render - 1])
        
    
    # create a legend as LineCollection doesn't have any by default
    l1_legend = Line2D([1, 0], [1, 0], color='#A325D9', linewidth=2)
    l2_legend = Line2D([1, 0], [1, 0], color='#04AD7B', linewidth=2)
    l3_legend = Line2D([1, 0], [1, 0], color='#D93425', linewidth=2)
    ax[0].legend([l1_legend, l2_legend, l3_legend], [r'$m_1 = {}$'.format(m1), r'$m_2 = {}$'.format(m2), r'$m_3 = {}$'.format(m3)],loc=(0.2, -0.08) ,fancybox=True, facecolor='white', edgecolor='black', frameon=True)
    #ax1.legend(loc=(0.2, -0.08), fancybox=True, facecolor='white', edgecolor='black', frameon=True)
    
    anim = FuncAnimation(fig, animate, frames=len(x_1) - n_points_to_render, interval=0)
    return anim

# if __name__ == '__main__':
#     pica_fade = plot2Dfade_solution(x_1, y_1, z_1, t, x_2, y_2, z_2, x_3, y_3, z_3, 
#              x_4, y_4, z_4, x_5, y_5, z_5, x_6, y_6, z_6, x_7, y_7, z_7, x_8, y_8, z_8,x_9, y_9, z_9)
    
#     pica_fade.save('/home/feradofogo/Relativistic 3b GW/testeFADING.mp4', writer = 'ffmpeg', fps = 30)
#     #anim.save('/home/feradofogo/Relativistic 3b GW/jesus.mp4', writer = 'ffmpeg', fps = 30)
#     #pica.save('pica.mp4', writer = 'ffmpeg', fps = 30)
