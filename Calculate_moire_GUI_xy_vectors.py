# -*- coding: utf-8 -*-
"""
Spyder Editor

Written by Evan Telford (ejt2133@columbia.edu)
"""
#%%
#load pertinent packages
import numpy as np
import sys
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib
import math
font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 3}
matplotlib.rc('font', **font)
matplotlib.use('Qt5Agg')
# import sys
from PyQt5 import QtCore
from PyQt5.QtWidgets import (
    QApplication,
    QLabel,
    QLineEdit,
    QMainWindow,
    QPushButton,
    QWidget,
    QGroupBox,
    QGridLayout
)
#%%
#creates a class to initialize the plotting canvas
class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=4, height=4, dpi=300):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        for axis in ['top','bottom','left','right']:
            self.axes.spines[axis].set_linewidth(0.5)
        self.axes.xaxis.set_tick_params(width=0.5)
        self.axes.yaxis.set_tick_params(width=0.5)
        super(MplCanvas, self).__init__(self.fig)
 
#%%
#main class for generating moire lattices
class App(QMainWindow):
    
    def __init__(self):
        super().__init__()
        self.title = 'Moire Calculator'
        self.left = 0
        self.top = 0
        self.width = 1500
        self.height = 500
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        
        self.dataplot1=MplCanvas(self,width=4,height=4,dpi=300)
        self.toolbar1=NavigationToolbar(self.dataplot1,self)
        self.dataplot2=MplCanvas(self,width=4,height=4,dpi=300)
        self.toolbar2=NavigationToolbar(self.dataplot2,self)
        self.dataplot3=MplCanvas(self,width=4,height=4,dpi=300)
        self.toolbar3=NavigationToolbar(self.dataplot3,self)
        
        self.layout=QGridLayout()
        self.box = QGroupBox()
        V1=QLabel("Enter the first vector set [x,y] (nm)")
        V2=QLabel("Enter the second vector set [x,y] (nm)")
        A=QLabel("Enter the angle (theta)")
        V1.setAlignment(QtCore.Qt.AlignCenter)
        V2.setAlignment(QtCore.Qt.AlignCenter)
        A.setAlignment(QtCore.Qt.AlignCenter)
        self.V1ax = QLineEdit()
        self.V1bx = QLineEdit()
        self.V2ax = QLineEdit()
        self.V2bx = QLineEdit()
        self.V1ay = QLineEdit()
        self.V1by = QLineEdit()
        self.V2ay = QLineEdit()
        self.V2by = QLineEdit()
        self.A = QLineEdit()
        l_button=QPushButton('Calculate Moire Vectors') 
        l_button.clicked.connect(lambda: self.calc_moire_vectors()) 
        l_button.clicked.connect(lambda: self.visualize_moire()) 
        self.MVa=QLabel('')    
        self.MVa.setAlignment(QtCore.Qt.AlignCenter)    
        self.MVb=QLabel('')    
        self.MVb.setAlignment(QtCore.Qt.AlignCenter)
        self.MVamag=QLabel('')    
        self.MVamag.setAlignment(QtCore.Qt.AlignCenter)    
        self.MVbmag=QLabel('')    
        self.MVbmag.setAlignment(QtCore.Qt.AlignCenter)
        
        self.layout.addWidget(V1, 0,0,1,4) 
        self.layout.addWidget(self.V1ax,1,0,1,1)
        self.layout.addWidget(self.V1ay,1,1,1,1)
        self.layout.addWidget(self.V1bx,1,2,1,1)
        self.layout.addWidget(self.V1by,1,3,1,1)
        self.layout.addWidget(V2,2,0,1,4)
        self.layout.addWidget(self.V2ax,3,0,1,1)
        self.layout.addWidget(self.V2ay,3,1,1,1)
        self.layout.addWidget(self.V2bx,3,2,1,1)
        self.layout.addWidget(self.V2by,3,3,1,1)
        self.layout.addWidget(A,4,0,1,4)
        self.layout.addWidget(self.A,5,0,1,4)
        self.layout.addWidget(l_button, 6,0,1,4)
        self.layout.addWidget(self.MVa,7,0,1,2)
        self.layout.addWidget(self.MVb,7,2,1,2)
        self.layout.addWidget(self.MVamag,8,0,1,2)
        self.layout.addWidget(self.MVbmag,8,2,1,2)
        
        self.layout.addWidget(self.toolbar1,0,4,1,4)
        self.layout.addWidget(self.toolbar2,0,8,1,4)
        self.layout.addWidget(self.toolbar3,0,12,1,4)
        
        self.layout.addWidget(self.dataplot1,1,4,7,4)
        self.layout.addWidget(self.dataplot2,1,8,7,4)
        self.layout.addWidget(self.dataplot3,1,12,7,4)
        
        self.box.setLayout(self.layout)
        
        self.widget=QWidget()
        self.widget.setLayout(self.layout)
        self.setCentralWidget(self.widget)
    
    #function to calculate the moire vectors
    def calc_moire_vectors(self):
        a=np.array([float(self.V1ax.text()),float(self.V1ay.text())])
        b=np.array([float(self.V1bx.text()),float(self.V1by.text())])
        c=np.array([float(self.V2ax.text()),float(self.V2ay.text())])
        d=np.array([float(self.V2bx.text()),float(self.V2by.text())])
        angle=float(self.A.text())*np.pi/180
        R=[[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]]
        
        rx,ry=self.moire_real_v3(a,b,np.dot(R,c),np.dot(R,d))
        self.MVa.setText('L1 = ['+str(np.round(rx[0],4))+' , '+str(np.round(rx[1],4))+'] (nm)')
        self.MVb.setText('L2 = ['+str(np.round(ry[0],4))+' , '+str(np.round(ry[1],4))+'] (nm)')
        self.MVamag.setText('|L1| = '+str(np.round(np.linalg.norm(rx),4))+' (nm)')
        self.MVbmag.setText('|L2| = '+str(np.round(np.linalg.norm(ry),4))+' (nm)')
        
    #caclulates reciprocal vector of two vectors
    def calculate_reciprocal(self,a,b):
        ax=a[0]
        ay=a[1]
        bx=b[0]
        by=b[1]
        a_vec=np.hstack([ax,ay])
        b_vec=np.hstack([bx,by])
        R=[[0,-1],[1,0]]
        b1=2*np.pi*(np.dot(R,b_vec))/np.dot(a,np.dot(R,b_vec))
        b2=2*np.pi*(np.dot(R,a_vec))/np.dot(b,np.dot(R,a_vec))
        
        return b1, b2

    #calculates the angle in degrees between two vectors
    def find_angle(self,a,b):
        angle=np.arccos(np.dot(a,b)/np.linalg.norm(a)/np.linalg.norm(b))
        demi=np.cross(a,b)
        if demi>=0:
            angle=angle
        else:
            angle=(np.pi-angle)+np.pi
        
        return angle*180/np.pi

    #determines the rotational symmetry between two vectors
    def find_symmetry(self,a,b):
        angle=self.find_angle(a,b)
        if angle==0:
            print("Angle between vectors is 0, can't determine the rotational symmetry!")
            return 0
        elif angle==360:
            print("Angle between vectors is 360, can't determine the rotational symmetry!")
            return 0
        else:
            symmetry=360/angle
            if np.abs(np.round(symmetry)-symmetry)<0.001: #determines if the symmetry is within 0.1% of the nearest integer.
                if np.abs(np.linalg.norm(a)-np.linalg.norm(b))<0.001: #determines if vectors lengths are within 0.1%
                    symmetry=np.round(symmetry)
                else: #accounts for breaking symmetry due to different vector lengths
                    symmetry=2
            else:
                symmetry=1
            
            return symmetry

    #correctly subtracts reciprocal vectors
    def sub_reciprocal_v3(self,a,b,c,d):
        kax,kay=self.calculate_reciprocal(a,b)
        kbx,kby=self.calculate_reciprocal(c,d)
        
        #find out what the lowest symmetry is between the two vectors.
        symmetry=self.find_symmetry(kax,kay)
        symmetry2=self.find_symmetry(kbx,kby)
        low_symm=np.min([symmetry,symmetry2])
        
        #determines the angle at which symmetry operations need to be performed and if necessary, performs the symmetry operation.
        theta_ref=self.find_angle(kax,kbx)
        rot_pass=360/low_symm/2
        if theta_ref>=rot_pass:
            offset=-rot_pass*2 - np.floor(((np.floor(theta_ref/rot_pass)-1)/2))*rot_pass*2
            R=[[np.cos(offset*np.pi/180),-np.sin(offset*np.pi/180)],[np.sin(offset*np.pi/180),np.cos(offset*np.pi/180)]]
            kbx=np.dot(R,kbx)
            kby=np.dot(R,kby)
        kx_diff=kax-kbx
        ky_diff=kay-kby
            
        return kx_diff, ky_diff

    #correctly determines the moire vectors by subtracting reciprocal vectors
    def moire_real_v3(self,a,b,c,d):
        kax,kay=self.calculate_reciprocal(a,b)
        kbx,kby=self.calculate_reciprocal(c,d)
        
        #find out what the lowest symmetry is between the two vectors.
        symmetry=self.find_symmetry(kax,kay)
        symmetry2=self.find_symmetry(kbx,kby)
        low_symm=np.min([symmetry,symmetry2])

        theta_ref=self.find_angle(kax,kbx)
        rot_pass=360/low_symm/2
        if theta_ref>=rot_pass:
            offset=-rot_pass*2 - np.floor(((np.floor(theta_ref/rot_pass)-1)/2))*rot_pass*2
            R=[[np.cos(offset*np.pi/180),-np.sin(offset*np.pi/180)],[np.sin(offset*np.pi/180),np.cos(offset*np.pi/180)]]
            kbx=np.dot(R,kbx)
            kby=np.dot(R,kby)
        kx_diff=kax-kbx
        ky_diff=kay-kby
            
        rx_diff,ry_diff=self.calculate_reciprocal(kx_diff,ky_diff)
        
        return rx_diff, ry_diff
    
    #visualizes the moire lattice by plotting the vectors, reciprocal vectors, and moire vectors.    
    def visualize_moire(self):
        a=np.array([float(self.V1ax.text()),float(self.V1ay.text())])
        b=np.array([float(self.V1bx.text()),float(self.V1by.text())])
        c=np.array([float(self.V2ax.text()),float(self.V2ay.text())])
        d=np.array([float(self.V2bx.text()),float(self.V2by.text())])
        angle=float(self.A.text())*np.pi/180
        
        #determine rotation matrix
        R=[[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]]
        
        #redefine vectors
        ax=a[0]
        ay=a[1]
        bx=b[0]
        by=b[1]
        cx=c[0]
        cy=c[1]
        dx=d[0]
        dy=d[1]
        a_vec=np.hstack([ax,ay])
        b_vec=np.hstack([bx,by])
        c_vec=np.hstack([cx,cy])
        d_vec=np.hstack([dx,dy])
        
        #rotate second vectors
        c_rot=np.dot(R,c_vec)
        d_rot=np.dot(R,d_vec)
        
        #clear plots
        self.dataplot1.axes.cla()
        self.dataplot2.axes.cla()
        self.dataplot3.axes.cla()
        
        #plot real-space starting vectors
        self.dataplot1.axes.quiver(0,0,a_vec[0],a_vec[1],angles='xy', scale_units='xy', scale=1,width=0.01,linestyle='-',edgecolor='k',linewidth=0.5,color=['r'],label='Vector 1a')
        self.dataplot1.axes.quiver(0,0,b_vec[0],b_vec[1],angles='xy', scale_units='xy', scale=1,width=0.01,linestyle='--',edgecolor='k',linewidth=0.5,color=['r'],label='Vector 1b')
        self.dataplot1.axes.quiver(0,0,c_rot[0],c_rot[1],angles='xy', scale_units='xy', scale=1,width=0.01,linestyle='-',edgecolor='k',linewidth=0.5,color=['b'],label='Vector 2a')
        self.dataplot1.axes.quiver(0,0,d_rot[0],d_rot[1],angles='xy', scale_units='xy', scale=1,width=0.01,linestyle='--',edgecolor='k',linewidth=0.5,color=['b'],label='Vector 2b')
        self.dataplot1.axes.legend(loc='lower left',fontsize=3)
        lim=np.max([np.linalg.norm(a_vec),np.linalg.norm(b_vec),np.linalg.norm(c_rot),np.linalg.norm(d_rot)])
        self.dataplot1.axes.set_xlim(-lim*1.1,lim*1.1)
        self.dataplot1.axes.set_ylim(-lim*1.1,lim*1.1)
        self.dataplot1.axes.tick_params(direction="in")
        self.dataplot1.axes.yaxis.set_ticks_position('both')
        self.dataplot1.axes.xaxis.set_ticks_position('both')
        self.dataplot1.axes.grid(linewidth=0.25,color='k',linestyle='--')
        self.dataplot1.axes.set_xlabel('$x (nm)$',labelpad=0.1)
        self.dataplot1.axes.set_ylabel('$y (nm)$',labelpad=0.1)
        self.dataplot1.fig.tight_layout()
        self.dataplot1.draw()
        
        #plot reciprocal space vectors
        a_recip,b_recip=self.calculate_reciprocal(a_vec, b_vec)
        c_recip,d_recip=self.calculate_reciprocal(c_rot, d_rot)
        self.dataplot2.axes.quiver(0,0, a_recip[0],a_recip[1],angles='xy', scale_units='xy', scale=1,width=0.01,linestyle='-',edgecolor='k',linewidth=0.5,color=['r'],label='Vector 1a')
        self.dataplot2.axes.quiver(0,0, b_recip[0],b_recip[1],angles='xy', scale_units='xy', scale=1,width=0.01,linestyle='--',edgecolor='k',linewidth=0.5,color=['r'],label='Vector 1b')
        self.dataplot2.axes.quiver(0,0, c_recip[0],c_recip[1],angles='xy', scale_units='xy', scale=1,width=0.01,linestyle='-',edgecolor='k',linewidth=0.5,color=['b'],label='Vector 2a')
        self.dataplot2.axes.quiver(0,0, d_recip[0],d_recip[1],angles='xy', scale_units='xy', scale=1,width=0.01,linestyle='--',edgecolor='k',linewidth=0.5,color=['b'],label='Vector 2b')
        lim=np.max([np.linalg.norm(a_recip),np.linalg.norm(b_recip),np.linalg.norm(c_recip),np.linalg.norm(d_recip)])
        self.dataplot2.axes.set_xlim(-lim*2.1,lim*2.1)
        self.dataplot2.axes.set_ylim(-lim*2.1,lim*2.1)
        #plot difference in reciprocal space
        kx,ky=self.sub_reciprocal_v3(a_vec, b_vec, c_rot, d_rot)
        self.dataplot2.axes.quiver(c_recip[0],c_recip[1], kx[0],kx[1],angles='xy', scale_units='xy', scale=1,width=0.01,linestyle='-',edgecolor='k',linewidth=0.5,color=['g'],label='Difference a')
        self.dataplot2.axes.quiver(d_recip[0],d_recip[1],ky[0],ky[1],angles='xy', scale_units='xy', scale=1,width=0.01,linestyle='--',edgecolor='k',linewidth=0.5,color=['g'],label='Difference b')
        self.dataplot2.axes.legend(loc='lower left',fontsize=3)
        #plot parameters
        self.dataplot2.axes.tick_params(direction="in")
        self.dataplot2.axes.yaxis.set_ticks_position('both')
        self.dataplot2.axes.xaxis.set_ticks_position('both')
        self.dataplot2.axes.grid(linewidth=0.25,color='k',linestyle='--')
        self.dataplot2.axes.set_xlabel('$x^{-1} (nm^{-1})$',labelpad=0.1)
        self.dataplot2.axes.set_ylabel('$y^{-1} (nm^{-1})$',labelpad=0.1)
        self.dataplot2.fig.tight_layout()
        self.dataplot2.draw()
        
        #plot moire in real space
        m_x,m_y=self.moire_real_v3(a_vec, b_vec, c_rot, d_rot)
        self.dataplot3.axes.quiver(0,0, m_x[0],m_x[1],angles='xy', scale_units='xy', scale=1,width=0.01,linestyle='-',edgecolor='k',linewidth=0.5,color=['b'],label='Moire vector 1, $\lambda_m$ = '+str(np.round(np.linalg.norm(m_x),3))+' nm')
        self.dataplot3.axes.quiver(0,0, m_y[0],m_y[1],angles='xy', scale_units='xy', scale=1,width=0.01,linestyle='--',edgecolor='k',linewidth=0.5,color=['b'],label='Moire vector 2, $\lambda_m$ = '+str(np.round(np.linalg.norm(m_y),3))+' nm')
        lim=np.max([np.linalg.norm(m_x),np.linalg.norm(m_y)])
        self.dataplot3.axes.set_xlim(-lim*1.1,lim*1.1)
        self.dataplot3.axes.set_ylim(-lim*1.1,lim*1.1)
        self.dataplot3.axes.legend(loc='lower left',fontsize=3)
        self.dataplot3.axes.tick_params(direction="in")
        self.dataplot3.axes.yaxis.set_ticks_position('both')
        self.dataplot3.axes.xaxis.set_ticks_position('both')
        self.dataplot3.axes.grid(linewidth=0.25,color='k',linestyle='--')
        self.dataplot3.axes.set_xlabel('$x (nm)$',labelpad=0.1)
        self.dataplot3.axes.set_ylabel('$y (nm)$',labelpad=0.1)
        self.dataplot3.fig.tight_layout()
        self.dataplot3.draw()
#%%    
app = QApplication(sys.argv)
window = App()
window.show()
app.exec()
