#################### FILE NAME: RadialStatorTool.py #######################
#==========================================================================
# author: Jozef Stuijt & Nitish Anand                                     |
# 	: Master Student,                                                     |
#	: Process and Energy Departmemt,                                      |
#	: TU Delft,                                                           |
#	: The Netherlands                                                     |
#                                                                         |
# email : joe.stuijt@gmail.com                                            | 								                                          |
# Description: Defines the main class for building a radial stator stage  |
# based on nozzle coordinates                                             |
#                                                                         |
# Functions computeUPandDOWNsection, cart2Radial and CurvedNozzle used    |
# from the original work of Nitish Anand (nitish.ug2@gmail.com)           |
#==========================================================================

### Importing Classes ###
import numpy as np
import math
from scipy.optimize import minimize
from scipy.interpolate import interp1d
from copy import deepcopy

try:                                       # check if import possible related to cluster
   import matplotlib.pyplot as plt         # for plotting routines
   imprtplt = True
except ImportError: 
   print('\nWarining: Plotting module import unsuccessful\n') 
   imprtplt = False

import pdb				                  # Debugging Module

from geomdl import NURBS as ns
from geomdl import BSpline as bs
from geomdl import utilities as utils
### End Importing Classes ###

### Start Main Class ###

class AxialSupersonicStator(object):
    
	def __init__(self, 
                  flowAngle,
                  radiusOut, 
                  ActualRin , 
                  ActualRout, 
                  RealRout, 
                  nBlades = 18,
                  outFileMoc = 'Nozzle_coords.out', 
		          kernelRadiusRatio = 50.0,
                  ScaleNoz = 1.0,
                  TEminT = 0.05):
		#, Throat = 1e-3):

		self._nBlades           = nBlades                                 # No of blades
		self._pitch             = 1.0 #360.0/nBlades                           # Pitch in degrees
		self._flowAngle         = flowAngle                               # Flow Angle
		self._radiusOut         = radiusOut                               # Nozzel wall outlet Radius
		self._outFileMoc        = outFileMoc                              # MoC File Name
		self._kernelRadiusRatio = kernelRadiusRatio                       # wall radius at Throat
		self._ActualRout        = ActualRout                              # Stator OUT
		self._ActualRin         = ActualRin                               # Rotor IN
		self._RealRout          = RealRout                                # 
		self._ScaleNOZ          = ScaleNoz
		self._TEminT            = TEminT
		self._Moc_xyz           = np.loadtxt(self._outFileMoc)
		self._halfThroat        = self._Moc_xyz[1,1]
		self._rho_t             = self._kernelRadiusRatio*self._halfThroat
		#self._Throat 			= Throat

		# Start at initialization for optimization purposes
		self._computeUPandDOWNsection()
		self._cart2Radial()
		self._AxialNozzle()
#-----------------------------------------------------------------------------------------------#
#
#
#
	def BladeBuild(self):
		self._NURBSctrlpts()
		self._connectNozzleUPnDOWN()
		plt.show()
		self._FilletTrailingEdge()
		self.Vane=np.append(self._nozzleUp, self._nozzleDown[-1:0:-1], axis=0)
		self.Vane=np.append(self.Vane, self.Vane[0:1],axis=0)
		self._i_LE = closest_point([self._LEpt[0],self._LEpt[1],0],self.Vane)
		self.LEcoords = self.Vane[self._i_LE]
#-----------------------------------------------------------------------------------------------#
#
#
#
	def _computeUPandDOWNsection(self):
		"Computes lower and upper section of the nozzle"
		self._nozzleUp   = deepcopy(self._Moc_xyz[1:-2:])
		self._nozzleDown = deepcopy(self._Moc_xyz[1:-2:])
		self._nozzleCenter=np.column_stack((self._Moc_xyz[1:-2:,0],np.zeros(len(self._Moc_xyz[1:-2:]))))
		self._nozzleCenter=np.column_stack((self._nozzleCenter,np.zeros(len(self._Moc_xyz[1:-2:]))))
		for i in range(len(self._nozzleDown)):self._nozzleDown[i,1] *=-1.0
#-----------------------------------------------------------------------------------------------#
#
#
#
	def _cart2Radial(self):
		'Rotate and scale'
		dr=max(self._nozzleDown[:,0])
		dr *= self._ScaleNOZ

		"Translate to  x = x - dr "
		for i in range(len(self._nozzleUp)):self._nozzleUp[i,0] -= dr
		for i in range(len(self._nozzleDown)):self._nozzleDown[i,0] -= dr
		for i in range(len(self._nozzleCenter)):self._nozzleCenter[i,0] -= dr
		 
		"Rotate of the rotate_angle"
		for i in range(len(self._nozzleUp)):   self._nozzleUp[i,:] = Rotate_about_center(self._nozzleUp[i],self._flowAngle,[-abs(dr),0])
		for i in range(len(self._nozzleDown)): self._nozzleDown[i,:] = Rotate_about_center(self._nozzleDown[i],self._flowAngle,[-abs(dr),0])
		for i in range(len(self._nozzleCenter)): self._nozzleCenter[i,:] = Rotate_about_center(self._nozzleCenter[i],self._flowAngle,[-abs(dr),0])

		"Scale down to the correct outlet radius"
		radiusOut = math.sqrt(self._nozzleDown[-1,0]**2+self._nozzleDown[-1,1]**2)
		UPx=[self._nozzleUp[-1,0]]
		UPy=[self._nozzleUp[-1,1]]
		DOWNx=[self._nozzleDown[-1,0]]
		DOWNy=[self._nozzleDown[-1,1]]
		Scale = (self._pitch)*math.cos(math.radians(self._flowAngle))/math.sqrt((UPx[0]-DOWNx[0])**2+(UPy[0]-DOWNy[0])**2)
		#Scale= self._radiusOut/radiusOut
		self._rho_t *= Scale
		self.scaledHalfThroat = Scale*self._halfThroat

		for i in range(len(self._nozzleDown)):
			self._nozzleDown[i,0] *= Scale 
			self._nozzleDown[i,1] *= Scale 
         
		for i in range(len(self._nozzleUp)):
			self._nozzleUp[i,0] *= Scale 
			self._nozzleUp[i,1] *= Scale     
        
		for i in range(len(self._nozzleCenter)):
			self._nozzleCenter[i,0] *= Scale 
			self._nozzleCenter[i,1] *= Scale    
            
		#plt.plot(self._nozzleUp[:,0], self._nozzleUp[:,1])
		#plt.plot(self._nozzleDown[:,0],self._nozzleDown[:,1])
		#plt.plot(self._nozzleCenter[:,0],self._nozzleCenter[:,1])
		#plt.axis('equal')
		#plt.show()
#-----------------------------------------------------------------------------------------------#
#
#
#    
	def _AxialNozzle(self):
		for i in range(len(self._nozzleDown)):
			self._nozzleDown[i,1]+=(self._pitch+self._TEminT);
		A=Intersection(self._nozzleDown[-1],self._nozzleUp[-1],[math.tan(math.pi/2+math.radians(self._flowAngle)),math.tan(math.radians(self._flowAngle))])
		center = [(A[0]+self._nozzleDown[-1,0])/2,(A[1]+self._nozzleDown[-1,1])/2,0]
		self._nozzleUp= np.append(self._nozzleUp,[[A[0],A[1],0]],axis=0)
		self._FilletTrailingEdge()
		plt.plot(self._nozzleUp[:,0],self._nozzleUp[:,1])
		plt.plot(self._nozzleDown[:,0],self._nozzleDown[:,1])
		plt.axis('equal')
		self.UP = self._nozzleUp
		self.DOWN = self._nozzleDown

#-----------------------------------------------------------------------------------------------#
#
#
#
	def _CurvedNozzle(self):
		"Find angles"  
		ang_1=math.atan(self._nozzleUp[-1,1]/self._nozzleUp[-1,0])
		ang_2=math.atan(self._nozzleDown[-1,1]/self._nozzleDown[-1,0])

		"1st point of the curved nozzle will be same as that of the straight nozzle"
		self._curveNozcenter=np.array([[self._nozzleCenter[0,0],self._nozzleCenter[0,1]]])

		"Transformation matrix"
		trans = np.array([[math.cos(math.radians(self._flowAngle)),-math.sin(math.radians(self._flowAngle))],[math.sin(math.radians(self._flowAngle)),math.cos(math.radians(self._flowAngle))]])

		d_curve=0
		UPx=[self._nozzleUp[0,0]]
		UPy=[self._nozzleUp[0,1]]
		DOWNx=[self._nozzleDown[0,0]]
		DOWNy=[self._nozzleDown[0,1]]
		flag=0
        
		while flag==0:
			"Vec from the axis of the stator"
			vec_r = np.array([[self._curveNozcenter[-1,0]],[self._curveNozcenter[-1,1]]])
			fac=0.001

			"Rotate the vec by 70 degrees at the local point"
			vec_ang=(vec_r)-(trans.dot(vec_r)*fac)

			self._curveNozcenter=np.append(self._curveNozcenter,vec_ang.transpose(),axis=0)
            
			"Calculate stream-wise distance in curved nozzle"
			d_curve=d_curve+math.sqrt((self._curveNozcenter[-1,0]-self._curveNozcenter[-2,0])**2+(self._curveNozcenter[-1,1]-self._curveNozcenter[-2,1])**2)
			j=0
            
			while (d_curve>(((self._nozzleUp[j,1]-self._nozzleCenter[0,1])**2+(self._nozzleUp[j,0]-self._nozzleCenter[0,0])**2)**0.5)):
				j+=1
				try:self._nozzleUp[j,1]
				except:
					flag=1
					break
					d_curve=0
                    
			"Vec normal to the nozzle stream"
			Nvec=trans.dot(vec_r)
			Nvec=UnitVec(Nvec[0],Nvec[1])

			if j==len(self._nozzleUp):
				break
            
			if j==0:
				j=1
				d_span_j=(((self._nozzleUp[j-1,1]-self._nozzleCenter[j-1,1])**2+(self._nozzleUp[j-1,0]-self._nozzleCenter[j-1,0])**2)**0.5)
				d_span_j_=(((self._nozzleUp[j,1]-self._nozzleCenter[j,1])**2+(self._nozzleUp[j,0]-self._nozzleCenter[j,0])**2)**0.5)
				d_stream_j=(((self._nozzleCenter[0,1]-self._nozzleCenter[j-1,1])**2+(self._nozzleCenter[0,0]-self._nozzleCenter[j-1,0])**2)**0.5)
				d_stream_j_=(((self._nozzleCenter[0,1]-self._nozzleCenter[j,1])**2+(self._nozzleCenter[0,0]-self._nozzleCenter[j,0])**2)**0.5)
				d_span=Inter(d_curve,d_stream_j,d_stream_j_,d_span_j,d_span_j_)
			else:
				"Interpolate to find the exact width of the nozzle"
				d_span_j=(((self._nozzleUp[j-1,1]-self._nozzleCenter[j-1,1])**2+(self._nozzleUp[j-1,0]-self._nozzleCenter[j-1,0])**2)**0.5)
				d_span_j_=(((self._nozzleUp[j,1]-self._nozzleCenter[j,1])**2+(self._nozzleUp[j,0]-self._nozzleCenter[j,0])**2)**0.5)
				d_stream_j=(((self._nozzleCenter[0,1]-self._nozzleCenter[j-1,1])**2+(self._nozzleCenter[0,0]-self._nozzleCenter[j-1,0])**2)**0.5)
				d_stream_j_=(((self._nozzleCenter[0,1]-self._nozzleCenter[j,1])**2+(self._nozzleCenter[0,0]-self._nozzleCenter[j,0])**2)**0.5)
				d_span=Inter(d_curve,d_stream_j,d_stream_j_,d_span_j,d_span_j_)
            
			UPvec=Rotate_about_center(Nvec,270,self._curveNozcenter[-1])
			UPvec=UnitVec(UPvec[0],UPvec[1])
			DOWNvec=[-UPvec[0],-UPvec[1]]

			"Defined width of the curved nozzle"
			Up=np.concatenate((UPvec[0]*d_span+vec_ang[0],UPvec[1]*d_span+vec_ang[1]),axis=0)
			Down=np.concatenate((DOWNvec[0]*d_span+vec_ang[0],DOWNvec[1]*d_span+vec_ang[1]),axis=0)

			UPx=UPx+[Up[0]]
			UPy=UPy+[Up[1]]
			DOWNx=DOWNx+[Down[0]]
			DOWNy=DOWNy+[Down[1]]
			
            
		#plt.plot(UPx,UPy,'or')
		#plt.plot(DOWNx,DOWNy,'og')
		#plt.axis('equal')
		#plt.show()
		"Writing the location of the curved nozzle"
		self.DOWN=np.column_stack((DOWNx,DOWNy))
		self.DOWN=np.column_stack((self.DOWN,np.zeros(len(self.DOWN))))
		self.UP=np.column_stack((UPx,UPy))
		self.UP=np.column_stack((self.UP,np.zeros(len(self.UP))))

		self.DOWN=np.delete(self.DOWN,0,axis=0)
		self.UP=np.delete(self.UP,0,axis=0)
		for i in range(len(self.DOWN)): self.DOWN[i,:] = Rotate_about_center(self.DOWN[i], -self._pitch,[0,0])
        
		#pdb.set_trace()
		#plt.plot(self.UP[:,0],self.UP[:,1],'*b')
		#plt.plot(self.DOWN[:,0],self.DOWN[:,1],'om')
		#plt.plot(self.UP[0,0],self.UP[0,1],'ob')
		#plt.show()

		"Extending longer section of the nozzle"
		d1=1
		prev_d=math.sqrt(self.DOWN[-1,0]**2+self.DOWN[-1,1]**2)
		while abs(ang_2-ang_1)<self._pitch:
			vec_r = np.array([[self.UP[-1,0]],[self.UP[-1,1]]])
			if d1<1e-1:fac=0.01
			else: fac=0.001
			vec_ang=(vec_r)-(trans.dot(vec_r)*fac)
			self.UP=np.append(self.UP,np.transpose(np.append(vec_ang,[[0.0]],axis=0)),axis=0)
			d1=(((self.DOWN[-1,0]-self.UP[-1,0])**2+(self.DOWN[-1,1]-self.UP[-1,1])**2)**0.5)
			if d1<prev_d:
				prev_d=d1
			else:break

		self.nozArea = ((self._nozzleUp[-1,0]-self._nozzleDown[-1,0])**2+(self._nozzleUp[-1,1]-self._nozzleDown[-1,1])**2)**0.5

		r_d  = (self.DOWN[-1,0]**2+self.DOWN[-1,1]**2)**0.5
		r_u  = (self.UP[-1,0]**2+self.UP[-1,1]**2)**0.5
		r_u2 = (self.UP[-2,0]**2+self.UP[-2,1]**2)**0.5

		d    = abs(((self.DOWN[-1,0]-self.UP[-1,0])**2+(self.DOWN[-1,1]-self.UP[-1,1])**2)**0.5-self._TEminT)
		r = abs(self._ActualRout-(self.DOWN[-1,0]**2+self.DOWN[-1,1]**2)**0.5)

		#self._lastPUP = self.UP[-1]

		self.thrtStrtUp = self.UP[0]
		self.thrtEndUp = self._curveNozcenter[0]
		self.thrtStrtD = self.DOWN[0]
		self.thrtEndD = Rotate_about_center(self.thrtEndUp,-self._pitch,[0,0])

		self.nozEndUp = self._curveNozcenter[-1]
		self.nozStrtUp = self.UP[closest_point([self.nozEndUp[0],self.nozEndUp[1],0],self.UP)]
		self.nozStrtD = self.DOWN[-1]
		self.nozEndD = Rotate_about_center(self.nozEndUp,-self._pitch,[0,0])

		#plt.plot(self._curveNozcenter[:,0],self._curveNozcenter[:,1])
		#plt.plot(self._nozzleCenter[:,0],self._nozzleCenter[:,1])

		#plt.plot(self.thrtStrtUp[0],self.thrtStrtUp[1],'*r')
		#plt.plot(self.thrtEndUp[0],self.thrtEndUp[1],'*r')
		#plt.plot(self.thrtStrtD[0],self.thrtStrtD[1],'*r')
		#plt.plot(self.thrtEndD[0],self.thrtEndD[1],'*r')

		#plt.plot(self.nozEndUp[0],self.nozEndUp[1],'*b')
		#plt.plot(self.nozStrtUp[0],self.nozStrtUp[1],'*b')
		#plt.plot(self.nozStrtD[0],self.nozStrtD[1],'*b')
		#plt.plot(self.nozEndD[0],self.nozEndD[1],'*b')

		#plt.plot(self.UP[:,0],self.UP[:,1])
		#plt.plot(self.DOWN[:,0],self.DOWN[:,1])
		#plt.axis('equal')
		#plt.show()
		#d_tU = ((self.thrtStrtUp[0]-self.thrtEndUp[0])**2+(self.thrtStrtUp[1]-self.thrtEndUp[1])**2)**0.5
		#d_tD = ((self.thrtStrtD[0]-self.thrtEndD[0])**2+(self.thrtStrtD[1]-self.thrtEndD[1])**2)**0.5
		#d_t = abs(d_tU+d_tD-self._Throat)

		p_i = LinePtsInter(self.DOWN[-1],self.DOWN[-2],self.UP[-1],self.UP[-2])
		r_pi = (p_i[0]**2+p_i[1]**2)**0.5

		#if (r_u>r_d): d = 1
		#if (r_u>r_d) or (self._radiusOut<0): d = 1
		if (r_pi>r_d) or (self._radiusOut<0): d = 1
		self.d = d+r
		self.TE_thk = ((self.DOWN[-1,0]-self.UP[-1,0])**2+(self.DOWN[-1,1]-self.UP[-1,1])**2)**0.5
		#if (self._radiusOut<self._RealRout): d_t = 1
		#self.d_t = d+d_t

#-----------------------------------------------------------------------------------------------#
#
#
#
	def _NURBSctrlpts(self):

		#plt.plot(self.UP[:,0],self.UP[:,1],'ok')
		#plt.plot(self.DOWN[:,0],self.DOWN[:,1],'or')
		#plt.plot('show')
		
		pn0 = [self.UP[0,0],self.UP[0,1]]
		pn16 = [self.DOWN[0,0],self.DOWN[0,1]]

		VecCenStart = [pn0[0]-self._nozzleCenter[0,0],pn0[1]-self._nozzleCenter[0,1]]
		VecCenStart = UnitVec(VecCenStart[0],VecCenStart[1])
		kernOriUp = [self._nozzleUp[0,0]+self._rho_t*VecCenStart[0],self._nozzleUp[0,1]+self._rho_t*VecCenStart[1]]

		eps = 1e-3
		for i, pt in enumerate(self.UP):
			r_pt = ((pt[0]-kernOriUp[0])**2+(pt[1]-kernOriUp[1])**2)**0.5
			if abs(r_pt-self._rho_t)<eps: i_r1 = i
		#"for robustness, could be improved?"
		#pdb.set_trace()
		try: pr1 = [self.UP[i_r1,0],self.UP[i_r1,1]]
		except NameError: raise ValueError('Kernel origin not found, possible issue with nozzle coordinates')

		#pdb.set_trace()
		
		vec_ou = [self.UP[0,0]-kernOriUp[0],self.UP[0,1]-kernOriUp[1]]
		vec_ou = UnitVec(vec_ou[0],vec_ou[1])

		vec_or1 =[self.UP[i_r1,0]-kernOriUp[0],self.UP[i_r1,1]-kernOriUp[1]]
		vec_or1 = UnitVec(vec_or1[0],vec_or1[1])

		ang_u = 5;#int(0.5*math.degrees(math.acos(vec_ou[0]*vec_or1[0]+vec_ou[1]*vec_or1[1])))

		#pdb.set_trace()

		'Enlarge existing part of nozzle using continuity circle'
		ContCircUp = np.zeros((ang_u,3))
		for i in range(ang_u): ContCircUp[i] = Rotate_about_center(pn0,-i,kernOriUp)

		'Find 2 control points using circle to assure C^2 continuity with NURBS'
		pdb.set_trace()
		pn1 = [ContCircUp[int(ang_u/2),0],ContCircUp[int(ang_u/2),1]]
		pn2 = [ContCircUp[-1,0],ContCircUp[-1,1]]
		pc1 = [ContCircUp[int(ang_u/4),0],ContCircUp[int(ang_u/4),1]]
		pc2 = [ContCircUp[int(3*ang_u/4),0],ContCircUp[int(3*ang_u/4),1]]

		'Find 2 control points for down nozzle mirroring circle'
		refnozctr0 = [self._nozzleCenter[0,0],self._nozzleCenter[0,1]]
		refnozctr1 = [self._nozzleCenter[1,0],self._nozzleCenter[1,1]]

		#pn14 = mirror_point([ContCircUp[int(ang_u/2),0],ContCircUp[int(ang_u/2),1]],refnozctr0,refnozctr1)
		pn14 = mirror_point(pn2,refnozctr0,refnozctr1)
		pn14 = Rotate_about_center(pn14,-self._pitch,[0,0])
		pn14 = [pn14[0],pn14[1]]

		#pn15 = mirror_point([ContCircUp[int(ang_u/4),0],ContCircUp[int(ang_u/8),1]],refnozctr0,refnozctr1)
		pn15 = mirror_point(pn1,refnozctr0,refnozctr1)
		pn15 = Rotate_about_center(pn15,-self._pitch,[0,0])
		pn15 = [pn15[0],pn15[1]]

		#pc3 = mirror_point([ContCircUp[int(ang_u/4),0],ContCircUp[int(ang_u/4),1]],refnozctr0,refnozctr1)
		pc3 = mirror_point(pc2,refnozctr0,refnozctr1)
		pc3 = Rotate_about_center(pc3,-self._pitch,[0,0])
		pc3 = [pc3[0], pc3[1]]

		#pc4 = mirror_point([ContCircUp[int(ang_u/16),0],ContCircUp[int(ang_u/16),1]],refnozctr0,refnozctr1)
		pc4 = mirror_point(pc1,refnozctr0,refnozctr1)
		pc4 = Rotate_about_center(pc4,-self._pitch,[0,0])
		pc4 = [pc4[0], pc4[1]]

		'Find point in nozzle up perpendicular to start nozzle down'
		i_r2 = closest_point(self.DOWN[0],self.UP)
		pr2 = [self.UP[i_r2,0],self.UP[i_r2,1]]

		"Find blade leading edge radius"
		pNC0 = self._nozzleCenter[0]
		pNC1 = self._nozzleCenter[1]
		pNC0adj = Rotate_about_center(pNC0,-self._pitch,[0,0])
		pNC1adj = Rotate_about_center(pNC1,-self._pitch,[0,0])
		m1,k1 = Line2Pts(pNC0,pNC1)
		m2,k2 = Line2Pts(pNC0adj,pNC1adj)
		m3 = -1/m1
		k3 = pNC0[0]*(m1-m3)+k1
		x_i = (k3-k2)/(m2-m3)
		y_i = m2*x_i+k2
		pn7 = [x_i,y_i]

		self._R_LE = (x_i**2+y_i**2)**0.5
		if self._R_LE>self._ActualRin: self._R_LE = self._ActualRin

		'Find first point and origin of convergent circle'
		[pn7_1,pn7_2] = InterLineCirc(self._nozzleCenter[0],pn0,self._R_LE)
		if (pn7_1[0]<pn0[0]): pn7 = pn7_1
		else: pn7 = pn7_2
		pn7 = [pn7[0],pn7[1]]
		CircOri = [0.5*(pn0[0]+pn7[0]),0.5*(pn0[1]+pn7[1])]

		'Add two more points for C^2 continuity'
		pn6 = Rotate_about_center(pn7,5,CircOri)
		pn8 = Rotate_about_center(pn7,-5,CircOri)
		pn6 = [pn6[0],pn6[1]]
		pn8 = [pn8[0],pn8[1]]

		'Find additional points in circle perpendicular to the diameter just found'
		pn4 = Rotate_about_center(pn0,-90,CircOri)
		pn4 = [pn4[0],pn4[1]]
		self._LEpt = pn4

		'Add two more points for C^2 continuity'
		pn3 = Rotate_about_center(pn4,5,CircOri)		
		pn5 = Rotate_about_center(pn4,-5,CircOri)
		pn3 = [pn3[0],pn3[1]]
		pn5 = [pn5[0],pn5[1]]

		'Find additional points between pn0, pn4, pn7'
		pa2 = Rotate_about_center(pn0,-45,CircOri)
		pa2 = [pa2[0],pa2[1]]
		pa3 = Rotate_about_center(pa2,5,CircOri)
		pa1 = Rotate_about_center(pa2,-5,CircOri)
		pa3 = [pa3[0],pa3[1]]
		pa1 = [pa1[0],pa1[1]]
		pa5 = Rotate_about_center(pn0,-135,CircOri)
		pa5 = [pa5[0],pa5[1]]
		pa4 = Rotate_about_center(pa5,5,CircOri)
		pa6 = Rotate_about_center(pa5,-5,CircOri)
		pa4 = [pa4[0],pa4[1]]
		pa6 = [pa6[0],pa6[1]]
		#self._LEpt = pn4

		'Find additional points in circle perpendicular to the diameter just found'
		pn4 = Rotate_about_center(pn0,-90,CircOri)
		pn4 = [pn4[0],pn4[1]]
		self._LEpt = pn4

		

		'Find the end point of convergent circle'
		# Mirror point 4
		#pr5 = mirror_point(pn4,CircOri,pn0)
		

		vec_n4 = [CircOri[0]-pn4[0],CircOri[1]-pn4[1]]
		vec_n4 = UnitVec(vec_n4[0],vec_n4[1])

		#vec_n4 = [pr5[0]-CircOri[0],pr5[1]-CircOri[1]]
		#vec_n4 = UnitVec(vec_n4[0],vec_n4[1])

		vec_pr1 = [pr1[0]-CircOri[0],pr1[1]-CircOri[1]]
		vec_pr1 = UnitVec(vec_pr1[0],vec_pr1[1])

		angCirc = 360-2*math.degrees(math.acos(vec_n4[0]*vec_pr1[0]+vec_n4[1]*vec_pr1[1]))

		#pn10 = Rotate_about_center(pn0,-angCirc,CircOri)
		pn10 = Rotate_about_center(pr1,-angCirc,CircOri)
		pn10 = [pn10[0],pn10[1]]

		'Add two more points for C^2 continuity'
		pn9 = Rotate_about_center(pn10,5,CircOri)
		pn11 = Rotate_about_center(pn10,-5,CircOri)
		pn9 = [pn9[0],pn9[1]]
		pn11 = [pn11[0],pn11[1]]

		'Find the last points connecting the circle and down section'
		UPshort = self.UP[i_r1:i_r2]

		i_r3 = FindXCurve(UPshort,0.45)+i_r1
		pr3 = [self.UP[i_r3,0],self.UP[i_r3,1]]

		ref_len = ((pn7[0]-pr1[0])**2+(pn7[1]-pr1[1])**2)**0.5

		if self.UP[i_r3+1,1]<self.UP[i_r3-1,1]: vecr3 = [self.UP[i_r3+1,0]-self.UP[i_r3-1,0],self.UP[i_r3+1,1]-self.UP[i_r3-1,1]]
		else: vecr3 = [self.UP[i_r3-1,0]-self.UP[i_r3+1,0],self.UP[i_r3-1,1]-self.UP[i_r3+1,1]]		
		vecr3 = UnitVec(vecr3[0],vecr3[1])
		vecr3 = [vecr3[1],-vecr3[0]]

		pn12 = [pr3[0]+0.25*ref_len*vecr3[0],pr3[1]+0.25*ref_len*vecr3[1]]

		i_r4 = FindXCurve(UPshort,0.55)+i_r1
		pr4 = [self.UP[i_r4,0],self.UP[i_r4,1]]

		if self.UP[i_r4+1,1]<self.UP[i_r4-1,1]: vecr4 = [self.UP[i_r4+1,0]-self.UP[i_r4-1,0],self.UP[i_r4+1,1]-self.UP[i_r4-1,1]]
		else: vecr4 = [self.UP[i_r4-1,0]-self.UP[i_r4+1,0],self.UP[i_r4-1,1]-self.UP[i_r4+1,1]]		
		vecr4 = UnitVec(vecr4[0],vecr4[1])
		vecr4 = [vecr4[1],-vecr4[0]]


		pn13 = [pr4[0]+0.25*ref_len*vecr4[0],pr4[1]+0.25*ref_len*vecr4[1]]
		
		#self._ctrlpts = [pn0,pn1,pn2,pn3,pn4,pn5,pn6,pn7,pn8,pn9,pn10,pn11,pn12,pn13,pn14,pn15,pn16]
		self._ctrlpts = [pn0,pn1,pc1,pn2,pc2,pa1,pa2,pa3,pn3,pn4,pn5,pa4,pa5,pa6,pn6,pn7,pn8,pn9,pn10,pn11,pn12,pn13,pn14,pc3,pn15,pc4,pn16]
		#self._ctrlpts = [pn0,pn1,pn2,pn3,pn4,pn5,pn6,pn7,pn8,pn12,pn13,pn14,pn15,pn16]

		self._distMshBndry = ((pn16[0]-pn14[0])**2+(pn16[1]-pn14[1])**2)**0.5
		#self._convRefPt = pn10

		plt.plot(pn0[0],pn0[1],'*r')
		plt.plot(pn1[0],pn1[1],'*b')
		plt.plot(pn2[0],pn2[1],'*b')	
		plt.plot(pn3[0],pn3[1],'*b')
		plt.plot(pn4[0],pn4[1],'*b')
		plt.plot(pn5[0],pn5[1],'*b')
		plt.plot(pn6[0],pn6[1],'*b')
		plt.plot(pn7[0],pn7[1],'*b')
		plt.plot(pn8[0],pn8[1],'*b')
		plt.plot(pn9[0],pn9[1],'*b')
		plt.plot(pn10[0],pn10[1],'*c')
		plt.plot(pn11[0],pn11[1],'*b')
		plt.plot(pn12[0],pn12[1],'*b')
		plt.plot(pn13[0],pn13[1],'*b')
		plt.plot(pn14[0],pn14[1],'*b')
		plt.plot(pn15[0],pn15[1],'*b')
		plt.plot(pn16[0],pn16[1],'*r')

		#plt.plot(pc1[0],pc1[1],'+b')
		#plt.plot(pc2[0],pc2[1],'+b')
		#plt.plot(pc3[0],pc3[1],'+b')
		#plt.plot(pc4[0],pc4[1],'+b')

		#plt.plot(pa1[0],pa1[1],'+c')
		#plt.plot(pa2[0],pa2[1],'+c')
		#plt.plot(pa3[0],pa3[1],'+c')
		#plt.plot(pa4[0],pa4[1],'+c')
		#plt.plot(pa5[0],pa5[1],'+c')
		#plt.plot(pa6[0],pa6[1],'+c')

		#plt.plot(pr1[0],pr1[1],'*g')
		#plt.plot(pr2[0],pr2[1],'*g')
		#plt.plot(pr3[0],pr3[1],'*g')
		#plt.plot(pr4[0],pr4[1],'*g')
		#plt.plot(pr5[0],pr5[1],'+g')

		#plt.plot(CircOri[0],CircOri[1],'ok')

#-----------------------------------------------------------------------------------------------#
#
#
#
	def _connectNozzleUPnDOWN(self):

		# Create a NURBS curve instance
		curve = ns.Curve()
		#curve = bs.Curve()
		# Set up the NURBS curve
		curve.degree  = len(self._ctrlpts)-1
		weights = [1,10,15,20,25,1,5,1,1,20,1,1,5,1,1,20,1,1,1,1,1,1,25,5,25,5,1]
		# Add weight (1 = no weight) to control points
		#curve.ctrlptsw = ctrlptsw
		for i, pt in enumerate(self._ctrlpts): self._ctrlpts[i] = [pt[0]*weights[i],pt[1]*weights[i],weights[i]]
		curve.ctrlpts = self._ctrlpts
		#for i, pt in enumerate(self._ctrlpts): self._ctrlpts[i] = [pt[0],pt[1]]
		#curve.ctrlpts = self._ctrlpts
		#weights = [1,1,1,1,1,1,1,1,1,1,1,1,1,1]
		#ctrlptsw = []
		#for i, pt in enumerate(self._ctrlpts): #self._ctrlpts[i] = [pt[0],pt[1],1]
			#curve.ctrlpts[i] = [self._ctrlpts[i][0], self._ctrlpts[i][1], weights[i]]
			#ctrlptsw.append([self._ctrlpts[i][0], self._ctrlpts[i][1], weights[i]])
		#curve.ctrlptsw = ctrlptsw

		# Auto-generate the knot vector
		#curve.knotvector = utils.generate_knot_vector(curve.degree, len(curve.ctrlpts))
		curve.knotvector = utils.generate_knot_vector(curve.degree, len(curve.ctrlpts))
		#curve.knotvector = utils.generate_knot_vector()
		# Calculate curve points
		curve.delta = 0.001
		curve.evaluate()
		# Arrange curve points for plotting
		curvepts_x = []
		curvepts_y = []

		for pt in curve.curvepts:
			curvepts_x.append(pt[0])
			curvepts_y.append(pt[1])
		curvepts_z = np.zeros(len(curve.curvepts))
		Section = np.column_stack((curvepts_x,curvepts_y,curvepts_z))

		self.DOWN = np.append(Section,self.DOWN,axis=0)

		self._nozzleUp   = self.UP
		self._nozzleDown = self.DOWN

		plt.plot(self._nozzleDown[:,0],self._nozzleDown[:,1],'*r')
		plt.plot(self._nozzleUp[:,0],self._nozzleUp[:,1],'*m')
		#plt.show()
#-----------------------------------------------------------------------------------------------#
#
#
#
		          
	def _FilletTrailingEdge(self):
		slope_Up = (self._nozzleUp[-2,1]-self._nozzleUp[-1,1])/(self._nozzleUp[-2,0]-self._nozzleUp[-1,0])
		rot=Rotate_about_center(self._nozzleDown[-2],90,self._nozzleDown[-1])

		slope_Down = (rot[1]-self._nozzleDown[-1,1])/(rot[0]-self._nozzleDown[-1,0])

		Up_Point=Intersection(self._nozzleDown[-1],self._nozzleUp[-1],[slope_Down,slope_Up])

		self._nozzleUp[-1]=[Up_Point[0],Up_Point[1],0.0]

		center = [(self._nozzleDown[-1,0]+Up_Point[0])/2,(self._nozzleDown[-1,1]+Up_Point[1])/2]
		self._TEctr = center

		angle=np.linspace(0,180,11)
		self.fillet=[[Up_Point[0],Up_Point[1],0.0]]
		self._nozzleUp=np.delete(self._nozzleUp, -2, axis=0)

		for i in angle:
			c=Rotate_about_center(Up_Point,i,center)
			self._nozzleUp=np.append(self._nozzleUp,[c],axis=0)

			self.fillet=np.append(self.fillet,[c],axis=0)

		self.TEpt = self.fillet[int(len(self.fillet)/2)]

		#plt.plot(self.fillet[0,0],self.fillet[0,1],'or')
		self._lastPUP = self.fillet[0]
		self.bldExitR = ((self._lastPUP[0])**2+(self._lastPUP[1])**2)**0.5
		#pdb.set_trace()

		self._nozzleUp=np.append(self._nozzleUp,[self._nozzleDown[-1]],axis=0)

		#self._lastPUP = self.UP[-1]
#-----------------------------------------------------------------------------------------------#
#
#
#
	# Meshing Boundary curves
	def MeshBoundariesBlade(self):

		i_TE = closest_point([self.TEpt[0],self.TEpt[1],0],self.Vane)

		self.LEcoords = self.Vane[self._i_LE]

		self.SUCTION = self.Vane[i_TE:self._i_LE+1]

		PRES1 = self.Vane[self._i_LE:-1]
		PRES2 = self.Vane[0:i_TE+1]

		self.PRESSURE = np.append(PRES1,PRES2,axis=0)

		#"Find Limit Points"
		p0 = self._curveNozcenter[0]
		p1 = self._curveNozcenter[-1]

		vec_LE_NozCtr = UnitVec(p0[0]-self.Vane[self._i_LE,0], p0[1]-self.Vane[self._i_LE,1])
		vec_ctr = UnitVec(p1[0]-p0[0],p1[1]-p0[1])
		ref_ang = math.degrees(math.acos(vec_LE_NozCtr[0]*vec_ctr[0]+vec_LE_NozCtr[1]*vec_ctr[1]))

		p2 = Rotate_about_center(self.Vane[self._i_LE],ref_ang/2,p0)

		p5 = [p0[0]-vec_ctr[0]*self._distMshBndry,p0[1]-vec_ctr[1]*self._distMshBndry]

		[p61,p62] = InterLineCirc(p2,p5,self._R_LE)
		if (p61[0]<p2[0]): p6 = p61
		else: p6 = p62

		[p31,p32] = InterLineCirc(p2,p5,self._ActualRin)
		if (p31[0]<p0[0]): p3 = p31
		else: p3 = p32

		pr4 = Rotate_about_center(self._lastPUP,self._pitch/2,[0,0])
		[p41,p42] = InterLineCirc([0,0],pr4,self._RealRout)
		if (p41[0]<pr4[0]): p4 = p41
		else: p4 = p42

		#plt.plot(p0[0],p0[1],'*g')
		#plt.plot(p1[0],p1[1],'*g')
		#plt.plot(p2[0],p2[1],'*g')
		#plt.plot(p3[0],p3[1],'*g')
		#plt.plot(p4[0],p4[1],'*g')
		#plt.plot(p5[0],p5[1],'+g')
		#plt.plot(p6[0],p6[1],'*g')

		self.PERIO1 = np.zeros(((7,3)))
		self.PERIO1[0] = [p3[0],p3[1],0]
		self.PERIO1[1] = [p6[0],p6[1],0]
		self.PERIO1[2] = [p2[0],p2[1],0]
		self.PERIO1[3] = [p5[0],p5[1],0]
		self.PERIO1[4] = [p0[0],p0[1],0]
		self.PERIO1[5] = [p1[0],p1[1],0]
		self.PERIO1[6] = [p4[0],p4[1],0]

		self.PERIO2 = np.zeros(((len(self.PERIO1),3)))
		for i in range(len(self.PERIO2)): self.PERIO2[i] = Rotate_about_center(self.PERIO1[i],-self._pitch,[0,0])

		self.INLET = np.zeros(((int(self._pitch),3)))
		self.OUTLET = np.zeros(((int(self._pitch),3)))
		for i in range(int(self._pitch)): 
			self.INLET[i] = Rotate_about_center(p3,-i,[0,0])
			self.OUTLET[i] = Rotate_about_center(p4,-i,[0,0])

		self.INLET[-1] = Rotate_about_center(p3,-self._pitch,[0,0])
		self.OUTLET[-1] = Rotate_about_center(p4,-self._pitch,[0,0])

		self.BladeRout = (self.TEpt[0]**2+self.TEpt[1]**2)**0.5
		self.pitchLen = 2*math.pi*(self._pitch/360)*self.BladeRout
		self.chordLen = ((self.LEcoords[0]-self.TEpt[0])**2+(self.LEcoords[1]-self.TEpt[1])**2)**0.5

#---------------------------------------------------------------------------------------------#
# HELP FUNCTIONS #

def UnitVec(u,v):
	mag = math.sqrt(u**2+v**2)
	u1 = u/mag
	v1 = v/mag
	return [u1,v1]
#-----------------------------------------------------------------------------------------------#
#
#
#
# Liner Interpolation function for Table
def Inter(X,x1,x2,y1,y2):
    return(y1+(X-x1)*((y2-y1)/(x2-x1)))
#-----------------------------------------------------------------------------------------------#
#
#
#
def Rotate_about_center(y_t,ang,rho_d):
	x = [y_t[0],y_t[1]]
	center = [rho_d[0],rho_d[1]]
	x=np.array(x)-np.array(center)
	ang = math.radians(ang)
	M = np.matrix([[math.cos(ang), -math.sin(ang)],[math.sin(ang),math.cos(ang)]])
	x = np.transpose(np.matrix(x))
	rot_x = np.matrix(M) * np.matrix(x)
	rot_x = np.transpose(rot_x)
	rot_x = np.mat(rot_x) + center    
	return rot_x[0,0],rot_x[0,1], 0.0

# Intersection(): Finds position of the intersection of two lines defined by x,y and slope
#-----------------------------------------------------------------------------------------------#
#
#
#
def Intersection(x1,x2,m):
	x3 = (m[0]*x1[0] - m[1]*x2[0] + x2[1] - x1[1])/(m[0]-m[1])
	y3 = (m[0]*(x3-x1[0]))+x1[1]
	return x3,y3
#-----------------------------------------------------------------------------------------------#
#
#
#
def LinePtsInter(p1,p2,p3,p4):
	m1,k1 = Line2Pts(p1,p2)
	m2,k2 = Line2Pts(p3,p4)
	x = (k2-k1)/(m1-m2)
	y = m1*x+k1
	return [x,y]
#-----------------------------------------------------------------------------------------------#
#
#
#
# Find the gradient and y-intercept of line between two points
def Line2Pts(p1,p2):
	if (p1[0]>p2[0]): m = (p1[1]-p2[1])/(p1[0]-p2[0])
	else: m = (p2[1]-p1[1])/(p2[0]-p1[0])
	k = p1[1]-m*p1[0]
	return m,k
#-----------------------------------------------------------------------------------------------#
#
#
#
# Find the intersection of line between two points and circle of radius r
def InterLineCirc(p1,p2,r):
	[m,k] = Line2Pts(p1,p2)
	A = m**2+1
	B = 2*m*k
	C = k**2-r**2
	x1 = (-B+(B**2-4*A*C)**0.5)/(2*A)
	y1 = x1*m+k
	x2 = (-B-(B**2-4*A*C)**0.5)/(2*A)
	y2 = x2*m+k
	p3 = [x1,y1]
	p4 = [x2,y2]
	return p3,p4

def FindMidCurve(array):
	length = []
	length.append(0)
	for i in range(len(array)):
		if (i==0): pass
		else:
			d = ((array[i,0]-array[i-1,0])**2+(array[i,1]-array[i-1,1])**2)**0.5
			length.append(d+length[i-1])
	return np.argmin(np.abs([length-(length[-1]/2)]))
#-----------------------------------------------------------------------------------------------#
#
#
#
# find point at x/L (0<x<1) length of curve
def FindXCurve(array,x):
	length = []
	length.append(0)
	for i in range(len(array)):
		if (i==0): pass
		else:
			d = ((array[i,0]-array[i-1,0])**2+(array[i,1]-array[i-1,1])**2)**0.5
			length.append(d+length[i-1])
	return np.argmin(np.abs([length-(length[-1]*x)]))

def closest_point(pt, array):
	deltas = array - pt
	dist_2 = np.einsum('ij,ij->i', deltas, deltas)
	return np.argmin(dist_2)
#-----------------------------------------------------------------------------------------------#
#
#
#
# https://stackoverflow.com/questions/8954326/how-to-calculate-the-mirror-point-along-a-line#8954454
def mirror_point(pt,P1,P2):
	A = P2[1]-P1[1]
	B = -(P2[0]-P1[0])
	C = -A*P1[0] - B*P1[1]

	M = math.sqrt(A*A+B*B)

	Ap = A/M
	Bp = B/M
	Cp = C/M

	D = Ap*pt[0]+Bp*pt[1]+Cp
	return [pt[0]-2*Ap*D,pt[1]-2*Bp*D]

### --- Optimization functions
#-----------------------------------------------------------------------------------------------#
#
#
#
def OptimizeScaling(args, flowAngle, ActualRin, ActualRout, RealRout, nBlades = 18, outFileMoc = 'Nozzle_coords.out', kernelRadiusRatio = 50.0, TEminT = 0.0005):
	opti = RadialSupersonicStator(flowAngle, args[0], ActualRin , ActualRout, RealRout, nBlades, outFileMoc, kernelRadiusRatio, args[1], TEminT)
	print("Tip Radius: %.10f Scaling Factor: %.10f RESIDUAL: %.10f"%(args[0],args[1],opti.d))
	return opti.d

#def OptimizeThroat(args, flowAngle, ActualRin, ActualRout, RealRout, nBlades = 18, outFileMoc = 'Nozzle_coords.out', kernelRadiusRatio = 50.0, 
#		TEminT = 0.0005, Throat = 1e-3):
#	opti = RadialSupersonicStator(flowAngle, args[0], ActualRin , ActualRout, RealRout, nBlades, outFileMoc, kernelRadiusRatio, args[1], TEminT, Throat)
#	print("Tip Radius: %.10f Scaling Factor: %.10f RESIDUAL: %.10f"%(args[0],args[1],opti.d_t))
#	return opti.d_t

#def OptimizeThroatRadius(args, flowAngle, ActualRin, ActualRout, RealRout, nBlades = 18, outFileMoc = 'Nozzle_coords.out', kernelRadiusRatio = 50.0, 
#		TEminT = 0.0005, Throat = 1e-3):
#	opti = RadialSupersonicStator(flowAngle, args[0], ActualRin , ActualRout, RealRout, nBlades, outFileMoc, kernelRadiusRatio, args[1], TEminT, Throat)
#	print("Tip Radius: %.10f Scaling Factor: %.10f RESIDUAL: %.10f"%(args[0],args[1],opti.d_t))
#	return opti.d_t

### END
#flowAngle, radiusOut, ActualRin , ActualRout, RealRout, nBlades = 18, outFileMoc = 'Nozzle_coords.out', 
		#kernelRadiusRatio = 50.0, ScaleNoz = 1.0, TEminT = 0.0005, Throat = 1e-3

Axial = AxialSupersonicStator(70,1,1,1,5)
Axial.BladeBuild()