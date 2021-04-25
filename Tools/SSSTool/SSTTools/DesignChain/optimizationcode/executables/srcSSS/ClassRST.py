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
import sys
from srcSSS.Functions import *
#from Tools.TestSSSTool.SST import flowAngle, outFileMoc

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

class SuperSonicStator(object):
	#-----------------------------------------------------------------------------------------------#
	#
	#
	#			
	def __init__(self,IN):
		"Automatic Initialization Function: Takes in the read file"
		self.IN = IN
		self.InitializeFromFile()
	#-----------------------------------------------------------------------------------------------#
	#
	#
	#					
	def InitializeFromFile(self):
		"Initializing values from the Read File"
		self._flowAngle 		= float(self.IN['flowAngle'])
		self._radiusOut 		= float(self.IN['radiusOut'])
		self._outFileMoc        = self.IN['outFileMoc']                              		# MoC File Name
		self._kernelRadiusRatio = float(self.IN['kernelRadiusRatio'])
		self._ActualRout        = float(self.IN['ActualRout'])                              # Stator OUT
		self._ActualRin         = float(self.IN['ActualRin'])                               # Stator IN
		self._RealRout          = float(self.IN['RealRout'])                                # Rotor IN
		self._ScaleNOZ          = float(self.IN['ScaleNoz'])
		self._TEminT            = float(self.IN['TEminT'])
		self._Moc_xyz           = np.loadtxt(self._outFileMoc)
		self._halfThroat        = self._Moc_xyz[1,1]
		self._rho_t             = self._kernelRadiusRatio*self._halfThroat
		self._AXIAL 			= False
		self._RADIAL 			= False
		"If you are trying to DEBUG put this to True"
		self._DEBUG				= False
		self._SecretScale		= 1
		if self.IN['STATOR_KIND'] == "AXIAL":
			self._AXIAL = True
			self._InitializeAxial()
		elif self.IN["STATOR_KIND"] == "RADIAL":
			self._RADIAL = True
			self._InitializeRadial()
		else:
			print("%s is invalid stator type... Exiting... "%(self.IN["STATOR_KIND"]))
			sys.exit()
			
	#-----------------------------------------------------------------------------------------------#
	#
	#
	#						
	def _InitializeAxial(self):
		"Initializing Axial Pitch for AXIAL Stator"
		try: 
			self._pitch = float(self.IN['AxialPitch'])
			print("\n\nAXIAL SUPERSONIC STATOR INITIALIZED:\n")
		except: 
			print("Error :: AxialPitch not defined in Configuration File. \n\nExiting.... \n\n")
			sys.exit()
	#-----------------------------------------------------------------------------------------------#
	#
	#
	#			
	def _InitializeRadial(self):
		"Initializing Radial Pitch for RADIAL Stator"
		try: 
			self._pitch = 360.0/float(self.IN['nBlades'])
			self._nBlades = int(self.IN['nBlades'])
			print("\n\nRADIAL SUPERSONIC STATOR INITIALIZED:\n")
		except: 
			print("Error :: nBlades not defined in Configuration File. \n\n Exiting.... \n\n")
			sys.exit()
	#-----------------------------------------------------------------------------------------------#
	#
	#
	#				
	def NewBladeBuild(self):
		"Main Function to build the stators"
		if self._RADIAL:
			self.LogSpiral()
			self._computeUPandDOWNsection()
			self._newcart2Radial()
			self._CurvedNozzle()
			if imprtplt:
				plt.plot(self.UP[:,0],self.UP[:,1])
				plt.plot(self.DOWN[:,0],self.DOWN[:,1])
				plt.show()
			dx = self.DOWN[-1,0] + self._ActualRout
			dy = self.DOWN[-1,1] 
			for i in range(len(self.DOWN)):
				self.DOWN[i,0] -=dx 
				self.DOWN[i,1] -=dy 
			for i in range(len(self.UP)):
				self.UP[i,0] -=dx
				self.UP[i,1] -=dy
			nLS = len(self._LogSpiral)
			self._MeshPropCalc()		
			for i in range(0,nLS):
				self.UP = np.append(self.UP, [self._LogSpiral[nLS-i-1]], axis=0)
			for i in range(len(self.DOWN)): self.DOWN[i,:] = Rotate_about_center(self.DOWN[i], -self._pitch,[0,0])
			self.BladeBuild()
		elif self._AXIAL:
				self._computeUPandDOWNsection()
				self._newcart2Radial()
				self._AxialNozzle()
				for i in range(len(self.DOWN)): self.DOWN[i,1] -=self._pitch
				self._MeshPropCalc()
				for i in range(len(self.DOWN)): self.DOWN[i,1] +=self._pitch
				self.BladeBuild()
		else:
			sys.exit()

	def _MeshPropCalc(self):

		# plt.plot(self.UP[:,0],self.UP[:,1])
		# plt.plot(self.DOWN[:,0],self.DOWN[:,1])
		# plt.axis('equal')
		# plt.show()

		self.nozArea = ((self._nozzleUp[-1,0]-self._nozzleDown[-1,0])**2+(self._nozzleUp[-1,1]-self._nozzleDown[-1,1])**2)**0.5
		
		r_d  = (self.DOWN[-1,0]**2+self.DOWN[-1,1]**2)**0.5
		r_u  = (self.UP[-1,0]**2+self.UP[-1,1]**2)**0.5
		r_u2 = (self.UP[-2,0]**2+self.UP[-2,1]**2)**0.5

		d    = abs(((self.DOWN[-1,0]-self.UP[-1,0])**2+(self.DOWN[-1,1]-self.UP[-1,1])**2)**0.5-self._TEminT)
		r = abs(self._ActualRout-(self.DOWN[-1,0]**2+self.DOWN[-1,1]**2)**0.5)

		#self._lastPUP = self.UP[-1]
		
		
		#self.thrtEndUp = self._curveNozcenter[0]
		self.thrtStrtUp = self.UP[0]
		self.thrtEndUp = [(self.UP[0,0]+self.DOWN[0,0])/2,(self.UP[0,1]+self.DOWN[0,1])/2]
		
		self.thrtStrtD = self.DOWN[0]
		self.thrtEndD = Rotate_about_center(self.thrtEndUp,-self._pitch,[0,0])

		self.nozEndUp = [(self.UP[-1,0]+self.DOWN[-1,0])/2,(self.UP[-1,1]+self.DOWN[-1,1])/2]
		self.nozStrtUp = self.UP[closest_point([self.nozEndUp[0],self.nozEndUp[1],0],self.UP)]
		self.nozStrtD = self.DOWN[-1]
		self.nozEndD = Rotate_about_center(self.nozEndUp,-self._pitch,[0,0])

	#-----------------------------------------------------------------------------------------------#
	#
	#
	#						
	def LogSpiral(self):
		"Function to generate the Logorithmic Spiral"
		b = Rotate_about_center([-self._ActualRout,0],-self._pitch,[0,0])
		vec = [b[0] + self._TEminT*math.cos(math.radians(-90+self._flowAngle)), b[1]+self._TEminT*math.sin(math.radians(-90+self._flowAngle))]
		if self._DEBUG:
			plt.plot(vec[0],vec[1],'.r')
		dAng = abs(math.degrees(math.atan(vec[1]/vec[0]))+self._pitch)
		A = math.sqrt(vec[0]**2+vec[1]**2)
		B = math.atan(math.radians(90-self._flowAngle))
		t = np.linspace(0, math.pi, 4000)
		flag=False
		i = 0
		LogSpiral = []
		while flag!=True:
			logxy = [-1*A*math.exp(B*t[i])*math.cos(t[i]),-1*A*math.exp(B*t[i])*math.sin(t[i])]
			d = Rotate_about_center(logxy,-self._pitch+(dAng),[0,0])
			if self._DEBUG:
				plt.plot(d[0],d[1],'.m')
				self._log_d=d
			slope = (0-d[1])/(-self._ActualRout-d[0])
			i+=1
			LogSpiral.append(d)
			#print(math.degrees(math.atan(slope)) , self._flowAngle - 90,math.degrees(math.atan(slope)) + 90 - self._flowAngle)
			if abs(math.degrees(math.atan(slope)) + 90 - self._flowAngle)<1:
				flag= True
		self._LogSpiral = np.array(LogSpiral)
		self._SecretScale = ((d[0]+self._ActualRout)**2+d[1]**2)**0.5
	
	def Build_Diverging(self):
		"DELETE"
		self._computeUPandDOWNsection()
		self._cart2Radial()
		self._CurvedNozzle()
		self._ExtendSuctionSide()

	#-----------------------------------------------------------------------------------------------#
	#
	#
	#			
	def _AxialNozzle(self):
		# self.DOWN = deepcopy(self._nozzleDown)
		dy = self._nozzleDown[-1,1]
		for i in range(len(self._nozzleDown)): self._nozzleDown[i,1] -=dy
		for i in range(len(self._nozzleUp)): self._nozzleUp[i,1] -=dy
		for i in range(len(self._nozzleCenter)): self._nozzleCenter[i,1] -=dy
		for i in range(len(self._nozzleDown)):
			self._nozzleDown[i,1]+=(self._pitch+self._TEminT);
		A=Intersection(self._nozzleDown[-1],self._nozzleUp[-1],[math.tan(math.pi/2+math.radians(self._flowAngle)),math.tan(math.radians(self._flowAngle))])
		center = [(A[0]+self._nozzleDown[-1,0])/2,(A[1]+self._nozzleDown[-1,1])/2,0]
		vec = [(self._nozzleUp[-1,0]-A[0])/50, (self._nozzleUp[-1,1]-A[1])/50]
		for i in range(0,50):
			self._nozzleUp= np.append(self._nozzleUp,[[self._nozzleUp[-1,0]-vec[0],self._nozzleUp[-1,1]-vec[1],0.0]],axis=0)	
		# self._nozzleUp= np.append(self._nozzleUp,[[A[0],A[1],0]],axis=0)
		if self._DEBUG:
			plt.plot(self._nozzleUp[:,0],self._nozzleUp[:,1])
			plt.plot(self._nozzleDown[:,0],self._nozzleDown[:,1])
			plt.axis('equal')
			plt.show()
		self.DOWN = deepcopy(self._nozzleDown)
		self.UP = self._nozzleUp
	#-----------------------------------------------------------------------------------------------#
	#
	#
	#
	def BladeBuild(self):
		self._NURBSctrlpts()
		self._connectNozzleUPnDOWN()
		self._FilletTrailingEdge()
		self.Vane=np.append(self._nozzleUp, self._nozzleDown[-1:0:-1], axis=0)
		self.Vane=np.append(self.Vane, self.Vane[0:1],axis=0)
		self._i_LE = closest_point([self._LEpt[0],self._LEpt[1],0],self.Vane)
		self.LEcoords = self.Vane[self._i_LE]
		
            
#-----------------------------------------------------------------------------------------------#
#
#
#
	def _CurvedNozzle(self):

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
			
		"Writing the location of the curved nozzle"
		self.DOWN=np.column_stack((DOWNx,DOWNy))
		self.UP=np.column_stack((UPx,UPy))
		self.DOWN=np.column_stack((self.DOWN,np.zeros(len(self.DOWN))))
		self.UP=np.column_stack((self.UP,np.zeros(len(self.UP))))
		self.DOWN=np.delete(self.DOWN,0,axis=0)
		self.UP=np.delete(self.UP,0,axis=0)
	
	def _ExtendSuctionSide(self):
		"Not in use -- do not delete"
		"Find angles"  
		ang_1=math.atan(self._nozzleUp[-1,1]/self._nozzleUp[-1,0])
		ang_2=math.atan(self._nozzleDown[-1,1]/self._nozzleDown[-1,0])
		
		trans = np.array([[math.cos(math.radians(self._flowAngle)),-math.sin(math.radians(self._flowAngle))],[math.sin(math.radians(self._flowAngle)),math.cos(math.radians(self._flowAngle))]])
		
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
			else: fac=0.0001
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

		# self.thrtStrtUp = self.UP[0]
		# self.thrtEndUp = self._curveNozcenter[0]
		# self.thrtStrtD = self.DOWN[0]
		# self.thrtEndD = Rotate_about_center(self.thrtEndUp,-self._pitch,[0,0])

		# self.nozEndUp = self._curveNozcenter[-1]
		# self.nozStrtUp = self.UP[closest_point([self.nozEndUp[0],self.nozEndUp[1],0],self.UP)]
		# self.nozStrtD = self.DOWN[-1]
		# self.nozEndD = Rotate_about_center(self.nozEndUp,-self._pitch,[0,0])

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
		print(self.TE_thk)
		#if (self._radiusOut<self._RealRout): d_t = 1
		#self.d_t = d+d_t

#-----------------------------------------------------------------------------------------------#
#
#
#
	def _NURBSctrlpts(self):

		# plt.plot(self.UP[:,0],self.UP[:,1],'.k')
		# plt.plot(self.DOWN[:,0],self.DOWN[:,1],'.r')
		# plt.axis('equal')
		# plt.show()
		
		pn0 = [self.UP[0,0],self.UP[0,1]]
		pn16 = [self.DOWN[0,0],self.DOWN[0,1]]

		#VecCenStart = [pn0[0]-self._nozzleCenter[0,0],pn0[1]-self._nozzleCenter[0,1]]
		#plt.plot(self._nozzleCenter[0,0],self._nozzleCenter[0,1],'oy')

		#pn16_adj = Rotate_about_center(pn16,self._pitch,[0,0])
		#pn16_adj = [pn16[0],pn16[1]]
		#newP = [0.5*(pn0[0]+pn16_adj[0]), 0.5*(pn0[1]+pn16_adj[1])]
		#newP = [0.5*(pn0[0]+pn16_adj[0]), 0.5*(pn0[1]+pn16_adj[1])]
		#VecCenStart = [pn0[0]-self._nozzleCenter[0,0],pn0[1]-self._nozzleCenter[0,1]]
		VecCenStart = [self.UP[0,0]-self.UP[1,0], self.UP[0,1]-self.UP[1,1]]
		VecCenStart_Pen = UnitVec(VecCenStart[1],-VecCenStart[0])

		#VecCenStart = [pn0[0]-newP[0],pn0[1]-newP[1]]
		#VecCenStart = [pn0[0],pn0[1]]
		#plt.plot(newP[0], newP[1],'*m')
		#plt.plot(pn16_adj[0], pn16_adj[1],'*c')
		#plt.plot(self._curveNozcenter[0,0],self._curveNozcenter[0,1],'og')
		#pdb.set_trace()
		#plt.plot(self._log_d[0],self._log_d[1],'ob')
		#VecCenStart = UnitVec(VecCenStart[0],VecCenStart[1])
		#kernOriUp = [self._nozzleUp[0,0]+self._rho_t*VecCenStart[0],self._nozzleUp[0,1]+self._rho_t*VecCenStart[1]]
		#kernOriUp = [pn0[0]+self._rho_t*VecCenStart[0],pn0[1]+self._rho_t*VecCenStart[1]]
		kernOriUp = [pn0[0]+self._rho_t*VecCenStart_Pen[0],pn0[1]+self._rho_t*VecCenStart_Pen[1]]

		if self._RADIAL:
			eps = 1e-3
		else:
			eps = 1e-2#5e-3
		for i, pt in enumerate(self.UP):
			r_pt = ((pt[0]-kernOriUp[0])**2+(pt[1]-kernOriUp[1])**2)**0.5
			if abs(r_pt-self._rho_t)<eps: i_r1 = i
		#"for robustness, could be improved?"
		try: pr1 = [self.UP[i_r1,0],self.UP[i_r1,1]]
		except NameError: raise ValueError('Kernel origin not found, possible issue with nozzle coordinates')

		#'New algorithm to find kernel radius'

		#pdb.set_trace()
		
		vec_ou = [self.UP[0,0]-kernOriUp[0],self.UP[0,1]-kernOriUp[1]]
		vec_ou = UnitVec(vec_ou[0],vec_ou[1])

		vec_or1 =[self.UP[i_r1,0]-kernOriUp[0],self.UP[i_r1,1]-kernOriUp[1]]
		vec_or1 = UnitVec(vec_or1[0],vec_or1[1])

		ang_u = int(0.5*math.degrees(math.acos(vec_ou[0]*vec_or1[0]+vec_ou[1]*vec_or1[1])))


		#pdb.set_trace()

		'Enlarge existing part of nozzle using continuity circle'
		ContCircUp = np.zeros((ang_u,3))
		# the i was -i earlier as per Jozefs Algorithm
		if self._AXIAL:
			fac = -1
		else:
			fac = -1
		#Extremely Small Bug Hence the above if statement TODO: Please Fix NITISH
		for i in range(ang_u): ContCircUp[i] = Rotate_about_center(pn0,fac*i,kernOriUp)
		#pdb.set_trace()
		'Find 2 control points using circle to assure C^2 continuity with NURBS'
		pn1 = [ContCircUp[int(ang_u/2),0],ContCircUp[int(ang_u/2),1]]
		pn2 = [ContCircUp[-1,0],ContCircUp[-1,1]]
		pc1 = [ContCircUp[int(ang_u/4),0],ContCircUp[int(ang_u/4),1]]
		pc2 = [ContCircUp[int(3*ang_u/4),0],ContCircUp[int(3*ang_u/4),1]]

		'Find 2 control points for down nozzle mirroring circle'
		'NITISH : Redefined the algo for refnozctr0 and 1'
		if self._RADIAL:
			temp_d1 = Rotate_about_center(self.DOWN[0],self._pitch,[0,0])
			temp_d2 = Rotate_about_center(self.DOWN[1],self._pitch,[0,0])
		else:
			temp_d1 = [self.DOWN[0,0], self.DOWN[0,1]-self._pitch]
			temp_d2 = [self.DOWN[1,0], self.DOWN[1,1]-self._pitch]
		
		refnozctr0 = [(self.UP[0,0]+temp_d1[0])/2,(self.UP[0,1]+temp_d1[1])/2]
		refnozctr1 = [(self.UP[1,0]+temp_d2[0])/2,(self.UP[1,1]+temp_d2[1])/2]
		
		#refnozctr0 = [self._nozzleCenter[0,0]+dx,self._nozzleCenter[0,1]+dy]
		#refnozctr1 = [self._nozzleCenter[1,0]+dx,self._nozzleCenter[1,1]+dy]

		#refnozctr0 = [self._nozzleCenter[0,0],self._nozzleCenter[0,1]]
		#refnozctr1 = [self._nozzleCenter[1,0],self._nozzleCenter[1,1]]

		pn14 = mirror_point(pn2,refnozctr0,refnozctr1)
		pn15 = mirror_point(pn1,refnozctr0,refnozctr1)
		pc3 = mirror_point(pc2,refnozctr0,refnozctr1)
		pc4 = mirror_point(pc1,refnozctr0,refnozctr1)

		#pn14 = mirror_point([ContCircUp[int(ang_u/2),0],ContCircUp[int(ang_u/2),1]],refnozctr0,refnozctr1)
		#pn14 = mirror_point(pn2,refnozctr0,refnozctr1)
		#pn15 = mirror_point(pn1,refnozctr0,refnozctr1)
		#pc3 = mirror_point(pc2,refnozctr0,refnozctr1)
		#pc4 = mirror_point(pc1,refnozctr0,refnozctr1)


		
		if self._RADIAL:		
			pn14 = Rotate_about_center(pn14,-self._pitch,[0,0])
			pn15 = Rotate_about_center(pn15,-self._pitch,[0,0])
			pc3 = Rotate_about_center(pc3,-self._pitch,[0,0])
			pc4 = Rotate_about_center(pc4,-self._pitch,[0,0])
		elif self._AXIAL:
			pn14[1]=pn14[1]+self._pitch
			pn15[1]=pn15[1]+self._pitch
			pc3[1]=pc3[1]+self._pitch
			pc4[1]=pc4[1]+self._pitch 
		else:
			sys.exit()
		pn14 = [pn14[0],pn14[1]]
		pn15 = [pn15[0],pn15[1]]
		pc3 = [pc3[0], pc3[1]]
		pc4 = [pc4[0], pc4[1]]
		#pn15 = mirror_point([ContCircUp[int(ang_u/4),0],ContCircUp[int(ang_u/8),1]],refnozctr0,refnozctr1)
		#pc3 = mirror_point([ContCircUp[int(ang_u/4),0],ContCircUp[int(ang_u/4),1]],refnozctr0,refnozctr1)
		#pc4 = mirror_point([ContCircUp[int(ang_u/16),0],ContCircUp[int(ang_u/16),1]],refnozctr0,refnozctr1)

		'Find point in nozzle up perpendicular to start nozzle down'
		i_r2 = closest_point(self.DOWN[0],self.UP)
		pr2 = [self.UP[i_r2,0],self.UP[i_r2,1]]
		
		
		
		"Find blade leading edge radius"
		'Jozef: Changed this to adapt to new blade'
		#pNC1 = [self.UP[1,0],self.UP[1,1]]
		pNC0 = [self.UP[0,0],self.UP[0,1]]
		pNC1 = kernOriUp
		if self._RADIAL:
			pNC0adj = [self.DOWN[0,0],self.DOWN[0,1]]
			pNC1adj = [self.DOWN[1,0],self.DOWN[1,1]]
		else:
			pNC0adj = [self.DOWN[0,0],self.DOWN[0,1]]
			pNC1adj = [self.DOWN[-1,0],self.DOWN[-1,1]]


		#pNC0 = refnozctr0 #self._nozzleCenter[0]
		#pNC1 = refnozctr1 #self._nozzleCenter[1]
		#VecLE_Down = [self.DOWN[0,0]-self.DOWN[1,0], self.DOWN[0,1]-self.DOWN[1,1]]
		#VecLE_Down = UnitVec(VecLE_Down[0], VecLE_Down[1])
		m1,k1 = Line2Pts(pNC0,pNC1)
		m2,k2 = Line2Pts(pNC0adj,pNC1adj)
		x_i = (k2-k1)/(m1-m2)
		y_i = m1*x_i+k1
		#m3 = -1/m1
		#k3 = pNC0[0]*(m1-m3)+k1
		#x_i = (k3-k2)/(m2-m3)
		#y_i = m2*x_i+k2
		pn7 = [x_i,y_i]

		#if (x_i**2+y_i**2)**0.5>self._ActualRin:
		#	[pn7_1,pn7_2] = InterLineCirc(pn0,kernOriUp,self._ActualRin)
		#	if (pn7_1[0]<pn0[0]): pn7 = pn7_1
		#	else: pn7 = pn7_2
		#	pn7 = [pn7[0],pn7[1]]



		#if self._RADIAL:
		#	pNC0adj = Rotate_about_center(pNC0,-self._pitch,[0,0])
		#	pNC1adj = Rotate_about_center(pNC1,-self._pitch,[0,0])
		#else:
		#	pNC0adj = [pNC0[0],pNC0[1] + self._pitch]
		#	pNC1adj = [pNC1[0],pNC1[1] + self._pitch]
		#if self._DEBUG:
		#	plt.plot([pNC0[0],pNC1[0]],[pNC0[1],pNC1[1]],'r')
		#	plt.plot([pNC0adj[0],pNC1adj[0]],[pNC0adj[1],pNC1adj[1]],'b')
		#m1,k1 = Line2Pts(pNC0,pNC1)
		#m2,k2 = Line2Pts(pNC0adj,pNC1adj)
		#m3 = -1/m1
		#k3 = pNC0[0]*(m1-m3)+k1
		#x_i = (k3-k2)/(m2-m3)
		#y_i = m2*x_i+k2
		#pn7 = [x_i,y_i]

		self._R_LE = (x_i**2+y_i**2)**0.5
		#if self._R_LE>self._ActualRin: self._R_LE = self._ActualRin
		
		#if self._RADIAL:
		#	'Find first point and origin of convergent circle'
		#	[pn7_1,pn7_2] = InterLineCirc(self._nozzleCenter[0],pn0,self._R_LE)
		#	if (pn7_1[0]<pn0[0]): pn7 = pn7_1
		#	else: pn7 = pn7_2
		#	pn7 = [pn7[0],pn7[1]]
		#if self._RADIAL:
		#	'Find first point and origin of convergent circle'
		#	[pn7_1,pn7_2] = InterLineCirc(newP,self.UP[0],self._R_LE)
		#	if (pn7_1[0]<pn0[0]): pn7 = pn7_1
		#	else: pn7 = pn7_2
		#	pn7 = [pn7[0],pn7[1]]

		#VecLE_UP = 
		#VecCenStart = [self.UP[0,0]-self.UP[1,0], self.UP[0,1]-self.UP[1,1]]
		#VecCenStart_Pen = UnitVec(VecCenStart[1],-VecCenStart[0])


		
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
		if self._DEBUG:
			plt.plot(pn0[0],pn0[1],'*r')
			plt.plot(pn1[0],pn1[1],'*b')
			plt.plot(pn2[0],pn2[1],'*b')	
			plt.plot(pn3[0],pn3[1],'*b')
			plt.plot(pn4[0],pn4[1],'*b')
			plt.plot(pn5[0],pn5[1],'*b')
			plt.plot(pn6[0],pn6[1],'*b')
			plt.plot(pn7[0],pn7[1],'*m')
			plt.plot(pn8[0],pn8[1],'*b')
			plt.plot(pn9[0],pn9[1],'*y')
			plt.plot(pn10[0],pn10[1],'*g')
			plt.plot(pn11[0],pn11[1],'*y')
			plt.plot(pn12[0],pn12[1],'*b')
			plt.plot(pn13[0],pn13[1],'*b')
			plt.plot(pn14[0],pn14[1],'*b')
			plt.plot(pn15[0],pn15[1],'*b')
			plt.plot(pn16[0],pn16[1],'*r')
	
			plt.plot(pc1[0],pc1[1],'+y')
			plt.plot(pc2[0],pc2[1],'+b')
			plt.plot(pc3[0],pc3[1],'+k')
			plt.plot(pc4[0],pc4[1],'+r')
	
			plt.plot(pa1[0],pa1[1],'+c')
			plt.plot(pa2[0],pa2[1],'+c')
			plt.plot(pa3[0],pa3[1],'+c')
			plt.plot(pa4[0],pa4[1],'+c')
			plt.plot(pa5[0],pa5[1],'+c')
			plt.plot(pa6[0],pa6[1],'+c')

			plt.plot(pr1[0],pr1[1],'oc')
			plt.plot(kernOriUp[0],kernOriUp[1],'ok')

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
	def _computeUPandDOWNsection(self):
		"Computes lower and upper section of the nozzle"
		self._nozzleUp   = deepcopy(self._Moc_xyz[1:-2:])
		self._nozzleDown = deepcopy(self._Moc_xyz[1:-2:])
		self._nozzleCenter=np.column_stack((self._Moc_xyz[1:-2:,0],np.zeros(len(self._Moc_xyz[1:-2:]))))
		self._nozzleCenter=np.column_stack((self._nozzleCenter,np.zeros(len(self._Moc_xyz[1:-2:]))))
		for i in range(len(self._nozzleDown)):self._nozzleDown[i,1] *=-1.0

	def _newcart2Radial(self):
		dist = self._Moc_xyz[-2,1]*2
		for i in range(len(self._nozzleUp)):   self._nozzleUp[i,:] = Rotate_about_center(self._nozzleUp[i],self._flowAngle,[0,0])
		for i in range(len(self._nozzleDown)): self._nozzleDown[i,:] = Rotate_about_center(self._nozzleDown[i],self._flowAngle,[0,0])
		for i in range(len(self._nozzleCenter)): self._nozzleCenter[i,:] = Rotate_about_center(self._nozzleCenter[i],self._flowAngle,[0,0])
		
		if self._RADIAL:
			Scale = self._SecretScale/dist
		elif self._AXIAL:
			UPx=[self._nozzleUp[-1,0]]
			UPy=[self._nozzleUp[-1,1]]
			DOWNx=[self._nozzleDown[-1,0]]
			DOWNy=[self._nozzleDown[-1,1]]
			dpitch=self._TEminT/math.cos(math.radians(self._flowAngle))
			Scale = (self._pitch-dpitch)*math.cos(math.radians(self._flowAngle))/math.sqrt((UPx[0]-DOWNx[0])**2+(UPy[0]-DOWNy[0])**2)
		
		for i in range(len(self._nozzleDown)):
			self._nozzleDown[i,0] *= Scale 
			self._nozzleDown[i,1] *= Scale 
         
		for i in range(len(self._nozzleUp)):
			self._nozzleUp[i,0] *= Scale
			self._nozzleUp[i,1] *= Scale
        
		for i in range(len(self._nozzleCenter)):
			self._nozzleCenter[i,0] *= Scale
			self._nozzleCenter[i,1] *= Scale

		self._rho_t *= Scale
		self.scaledHalfThroat = Scale

		if self._RADIAL:		
			'Translate the Nozzle'
			dx = self._nozzleDown[-1,0] + self._ActualRout + self._TEminT
			dy = self._nozzleDown[-1,1]
		elif self._AXIAL:
			dx = self._nozzleDown[-1,0] + self._ActualRout
			dy = 0
			
		for i in range(len(self._nozzleDown)):
			self._nozzleDown[i,0] -=dx 
			self._nozzleDown[i,1] -=dy 
         
		for i in range(len(self._nozzleUp)):
			self._nozzleUp[i,0] -=dx
			self._nozzleUp[i,1] -=dy
        
		for i in range(len(self._nozzleCenter)):
			self._nozzleCenter[i,0] -=dx
			self._nozzleCenter[i,1] -=dy
		
		if self._DEBUG:
			plt.plot(self._nozzleUp[:,0], self._nozzleUp[:,1])
			plt.plot(self._nozzleDown[:,0],self._nozzleDown[:,1])
			plt.plot(self._nozzleCenter[:,0],self._nozzleCenter[:,1])
			plt.plot(self._nozzleCenter[0,0],self._nozzleCenter[0,1],'*r')
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
		self.TE_thk =  2*math.sqrt((Up_Point[0]-center[0])**2+(Up_Point[1]-center[1])**2)
		print('Obtained TE Thickenss :: %f'%(self.TE_thk))
		#self._lastPUP = self.UP[-1]
	#-----------------------------------------------------------------------------------------------#
	#
	#
	#
	# Meshing Boundary curves
	def MeshBoundariesBlade(self):

		"Jozef: added a delta R for the spacing of axial boundary limits"
		if not self._RADIAL:
			#pdb.set_trace()
			fac = 0.1
			dR = fac*self._R_LE

		i_TE = closest_point([self.TEpt[0],self.TEpt[1],0],self.Vane)

		self.LEcoords = self.Vane[self._i_LE]

		self.SUCTION = self.Vane[i_TE:self._i_LE+1]

		PRES1 = self.Vane[self._i_LE:-1]
		PRES2 = self.Vane[0:i_TE+1]

		self.PRESSURE = np.append(PRES1,PRES2,axis=0)

		#"Find Limit Points"
		p0 = self.thrtEndUp
		p1 = self.nozEndUp

		vec_LE_NozCtr = UnitVec(p0[0]-self.Vane[self._i_LE,0], p0[1]-self.Vane[self._i_LE,1])
		vec_ctr = UnitVec(p1[0]-p0[0],p1[1]-p0[1])
		ref_ang = math.degrees(math.acos(vec_LE_NozCtr[0]*vec_ctr[0]+vec_LE_NozCtr[1]*vec_ctr[1]))

		#if self._RADIAL:
		p2 = Rotate_about_center(self.Vane[self._i_LE],ref_ang/2,p0)
		#else:


		p5 = [p0[0]-vec_ctr[0]*self._distMshBndry,p0[1]-vec_ctr[1]*self._distMshBndry]

		if self._RADIAL:
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
		else:
			p6 = [-self._R_LE,p2[1]]
			p3 = [p6[0]-dR,p6[1]]
			p4 = [p1[0]+dR,p1[1]]


		if imprtplt:
			plt.plot(p0[0],p0[1],'*r')
			plt.plot(p1[0],p1[1],'*g')
			plt.plot(p2[0],p2[1],'*b')
			plt.plot(p3[0],p3[1],'*y')
			plt.plot(p4[0],p4[1],'*k')
			plt.plot(p5[0],p5[1],'*m')
			plt.plot(p6[0],p6[1],'*c')

		self.PERIO1 = np.zeros(((7,3)))
		self.PERIO1[0] = [p3[0],p3[1],0]
		self.PERIO1[1] = [p6[0],p6[1],0]
		self.PERIO1[2] = [p2[0],p2[1],0]
		self.PERIO1[3] = [p5[0],p5[1],0]
		self.PERIO1[4] = [p0[0],p0[1],0]
		self.PERIO1[5] = [p1[0],p1[1],0]
		self.PERIO1[6] = [p4[0],p4[1],0]
		
		if imprtplt:
			self.plotStatorStage()
			plt.axis('equal')

		if self._AXIAL:
			for i in range(len(self.PRESSURE)): self.PRESSURE[i,1] -=self._pitch
			for i in range(len(self.SUCTION)): self.SUCTION[i,1] -=self._pitch

		self.PERIO2 = np.zeros(((len(self.PERIO1),3)))
		if self._RADIAL:
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
		else:
			for i in range(len(self.PERIO2)): self.PERIO2[i] = [self.PERIO1[i,0], self.PERIO1[i,1] - (self._pitch), 0.0]
			self.INLET = np.zeros(((11,3)))
			self.OUTLET = np.zeros(((11,3)))
			for i in range(int(11)): 
				self.INLET[i] = [p3[0], p3[1]- i*(self._pitch/10), 0.0]
				self.OUTLET[i] = [p4[0], p4[1]- i*(self._pitch/10), 0.0]
			self.BladeRout = (self.TEpt[0]**2+self.TEpt[1]**2)**0.5
			self.pitchLen = 2*math.pi*(self._pitch/360)*self.BladeRout
			self.chordLen = ((self.LEcoords[0]-self.TEpt[0])**2+(self.LEcoords[1]-self.TEpt[1])**2)**0.5
		
	#---------------------------------------------------------------------------------------------#
	#
	# Plots the complete stator stage including all the vanes
	#
	def plotStatorStage(self):
		#nBlades *= nBlades
		nBlades = int(self.IN['nBldPlt'])
		pltVane = np.empty((len(self.Vane),3))
		#nBlades = 5
		if (nBlades%2 == 0): rng = np.linspace(-nBlades/2+1,nBlades/2,nBlades)	
		else: rng = np.linspace(math.ceil(-nBlades/2),math.floor(nBlades/2),nBlades)
		#pdb.set_trace()
		Info = ("|| Flow_Angle = %0.1f || pitch = %0.1f || TE_Thick = %0.1e ||"%(self._flowAngle,self._pitch, self.TE_thk))
		for i in rng:
			if self._RADIAL:
				for j in range(len(self.Vane)): pltVane[j,:] = Rotate_about_center(self.Vane[j,:],i*self._pitch,[0,0])
			else:
				for j in range(len(self.Vane)): pltVane[j,:] = [self.Vane[j,0],self.Vane[j,1]-i*self._pitch,0.0]
			plt.plot(pltVane[:,0],pltVane[:,1],'k')
		plt.title(Info)
		plt.grid('on')
		plt.xlabel('x-Axis')
		plt.ylabel('y-Axis')