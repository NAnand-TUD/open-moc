####################FILE NAME: CLASSES.py#######################
#================================================================
# author: Nitish Anand						|
# 	: Master Student, 					|
#	:Process and Energy Departmemt,				|
#	:TU Delft,						|
#	:The Netherlands					|
#								|
# email : nitish.ug2@gmail.com					|
# 								|
# Description: Defines the class variable for the MAIN code. 	|
#	: All the property variable defined in the MOC is a 	|
#	: object of the CLASS mentioned in this module.		|
# CLASS Function:						|
#	:							|
#================================================================

import time
import sys
import os        
import math as calc
from sympy import Symbol, nsolve, log, ln
from sympy.solvers import solve
#from scipy.optimize import minimize
import pdb
import numpy as np
from src.IO import *
import subprocess
#Jozef Add import for CoolProp
#import fluidprop_py as fp

#
# CLASS to store only points
#

class Points_mini:
	# Dummy Class for Data Passing Only
	def __init__(self,x=0,y=0,u=0,v=0):
		#if Input == None:
		self.x=x
		self.y=y
		self.u=u
		self.v=v
#
# CLASS-> MAIN CLASS, contains 1.ThermoProp()
#			       2. Sauer_Line(real,perfect)
#			       3. Reflex_Zone()
#			       4. OtherSupportingFunction
		
class Points(Points_mini):
	# Function:
	#	a. SauerRealGas()
	#	b. SauerMain() -> ThermoProp Not Needed
	#	c. ThermoProp() -> PRFT , EoS , EoS_TAB, RefProp, CoolProp
	#	d. ReflexLine()
	#	e. PointInPoly()
	# Self Initialized Variables

	# Argument Enabled Initialization <1st argument name will be considered as the InputFile>
	try:INFile = sys.argv[1]
	except:INFile = 'MOC_Config.in'
	IN=ReadUserInput(INFile)
	To = float(IN['To'])
	Po = float(IN['Po'])
	gamma = float(IN['gamma'])
	R = float(IN['R'])/float(IN['M'])
	y_t = float(IN['y_t'])
	n = float(IN['n'])
	rho_t = float(IN['rho_t'])
	delta = 0
	Tc = float(IN['Tc'])
	Pc = float(IN['Pc'])
	Vc = float(IN['Vc'])/float(IN['M'])
	GasEqu = IN['GAS_EQU']
	Ho = float(IN['Ho'])
	ao = calc.sqrt(gamma*R*To)
	Omega=0
	IT_rho=0
	IT_P=0
	IT_c=0
	IT_v=0
	IT_h=0
	IT_T=0
	
	#Initlaizating Fluid for RefProp ONLY
	if (GasEqu=='RefProp'or GasEqu=='RefProp_TAB'):
		FLD = IN['FLDNAME']
		s = fp.fluidprop_py(FLD,'PT',float(IN['Po'])/1e5,float(IN['To'])-273,'entropy')
		Ho = fp.fluidprop_py(FLD,'PT',float(IN['Po'])/1e5,float(IN['To'])-273,'enthalpy')
#---------------------------------------------------------------------------------------------#
	# Printing the entire Class Variable using str(ObjectName)
	def __str__(self):
		print("x = {} \ny = {} \nu = {} \nv = {} \nM = {} \na = {} \nP = {} \nT = {} \nrho = {} \nV = {} \nh = {} \n".format(self.x,self.y,self.u,self.v,self.M,self.a,self.Press,self.Temp,self.rho,self.V,self.h))

#---------------------------------------------------------------------------------------------#
	# Initializing the variable (Optional)[Default=0]
	def __init__(self,x=0,y=0,u=0,v=0):
		Points_mini.__init__(self,x,y,u,v)
		self.M=0
		self.a=0
		self.V=0
		self.Press=0
		self.Temp=0
		self.rho=0
		self.h=0
		self.r = 0

#---------------------------------------------------------------------------------------------#

	# Sauer Analysis with RefProp (Ref: Under Inspection)
	def Sauer_RefTab(self):
		V_flag =95
		f=1
		while abs(self.M-1)>1e-5:
			self.u = V_flag
			self.ThermoProp(self.IN)
			f=self.a**2/2000 + self.h - self.Ho
			print(f)
			print('M',self.M)
			V_flag = self.a - (0.1*f*(abs(f)/1000)**0.5/abs(f))
			print('V_flag',V_flag)

#---------------------------------------------------------------------------------------------#
	# Sauer Analysis with RefProp (Ref: Alberto Code, Contact Matteo Pini for more Info.)[Special CASE for TOL]
	# Gustavo: Do not Review this function (similar function is at the end of the code)
	def Sauer_RealGas1(self):
		IT_rho=fp.fluidprop_py(self.FLD,'hs',self.Ho,self.s,'density')
		IT_T=fp.fluidprop_py(self.FLD,'vs',1/IT_rho,self.s,'temperature')
		for i in range(0,10000):
			IT_T=fp.fluidprop_py(self.FLD,'vs',1/IT_rho,self.s,'temperature')
			IT_c=fp.fluidprop_py(self.FLD,'Tv',IT_T,1/IT_rho,'soundspeed')			
			IT_G=fp.fluidprop_py(self.FLD,'Tv',IT_T,1/IT_rho,'gamma')
			IT_h=fp.fluidprop_py(self.FLD,'Tv',IT_T,1/IT_rho,'enthalpy')
			IT_P=fp.fluidprop_py(self.FLD,'Tv',IT_T,1/IT_rho,'pressure')
			f=IT_c**2/2000 +IT_h-self.Ho
			df=-IT_G*(IT_c**2*IT_rho)
			IT_v=1/IT_rho - f/df
			IT_rho=1/IT_v
			print(f)
			if 1:#(abs(f)<1e-3):
				Throat=np.loadtxt('ThroatProp.out')
				self.a=Throat[0]#140.455#IT_c
				self.Press=Throat[1]#13.9751#IT_P
				self.Temp=Throat[2]#288.327#IT_T
				self.rho=Throat[3]#59.773443394897484#IT_rho
				self.gamma=Throat[4]#0.829105#IT_G
				self.h=Throat[5]#538.028#IT_h
				#print(IT_c,IT_P,IT_T,IT_rho,IT_G,IT_h)
				#pdb.set_trace()
				break
		alpha = calc.sqrt((1+self.delta)/(2*(self.gamma*self.rho_t*self.y_t)))
		epsilon = -alpha*self.y_t*self.gamma/(3+self.delta)
		self.x = -self.gamma*alpha*self.y**2/(3+self.delta)
		self.u = self.a*(1+(alpha*self.x+self.gamma*(alpha*self.y)**2/(1+self.delta)))
		self.v=self.a*((2*self.gamma*alpha**2*self.x*self.y/(1+self.delta))+(2*self.gamma**2*(alpha*self.y)**3)/((1+self.delta)*(3+self.delta)))
		self.x=self.x-epsilon

#---------------------------------------------------------------------------------------------#
	# Main Function for Sauer Analysis
	def SauerMain(self):
		#if (self.GasEqu=='RefProp'or self.GasEqu=='RefProp_TAB'):
		#	self.Sauer_RealGas()
			#self.Sauer_RefTab()
			#print('GAMMA',self.gamma)
			#pdb.set_trace()
		if 1:
			ao = calc.sqrt(self.gamma*self.R*self.To)
			epsilon = -1*(self.y_t/(2*(3+self.delta)))*(((self.gamma+1)*(1+self.delta)/(self.rho_t/self.y_t))**0.5);
			a_star = ((2*self.gamma*self.R*self.To)/(self.gamma+1))**0.5;
			alpha = ((1+self.delta)/((self.gamma+1)*self.rho_t*self.y_t))**0.5
			flag = 0		
			x_ = (-(self.gamma+1)*alpha*self.y**2/(2*(flag+1+self.delta)));
			flag = 2
			self.x = (-(self.gamma+1)*alpha*self.y**2/(2*(flag+1+self.delta)))
			u_dash = (alpha*self.x)+((((self.gamma+1)*alpha**2*self.y**2))/(2*(self.delta+1)))
			v_dash = ((((self.gamma+1)*alpha**2*self.y*self.x))/((self.delta+1)))+(((self.gamma+1)**2*alpha**3*self.y**3)/(2*(1+self.delta)*(3+self.delta)))
			self.u = a_star*(1+u_dash)
			v_tilde = a_star*v_dash	
			if (self.GasEqu!='RefProp' and self.GasEqu!='RefProp_TAB'):
				self.a = calc.sqrt(ao**2 - ((self.gamma-1)*self.u**2/(2)))
				self.M = self.u / self.a
				self.Press = self.Po/(1+((self.gamma-1)*self.M**2/(2)))**((self.gamma)/(self.gamma-1))
				self.Temp = self.To/(1+((self.gamma-1)*self.M**2/(2)))
				self.rho = self.Press/(self.R*self.Temp)
				self.x=self.x-epsilon
		self.M = calc.sqrt(self.u**2+self.v**2) / self.a

#---------------------------------------------------------------------------------------------#
	# Thermodynamics Properties Calculator (PRFT, EoS, EoS_TAB,RefProp)
	def ThermoProp(self,IN):
		self.V = calc.sqrt((self.u**2 + self.v**2))
		#PRFT: Solve using PERFECT gas equations [gamma necessary]
		if self.GasEqu == 'PRFT':
			Cp = self.R*self.gamma/(self.gamma-1)
			self.Temp = self.To - (self.V**2/(2*Cp))
			self.Press = self.Po*(self.Temp/self.To)**(self.gamma/(self.gamma-1))
			self.rho = self.Press / (self.R*self.Temp)
			self.a = calc.sqrt((self.ao**2)-((self.gamma-1)*self.V**2/2))
			self.M = self.V/self.a
		#EoS: Solve by directly solving Van der Waals EoS
		elif self.GasEqu == 'EoS':
			eqP = Symbol('eqP')
			eqT = Symbol('eqT')
			Vrho = Symbol('Vrho')
			b1 = self.R*self.Tc/(8*self.Pc)
			a1 = 27*(self.R*self.Tc)**2/(64*self.Pc)      
			Cp = self.gamma*self.R/(self.gamma-1)
			hig = self.Ho - (Cp)*(self.To-eqT)
			h = self.Ho - (0.5*self.V**2)
			VDW_h = - (2*a1/Vrho) + (b1*self.R*eqT/(Vrho-b1)) + hig - h
			VDW = (eqP*Vrho**3) -(eqP*b1+self.R*eqT)*Vrho**2 + a1*Vrho - a1*b1
			VDW_s = (self.gamma*log(eqT/self.To))/(self.gamma-1) - log(eqP/self.Po) + log(Vrho-b1) - log(Vrho) + log(eqP*Vrho/self.R/eqT)      
			try:[Vrho,T,P] = nsolve((VDW,VDW_h,VDW_s),(Vrho,eqT,eqP),(float(IN['VrhoGuess']),float(IN['TGuess']),float(IN['PGuess'])))
			except:pdb.set_trace()
			try:Vrho = float(Vrho.real)
			except:pdb.set_trace()
			try:self.Temp = float(T)
			except:self.Temp = float(T.real)
			try:self.Press = float(P)
			except:self.Press = float(P.real)
			self.rho = 1/float(Vrho)
			D1 = (2-self.gamma)*(Vrho-b1)/(self.gamma-1)
			D2 = self.R*self.Temp/self.Press
			N1 = self.R*T*((1/(Vrho-b1))+(1/Vrho))
			N2 = (P-(a1/Vrho**2) + (2*a1*b1/Vrho**3))*((2-self.gamma)/(self.gamma-1))
			try:a = calc.sqrt(Vrho**2*(N1.real+N2.real)/(D1.real+D2.real))
			except: pdb.set_trace()
			self.a = a
			self.M = self.V/self.a
			IN['VrhoGuess'] = float(Vrho.real)
			IN['TGuess'] = float(T.real)
			IN['PGuess'] = float(P.real)
		# EoS_TAB = Make and solve using van der Waals EoS table
		elif self.GasEqu == 'EoS_TAB':
			eqP = Symbol('eqP')
			eqT = Symbol('eqT')
			eqVrho = Symbol('eqVrho')
			b1 = self.R*self.Tc/(8*self.Pc)
			a1 = 27*(self.R*self.Tc)**2/(64*self.Pc)
			Cp = self.gamma*self.R/(self.gamma-1)
			hig = self.Ho - (Cp)*(self.To-eqT)
			# CREATE TABLE
			if IN['Table_Oper'] == 'Create':
				try:os.remove(os.getcwd()+'/Table.TAB')
				except: print ("Continue: No File to Delete <Table> \n Location:",os.getcwd())
				minV = int(IN['LLim'])
				maxV = int(IN['ULim'])
				dV = float(IN['dV'])
				V1= np.arange(minV,maxV,dV)
				for V in V1:
					self.V = V
					h = self.Ho - (0.5*self.V**2)
					VDW_h = - (2*a1/eqVrho) + (b1*self.R*eqT/(eqVrho-b1)) + hig - h
					VDW = (eqP*eqVrho**3) -(eqP*b1+self.R*eqT)*eqVrho**2 + a1*eqVrho - a1*b1
					VDW_s = (self.gamma*log(eqT/self.To))/(self.gamma-1) - log(eqP/self.Po) + log(eqVrho-b1) - log(eqVrho) + log(eqP*eqVrho/self.R/eqT)      
					try:[Vrho,T,P] = nsolve((VDW,VDW_h,VDW_s),(eqVrho,eqT,eqP),(float(IN['VrhoGuess']),float(IN['TGuess']),float(IN['PGuess'])))
					except:pdb.set_trace()
					try:Vrho = float(Vrho.real)
					except:pdb.set_trace()
					try:self.Temp = float(T)
					except:self.Temp = float(T.real)
					try:self.Press = float(P)    
					except:self.Press = float(P.real)
					self.rho = 1/float(Vrho)
            				###Refer Thesis for derivation###
					D1 = (2-self.gamma)*(Vrho-b1)/(self.gamma-1)
					D2 = self.R*self.Temp/self.Press
					N1 = self.R*T*((1/(Vrho-b1))+(1/Vrho))    
					N2 = (self.Press-(a1/Vrho**2) + (2*a1*b1/Vrho**3))*((2-self.gamma)/(self.gamma-1))        
					a = calc.sqrt(Vrho**2*(N1.real+N2.real)/(D1.real+D2.real))
					self.a = a
					self.M = self.V/self.a
					IN['VrhoGuess'] = float(Vrho.real)
					IN['TGuess'] = float(T.real)
					IN['PGuess'] = float(P.real)
					WriteDataFile('a',self,'Table.TAB','TAB')
				print('***Table Create***')
				sys.exit(-1)
			# CALCULATE FROM TABLE
			if IN['Table_Oper'] == 'Calc':
				TABProp = np.loadtxt('Table.TAB')
				for i in range(len(TABProp)):
					if TABProp[i][0]>self.V:
						self.M =  self.Inter(TABProp[i][0],TABProp[i-1][0],TABProp[i][1],TABProp[i-1][1])
						self.rho=  self.Inter(TABProp[i][0],TABProp[i-1][0],TABProp[i][2],TABProp[i-1][2])
						self.a=  self.Inter(TABProp[i][0],TABProp[i-1][0],TABProp[i][3],TABProp[i-1][3])
						self.Temp=  self.Inter(TABProp[i][0],TABProp[i-1][0],TABProp[i][4],TABProp[i-1][4])
						self.Press=  self.Inter(TABProp[i][0],TABProp[i-1][0],TABProp[i][5],TABProp[i-1][5])
						break
		# RefProp_TAB = Table generated from RefProp
		elif self.GasEqu == 'RefProp_TAB':
			if IN['Table_Oper']=='Create':
				try:os.remove(os.getcwd()+'/Table.TAB')
				except: print ("Continue: No File to Delete <Table> \n Location:",os.getcwd())
				minV = int(IN['LLim'])
				maxV = int(IN['ULim'])
				dV = float(IN['dV'])
				V1= np.arange(minV,maxV,dV)
				print('***Creating File***')
				for V in V1:
					self.V = V
					self.h = self.Ho - (0.5*self.V**2/1000)
					self.Press=fp.fluidprop_py(self.FLD,'hs',self.h,self.s,'pressure')
					self.Temp=fp.fluidprop_py(self.FLD,'hs',self.h,self.s,'temperature')
					self.rho=fp.fluidprop_py(self.FLD,'hs',self.h,self.s,'density')
					self.a=fp.fluidprop_py(self.FLD,'hs',self.h,self.s,'soundspeed')
					self.M=self.V/self.a
					self.gamma =fp.fluidprop_py(self.FLD,'hs',self.h,self.s,'gamma')
					WriteDataFile('a',self,'Table.TAB','TAB')
					#pdb.set_trace()
				print('***Table Create***')
				sys.exit(-1)
			if IN['Table_Oper']=='Calc':
				self.h = self.Ho - (0.5*self.V**2/1000)
				TABProp = np.loadtxt('Table.TAB') #[HARDCODED]
				for i in range(len(TABProp)):
					if TABProp[i][0]>self.V:
						#print(TABProp[i][0],self.V)
						self.M =  self.Inter(TABProp[i][0],TABProp[i-1][0],TABProp[i][1],TABProp[i-1][1])
						self.rho=  self.Inter(TABProp[i][0],TABProp[i-1][0],TABProp[i][2],TABProp[i-1][2])
						self.a=  self.Inter(TABProp[i][0],TABProp[i-1][0],TABProp[i][3],TABProp[i-1][3])
						self.Temp=  self.Inter(TABProp[i][0],TABProp[i-1][0],TABProp[i][4],TABProp[i-1][4])
						self.Press=  self.Inter(TABProp[i][0],TABProp[i-1][0],TABProp[i][5],TABProp[i-1][5])
						self.gamma=  self.Inter(TABProp[i][0],TABProp[i-1][0],TABProp[i][6],TABProp[i-1][6])
						#print(self.Temp,self.rho,self.a,self.gamma,self.M)
						break
		# RefProp = RefProp interpretor via python [Module fluidprop.py]
		elif self.GasEqu == 'RefProp':
			self.h = self.Ho - (0.5*self.V**2/1000)
			self.Press=fp.fluidprop_py(self.FLD,'hs',self.h,self.s,'pressure')
			self.Temp=fp.fluidprop_py(self.FLD,'hs',self.h,self.s,'temperature')
			self.rho=fp.fluidprop_py(self.FLD,'hs',self.h,self.s,'density')
			self.a=fp.fluidprop_py(self.FLD,'hs',self.h,self.s,'soundspeed')
			self.M=self.V/self.a
			#print(self.M)
			self.gamma =fp.fluidprop_py(self.FLD,'hs',self.h,self.s,'gamma')
            
		# CoolProp = Import fluid Properties using CoolProp
		elif self.GasEqu == 'CoolProp':
			#Jozef Add Here#
			print("Jozef will add here")
            		
		elif self.GasEqu == 'LUT':
        	    	#Aryamann Add Here#
			print("Aryamann will add here")
            	
		# ERROR Message
		else:
			print("\n\n\n\t*********INVALID INPUT DATA**********\n\t *GAS_EQU : ise -> isentropic case \n\t *GAS_EQS : EoS -> Equations of State\n\t***************************************\n\n\n")
			sys.exit(1)

#---------------------------------------------------------------------------------------------#
	# Reflex Line position calculator (ReflexZone)
	def ReflexLine(self,n,l):
		RetObj = Points()
		RetObj.x = self.x + 1*n*l*calc.cos(((1*calc.asin(1/self.M))+(calc.atan(self.v/self.u))))
		RetObj.y = self.y + 1*n*l*calc.sin(((1*calc.asin(1/self.M))+(calc.atan(self.v/self.u))))
		RetObj.u = self.u
		RetObj.v = 0
		return(RetObj)

#---------------------------------------------------------------------------------------------#
	# Check if the characteristic point is inside the nozzle [Exclusive for Reflex Zone]
	def PointInPoly(self,poly):
		n = len(poly)
		inside = False
		p1x = poly[0].x
		p1y = poly[0].y
		for i in range(n+1):
			p2x = poly[i % n].x
			p2y = poly[i % n].y
			if self.y > min(p1y,p2y):
				if self.y <= max(p1y,p2y):
					if self.x <= max(p1x,p2x):
						if p1y != p2y:
							xints = (self.y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
						if p1x == p2x or self.x <= xints:
							inside = not inside
			p1x,p1y = p2x,p2y
		return inside

#---------------------------------------------------------------------------------------------#
	# Liner Interpolation function for Table
	def Inter(self,x1,x2,y1,y2):
		if self.GasEqu == 'EoS_TAB':X = self.V
		else: X = self.h
		k=y1+(X-x1)*((y2-y1)/(x2-x1));
		return(k);
#---------------------------------------------------------------------------------------------#
# Sauer Analysis with RefProp (Ref: Alberto Code) [General Case]
# Guardone 2012, supersonic nozzle design code [point of contact: Matteo Pini]
def Sauer_RealGas(P):
	IT_rho=fp.fluidprop_py(P.FLD,'hs',P.Ho,P.s,'density')
	#pdb.set_trace()
	IT_T=fp.fluidprop_py(P.FLD,'vs',1/IT_rho,P.s,'temperature')
	#pdb.set_trace()
	for i in range(0,10000):
		IT_T=fp.fluidprop_py(P.FLD,'vs',1/IT_rho,P.s,'temperature')
		IT_c=fp.fluidprop_py(P.FLD,'Tv',IT_T,1/IT_rho,'soundspeed')			
		IT_G=fp.fluidprop_py(P.FLD,'Tv',IT_T,1/IT_rho,'gamma')
		IT_h=fp.fluidprop_py(P.FLD,'Tv',IT_T,1/IT_rho,'enthalpy')
		IT_P=fp.fluidprop_py(P.FLD,'Tv',IT_T,1/IT_rho,'pressure')
		f=IT_c**2/2000 +IT_h-P.Ho
		if i==0: f_tot=f
		df=-IT_G*(IT_c**2*IT_rho)
		#pdb.set_trace()
		IT_v=1/IT_rho - f/df
		IT_rho=1/IT_v
		#print(f)
		printProgress(f_tot-f,f_tot,'REAL GAS: Sauer Iteration')
		if (abs(f)<1e-3):
			P.a=IT_c
			P.Press=IT_P
			P.Temp=IT_T
			P.rho=IT_rho
			P.gamma=IT_G
			P.h=IT_h
			print("\n\n*****Iteation Summary****\nsos\t:%f\nPress\t:%f\nTemp\t:%f\nRho\t:%f\ngamma\t:%f\nh\t:%f\n*****Summary End*****\n"%(IT_c,IT_P,IT_T,IT_rho,IT_G,IT_h))
			print("Writing Throat Prop. File\n\n")
			f=open('ThroatProp.out','w')
			f.write("%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t"%(IT_c,IT_P,IT_T,IT_rho,IT_G,IT_h))
			f.close()
			break
	alpha = calc.sqrt((1+P.delta)/(2*(P.gamma*P.rho_t*P.y_t)))
	epsilon = -alpha*P.y_t*P.gamma/(3+P.delta)
	P.x = -P.gamma*alpha*P.y**2/(3+P.delta)
	P.u = P.a*(1+(alpha*P.x+P.gamma*(alpha*P.y)**2/(1+P.delta)))
	P.v=P.a*((2*P.gamma*alpha**2*P.x*P.y/(1+P.delta))+(2*P.gamma**2*(alpha*P.y)**3)/((1+P.delta)*(3+P.delta)))
	P.x=P.x-epsilon
	return P
	#pdb.set_trace()

##
## END
##

=======
####################FILE NAME: CLASSES.py#######################
#================================================================
# author: Nitish Anand	& Jozef Stuijt							|
# 	:Master Student, 											|
#	:Process and Energy Departmemt,								|
#	:TU Delft,													|
#	:The Netherlands											|
#																|
# email : nitish.ug2@gmail.com									|
# 																|
# Description: Defines the class variable for the MAIN code. 	|
#	: All the property variable defined in the MOC is a 		|
#	: object of the CLASS mentioned in this module.				|
# CLASS Function:												|
#	:															|
#================================================================

# Imports
import time
import sys
import os        
import math as calc
from sympy import Symbol, nsolve, log, ln
from sympy.solvers import solve
import pdb
import numpy as np
from srcMOC.IO import *
import subprocess
import CoolProp.CoolProp as cp

#
# CLASS to store only points
#

class Points_mini:
	# Dummy Class for Data Passing Only
	def __init__(self,x=0,y=0,u=0,v=0):
		#if Input == None:
		self.x=x
		self.y=y
		self.u=u
		self.v=v
#
# CLASS-> MAIN CLASS, contains 
#				   1.ThermoProp()
#			       2. Sauer_Line(real,perfect)
#			       3. Reflex_Zone()
#			       4. OtherSupportingFunction
		
class Points(Points_mini):
	# Function:
	#	a. SauerRealGas()
	#	b. SauerMain() -> ThermoProp Not Needed
	#	c. ThermoProp() -> PRFT , EoS , EoS_TAB, RefProp, CoolProp
	#	d. ReflexLine()
	#	e. PointInPoly()
	# Self Initialized Variables

	# Argument Enabled Initialization <1st argument name will be considered as the InputFile>
	DIR = os.getcwd()+'/'
	try:INFile = DIR+sys.argv[1]
	except:INFile = DIR+'MOC_Config.in'
	IN     = ReadUserInput(INFile)
	To     = float(IN['To'])
	Po     = float(IN['Po'])
	gamma  = float(IN['gamma'])
	R      = float(IN['R'])/float(IN['M'])
	y_t    = float(IN['y_t'])
	n      = float(IN['n'])
	rho_t  = float(IN['rho_t'])
	delta  = 0
	Tc     = float(IN['Tc'])
	Pc     = float(IN['Pc'])
	Vc     = float(IN['Vc'])/float(IN['M'])
	GasEqu = IN['GAS_EQU']
	Ho     = float(IN['Ho'])
	ao     = calc.sqrt(gamma*R*To)
	Omega  = 0
	IT_rho = 0
	IT_P   = 0
	IT_c   = 0
	IT_v   = 0
	IT_h   = 0
	IT_T   = 0

	# Added by Jozef to check throat file
	THROAT = IN['THROAT_PROP']
	if THROAT=='ITER': pass
	else:
		try: THROAT_FILE = np.loadtxt(THROAT)
		except: 
			raise IOError('Invalid Input: THROAT_PROP')
			sys.exit()
	
	#Initlaizating Fluid for RefProp or CoolProp ONLY
	if (GasEqu=='RefProp' or GasEqu=='RefProp_TAB'):
		FLD = IN['FLDNAME']
		s   = fp.fluidprop_py(FLD,'PT',float(IN['Po'])/1e5,float(IN['To'])-273,'entropy')
		Ho  = fp.fluidprop_py(FLD,'PT',float(IN['Po'])/1e5,float(IN['To'])-273,'enthalpy')
	elif (GasEqu == 'CoolProp'):
		FLD = IN['FLDNAME']
		s   = cp.PropsSI('S','P',float(IN['Po']),'T',float(IN['To']),FLD)/1000
		Ho  = cp.PropsSI('H','P',float(IN['Po']),'T',float(IN['To']),FLD)/1000
	#else:
	#	raise IOError('Invalid Input: GasEqu')
	#	sys.exit()
#---------------------------------------------------------------------------------------------#
	# Printing the entire Class Variable using str(ObjectName)
	def __str__(self):
		print("x = {} \ny = {} \nu = {} \nv = {} \nM = {} \na = {} \nP = {} \nT = {} \nrho = {} \nV = {} \nh = {} \n".format(self.x,self.y,self.u,self.v,self.M,self.a,self.Press,self.Temp,self.rho,self.V,self.h))

#---------------------------------------------------------------------------------------------#
	# Initializing the variable (Optional)[Default=0]
	def __init__(self,x=0,y=0,u=0,v=0):
		Points_mini.__init__(self,x,y,u,v)
		self.M     = 0
		self.a     = 0
		self.V     = 0
		self.Press = 0
		self.Temp  = 0
		self.rho   = 0
		self.h     = 0
		self.r     = 0

#---------------------------------------------------------------------------------------------#

	# Sauer Analysis with RefProp (Ref: Under Inspection)
	def Sauer_RefTab(self):
		V_flag =95
		f=1
		while abs(self.M-1)>1e-5:
			self.u = V_flag
			self.ThermoProp(self.IN)
			f=self.a**2/2000 + self.h - self.Ho
			print(f)
			print('M',self.M)
			V_flag = self.a - (0.1*f*(abs(f)/1000)**0.5/abs(f))
			print('V_flag',V_flag)

#---------------------------------------------------------------------------------------------#
	# Sauer Analysis with RefProp (Ref: Alberto Code, Contact Matteo Pini for more Info.)[Special CASE for TOL]
	# Gustavo: Do not Review this function (similar function is at the end of the code)
	# Jozef: Added this code to the Sauer function
	'''def Sauer_RealGas1(self):
		if (self.GasEqu=='RefProp' or self.GasEqu=='RefProp_TAB'):
			IT_rho=fp.fluidprop_py(self.FLD,'hs',self.Ho,self.s,'density')
			IT_T=fp.fluidprop_py(self.FLD,'vs',1/IT_rho,self.s,'temperature')
			for i in range(0,10000):
				IT_T=fp.fluidprop_py(self.FLD,'vs',1/IT_rho,self.s,'temperature')
				IT_c=fp.fluidprop_py(self.FLD,'Tv',IT_T,1/IT_rho,'soundspeed')			
				IT_G=fp.fluidprop_py(self.FLD,'Tv',IT_T,1/IT_rho,'gamma')
				IT_h=fp.fluidprop_py(self.FLD,'Tv',IT_T,1/IT_rho,'enthalpy')
				IT_P=fp.fluidprop_py(self.FLD,'Tv',IT_T,1/IT_rho,'pressure')
				f=IT_c**2/2000 +IT_h-self.Ho
				df=-IT_G*(IT_c**2*IT_rho)
				IT_v=1/IT_rho - f/df
				IT_rho=1/IT_v
				print(f)
				if 1:#(abs(f)<1e-3):
					Throat=np.loadtxt('ThroatProp.out')
					self.a=Throat[0]#140.455#IT_c
					self.Press=Throat[1]#13.9751#IT_P
					self.Temp=Throat[2]#288.327#IT_T
					self.rho=Throat[3]#59.773443394897484#IT_rho
					self.gamma=Throat[4]#0.829105#IT_G
					self.h=Throat[5]#538.028#IT_h
					#print(IT_c,IT_P,IT_T,IT_rho,IT_G,IT_h)
					#pdb.set_trace()
					break
			alpha = calc.sqrt((1+self.delta)/(2*(self.gamma*self.rho_t*self.y_t)))
			epsilon = -alpha*self.y_t*self.gamma/(3+self.delta)
			self.x = -self.gamma*alpha*self.y**2/(3+self.delta)
			self.u = self.a*(1+(alpha*self.x+self.gamma*(alpha*self.y)**2/(1+self.delta)))
			self.v=self.a*((2*self.gamma*alpha**2*self.x*self.y/(1+self.delta))+(2*self.gamma**2*(alpha*self.y)**3)/((1+self.delta)*(3+self.delta)))
			self.x=self.x-epsilon
		elif (self.GasEqu == 'CoolProp'):
			#print('Sauer_RealGas1 with CoolProp')
			IT_rho = cp.PropsSI('D','H',self.Ho*1000,'S',self.s*1000,self.FLD)
			IT_T = cp.PropsSI('T','H',self.Ho*1000,'S',self.s*1000,self.FLD)-273
			for i in range(0,10000):
				IT_T = cp.PropsSI('T','D',IT_rho,'S',self.s*1000,self.FLD)-273
				IT_c = cp.PropsSI('A','D',IT_rho,'T',IT_T-273,self.FLD)
				IT_G = (cp.PropsSI('CPMASS','D',IT_rho,'T',IT_T-273,self.FLD)/cp.PropsSI('CVMASS','D',IT_rho,'T',IT_T-273,self.FLD))
				IT_h = cp.PropsSI('H','D',IT_rho,'T',IT_T-273,self.FLD)/1000
				IT_P = cp.PropsSI('P','D',IT_rho,'T',IT_T-273,self.FLD)/1e5
				f=IT_c**2/2000 +IT_h-self.Ho
				df=-IT_G*(IT_c**2*IT_rho)
				IT_v=1/IT_rho - f/df
				IT_rho=1/IT_v
				print(f)
				if 1:#(abs(f)<1e-3):
					Throat=np.loadtxt('ThroatProp.out')
					self.a=Throat[0]#140.455#IT_c
					self.Press=Throat[1]#13.9751#IT_P
					self.Temp=Throat[2]#288.327#IT_T
					self.rho=Throat[3]#59.773443394897484#IT_rho
					self.gamma=Throat[4]#0.829105#IT_G
					self.h=Throat[5]#538.028#IT_h
					#print(IT_c,IT_P,IT_T,IT_rho,IT_G,IT_h)
					#pdb.set_trace()
					break
			alpha = calc.sqrt((1+self.delta)/(2*(self.gamma*self.rho_t*self.y_t)))
			epsilon = -alpha*self.y_t*self.gamma/(3+self.delta)
			self.x = -self.gamma*alpha*self.y**2/(3+self.delta)
			self.u = self.a*(1+(alpha*self.x+self.gamma*(alpha*self.y)**2/(1+self.delta)))
			self.v=self.a*((2*self.gamma*alpha**2*self.x*self.y/(1+self.delta))+(2*self.gamma**2*(alpha*self.y)**3)/((1+self.delta)*(3+self.delta)))
			self.x=self.x-epsilon
		else:
			print("Invalid Input: GasEqu")'''

#---------------------------------------------------------------------------------------------#
	# Main Function for Sauer Analysis
	def SauerMain(self):
		#if (self.GasEqu=='RefProp'or self.GasEqu=='RefProp_TAB'):
		#	self.Sauer_RealGas()
			#self.Sauer_RefTab()
			#print('GAMMA',self.gamma)
			#pdb.set_trace()
		if 1:
			ao = calc.sqrt(self.gamma*self.R*self.To)
			epsilon = -1*(self.y_t/(2*(3+self.delta)))*(((self.gamma+1)*(1+self.delta)/(self.rho_t/self.y_t))**0.5);
			a_star = ((2*self.gamma*self.R*self.To)/(self.gamma+1))**0.5;
			alpha = ((1+self.delta)/((self.gamma+1)*self.rho_t*self.y_t))**0.5
			flag = 0		
			x_ = (-(self.gamma+1)*alpha*self.y**2/(2*(flag+1+self.delta)));
			flag = 2
			self.x = (-(self.gamma+1)*alpha*self.y**2/(2*(flag+1+self.delta)))
			u_dash = (alpha*self.x)+((((self.gamma+1)*alpha**2*self.y**2))/(2*(self.delta+1)))
			v_dash = ((((self.gamma+1)*alpha**2*self.y*self.x))/((self.delta+1)))+(((self.gamma+1)**2*alpha**3*self.y**3)/(2*(1+self.delta)*(3+self.delta)))
			self.u = a_star*(1+u_dash)
			v_tilde = a_star*v_dash	
			if (self.GasEqu!='RefProp' and self.GasEqu!='RefProp_TAB'):
				self.a = calc.sqrt(ao**2 - ((self.gamma-1)*self.u**2/(2)))
				self.M = self.u / self.a
				self.Press = self.Po/(1+((self.gamma-1)*self.M**2/(2)))**((self.gamma)/(self.gamma-1))
				self.Temp = self.To/(1+((self.gamma-1)*self.M**2/(2)))
				self.rho = self.Press/(self.R*self.Temp)
				self.x=self.x-epsilon
		self.M = calc.sqrt(self.u**2+self.v**2) / self.a

#---------------------------------------------------------------------------------------------#
	# Thermodynamics Properties Calculator (PRFT, EoS, EoS_TAB,RefProp)
	def ThermoProp(self,IN):
		self.V = calc.sqrt((self.u**2 + self.v**2))
		#PRFT: Solve using PERFECT gas equations [gamma necessary]
		if self.GasEqu == 'PRFT':
			Cp = self.R*self.gamma/(self.gamma-1)
			self.Temp = self.To - (self.V**2/(2*Cp))
			self.Press = self.Po*(self.Temp/self.To)**(self.gamma/(self.gamma-1))
			self.rho = self.Press / (self.R*self.Temp)
			self.a = calc.sqrt((self.ao**2)-((self.gamma-1)*self.V**2/2))
			self.M = self.V/self.a
		#EoS: Solve by directly solving Van der Waals EoS
		elif self.GasEqu == 'EoS':
			eqP = Symbol('eqP')
			eqT = Symbol('eqT')
			Vrho = Symbol('Vrho')
			b1 = self.R*self.Tc/(8*self.Pc)
			a1 = 27*(self.R*self.Tc)**2/(64*self.Pc)      
			Cp = self.gamma*self.R/(self.gamma-1)
			hig = self.Ho - (Cp)*(self.To-eqT)
			h = self.Ho - (0.5*self.V**2)
			VDW_h = - (2*a1/Vrho) + (b1*self.R*eqT/(Vrho-b1)) + hig - h
			VDW = (eqP*Vrho**3) -(eqP*b1+self.R*eqT)*Vrho**2 + a1*Vrho - a1*b1
			VDW_s = (self.gamma*log(eqT/self.To))/(self.gamma-1) - log(eqP/self.Po) + log(Vrho-b1) - log(Vrho) + log(eqP*Vrho/self.R/eqT)      
			try:[Vrho,T,P] = nsolve((VDW,VDW_h,VDW_s),(Vrho,eqT,eqP),(float(IN['VrhoGuess']),float(IN['TGuess']),float(IN['PGuess'])))
			except:pdb.set_trace()
			try:Vrho = float(Vrho.real)
			except:pdb.set_trace()
			try:self.Temp = float(T)
			except:self.Temp = float(T.real)
			try:self.Press = float(P)
			except:self.Press = float(P.real)
			self.rho = 1/float(Vrho)
			D1 = (2-self.gamma)*(Vrho-b1)/(self.gamma-1)
			D2 = self.R*self.Temp/self.Press
			N1 = self.R*T*((1/(Vrho-b1))+(1/Vrho))
			N2 = (P-(a1/Vrho**2) + (2*a1*b1/Vrho**3))*((2-self.gamma)/(self.gamma-1))
			try:a = calc.sqrt(Vrho**2*(N1.real+N2.real)/(D1.real+D2.real))
			except: pdb.set_trace()
			self.a = a
			self.M = self.V/self.a
			IN['VrhoGuess'] = float(Vrho.real)
			IN['TGuess'] = float(T.real)
			IN['PGuess'] = float(P.real)
		# EoS_TAB = Make and solve using van der Waals EoS table
		elif self.GasEqu == 'EoS_TAB':
			eqP = Symbol('eqP')
			eqT = Symbol('eqT')
			eqVrho = Symbol('eqVrho')
			b1 = self.R*self.Tc/(8*self.Pc)
			a1 = 27*(self.R*self.Tc)**2/(64*self.Pc)
			Cp = self.gamma*self.R/(self.gamma-1)
			hig = self.Ho - (Cp)*(self.To-eqT)
			# CREATE TABLE
			if IN['Table_Oper'] == 'Create':
				try:os.remove(os.getcwd()+'/Table.TAB')
				except: print ("Continue: No File to Delete <Table> \n Location:",os.getcwd())
				minV = int(IN['LLim'])
				maxV = int(IN['ULim'])
				dV = float(IN['dV'])
				V1= np.arange(minV,maxV,dV)
				for V in V1:
					self.V = V
					h = self.Ho - (0.5*self.V**2)
					VDW_h = - (2*a1/eqVrho) + (b1*self.R*eqT/(eqVrho-b1)) + hig - h
					VDW = (eqP*eqVrho**3) -(eqP*b1+self.R*eqT)*eqVrho**2 + a1*eqVrho - a1*b1
					VDW_s = (self.gamma*log(eqT/self.To))/(self.gamma-1) - log(eqP/self.Po) + log(eqVrho-b1) - log(eqVrho) + log(eqP*eqVrho/self.R/eqT)      
					try:[Vrho,T,P] = nsolve((VDW,VDW_h,VDW_s),(eqVrho,eqT,eqP),(float(IN['VrhoGuess']),float(IN['TGuess']),float(IN['PGuess'])))
					except:pdb.set_trace()
					try:Vrho = float(Vrho.real)
					except:pdb.set_trace()
					try:self.Temp = float(T)
					except:self.Temp = float(T.real)
					try:self.Press = float(P)    
					except:self.Press = float(P.real)
					self.rho = 1/float(Vrho)
            				###Refer Thesis for derivation###
					D1 = (2-self.gamma)*(Vrho-b1)/(self.gamma-1)
					D2 = self.R*self.Temp/self.Press
					N1 = self.R*T*((1/(Vrho-b1))+(1/Vrho))    
					N2 = (self.Press-(a1/Vrho**2) + (2*a1*b1/Vrho**3))*((2-self.gamma)/(self.gamma-1))        
					a = calc.sqrt(Vrho**2*(N1.real+N2.real)/(D1.real+D2.real))
					self.a = a
					self.M = self.V/self.a
					IN['VrhoGuess'] = float(Vrho.real)
					IN['TGuess'] = float(T.real)
					IN['PGuess'] = float(P.real)
					WriteDataFile('a',self,'Table.TAB','TAB')
				print('***Table Create***')
				sys.exit(-1)
			# CALCULATE FROM TABLE
			if IN['Table_Oper'] == 'Calc':
				TABProp = np.loadtxt('Table.TAB')
				for i in range(len(TABProp)):
					if TABProp[i][0]>self.V:
						self.M =  self.Inter(TABProp[i][0],TABProp[i-1][0],TABProp[i][1],TABProp[i-1][1])
						self.rho=  self.Inter(TABProp[i][0],TABProp[i-1][0],TABProp[i][2],TABProp[i-1][2])
						self.a=  self.Inter(TABProp[i][0],TABProp[i-1][0],TABProp[i][3],TABProp[i-1][3])
						self.Temp=  self.Inter(TABProp[i][0],TABProp[i-1][0],TABProp[i][4],TABProp[i-1][4])
						self.Press=  self.Inter(TABProp[i][0],TABProp[i-1][0],TABProp[i][5],TABProp[i-1][5])
						break
		# RefProp_TAB = Table generated from RefProp
		elif self.GasEqu == 'RefProp_TAB':
			if IN['Table_Oper']=='Create':
				try:os.remove(os.getcwd()+'/Table.TAB')
				except: print ("Continue: No File to Delete <Table> \n Location:",os.getcwd())
				minV = int(IN['LLim'])
				maxV = int(IN['ULim'])
				dV = float(IN['dV'])
				V1= np.arange(minV,maxV,dV)
				print('***Creating File***')
				for V in V1:
					self.V = V
					self.h = self.Ho - (0.5*self.V**2/1000)
					self.Press=fp.fluidprop_py(self.FLD,'hs',self.h,self.s,'pressure')
					self.Temp=fp.fluidprop_py(self.FLD,'hs',self.h,self.s,'temperature')
					self.rho=fp.fluidprop_py(self.FLD,'hs',self.h,self.s,'density')
					self.a=fp.fluidprop_py(self.FLD,'hs',self.h,self.s,'soundspeed')
					self.M=self.V/self.a
					self.gamma =fp.fluidprop_py(self.FLD,'hs',self.h,self.s,'gamma')
					WriteDataFile('a',self,'Table.TAB','TAB')
					#pdb.set_trace()
				print('***Table Create***')
				sys.exit(-1)
			if IN['Table_Oper']=='Calc':
				self.h = self.Ho - (0.5*self.V**2/1000)
				TABProp = np.loadtxt('Table.TAB') #[HARDCODED]
				for i in range(len(TABProp)):
					if TABProp[i][0]>self.V:
						#print(TABProp[i][0],self.V)
						self.M =  self.Inter(TABProp[i][0],TABProp[i-1][0],TABProp[i][1],TABProp[i-1][1])
						self.rho=  self.Inter(TABProp[i][0],TABProp[i-1][0],TABProp[i][2],TABProp[i-1][2])
						self.a=  self.Inter(TABProp[i][0],TABProp[i-1][0],TABProp[i][3],TABProp[i-1][3])
						self.Temp=  self.Inter(TABProp[i][0],TABProp[i-1][0],TABProp[i][4],TABProp[i-1][4])
						self.Press=  self.Inter(TABProp[i][0],TABProp[i-1][0],TABProp[i][5],TABProp[i-1][5])
						self.gamma=  self.Inter(TABProp[i][0],TABProp[i-1][0],TABProp[i][6],TABProp[i-1][6])
						#print(self.Temp,self.rho,self.a,self.gamma,self.M)
						break
		# RefProp = RefProp interpretor via python [Module fluidprop.py]
		elif self.GasEqu == 'RefProp':
			self.h = self.Ho - (0.5*self.V**2/1000)
			self.Press=fp.fluidprop_py(self.FLD,'hs',self.h,self.s,'pressure')
			self.Temp=fp.fluidprop_py(self.FLD,'hs',self.h,self.s,'temperature')
			self.rho=fp.fluidprop_py(self.FLD,'hs',self.h,self.s,'density')
			self.a=fp.fluidprop_py(self.FLD,'hs',self.h,self.s,'soundspeed')
			self.M=self.V/self.a
			#print(self.M)
			self.gamma =fp.fluidprop_py(self.FLD,'hs',self.h,self.s,'gamma')
		# CoolProp = Import fluid Properties using CoolProp
		elif self.GasEqu=='CoolProp':
			self.h     = self.Ho - (0.5*self.V**2/1000)
			self.Press = cp.PropsSI('P','H',self.h*1000,'S',self.s*1000,self.FLD)/1e5
			self.Temp  = cp.PropsSI('T','H',self.h*1000,'S',self.s*1000,self.FLD)-273
			self.rho   = cp.PropsSI('D','H',self.h*1000,'S',self.s*1000,self.FLD)
			self.a     = cp.PropsSI('A','H',self.h*1000,'S',self.s*1000,self.FLD)
			self.M     = self.V/self.a
			self.gamma = (cp.PropsSI('CPMASS','H',self.h*1000,'S',self.s*1000,self.FLD)/cp.PropsSI('CVMASS','H',self.h*1000,'S',self.s*1000,self.FLD))

		# ERROR Message
		else:
			print("\n\n\n\t*********INVALID INPUT DATA**********\n\t *GAS_EQU : PRFT -> Isentropic case \n\t *GAS_EQU : EoS -> Equations of State \n\t *GAS_EQU : EoS_TAB -> Equations of State using table \n\t *GAS_EQU : CoolProp -> Calcuate properties using CoolProp \n\t***************************************\n\n\n")
			sys.exit(1)

#---------------------------------------------------------------------------------------------#
	# Reflex Line position calculator (ReflexZone)
	def ReflexLine(self,n,l):
		RetObj   = Points()
		RetObj.x = self.x + 1*n*l*calc.cos(((1*calc.asin(1/self.M))+(calc.atan(self.v/self.u))))
		RetObj.y = self.y + 1*n*l*calc.sin(((1*calc.asin(1/self.M))+(calc.atan(self.v/self.u))))
		RetObj.u = self.u
		RetObj.v = 0
		return(RetObj)

#---------------------------------------------------------------------------------------------#
	# Check if the characteristic point is inside the nozzle [Exclusive for Reflex Zone]
	def PointInPoly(self,poly):
		n = len(poly)
		inside = False
		p1x = poly[0].x
		p1y = poly[0].y
		for i in range(n+1):
			p2x = poly[i % n].x
			p2y = poly[i % n].y
			if self.y > min(p1y,p2y):
				if self.y <= max(p1y,p2y):
					if self.x <= max(p1x,p2x):
						if p1y != p2y:
							xints = (self.y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
						if p1x == p2x or self.x <= xints:
							inside = not inside
			p1x,p1y = p2x,p2y
		return inside

#---------------------------------------------------------------------------------------------#
	# Liner Interpolation function for Table
	def Inter(self,x1,x2,y1,y2):
		if self.GasEqu == 'EoS_TAB':X = self.V
		else: X = self.h
		k=y1+(X-x1)*((y2-y1)/(x2-x1));
		return(k);
#---------------------------------------------------------------------------------------------#
# Sauer Analysis with RefProp (Ref: Alberto Code) [General Case]
# Guardone 2012, supersonic nozzle design code [point of contact: Matteo Pini]
def Sauer_RealGas(P):
	if (P.THROAT=='ITER'):
		if (P.GasEqu=='RefProp' or P.GasEqu=='RefProp_TAB'):
			IT_rho = fp.fluidprop_py(P.FLD,'hs',P.Ho,P.s,'density')
			IT_T   = fp.fluidprop_py(P.FLD,'vs',1/IT_rho,P.s,'temperature')
			for i in range(0,10000):
				IT_T   = fp.fluidprop_py(P.FLD,'vs',1/IT_rho,P.s,'temperature')
				IT_c   = fp.fluidprop_py(P.FLD,'Tv',IT_T,1/IT_rho,'soundspeed')			
				IT_G   = fp.fluidprop_py(P.FLD,'Tv',IT_T,1/IT_rho,'gamma')
				IT_h   = fp.fluidprop_py(P.FLD,'Tv',IT_T,1/IT_rho,'enthalpy')
				IT_P   = fp.fluidprop_py(P.FLD,'Tv',IT_T,1/IT_rho,'pressure')
				f      = IT_c**2/2000 +IT_h-P.Ho
				if i==0: f_tot=f
				df     = -IT_G*(IT_c**2*IT_rho)
				IT_v   = 1/IT_rho - f/df
				IT_rho = 1/IT_v
				printProgress(f_tot-f,f_tot,'REAL GAS: Sauer Iteration')
				if (abs(f)<1e-3):
					P.a     = IT_c
					P.Press = IT_P
					P.Temp  = IT_T
					P.rho   = IT_rho
					P.gamma = IT_G
					P.h     = IT_h
					print("\n\n*****Iteration Summary****\nsos\t:%f\nPress\t:%f\nTemp\t:%f\nRho\t:%f\ngamma\t:%f\nh\t:%f\n*****Summary End*****\n"%(IT_c,IT_P,IT_T,IT_rho,IT_G,IT_h))
					print("Writing Throat Prop. File\n\n")
					f = open('ThroatProp.out','w')
					f.write("%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t"%(IT_c,IT_P,IT_T,IT_rho,IT_G,IT_h,mu))
					f.close()
					break
		elif (P.GasEqu == 'CoolProp'):
			#print('DOING THIS WITH COOLPROP')
			IT_rho = cp.PropsSI('D','H',P.Ho*1000,'S',P.s*1000,P.FLD)
			IT_T   = cp.PropsSI('T','H',P.Ho*1000,'S',P.s*1000,P.FLD)-273
			#for i in range(0,10000):
			f = 1
			i = 0
			while (abs(f)>1e-3):
				IT_T   = cp.PropsSI('T','D',IT_rho,'S',P.s*1000,P.FLD)-273
				IT_c   = cp.PropsSI('A','D',IT_rho,'T',IT_T+273,P.FLD)
				IT_G   = (cp.PropsSI('CPMASS','D',IT_rho,'T',IT_T+273,P.FLD)/cp.PropsSI('CVMASS','D',IT_rho,'T',IT_T+273,P.FLD))
				IT_h   = cp.PropsSI('H','D',IT_rho,'T',IT_T+273,P.FLD)/1000
				IT_P   = cp.PropsSI('P','D',IT_rho,'T',IT_T+273,P.FLD)/1e5
				f      = IT_c**2/2000 +IT_h-P.Ho
				if i==0: f_tot=f
				df     = -IT_G*(IT_c**2*IT_rho)
				IT_v   = 1/IT_rho - f/df
				IT_rho = 1/IT_v
				i = i+1
				printProgress(f_tot-f,f_tot,'REAL GAS: Sauer Iteration')

			#print('Sauer finished after '+str(i)+' iterations')
			
			P.a     = IT_c
			P.Press = IT_P
			P.Temp  = IT_T
			P.rho   = IT_rho
			P.gamma = IT_G
			P.h     = IT_h
			mu = cp.PropsSI('V','D',IT_rho,'T',IT_T+273,P.FLD)
			print("\n\n*****Iteration Summary****\nsos\t:%f\nPress\t:%f\nTemp\t:%f\nRho\t:%f\ngamma\t:%f\nh\t:%f\n*****Summary End*****\n"%(IT_c,IT_P,IT_T,IT_rho,IT_G,IT_h))
			print("Writing Throat Prop. File\n\n")
			f=open('ThroatProp.out','w')
			f.write('Title = MoC Throat Propery File\nProperties =\nc\tP\tT\trho\tG\th\tmu\n')
			f.write("%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t"%(IT_c,IT_P,IT_T,IT_rho,IT_G,IT_h,mu))
			f.close()	
	else:
		throat_props = np.loadtxt(P.THROAT)
		P.a     = throat_props[0]
		P.Press = throat_props[1]
		P.Temp  = throat_props[2]
		P.rho   = throat_props[3]
		P.gamma = throat_props[4]
		P.h     = throat_props[5]

	alpha   = calc.sqrt((1+P.delta)/(2*(P.gamma*P.rho_t*P.y_t)))
	epsilon = -alpha*P.y_t*P.gamma/(3+P.delta)
	P.x     = -P.gamma*alpha*P.y**2/(3+P.delta)
	P.u     = P.a*(1+(alpha*P.x+P.gamma*(alpha*P.y)**2/(1+P.delta)))
	P.v     = P.a*((2*P.gamma*alpha**2*P.x*P.y/(1+P.delta))+(2*P.gamma**2*(alpha*P.y)**3)/((1+P.delta)*(3+P.delta)))
	P.x     = P.x-epsilon
	#print('Epsilon = '+str(epsilon)+' & P.x = '+str(P.x))
	return P

##
## END
##
