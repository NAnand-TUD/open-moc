#!/usr/bin/env python3
#!/bin/sh
# -*- coding: utf-8 -*-

### IMPORT PACKAGES ###
#import subprocess
from subprocess import Popen, PIPE
import os
import pdb
import sys
import re
import time
from scipy.optimize import minimize
import ast
import shlex
import math
import time

### -------------------------------------------------------------------------------------------------------- ###

def modify_spacingctrl_specs(TE_coords, rad, pitch = 20, name = 'stator'):
	filename = 'Db/spacingcontrol.'+name
	line_pi = None
	line_p4 = None
	line_c = None

	"First read the ratios of the template file"
	with open(filename,'r') as file:
		for num, line in enumerate(file):
			
			if 'PITCH' in line: line_pi = num+1
			if ('PERIO' in line) and  ('4' in line): line_p4 = num+1 
			if 'RADIUS' in line: line_c = num+1

			if num == line_p4: h_p4 = list(filter(None,re.split(r' |\n',line)))

			if num == line_c: h_c = list(filter(None,re.split(r' |\n',line)))

	lines = []
	with open(filename,'r') as file:
		for num, line in enumerate(file):
			if (num==line_pi): 		
				vals = re.split(r' |\n',line)
				vals = list(filter(None,vals))
				line = line.replace(vals[0],"{:.1f}".format(pitch))
			if (num==line_c):
				line = line.replace(h_c[0],"{:.4f}".format(rad))
				line = line.replace(h_c[1],"{:.4f}".format(TE_coords[0]))
				line = line.replace(h_c[2],"{:.4f}".format(TE_coords[1]))
			lines.append(line)


		with open(filename, 'w') as outfile:
			for line in lines: outfile.write(line)

	return float(h_p4[0])


def spcingctrl_from_template(ref_el_size, name = 'stator', sf = 3, opti = False):
	
	filename = 'Db/spacingcontrol.'+name

	line_i = None
	line_o = None
	line_p4 = None
	line_p5 = None
	line_b = None
	line_c = None
	
	"First read the ratios of the template file"
	with open(filename,'r') as file:
		for num, line in enumerate(file):
			
			if 'INFLOW' in line: line_i = num+1
			if 'OUTFLOW' in line: line_o = num+1
			if ('PERIO' in line) and  ('5' in line): line_p5 = num+1
			if ('PERIO' in line) and  ('4' in line): line_p4 = num+1 
			if 'BLADE' in line: line_b = num+1
			if 'RADIUS' in line: line_c = num+1

			if num == line_i: h_i = list(filter(None,re.split(r' |\n',line)))

			if num == line_o: h_o = list(filter(None,re.split(r' |\n',line)))

			if num == line_p5: h_p5 = list(filter(None,re.split(r' |\n',line)))

			if num == line_p4: h_p4 = list(filter(None,re.split(r' |\n',line)))

			if num == line_b: h_b = list(filter(None,re.split(r' |\n',line)))

			if num == line_c: h_c = list(filter(None,re.split(r' |\n',line)))

	ref = float(h_p4[0])

	fac_i = [float(h_i[0])/ref,float(h_i[1])/ref]
	fac_o = [float(h_o[0])/ref,float(h_o[1])/ref]
	fac_p5 = [float(h_p5[0])/ref,float(h_p5[1])/ref]
	fac_p4 = float(h_p4[1])/ref
	fac_b = [float(h_b[0])/ref,float(h_b[1])/ref]
	fac_c = float(h_c[3])/ref

	in_min = ref_el_size*fac_i[0]
	in_max = ref_el_size*fac_i[1]
	out_min = ref_el_size*fac_o[0]
	out_max = ref_el_size*fac_o[1]
	per5_min = ref_el_size*fac_p5[0]
	per5_max = ref_el_size*fac_p5[1]
	per4_min = ref_el_size
	per4_max = ref_el_size*fac_p4
	bld_min = ref_el_size*fac_b[0]
	bld_max = ref_el_size*fac_b[1]
	te = ref_el_size*fac_c

	if opti:
		in_min = (re.sub('[\[ \]]', '', str(round_sf(in_min,sf))))
		in_max = (re.sub('[\[ \]]', '', str(round_sf(in_max,sf))))
		out_min = (re.sub('[\[ \]]', '', str(round_sf(out_min,sf))))
		out_max = (re.sub('[\[ \]]', '', str(round_sf(out_max,sf))))
		per5_min = (re.sub('[\[ \]]', '', str(round_sf(per5_min,sf))))
		per5_max = (re.sub('[\[ \]]', '', str(round_sf(per5_max,sf))))
		per4_min = (re.sub('[\[ \]]', '', str(round_sf(per4_min,sf))))
		per4_max = (re.sub('[\[ \]]', '', str(round_sf(per4_max,sf))))
		bld_min = (re.sub('[\[ \]]', '', str(round_sf(bld_min,sf))))
		bld_max = (re.sub('[\[ \]]', '', str(round_sf(bld_max,sf))))
		te = (re.sub('[\[ \]]', '', str(round_sf(te,sf))))

	else:
		in_min = str(round_sf(in_min,sf))
		in_max = str(round_sf(in_max,sf))
		out_min = str(round_sf(out_min,sf))
		out_max = str(round_sf(out_max,sf))
		per5_min = str(round_sf(per5_min,sf))
		per5_max = str(round_sf(per5_max,sf))
		per4_min = str(round_sf(per4_min,sf))
		per4_max = str(round_sf(per4_max,sf))
		bld_min = str(round_sf(bld_min,sf))
		bld_max = str(round_sf(bld_max,sf))
		te = str(round_sf(te,sf))

	lines = []
	with open(filename,'r') as file:
		for num, line in enumerate(file):
			if (num==line_i):		
				line = line.replace(h_i[0],in_min)
				line = line.replace(h_i[1],in_max)
			if (num==line_o):		
				line = line.replace(h_o[0],out_min)
				line = line.replace(h_o[1],out_max)
			if (num==line_p5):		
				line = line.replace(h_p5[0],per5_min)
				line = line.replace(h_p5[1],per5_max)
			if (num==line_p4):		
				line = line.replace(h_p4[0],per4_min)
				line = line.replace(h_p4[1],per4_max)
			if (num==line_b): 		
				line = line.replace(h_b[0],bld_min)
				line = line.replace(h_b[1],bld_max)
			if (num==line_c):
				line = line.replace(h_c[3],te)
			lines.append(line)

		with open(filename, 'w') as outfile:
			for line in lines: outfile.write(line)

def round_sf(x, sig=1):
	return round(x, sig-int(math.floor(math.log10(abs(x))))-1)

def modify_options(BL_thickness, name = 'stator', per_tol = '1.e-6'):
	lines = []
	with open('Db/options','r') as infile:
		for num, line in enumerate(infile):
			if 'name' in line: line_name = num+1
			try: 
				if (num==line_name): line = '\'grd\'    \'{}\'\n'.format(name)
			except NameError: pass
			if 'BL thickness' in line: line_num = num+1
			try: 
				#if (num==line_num): line = '0                            {}\n'.format(eformat(BL_thickness, 1,1))
				if (num==line_num): line = '0                            {:1.1e}\n'.format(BL_thickness)
			except NameError: pass
			if 'Periodic geometry' in line: line_tol = num+1
			try: 
				#if (num==line_num): line = '0                            {}\n'.format(eformat(BL_thickness, 1,1))
				if (num==line_tol): line = '.true.                  {}\n'.format(per_tol)
			except NameError: pass
			lines.append(line)

	with open('Db/options', 'w') as outfile:
			for line in lines: outfile.write(line)

def opti_mesh_refel(args):
	spcingctrl_from_template(args[0], name, opti = True)
	ec, out, err = run_cmd('UMG2.sh')
	success = umg2_success(out)
	if ec == 0 and success[0]:
		try:
			NELEM = read_su2_eln()
			os.remove('su2mesh.su2')
			res = abs((NELEM-MIN_NELEM)/MIN_NELEM)
			if NELEM<MIN_NELEM: res = 1
		except: 
			res = 1
			NELEM = 0
	else:
		if success[1] == 1:
			print('\nError with mesh.. Attempting to increase periodic tolerance...')
			increase_per_tol()
			NELEM = read_su2_eln()
			res = abs((NELEM-MIN_NELEM)/MIN_NELEM)
		elif success[1] == 2:
			print('\nError with mesh.. Attempting to reduce mesh displacement...')
			reduce_bl_thk()
			NELEM = read_su2_eln()
			res = abs((NELEM-MIN_NELEM)/MIN_NELEM)
		else:
			res = 1
			NELEM = 0
	print("REF ELEMENT SIZE: %.6f # ELEMENTS: %.0f RESIDUAL: %.3f "%(args[0], NELEM, res))
	return res

def read_su2_eln(name='su2mesh'):
	filename = name+'.su2'
	lines = []
	with open(filename,'r') as file:
		for line in file:
			if not any(line.startswith(s) for s in ['  ',' ','        ']):
				lines.append(line)
	
	for line in lines:
		words = re.split(r'=| |%|\n|#',line)
		if not any(words[0] in s for s in ['\n','%',' ','#']):
			words = list(filter(None,words))
			if (words[0]== 'NELEM'): NELEM = int(words[1])
	return NELEM

def ReadInputFile(name):
	IN = {}
	infile = open(name,'r')
	for line in infile:
		words = re.split(r'=| |%|\n|#',line)
		if not any(words[0] in s for s in ['\n','%',' ','#']):
			words = list(filter(None,words))
			IN[words[0]] = words[1]
	return IN

def run_cmd(cmd):
	args = shlex.split(cmd)
	proc = Popen(args, stdout=PIPE, stderr=PIPE)
	out, err = proc.communicate()
	exitcode = proc.returncode

	return exitcode, out, err

def round_sf(x, sig=1):
	return round(x, sig-int(math.floor(math.log10(abs(x))))-1)

def umg2_success(umg2_out):
	out_char = umg2_out.decode("utf-8")
	#pdb.set_trace()
	# Add errors here
	if 'Error' in out_char:
		if 'periodic' in out_char:
			return [False, 1]
		elif 'Decrease the displacement!' in out_char:
			return [False, 2]
		else:
		#pdb.set_trace()
			return [False, 0]
	elif 'error' in out_char:
		#pdb.set_trace()
		return [False, 0]
	else:
		return [True, 0]

def get_quality(umg2_out):
	out_char = umg2_out.decode("utf-8")
	#try:
	worst = float('{}.{}{}'.format(out_char[-41],out_char[-39],out_char[-38]))
	#except ValueError:
	#	pdb.set_trace()

	return worst

def umg2_quality(args):
	modify_options(BL_h, name)
	#umg2_success = False
	if args[0]<=0:
		print('FACTOR: NA, WORST ELEMENT: NA') 
		return 1
	rad = (BL_h+0.5*TE_thk)/(10*args[0])
	modify_spacingctrl_specs(TE_coords, rad, name = 'stator')
	ec, out, err = run_cmd('UMG2.sh')
	success = umg2_success(out)
	#pdb.set_trace()
	if ec == 0 and success[0]:
	#if ec == 0 and umg2_success(out):
		print('REFINEMENT RADIUS: {:.3f}, WORST ELEMENT: {:1.2f}'.format(rad,get_quality(out))) 
		return get_quality(out)/10
	else:
		if success[1] == 1:
			print('\nError with mesh.. Attempting to increase periodic tolerance...')
			increase_per_tol()
			NELEM = read_su2_eln()
			res = abs((NELEM-MIN_NELEM)/MIN_NELEM)
			print('REFINEMENT RADIUS: {:.3f}, WORST ELEMENT: {:1.2f}'.format(rad,get_quality(out))) 
			return get_quality(out)/10
		elif success[1] == 2:
			print('\nError with mesh.. Attempting to reduce mesh displacement...')
			reduce_bl_thk()
			ec, out, err = run_cmd('UMG2.sh')
			#pdb.set_trace()
			print('REFINEMENT RADIUS: {:.3f}, WORST ELEMENT: {:1.2f}'.format(rad,get_quality(out))) 
			return get_quality(out)/10
		else:
			print('FACTOR: NA, WORST ELEMENT: NA') 
			return 1

def reduce_bl_thk(stp = 1e-4):
	mesh = False
	#lines = []
	line_bl = -1
	with open('Db/options','r') as fob:
		for num, line in enumerate(fob):
			if 'BL thickness' in line: 
				line_bl = num+1
			if num == line_bl:
				BL_thk = float(line.split(' ')[-1])
	#pdb.set_trace()
	while not mesh:
		BL_thk -= stp
		lines = []
		with open('Db/options','r') as fob:
			for num, line in enumerate(fob):
			#if 'BL thickness' in line: 
			#	line_bl = num+1
				if num == line_bl:
					line = '0                            {:1.1e}\n'.format(BL_thk)

				#try: 
				#if (num==line_name): line = '\'grd\'    \'{}\'\n'.format(name)
				#except NameError: pass
				#if 'BL thickness' in line: line_num = num+1
				#try: 
				#	#if (num==line_num): line = '0                            {}\n'.format(eformat(BL_thickness, 1,1))
				#	if (num==line_num): line = '0                            {:1.1e}\n'.format(BL_thickness)
				#except NameError: pass
				lines.append(line)
		#pdb.set_trace()

		with open('Db/options', 'w') as outfile:
				for line in lines: outfile.write(line)
		print('Reducing displacement to {:1.1e}...'.format(BL_thk))
		ec, out, err = run_cmd('hyb2d.exe')
		success = umg2_success(out)
		#ec, out, err = run_cmd('UMG2.sh')
		if ec == 0 and success[0]:
			mesh = True
		if BL_thk<0:
			raise ValueError('Not possible to build original mesh!')
	
'''
def increase_per_tol(stp = 5e-7):
	mesh = False
	#lines = []
	line_tol = -1
	with open('Db/options','r') as fob:
		for num, line in enumerate(fob):
			if 'tol perio' in line: 
				line_tol = num+1
			if num == line_tol:
				tol_per = float(line.split(' ')[-1])
	while not mesh:
		tol_per += stp
		lines = []
		with open('Db/options','r') as fob:
			for num, line in enumerate(fob):
				#if 'BL thickness' in line: 
				#	line_bl = num+1
				if num == line_tol:
					line = '.true.                  {:1.1e}\n'.format(tol_per)

				#try: 
				#if (num==line_name): line = '\'grd\'    \'{}\'\n'.format(name)
				#except NameError: pass
				#if 'BL thickness' in line: line_num = num+1
				#try: 
				#	#if (num==line_num): line = '0                            {}\n'.format(eformat(BL_thickness, 1,1))
				#	if (num==line_num): line = '0                            {:1.1e}\n'.format(BL_thickness)
				#except NameError: pass
				lines.append(line)

		with open('Db/options', 'w') as outfile:
				for line in lines: outfile.write(line)
		print('Reducing periodic tolerance to {:1.1e}...'.format(tol_per))
		ec, out, err = run_cmd('UMG2.sh')
		success = umg2_success(out)
		if ec == 0 and success[0]:
			mesh = True
		if tol_per>1e-3:
			raise ValueError('Not possible to build original mesh!')
'''

def increase_per_tol(stp = 10**0.5):
	mesh = False
	line_tol = -1
	with open('Db/options','r') as fob:
		for num, line in enumerate(fob):
			if 'tol perio' in line: 
				line_tol = num+1
			if num == line_tol:
				tol_per = float(line.split(' ')[-1])
	while not mesh:
		tol_per *= stp
		lines = []
		with open('Db/options','r') as fob:
			for num, line in enumerate(fob):
				#if 'BL thickness' in line: 
				#	line_bl = num+1
				if num == line_tol:
					line = '.true.                  {:1.1e}\n'.format(tol_per)

				#try: 
				#if (num==line_name): line = '\'grd\'    \'{}\'\n'.format(name)
				#except NameError: pass
				#if 'BL thickness' in line: line_num = num+1
				#try: 
				#	#if (num==line_num): line = '0                            {}\n'.format(eformat(BL_thickness, 1,1))
				#	if (num==line_num): line = '0                            {:1.1e}\n'.format(BL_thickness)
				#except NameError: pass
				lines.append(line)

		with open('Db/options', 'w') as outfile:
				for line in lines: outfile.write(line)
		print('Reducing periodic tolerance to {:1.1e}...'.format(tol_per))
		ec, out, err = run_cmd('UMG2.sh')
		success = umg2_success(out)
		if ec == 0 and success[0]:
			mesh = True
		if tol_per>1e-3:
			raise ValueError('Not possible to build original mesh!')

def opti_nelem_mesh():
	modify_options(BL_h, name)
	ref_el = modify_spacingctrl_specs(TE_coords, radius, name = 'stator')
	ec, out, err = run_cmd('UMG2.sh')
	success = umg2_success(out)
	if ec != 0 or not success[0]: 
		#pdb.set_trace()
		if success[1] == 1:
			print('\nError with mesh.. Attempting to reduce periodic tolerance...')
			increase_per_tol()
			#raise ValueError('Periodic error in original mesh... Decrease tolerance!')
		elif success[1] == 2:
			print('\nError with mesh.. Attempting to reduce mesh displacement...')
			reduce_bl_thk()
			#raise ValueError('Displacement error in original mesh... Decrease BL!')
		else:
			raise ValueError('Not possible to build original mesh!')
	NELEM = read_su2_eln()
	if NELEM<MIN_NELEM:
		while NELEM<MIN_NELEM:
			modify_options(BL_h, name)
			ref_el = ref_el*(1-abs(NELEM-MIN_NELEM)/MIN_NELEM)
			res = minimize(opti_mesh_refel, [ref_el], method='Nelder-Mead', tol=0.025)
			ref_el = res.x[0]
			spcingctrl_from_template(ref_el, name)
			ec, out, err = run_cmd('UMG2.sh')
			success = umg2_success(out)
			if ec == 0 and success[0]: 
				NELEM = read_su2_eln()
			else:
				NELEM = 0.5*MIN_NELEM
	else:
		res = minimize(opti_mesh_refel, [ref_el], method='Nelder-Mead', tol=0.025)
	ref_el = res.x[0]
	spcingctrl_from_template(ref_el, name)
	run_cmd('UMG2.sh')
	return read_su2_eln(name='su2mesh')

def opti_mesh_quality():
	res = minimize(umg2_quality, [1/(10*rad_fac)], method='Nelder-Mead')
	rad = (BL_h+0.5*TE_thk)/(10*res.x[0])
	modify_spacingctrl_specs(TE_coords, rad, name = 'stator')
	ec, out, err = run_cmd('UMG2.sh')
	return get_quality(out)

#### ----------------------------------------------------------###

if __name__ == "__main__":
	try:    INFile  = sys.argv[1]
	except: INFile  = 'Mesh_Config.in'

	IN              = ReadInputFile(INFile)
	name = IN['FILE_NAME']
	rad_fac = float(IN['RAD_FAC'])
	MIN_NELEM = float(IN['MIN_NELEM'])

	try: scale = float(IN['scale'])
	except KeyError: scale = 1

	if IN['SPECS_FILE'] == 'NO':
		TE_coords = ast.literal_eval(IN['TE_coords'])
		y_t = float(IN['Throat'])
		TE_thk = float(IN['TE_thk'])
	else:
		is_specs = False
		t_strt   = time.time()
		print('\nSearching for blade specification file...')
		while not is_specs:
			try:
				specs = ReadInputFile(IN['SPECS_FILE'])
				is_specs = True
				print('\n\nReading specification file...')
			except:
				time.sleep(5)
			wait_time = time.time()-t_strt
			#print(wait_time)
			if wait_time>120:
				raise KeyError('No blade specification files found')	
		TE_coords = ast.literal_eval(specs['TE_coords'])
		y_t = float(specs['Throat'])
		TE_thk = float(specs['TE_thk'])

	if IN['BL_IN'] == 'FAC':
		BL_FAC = float(IN['BL_FAC'])
		BL_h = y_t*scale*BL_FAC
	if IN['BL_IN'] == 'LEN':
		BL_h = float(IN['BL_H'])

	radius = (0.5*TE_thk+BL_h)*rad_fac

	print('\n\nStarting element number optimization...\n')
	final_nelem = opti_nelem_mesh()
	print('\n\nStarting refinement radius optimization...\n')
	if rad_fac != 0:
		final_qual = opti_mesh_quality()
	else:
		ec, out, err = run_cmd('UMG2.sh')
	final_qual = get_quality(out)


	print('\n\nFinal Mesh Specifications: Mesh Elements = {:1.1e}, Worst Element Quality = {:.2f}'.format(final_nelem,final_qual))
	print('\n\n\t\t Succesfully generated SU2 mesh with UMG2 \n\n')


