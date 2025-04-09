#!/usr/bin/env python

import math
import numpy as np
import os
import shutil
import python_subs as my
import sys


def PyRun(index,do_create_planet,do_put_in_core,do_relax,do_relax_model,do_evolve_planet,createmodel,coremodel,relaxedmodel,evolvemodel,mp_wo_core,rp,y,z,mcore,rhocore,irrad_col,flux_dayside,maxage,mesh_delta_coeff,time_delta_coeff,v_avg,lmin,lmax,p_atm, p_dyn,scaling_b,fohm,fdip,ohm_on,ohm_full,f_avg,w_age_max,pcmin,s0):

	# Prompt user before executing run.
	print("Compiling...")
	print()
	os.system("./mk")
	print()

	# Some additional output files.
	f = open('logfile_'+str(index)+'.txt', 'w')

	if do_create_planet:
		inlist1 = "inlist_create_" + str(index)
		run_time = my.create_planet(mp_wo_core,rp,y,z,inlist1,createmodel)

	success = True
	if (os.path.exists('LOGS/history.data')):
		k=open('LOGS/history.data','r')
		for line in k.readlines():
			pass
			last_temp=line
		k.close()
		last=last_temp.split()
		print("final model number in create=",last[0])
		print("last[0]==1000",last[0]=="1000")
		if last[0]=="1000":
			success=False
			outstring = '%6.3f\t%6.3f\t%6.3f\t%s\n' % (mp_wo_core, rp, run_time, success )
			f.write(outstring)
		if not success:
			print("Stopped at creating model")
			return

	if (mcore == 0.0):
		input_relax = createmodel
	else:
		if (do_put_in_core):
			print("INPUT CORE: ",createmodel)
			if not os.path.exists(createmodel):
				print("ERROR: createmodel not found",createmodel)
				return
			inlist2 = "inlist_core_" + str(index)
			run_time = my.put_core_in_planet(mcore,rhocore,inlist2,createmodel,coremodel)
			print("Core done")
		input_relax = coremodel

	if do_relax:
		print("INPUT RELAX: ",input_relax)
		if not os.path.exists(input_relax):
			print("ERROR: input_relax model not found",input_relax)
			return

		if do_relax_model:
			inlist3 = "inlist_relax_" + str(index)
			run_time = my.relax_model(irrad_col, flux_dayside, inlist3, input_relax, relaxedmodel)
		if not os.path.exists(relaxedmodel):
			print("ERROR: relaxedmodel not found",relaxedmodel)
			return
		input_evolve = relaxedmodel
	else:
		input_evolve = input_relax

	if do_evolve_planet:
		print("INPUT EVOLVE: ",input_evolve)
		if not os.path.exists(input_evolve):
			print("ERROR: relaxedmodel not found",input_evolve)
		inlist4 = "inlist_evolve_" + str(index)
		run_time = my.evolve_planet(irrad_col, flux_dayside, maxage, mesh_delta_coeff, time_delta_coeff, inlist4, input_evolve, evolvemodel, v_avg, lmin, lmax, p_atm, p_dyn, scaling_b, fohm, fdip, ohm_on, ohm_full, f_avg, w_age_max, pcmin, s0)

	if not os.path.exists(evolvemodel):
		print("ERROR: evolvemodel not found",evolvemodel)
		return
	
	f.close()
	

# Function for cleaning all contents of a directory (both files and subfolders).
def empty_directory(directory):
		for item in os.listdir(directory):
			item_path = os.path.join(directory, item)
			try:
				if os.path.isfile(item_path):
					os.remove(item_path)
				elif os.path.isdir(item_path):
					shutil.rmtree(item_path)
			except Exception as e:
				print(f"Failed to remove {item_path}. Reason: {e}")


if __name__=="__main__":

	import sys
	mesa_dir = "/home/daniele/codes/MESA/v24.08.1/"
	ctrl_path = mesa_dir+"mesa-24.08.1/star/defaults/"
	curr_path = mesa_dir+"make_planets/"
	model_path = curr_path+"models/"

	if not os.path.exists(model_path):
		os.mkdir(model_path)
			
	mjup = 1.8986e30
	mearth = 5.97e27
	sigma_sb = 5.670374419e-5		# Stefan-Boltzmann constant [erg cm^-2 s^-1 K^-4].
	
	rpin = 4.0		# inital planet radius in rjup
	y = 0.24		# helium fraction of both planet and (initial) star (default 0.24)

	# parameters having to do with irradiation
	irrad_col = 250.0		# column depth for depositing stellar radiation as heat (default 250)

	# Time parameters
	maxage = 1.e10			# ending age
	mesh_delta_coeff = 0.5	# grid step (default 1)
	time_delta_coeff = 0.5	# time step (default 1)

	# Ohmic model parameters
	lmin = 2
	lmax = 2
	p_atm = 10
	p_dyn = 1e6
	scaling_b = 'C'
	fohm = 0.5
	fdip = 0.1
	f_avg = 0.2
	w_age_max = 0.5
	ohm_on = 1e6
	ohm_full = 2e7
	pcmin = 0
	s0 = 0.25


	# Output directory
	out_path = mesa_dir+"RESULTS/2025.04.08-exploration_tohm20Myr_scale6_fohm05_fdip01_favg02_wavg05/"
	if not os.path.exists(out_path):
		os.mkdir(out_path)
		print("Creating the output directory: ",out_path)
	else:
		print("The output directory exists: ",out_path)
		user_input = input("Do you want to empty it? 'n' not to and exit, any other key to clear it and run the program: ")
		if user_input.lower() == "n":
			print("Exiting")
			sys.exit()
		else:
			empty_directory(out_path)
			print("Directory cleared")

	# ask = True for prompting messages in case of need of creating new models
	ask = False

	# Main planetary parameters to sweep
	zbase = 0.02		# metallicity of both planet and star (default 0.02)

	# Span M, Teq, vavg
	mass_arr = [0.7,1.,2.,4.]	# Planetary mass [Jupiter mass].
	mcore_base = 10.
	Teq_arr = [1500,1750,2000,2250]
	v_avg_arr = np.linspace(0,1.5e5,4)

	scale_inv_mass = True
	scale_inv_Teq = True

	# flags to do new models (True), or use the existing ones if they are there (False)
	# These flags don't change in the scripts
	new_createmodel = False
	new_coremodel = False

	# flags to do relaxation (True), and to use (False) or not (True) existing models, if they are there
	do_relax = True
	new_relaxedmodel = False

	# This is to evolve the planet, and should be always true
	do_evolve_planet = True

	# Creates the LOGS directory, needed to run
	if not os.path.exists('LOGS'):
		os.makedirs('LOGS')

	# Files to keep track of the run models, and the ones that crashed
	file_model = open(out_path+'model_list.txt', 'w')
	file_notstart = open(out_path+'model_notstarted_list.txt', 'w')
	file_crash = open(out_path+'model_crashed_list.txt', 'w')
	print("index      M[Mj]   Mc[Me] rhoc[gcc]  Rp[Rj] Teq[K] S*[cm2/g]      Y       Z     v_avg p_atm p_dyn sc.law fohm    fdip  pcmin  s0  ohm_on/full[yr]  f_avg  w_age_max  lmin-lmax",file=file_model)
	print("index (see model_list.txt for the parameters)",file=file_notstart)
	print("index (see model_list.txt for the parameters)",file=file_crash)
	file_model.close()
	file_notstart.close()
	file_crash.close()

	index = 0

	for im in range(len(mass_arr)):
		mass = mass_arr[im]
		mcore = mcore_base*np.sqrt(mass_arr[im])
		# z = zbase*np.sqrt(mass_arr[im])
		z = zbase
		mp_wo_core = mass - mcore*mearth/mjup
		rhocore = 10.*(mass_arr[im])		# Scaling relation, in order to have rhocore always larger than the bottom density of the envelope at old ages

		# For small masses, we need a smaller initial guess for the radius
		if (mass < 0.55):
			rp = 2.
		else:
			rp = rpin

		# Names for the models
		suffmc = str(mp_wo_core)[0:6] + "MJ"
		suffm = str(mass)[0:6] + "MJ"
		suffc = str(mcore)[0:6] + "Me"
		createmodel = model_path + "planet_create_" + suffmc +".mod"
		coremodel = model_path + "planet_core_" + suffm + suffc + ".mod"

		do_create_planet = new_createmodel
		# Create model in case the createmodel file doesn't exist
		if ((not new_createmodel) and (not os.path.exists(createmodel))):
			print("WARNING: createmodel doesn't exist: creating it",createmodel)
			print("do create model,new_createmodel:",do_create_planet,new_createmodel,"im",im)
			if (ask):
				user_input = input("'n' to exit, any other key to go on: ")
				if user_input.lower() == "n":
					print("Exiting")
					sys.exit()
			do_create_planet = True

		do_put_in_core = new_coremodel
		# Put core in case the coremodel file doesn't exist
		if ((not new_coremodel) and (not os.path.exists(coremodel))):
			print("WARNING: coremodel doesn't exist: doing it",coremodel)
			print("do put in core,new_coremodel:",do_put_in_core,new_coremodel,"im:",im)
			if (ask):
				user_input = input("'n' to exit, any other key to go on: ")
				if user_input.lower() == "n":
					print("Exiting")
					sys.exit()
			do_put_in_core = True

		for it in range(len(Teq_arr)):
			Teq = Teq_arr[it]
			flux_dayside = 4.0 * sigma_sb * Teq**4

			# Name for the relaxed model
			sufft = str(Teq)[0:6] + "K"
			relaxedmodel = model_path + "planet_relaxed_" + suffm + suffc + sufft + ".mod"

			# Relax model in case relaxation is required (do_relax),
			# an old relaxedmodel is required (new_relax_model=False) but the file doesn't exist
			do_relax_model = new_relaxedmodel
			
			if ((do_relax) and (not new_relaxedmodel) and (not os.path.exists(relaxedmodel))):
				do_relax_model = True
				new_relaxedmodel = True
				print("WARNING: relaxedmodel doesn't exist: relaxing it",relaxedmodel)
				print("do relax model,new_relaxedmodel:",do_relax_model,new_relaxedmodel,"im,it:",im,it)
				if ask:
					user_input = input("'n' to exit, any other key to go on: ")
					if user_input.lower() == "n":
						print("Exiting")
						sys.exit()

			for iv in range(len(v_avg_arr)):
				v_avg = v_avg_arr[iv]
				if (scale_inv_mass):
					v_avg = v_avg/np.sqrt(mass)
				if (scale_inv_Teq and Teq > 1500.):
					v_avg = v_avg/(Teq/1500.)**6

				# Create and put core only for the first temperature and first vavg case, for a given Teq
				# Relax model only for the first vavg case, for a given M and Teq, if the relaxed model doesn't exist
				if (iv > 0 or (not new_relaxedmodel)):
					if (it > 0):
						do_create_planet = False
						do_put_in_core = False
					do_relax_model = False
				else:
					do_relax_model = do_relax

				suffv = str(v_avg)[0:6] + "cms"
				evolvemodel = model_path + "planet_evolve_" + suffm + suffc + sufft + suffv + ".mod"

				print("Index,im,it,iv,create,core,relax,relax_model,evolve:",index,im,it,iv,do_create_planet,do_put_in_core,do_relax,do_relax_model,do_evolve_planet)

				# Print in the list of models the input parameters
				index = index + 1
				file_model = open(out_path+'model_list.txt', 'a')
				label="{:5d}".format(index),"{:8.1f}".format(mass),"{:8.1f}".format(math.floor(mcore)),"{:8d}".format(math.floor(rhocore)),"{:8d}".format(math.floor(rp)),"{:8d}".format(math.floor(Teq)),"{:8d}".format(math.floor(irrad_col)),"{:8.2f}".format(y),"{:8.3f}".format(z),"{:8d}".format(math.floor(v_avg)),"{:8.0e}".format(p_atm).replace("+0",""),"{:8.0e}".format(p_dyn).replace("+0",""),scaling_b,"{:8.2f}".format(fohm),"{:8.2f}".format(fdip),"{:8.0e}".format(pcmin).replace("+0",""),"{:8.2f}".format(s0),"{:8.0e}".format(ohm_on).replace("+0",""),"{:8.0e}".format(ohm_full).replace("+0",""),"{:8.2f}".format(f_avg),"{:8.2f}".format(w_age_max),"{:5d}".format(lmin),"{:5d}".format(lmax)
				   
				print("{:5d}".format(index),"{:8.1f}".format(mass),"{:8.1f}".format(math.floor(mcore)),"{:8d}".format(math.floor(rhocore)),"{:8d}".format(math.floor(rp)),"{:8d}".format(math.floor(Teq)),"{:8d}".format(math.floor(irrad_col)),"{:8.2f}".format(y),"{:8.3f}".format(z),"{:8d}".format(math.floor(v_avg)),"{:8.0e}".format(p_atm).replace("+0",""),"{:8.0e}".format(p_dyn).replace("+0",""),scaling_b,"{:8.2f}".format(fohm),"{:8.2f}".format(fdip),"{:8.0e}".format(pcmin).replace("+0",""),"{:8.2f}".format(s0),"{:8.0e}".format(ohm_on).replace("+0",""),"{:8.0e}".format(ohm_full).replace("+0",""),"{:8.2f}".format(f_avg),"{:8.2f}".format(w_age_max),"{:5d}".format(lmin),"{:5d}".format(lmax),file=file_model)
				file_model.close()

			    # RUN SIMULATIONS, WITH ALL NEEDED STEPS (create, core, relax, evolve)
				PyRun(index,do_create_planet,do_put_in_core,do_relax,do_relax_model,do_evolve_planet,createmodel,coremodel,relaxedmodel,evolvemodel,mp_wo_core,rp,y,z,mcore,rhocore, irrad_col,flux_dayside,maxage,mesh_delta_coeff,time_delta_coeff, v_avg,lmin,lmax,p_atm,p_dyn,scaling_b,fohm,fdip,ohm_on,ohm_full,f_avg,w_age_max,pcmin,s0)

				name = str(math.floor(index))
				print("="*20)
				print("Planet: ",label)
				print("="*20)
				   
				if os.path.isdir(curr_path+"LOGS"):

					shutil.copy(ctrl_path+"controls.defaults",out_path+"controls.defaults")
					shutil.copy(curr_path+"runscript_loops.py",out_path+"runscript_loops.py")
					shutil.copy(curr_path+"src/run_star_extras.f90",out_path+"run_star_extra.f90")
					new_path = out_path+"LOGS_"+name
					shutil.move(curr_path+"LOGS",new_path)
					os.makedirs(curr_path+"LOGS")

					if (os.path.exists(new_path+"/evolution.txt")):
						ag = np.loadtxt(new_path+"/evolution.txt",usecols=(0),unpack=True)
						# Add to crashed simulations if the age is too short:
						if (max(ag) < 0.9*maxage):
							print("Short evolution for model",name,":",max(ag)) 
							file_crash = open(out_path+'model_crashed_list.txt', 'a')
							print("Model index="+"{:3d}".format(index)+" at age (Gyr): "+"{:7.3f}".format(max(ag/1e9)),file=file_crash)
							file_crash.close()
						shutil.copy(new_path+"/evolution.txt",out_path+"evolution_"+name+".txt")

					if (os.path.exists("logfile_"+str(index)+".txt")):
						shutil.move("logfile_"+str(index)+".txt",out_path+"logfile_"+str(index)+".txt")
					if (os.path.exists("inlist_create_"+str(index))):
						shutil.move("inlist_create_"+str(index),out_path+"inlist_create_"+str(index))
					if (os.path.exists("inlist_core_"+str(index))):
						shutil.move("inlist_core_"+str(index),out_path+"inlist_core_"+str(index))
					if (os.path.exists("inlist_relax_"+str(index))):
						shutil.move("inlist_relax_"+str(index),out_path+"inlist_relax_"+str(index))
					if (os.path.exists("inlist_evolve_"+str(index))):
						shutil.move("inlist_evolve_"+str(index),out_path+"inlist_evolve_"+str(index))
					if (os.path.exists("inlist_ohmic")):
						shutil.move("inlist_ohmic",out_path+"inlist_ohmic_"+str(index))
					if (os.path.exists("inlist_create_"+str(index))):
						shutil.move("inlist_create_"+str(index),out_path+"inlist_create_"+str(index))
					print("="*20)
					print("The output folder moved to "+new_path)

					
				else:
					file_notstart = open(out_path+'model_notstarted_list.txt', 'a')
					print("Simulation not started\n")
					print("Model index="+"{:3d}".format(index),file=file_notstart)
					file_notstart.close()


	shutil.copytree(curr_path+"models",out_path+"models")

	print("Python Script: All done!")

