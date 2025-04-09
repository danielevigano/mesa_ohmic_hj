#!/usr/bin/env python
import numpy as np
import os
import shutil
import time
import sys as s

#---------------------------------------------------------------------------
# Subroutines needed by the main python_run.py script.
#---------------------------------------------------------------------------

# Constants.
msun = 1.9884099e33
mjup = 1.8981246e30
rjup = 7.1492e9
mearth = 5.9721679e27

#---------------------------------------------------------------------------
# Create the initial planet without the core.
#---------------------------------------------------------------------------
def create_planet(mp_wo_core,rp,y,z,inlist1,createmodel):
    # Keep track of run time (see further below).
    start_time = time.time()

    print()
    print("Create initial planet...")

    # Read and write inlist.
    f = open('inlist_create','r')
    g = f.read()
    f.close()

    g = g.replace("<<m>>",str(mp_wo_core*mjup))
    g = g.replace("<<r>>",str(rp*rjup))
    g = g.replace("<<z>>",str(z))
    g = g.replace("<<y>>",str(y))
    g = g.replace("<<ritefile>>",'"' + createmodel + '"')

    h = open(inlist1,'w')
    h.write(g)
    h.close()
    shutil.copyfile(inlist1,"inlist")

    # Run.
    os.system('./star')
    print("SHOULD HAVE CREATED A MODEL:",createmodel)

    # Keep track of run time (see above).
    run_time = time.time() - start_time
    print("Run time for create_planets in seconds:",run_time)
    return run_time

#---------------------------------------------------------------------------
# Put in the core.
#---------------------------------------------------------------------------
def put_core_in_planet(mcore,rhocore,inlist2,createmodel,inputmodel):
    start_time = time.time()

    print()
    print("Put in core...")

    # Read and write inlist.
    f = open('inlist_core','r')
    g = f.read()
    f.close()

    g = g.replace("<<new_core_mass>>",str(mcore*mearth/msun))
    g = g.replace("<<core_avg_rho>>",str(rhocore))
    g = g.replace("<<loadfile>>",'"' + createmodel + '"')
    g = g.replace("<<ritefile>>",'"' + inputmodel + '"')

    h = open(inlist2,'w')
    h.write(g)
    h.close()
    shutil.copyfile(inlist2,"inlist")

    # Run.
    # print("*** Core: Uncomment the execution!")
    os.system('./star')

    run_time = time.time() - start_time
    print("Run time to put in core in seconds:",run_time)
    return run_time

#---------------------------------------------------------------------------
# Relax the model with the core.
#---------------------------------------------------------------------------
def relax_model(irrad_col, flux_dayside, inlist3, coremodel, inputmodel):
    start_time = time.time()

    # Read and write inlist.
    f = open('inlist_relax','r')
    g = f.read()
    f.close()

    g = g.replace("<<irrad_col>>",str(irrad_col))
    g = g.replace("<<flux_dayside>>",str(flux_dayside))
    g = g.replace("<<loadfile>>",'"' + coremodel + '"')
    g = g.replace("<<ritefile>>",'"' + inputmodel + '"')

    h = open(inlist3,'w')
    h.write(g)
    h.close()
    shutil.copyfile(inlist3,"inlist")

    # Run.
    # print("*** Core: Uncomment the execution!")
    os.system('./star')

    run_time = time.time() - start_time
    print("Run time to put in core in seconds:",run_time)
    return run_time

#---------------------------------------------------------------------------
# Evolve the planet.
#---------------------------------------------------------------------------
def evolve_planet(irrad_col, flux_dayside, maxage, mesh_delta_coeff, time_delta_coeff, inlist4, inputmodel, evolvemodel, v_avg, lmin, lmax, p_atm, p_dyn, scaling_b, fohm, fdip, ohm_on, ohm_full, f_avg, w_age_max, pcmin, s0):

    start_time = time.time()

    print()
    print("Evolve planet...",inputmodel)
	
    # Read and write inlist.
    f = open('inlist_evolve','r')
    g = f.read()
    f.close()

    g = g.replace("<<irrad_col>>",str(irrad_col))
    g = g.replace("<<flux_dayside>>",str(flux_dayside))
    g = g.replace("<<maxage>>",str(maxage))
    g = g.replace("<<mesh_delta_coeff>>",str(mesh_delta_coeff))
    g = g.replace("<<time_delta_coeff>>",str(time_delta_coeff))
    g = g.replace("<<loadfile>>",'"' + inputmodel + '"')
    g = g.replace("<<ritefile>>",'"' + evolvemodel + '"')

    h = open(inlist4,'w')
    h.write(g)
    h.close()
    shutil.copyfile(inlist4,"inlist")

    # Generate the Ohmic model input from the template.
    fadd = open('inlist_ohmic_template','r')
    gadd = fadd.read()
    fadd.close()

    gadd = gadd.replace('<<v_avg>>',str(v_avg))
    gadd = gadd.replace('<<lmin>>',str(lmin))
    gadd = gadd.replace('<<lmax>>',str(lmax))
    gadd = gadd.replace('<<p_atm>>',str(p_atm))
    gadd = gadd.replace('<<p_dyn>>',str(p_dyn))
    gadd = gadd.replace('<<scaling_b>>',scaling_b)
    gadd = gadd.replace('<<fohm>>',str(fohm))
    gadd = gadd.replace('<<fdip>>',str(fdip))
    gadd = gadd.replace('<<age_ohm_on>>',str(ohm_on))
    gadd = gadd.replace('<<age_ohm_full>>',str(ohm_full))
    gadd = gadd.replace('<<f_avg>>',str(f_avg))
    gadd = gadd.replace('<<w_age_max>>',str(w_age_max))
    gadd = gadd.replace('<<pcmin>>',str(pcmin))
    gadd = gadd.replace('<<s0>>',str(s0))

    hadd = open('inlist_ohmic','w')
    hadd.write(gadd)
    hadd.close()

    # Run.
    os.system('./star')

    run_time = time.time() - start_time
    print("Run time to evolve in seconds:",run_time)
    return run_time

#---------------------------------------------------------------------------
# Print utility.
#---------------------------------------------------------------------------
def print_parameters(mass,rp,mcore,rhocore,irrad_col,flux_dayside,Teq,y,z,maxage):
    print('######################################################')
    print('Parameters:')
    print('mp/mj =',mass)
    print('rp/rj =',rp)
    print('mcore/me =',mcore)
    print('rhocore (g cm^-3) =',rhocore)
    print('irrad_col (cm^2 g^-1) =',irrad_col)
    print('flux_dayside (Gerg cm^-2 s^-1) =',flux_dayside/1.e9)
    print('Teq (K)=',Teq)
    print('Y =',y)
    print('Z =',z)
    print('evolve to age/Gyr =',maxage/1.e9)
    print('######################################################')
    return
