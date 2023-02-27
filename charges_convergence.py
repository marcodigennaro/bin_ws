#!/home/mdi0316/anaconda3/bin/python

import os, sys, re
import shutil
import numpy as np
import pandas as pd

user = 'mdi0316'
scripts_dir = '/home/{}/FUNCTIONS'.format(user)
classes_dir = '/home/{}/CLASSES'.format(user)
work_dir = '/data/{}/WORK'.format(user)

sys.path.insert(0, scripts_dir)
sys.path.insert(0, classes_dir)
sys.path.insert(0, work_dir)

import IONIC_LIQUID as IL
import GAMESS

import dimers
from dimers import read_object, read_zmat, get_gms_object, print_dft_results, print_mp2_results, read_dimer
from Functions import print_tab, running_jobs, compose_zmatrices, running_label, center_of_charge, center_of_mass, Coulomb_Energy, angle_between, Ang2Bohr, read_input_file


def main():

    head_dir = '/data/mdi0316/WORK/CONVERGENCE/T_90_P_90_R_3.5'
    
    templ_dir = '/data/mdi0316/WORK/DIMERS/EMIM_BF4/RUNS/SCAN/OPT/DFT/N311/B3LYP/T_90/P_90/R_3.5'
    templ_inp = 'gms_SCAN_emim_bf4_OPT_DFT_N311_B3LYP_T_90_P_90_R_3.5.inp'
    templ_log = 'log.gms_SCAN_emim_bf4_OPT_DFT_N311_B3LYP_T_90_P_90_R_3.5_1204133'
    templ_obj = GAMESS.GAMESS( inp_name = templ_inp, run_dir = templ_dir, out_name = 'templ_log' )
    templ_out_zmat = read_zmat(templ_obj.zmat_file)
    
    dft_csv = os.path.join( head_dir, 'basis_convergence_dft.csv' )
    mp2_csv = os.path.join( head_dir, 'basis_convergence_mp2.csv' )
    
    dft_df = pd.DataFrame()
    mp2_df = pd.DataFrame()

    for basis in GAMESS.full_gbasis_list:
      for funct in GAMESS.full_functionals_list:
        try:
           opt_dimer, cat_nat, ani_nat, zero_dft_en, zero_mp2_en = read_dimer( 'EMIM_BF4', basis, funct )

           run_lab = 'cc_{}_{}'.format(basis, funct)
           run_dir = os.path.join( head_dir, 'RUNS', basis, funct )
           run_inp = '{}.inp'.format(run_lab)
           os.makedirs(run_dir, exist_ok = True)
           gms_obj = GAMESS.GAMESS( inp_name = run_inp, run_dir = run_dir, run_type = 'OPTIMIZE', post_scf = 'DFTTYP',
                                     natoms = 24, nat_cat = 19, nat_ani = 5,
                                     basis = basis, functional = funct, ifreeze = '52,53,54' )

           gms_exec, gms_exec_err, gms_err, scf, geom, time = read_object( gms_obj, read_dict = templ_out_zmat,
                                                                           read_msg = 'charge.convergence' )
       

           if [ gms_exec, gms_exec_err, gms_err, scf, geom ] == [ 'TERMINATED.NORMALLY', False, False, 'CONVERGED', 'LOCATED' ]:
       
              gms_out_zmat = read_zmat(gms_obj.zmat_file)
       
              opt_dir = os.path.join( run_dir, 'OPT' )
              opt_inp = '{}_OPT.inp'.format(run_lab)
              os.makedirs(opt_dir, exist_ok = True)
              opt_obj = GAMESS.GAMESS( inp_name = opt_inp, run_dir = opt_dir, run_type = 'OPTIMIZE', post_scf = 'DFTTYP',
                                       natoms = 24, nat_cat = 19, nat_ani = 5,
                                       basis = basis, functional = funct )
              opt_exec, opt_exec_err, opt_err, opt_scf, opt_geom, opt_time = read_object( opt_obj, read_dict = gms_out_zmat,
                                                                                          read_msg = 'OPT.{}.{}'.format(basis,funct))
              
              if [ opt_exec, opt_exec_err, opt_err, opt_scf, opt_geom ] == [ 'TERMINATED.NORMALLY', False, False, 'CONVERGED', 'LOCATED' ]:
                 opt_out_dict = opt_obj.get_job_results()[1]
                 rlx_R = opt_out_dict['FINAL']['ZMAT'][19]['STR']['val'] 
                 df_line = pd.DataFrame( [{ 'RLX.R' :  rlx_R, 'BASIS' : basis, 'FUNCT' : funct }] )
              
                 opt_zmat = read_zmat(opt_obj.zmat_file)
                 for post_scf in ['DFTTYP', 'MP2']:
                     post_dir = os.path.join( run_dir, post_scf )
                     post_inp = '{}_{}.inp'.format(run_lab, post_scf)
                     os.makedirs( post_dir, exist_ok = True )
                     post_obj = GAMESS.GAMESS( inp_name = post_inp, run_dir = post_dir, run_type = 'ENERGY', post_scf = post_scf,
                                               natoms = 24, nat_cat = 19, nat_ani = 5,
                                               basis = basis, functional = funct )
                     post_exec, post_exec_err, post_err = read_object( post_obj, read_dict = opt_zmat, read_msg = '{}.{}.{}'.format(post_scf,basis,funct))[:3]
                     if [ post_exec, post_exec_err, post_err ] == [ 'TERMINATED.NORMALLY', False, False ]:
                        print( post_obj.get_job_results()[1].keys() )
                        post_out_dict = post_obj.get_job_results()[1]
                        if post_scf == 'DFTTYP':
                           dft_line = print_dft_results( post_obj, post_out_dict, dimer=True, distance=float(rlx_R), zero_dft_ener=zero_dft_en )
                           dft_line = pd.concat( [ df_line, dft_line ], axis=1, join='inner' )
                        elif post_scf == 'MP2':
                           mp2_line = print_mp2_results( post_obj, post_out_dict, dimer=True, distance=float(rlx_R), zero_mp2_ener=zero_mp2_en )
                           mp2_line = pd.concat( [ df_line, mp2_line ], axis=1, sort=False )

                 dft_df = dft_df.append( dft_line, ignore_index=True )
                 mp2_df = mp2_df.append( mp2_line, ignore_index=True )

        except(FileNotFoundError):
           print( 'missing monomer' ) 
   
    dft_df.to_csv( dft_csv )
    mp2_df.to_csv( mp2_csv )

if __name__ == '__main__':
  main()
