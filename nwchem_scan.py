#!/home/mdi0316/anaconda3/bin/python

import os, sys, re
import numpy as np
import shutil
import pandas as pd
import subprocess as sp

scripts_dir = '/home/mdi0316/FUNCTIONS'
classes_dir = '/home/mdi0316/CLASSES'
work_dir = '/data/mdi0316/WORK'

sys.path.insert(0, scripts_dir)
sys.path.insert(0, classes_dir)
sys.path.insert(0, work_dir)

import dimers
from dimers import print_ionic_charges
from Functions import *
import NWCHEM

head_dir='/data/mdi0316/WORK/NWCHEM'
templ_dir='/home/mdi0316/Inputfiles/NWCHEM'
slurm_sh='/home/mdi0316/SCRIPTS/submit_nwchem.sh'

bf4_template_lines = open( os.path.join( templ_dir, 'bf4.nw' ), 'r' ).readlines()
emim_template_lines = open( os.path.join( templ_dir, 'emim.nw' ), 'r' ).readlines()
emim_bf4_template_lines = open( os.path.join( templ_dir, 'emim_bf4.nw' ), 'r' ).readlines()

nwc_opt_df = pd.DataFrame()
nwc_ene_df = pd.DataFrame()

bf4_nat = 5
bf4_dir = os.path.join( head_dir, 'MONOMERS', 'BF4', 'OPT' )
bf4_log = find_last_log(bf4_dir)
bf4_obj = NWCHEM.NWCHEM( 'bf4.nw', bf4_dir, bf4_nat, out_name=bf4_log )
bf4_ene = bf4_obj.read_out()[2]

emim_nat = 19
emim_dir = os.path.join( head_dir, 'MONOMERS', 'EMIM', 'OPT' )
emim_log = find_last_log(emim_dir)
emim_obj = NWCHEM.NWCHEM( 'emim.nw', emim_dir, emim_nat, out_name=emim_log )
emim_ene = emim_obj.read_out()[2]

emim_bf4_nat = 24
basis = '6-311g'
funct = 'B3LYP'

for T in [ 90 ]:
  for P in [ 0, 90, 180, 270 ]:
  #for P in [ 90 ]:
    for R in list(np.arange( 2, 10, 0.5)) + list(np.arange( 10, 30, 5.)):
      print( 'TPR= {}, {}, {}'.format( T, P, R ) )
      opt_dir = os.path.join(head_dir, 'DIMERS', 'EMIM_BF4', 'SCAN', basis, funct, 'T_{}'.format(T), 'P_{}'.format(P), 'R_{}'.format(R), 'OPT' )
      opt_label = 'emim_bf4_opt_T_{}_P_{}_R_{}'.format(T,P,R) 
      opt_input = '{}.nw'.format(opt_label)
      opt_output = find_last_log(opt_dir)
  
      ene_dir = os.path.join(head_dir, 'DIMERS', 'EMIM_BF4', 'SCAN', basis, funct, 'T_{}'.format(T), 'P_{}'.format(P), 'R_{}'.format(R), 'ENE' )
      ene_label = 'emim_bf4_ene_T_{}_P_{}_R_{}'.format(T,P,R) 
      ene_input = '{}.nw'.format(ene_label)
      ene_output = find_last_log(ene_dir)

      if opt_output: 
         opt_obj = NWCHEM.NWCHEM( opt_input, opt_dir, emim_bf4_nat, out_name=opt_output )
         opt_status, opt_error = opt_obj.read_status()
         if [opt_status, opt_error] == ['Complete', False]:
            opt_geom_line, opt_step, opt_energy = opt_obj.read_out()[:3]
            opt_zmat = opt_obj.read_zmat(opt_geom_line)
            nwc_opt_df = nwc_opt_df.append( [ { 'Radius' : R, 'Theta': T, 'Phi' : P, 'INT.ENE' : opt_energy-bf4_ene-emim_ene} ], ignore_index = True ) 
         else:
            print( 'opt_stauts, opt_error = {}, {} in {}'.format( opt_status, opt_error, opt_dir ) )
      else:
         print( 'opt_output = False in ', opt_dir )
         opt_obj.write_input( { 'BC19' : R, 'dih19': P } )

      if opt_output: 
         if [opt_status, opt_error] == ['Complete', False]:
           ### find gamess result for charge comparison 
           #for basis in ['N311', 'STO', 'APCseg-1', 'CCQ']:
           #   gms_dir = '/data/mdi0316/WORK/DIMERS/{}/RUNS/SCAN/ENE/OPT_from_DFT_N311_B3LYP/DFT/{}/{}/T_{}/P_{}/R_{}'.format( 'EMIM_BF4', basis, 'B3LYP', T, P, R ) 
           #   if os.path.exists( gms_dir ):
           #      gms_logs = [ f for f in os.listdir(gms_dir) if f.startswith('log') ]
           #      if len(gms_logs) == 0:
           #         print( 'NO logs' )
           #      elif len(gms_logs) == 1:
           #         gms_log = gms_logs[0]
           #         gms_inp = [ f for f in os.listdir(gms_dir) if f.endswith('.inp') ][0]
           #         gms_obj = GAMESS.GAMESS( inp_name = gms_inp, run_dir = gms_dir, out_name = gms_log, natoms = 24 )
           #         gms_out_dict = gms_obj.get_job_results()[1]
           #         for gms_k, gms_v in gms_out_dict['CHARGE.ANALYSIS'].items():
           #             gms_line = { 'Radius': R, 'Theta': str(T), 'Phi': str(P), 
           #                          'Mull.charge' : gms_v['Mull.charge'], 'Lowd.charge' : gms_v['Lowd.charge'], 
           #                          'index' : '{}{}'.format( gms_v['elem.'], gms_v['idx.']), 'basis' : basis, 'funct' : 'B3LYP' }
           #             gms_charges_df = gms_charges_df.append( gms_line, ignore_index = True )
           #      else:
           #         print( 'more than 1 logs' )
           #   else:
           #      print( 'missing', gms_dir)
           ###

           ene_obj = NWCHEM.NWCHEM( ene_input, ene_dir, emim_bf4_nat, out_name=ene_output, templ_name = 'emim_bf4.nw' )
           if ene_output: 
              ene_status, ene_error = ene_obj.read_status()
              if [ene_status, ene_error] == ['Complete', False]:
                 ## energy 
                 ene_geom_line, ene_step, ene_energy, ene_mull_line, ene_lodw_line = ene_obj.read_out()
                 ## charges
                 cat_charge, ani_charge = ene_obj.read_charges(ene_mull_line)
                 nwc_line = { 'Radius' : R, 'Theta' : T, 'Phi' : P, 'INT.ENE' : round(ene_energy-bf4_ene-emim_ene,8),
                              'MULL.CAT.CHARGE' : round(cat_charge,6), 'MULL.ANI.CHARGE' : round(ani_charge,6) }
                 nwc_ene_df = nwc_ene_df.append( nwc_line, ignore_index = True )
              else:
                 print( 'ene_stauts, ene_error = {}, {} in {}'.format( ene_status, ene_error, ene_dir ) )
           else:
              print( 'ene_output = False in ', ene_dir )
              ene_obj.write_input( { 'BC19' : R, 'dih19': P }, add_lines_dict = {'dft': ' mulliken \n print mulliken ao\n'}, 
                                                               modify_lines_dict = opt_zmat )

nwc_opt_df.to_csv( 'nwchem_opt_energy.csv' )
nwc_ene_df.to_csv( 'nwchem_dft_energy.csv' )
