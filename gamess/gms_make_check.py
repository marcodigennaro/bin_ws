#!/home/mdi0316/anaconda3/bin/python

import json, os, sys
import numpy as np
import shutil
import math
import subprocess as sp

scripts_dir = '/home/mdi0316/FUNCTIONS'
classes_dir = '/home/mdi0316/CLASSES'
sys.path.insert(0, scripts_dir)
sys.path.insert(0, classes_dir)

import GAMESS
import SLURM
from Functions import print_tab, read_input_file, running_label

print( '>>> gms_make_check.py starts <<<' )
run_dir = os.getcwd()
print( run_dir )
check_dir = os.path.join( run_dir, 'CHECK' )
inp_file_name = [ f for f in os.listdir() if f.endswith('inp') ][0]
read_run = running_label( inp_file_name )
if read_run:
   print( 'Job is already running' )
else:
   try:
     log_file_name = [ f for f in os.listdir() if f.startswith('log') ][0]
   except(IndexError):
     print( 'missing output in {}'.format( run_dir ))
     sp.call( 'gms_submit_input.sh', shell=True )
     exit()
   
   print_tab(1, "run dir = {}".format( run_dir ))
   print_tab(1, "input file = {}".format( inp_file_name ))
   print_tab(1, "log file = {}".format( log_file_name ))
   
   orig_inp_file = os.path.join( run_dir, inp_file_name )
   orig_inp_dict = read_input_file( orig_inp_file, verbose=1 )
   orig_charge = orig_inp_dict['CONTRL']['ICHARG']
   orig_mult   = orig_inp_dict['CONTRL']['MULT']
   orig_runtyp = orig_inp_dict['CONTRL']['RUNTYP'] 
   orig_basis  = orig_inp_dict['BASIS']['GBASIS']
   orig_mwords = orig_inp_dict['SYSTEM']['MWORDS'] 
   orig_memddi = orig_inp_dict['SYSTEM']['MEMDDI'] 
   
   if 'CCTYP' in orig_inp_dict['CONTRL']:
      orig_scf = 'CCSDT'
      orig_funct = False
   elif 'MPLEVL' in orig_inp_dict['CONTRL']:
      orig_scf = 'MP2'
      orig_funct = False
   elif 'DFTTYP' in orig_inp_dict['CONTRL']: 
      orig_scf = 'DFT'
      orig_funct = orig_inp_dict['CONTRL']['DFTTYP']
   else:
      print( unknown_value )
   
   natoms_dict = { 'BF4' : 5, 'PF6' : 7, 'EMIM' : 19, 'EMIM_BF4' : 24, 'EMIM_PF6' : 26 }
   if 'CONVERGENCE' in run_dir or inp_file_name.startswith('cc_'):
      system_label = 'EMIM_BF4'
   else:
      system_label = [ item for item in run_dir.split('/') if item in natoms_dict.keys() ][0]
   
   orig_log_file = os.path.join( run_dir, log_file_name )
   log_lines = open( orig_log_file, 'r' ).readlines()
   log_lines.reverse() 
   for line in log_lines:
     if 'EXECUTION OF GAMESS TERMINATED NORMALLY' in line:
       print( 'job has terminated normally' )
       exit()

   if os.path.exists( check_dir ):
     #try:
        print_tab( 1, 'create dictionary of previous essais' )
        previous_dict = {}
        previous_idxs = os.listdir( check_dir ) 
        previous_idxs.sort()
        for check_idx in previous_idxs:
            print_tab( 3, '>> {}'.format( check_idx ) )
            tmp_inp_file = os.path.join( check_dir, check_idx, inp_file_name )
            tmp_inp_dict = read_input_file( tmp_inp_file, verbose=1 )
            previous_dict[check_idx] = {} 
            previous_dict[check_idx]['MWORDS'] = tmp_inp_dict['SYSTEM']['MWORDS']
            previous_dict[check_idx]['MEMDDI'] = tmp_inp_dict['SYSTEM']['MEMDDI']

        # check last tentative 
        print_tab( 1, 'check tentative: {}'.format(check_idx) )
        last_dir = os.path.join( check_dir, check_idx )
        #last_inp_name = [ f for f in os.listdir( last_dir ) if f.endswith('inp') ][0]
        last_log_name = [ f for f in os.listdir( last_dir ) if f.startswith('log') ][0]
        last_inp_file_name = os.path.join( last_dir, inp_file_name ) 
        last_log_file = os.path.join( last_dir, last_log_name ) 
        print_tab( 2, 'last_dir: {}'.format(last_dir) )

        last_obj = GAMESS.GAMESS( inp_name = inp_file_name, out_name = last_log_name, run_dir = last_dir,
                                  icharge = orig_charge,
                                  run_type = orig_runtyp, post_scf = orig_scf,
                                  functional = orig_funct, basis = orig_basis,
                                  mult = orig_mult, natoms = natoms_dict[system_label] )
   
        last_inp_dict = read_input_file( last_inp_file_name, verbose=1 )
        old_mwords = last_inp_dict['SYSTEM']['MWORDS']
        old_memddi = last_inp_dict['SYSTEM']['MEMDDI']
        print_tab( 2, 'old_mwords = {}'.format(old_mwords) )
        print_tab( 2, 'old_memddi = {}'.format(old_memddi) )
        new_mwords = old_mwords
        new_memddi = old_memddi

        ## read error
        print_tab( 2, [ 'Reading last log file', last_obj.out_file ] )
        last_exec, last_exec_err = last_obj.get_job_exec()
        print_tab( 2, '{}, {}'.format(last_exec, last_exec_err) ) 
        if [last_exec, last_exec_err] == [ 'TERMINATED.NORMALLY', False ]:
            print_tab( 3, 'check complete !!!' )
            os.chdir( run_dir )
            if new_mwords == orig_mwords and new_memddi == orig_memddi:
               print_tab( 4, [ 'WARNING:', 'new_mwords = orig_mwords = {}'.format(new_mwords),  'new_memddi = orig_memddi = {}'.format(new_memddi) ] ) 
               if int(new_mwords) < 10:
                  new_mwords = int(new_mwords) + 1
               else:
                  new_mwords = 1 + int(math.ceil( int(new_mwords) / 10.0 )) * 10   #round up to next 10
               new_memddi = int(math.ceil( int(new_memddi) / 100.0)) * 100  #round up to next 100
               print_tab( 4, [ 'new_mwords --> {}'.format(new_mwords),  'new_memddi --> {}'.format(new_memddi) ] ) 


            print_tab( 3, [ 'modifing input file:', orig_inp_file ] )
            old_lines = open( orig_inp_file, 'r' ).readlines()
            with open( orig_inp_file, 'w+' ) as new_inp_file:
              for line in old_lines:
                new_line = line
                if line.strip().startswith( 'MWORDS' ):
                   new_line = '  MWORDS={}\n'.format( new_mwords )
                elif line.strip().startswith( 'MEMDDI' ):
                   new_line = '  MEMDDI={}\n'.format( new_memddi )
                new_inp_file.write( new_line )
            print_tab( 3, 'Resubmitting with latest memory settings' )
            sp.call( 'sh /home/mdi0316/bin/gms_submit_input.sh', shell=True )
            print( 'Exiting...' )
            exit()    ## DON'T REMOVE
        else:
          last_log_lines = open( last_obj.out_file, 'r' ).readlines()
          last_log_lines.reverse()
          for last_count, last_line in enumerate(last_log_lines):
            if last_line.strip().startswith( 'PROCESS NO.' ) and all([ x in last_line for x in [ 'PROCESS', 'NO.', 'WORDS', 'REQUIRED=', 'AVAILABLE='] ]):
              print_tab( 2, last_line.strip() )
              new_mwords = int(last_line.split()[5]) // 1000000 + 1 
              break
            elif last_line.strip().startswith( 'THE DISTRIBUTED MEMORY REQUIRED FOR THIS STEP IS MEMDDI' ):
              print_tab( 2, last_line.strip() )
              new_memddi = last_line.split()[-2] 
              break
            elif last_line.strip().startswith('ADD') and all([ x in last_line for x in ['MORE', 'WORDS']]):
              print_tab( 2, last_line.strip() )
              new_mwords = int(last_line.split()[1]) // 1000000 + 1 + int(old_mwords)
              break
            elif all([ x in last_line for x in ['MEMDDI =', 'BUT NEEDS TO BE =']]):
              print_tab( 2, last_line.strip() )
              new_memddi = last_line.split()[-1]
              break
            elif last_line.strip() == '* INSUFFICIENT DISTRIBUTED MEMORY REQUESTED *':
              print_tab( 2, last_line.strip() )
              new_memddi = last_log_lines[last_count+2].split()[1]
              break
            elif last_line.strip() == '*** INSUFFICIENT DISTRIBUTED MEMORY REQUESTED ***':
              print_tab( 2, last_line.strip() )
              new_memddi = last_log_lines[last_count-1].split()[6]
              break
            elif last_line.strip() == 'PLEASE INCREASE -MWORDS- IN $SYSTEM APPROPRIATELY' and \
                 all([ x in last_log_lines[last_count+1] for x in ['NEED=', 'AVAILABLE=']]):
                 print_tab( 2, last_line.strip() )
                 print_tab( 2, last_log_lines[last_count+1].strip() )
                 new_mwords = int(last_log_lines[last_count+1].strip().split()[1]) // 1000000 + 1 
                 new_memddi = 10  # lower replicated memory
                 break
            else:
              pass #print( last_line )
   
          if int(new_mwords) % 10 != 0:
             new_mwords = int(math.ceil( int(new_mwords) / 10.0 )) * 10   #round up to next 10
          if int(new_memddi) % 100 != 0:
             new_memddi = int(math.ceil( int(new_memddi) / 100.0)) * 100  #round up to next 100
          print( 'new_mwords = {}'.format(new_mwords) )
          print( 'new_memddi = {}'.format(new_memddi) )
          for k, v in previous_dict.items(): 
            if int(new_mwords) == int(v['MWORDS']) and int(new_memddi) == int(v['MEMDDI']):
               print( '  same values as in "{:02d}". Skipping... '.format(int(k)) )
               #exit()
               new_mwords = input( 'insert new_mwords' )
               new_memddi = input( 'insert new_memddi' )

          ## make next input dict
          next_dir = os.path.join( check_dir, '{:02d}'.format( int(check_idx) + 1 ))
          next_inp_dict = dict(last_inp_dict)
          next_inp_dict['SYSTEM']['MWORDS'] = new_mwords
          next_inp_dict['SYSTEM']['MEMDDI'] = new_memddi
   
          ## make next object
          next_obj = GAMESS.GAMESS( inp_name = inp_file_name, out_name = last_log_name, run_dir = next_dir,
                                    icharge = orig_charge,
                                    run_type = orig_runtyp, post_scf = orig_scf,
                                    functional = orig_funct, basis = orig_basis,
                                    mult = orig_mult, natoms = natoms_dict[system_label] )
          next_obj.write_input_file_DICT( next_inp_dict )
   
     #except(ValueError):
     #  shutil.rmtree( check_dir )	
   else:
     # make first check with default settings
     next_dir = os.path.join( check_dir, '01' )
     os.makedirs( next_dir )
     sp.call( 'sed "s/EXETYP=RUN/EXETYP=CHECK/g" {} > {}/{}'.format( inp_file_name, next_dir, inp_file_name ), shell=True)
   
   
   ## submit new job
   slurm_obj = SLURM.SLURM( next_dir , 'GAMESS', job_name = inp_file_name )
   slurm_obj.write_batch()
   slurm_obj.submit_batch()
print( '--- gms_make_check.py ends ---' )
