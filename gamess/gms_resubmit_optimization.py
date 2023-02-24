#!/home/mdi0316/anaconda3/bin/python

import json, os, sys
import shutil
import subprocess as sp

scripts_dir = '/home/mdi0316/FUNCTIONS'
classes_dir = '/home/mdi0316/CLASSES'
sys.path.insert(0, scripts_dir)
sys.path.insert(0, classes_dir)

import GAMESS
from Functions import create_missing_zmat_labels, print_tab, read_input_file, running_label

work_dir = '/data/mdi0316/WORK/IONIC_LIQUIDS'
mono_json = os.path.join( work_dir, 'monomers_mdi0316.json' )
with open(mono_json,'r') as json_file:
   mono_dict = json.load(json_file)

def main():
    cwd_list = os.getcwd().split('/')
    cwd = os.getcwd()
    print( 'checking ', cwd )
    inp_file = [ f for f in os.listdir() if f.endswith('inp') ][0]
    log_file = [ f for f in os.listdir() if f.startswith('log') ][0]
    print( 'input file:', inp_file )
    print( 'log file:', log_file )

    if 'EMIM_BF4' in cwd:
       nats = 24
    elif 'EMIM_PF6' in cwd:
       nats = 26
    else:
       if 'MONOMERS' in cwd:
          pass
       elif 'CONVERGENCE' in cwd:
          nats = 24
       else:
          print( unknown_nats )

    proceed = True
    read_run = running_label( inp_file )
    if read_run:
       print( 'Job is already running' )
       proceed = False
    else:
       log_lines = open( log_file, 'r' ).readlines()
       log_lines.reverse()
       for line in log_lines:
         if 'EXECUTION OF GAMESS TERMINATED NORMALLY' in line:
           print_tab( 1, line )
           proceed = False
         if 'THE GEOMETRY SEARCH IS NOT CONVERGED' in line:
           print_tab( 1, line )
           proceed = True

    if proceed:
       prev_zmat_msg_list = [ l for l in open(inp_file).readlines() if l.startswith('zmat') ]
       if len(prev_zmat_msg_list) == 0:
          prev_zmat_msg = 'zmat.default'
       else:
          prev_zmat_msg = prev_zmat_msg_list[0].strip()
       print( 'previous zmat msg: ', prev_zmat_msg )

       if cwd_list[7] == 'MONOMERS':

          print( 'monomers optimization' )
          mono_label = cwd_list[9]
          mono_v = mono_dict[mono_label]
          mono_natom  = mono_v['nat']
          mono_charge = mono_v['charge']
          mono_mult   = mono_v['mult']
          gms_basis = cwd_list[10]
          gms_funct = cwd_list[11]
           
          gms_obj = GAMESS.GAMESS( inp_name = inp_file, out_name = log_file, run_dir = cwd,
                                   icharge = mono_charge,
                                   run_type = 'OPTIMIZE', post_scf = 'DFTTYP', 
                                   functional = gms_funct, basis = gms_basis,
                                   mult = mono_mult, 
                                   natoms = mono_natom ) 

       elif cwd_list[7] == 'DIMERS':
          print( 'dimers optimization' )
          gms_basis = cwd_list[13]
          gms_funct = cwd_list[14]

          gms_obj = GAMESS.GAMESS( inp_name = inp_file, out_name = log_file, run_dir = cwd,
                                   run_type = 'OPTIMIZE', post_scf = 'OPT', 
                                   functional = gms_funct, basis = gms_basis,
                                   natoms = nats ) 

          
       elif cwd_list[7] == 'CONVERGENCE': 
          print( 'convergence' )
          gms_basis = cwd_list[10]
          gms_funct = cwd_list[11]

          gms_obj = GAMESS.GAMESS( inp_name = inp_file, out_name = log_file, run_dir = cwd,
                                   run_type = 'OPTIMIZE', post_scf = 'DFTTYP', 
                                   functional = gms_funct, basis = gms_basis,
                                   natoms = nats ) 

       if len(sys.argv) == 1:
          print( 'will read last iteration from\n', log_file )
          gms_out_dict = gms_obj.get_out_dict()
          last_iteration = sorted(list(gms_out_dict['OPTIMIZE'].keys()))[-1]
          last_zmat = gms_out_dict['OPTIMIZE'][last_iteration]['ZMAT'] 
          last_msg = 'zmat.from.{}'.format(log_file) 
       else:   
          last_iter_log_file = sys.argv[1]
          read_dir = '/'.join(last_iter_log_file.split('/')[:-1])
          print( 'will read last iteration from\n', last_iter_log_file )
          print( 'read_dir:', read_dir )
          read_inp = [ f for f in os.listdir(read_dir) if f.endswith('inp') ][0]
          read_log = [ f for f in os.listdir(read_dir) if f.startswith('log') ][0]
          print_tab( 1, 'dir: {} '.format(read_dir)) 
          print_tab( 1, 'inp: {} '.format(read_inp)) 
          print_tab( 1, 'log: {} '.format(read_log))
          
          read_inp_file = os.path.join( read_dir, read_inp )
          read_inp_dict = read_input_file( read_inp_file )
          read_funct = read_inp_dict['CONTRL']['DFTTYP'] 
          read_basis = read_inp_dict['BASIS']['GBASIS'] 

          if cwd_list[7] == 'MONOMERS':
             read_obj = GAMESS.GAMESS( inp_name = read_inp, out_name = read_log, run_dir = read_dir,
                                       icharge = mono_charge,
                                       run_type = 'OPTIMIZE', post_scf = 'DFTTYP', 
                                       functional = read_funct, basis = read_basis,
                                       mult = mono_mult, natoms = mono_natom ) 
          else:
             read_obj = GAMESS.GAMESS( inp_name = read_inp, out_name = read_log, run_dir = read_dir,
                                       run_type = 'OPTIMIZE', post_scf = 'DFTTYP', 
                                       functional = read_funct, basis = read_basis, natoms = nats )

          read_out_dict = read_obj.get_out_dict()
          last_iteration = sorted(list(read_out_dict['OPTIMIZE'].keys()))[-1]
          last_zmat = read_out_dict['OPTIMIZE'][last_iteration]['ZMAT'] 
          last_msg = 'zmat.from.{}'.format(read_log) 

       gms_inp_dict = gms_obj.read_input_file()
       gms_inp_dict.pop( 'DATA' )
       gms_inp_dict.pop( 'ZMAT' )

       try:
         os.makedirs( 'OPT.STOPPED/{}'.format(prev_zmat_msg), exist_ok=True )
         shutil.move( inp_file , 'OPT.STOPPED/{}'.format(prev_zmat_msg) )
         shutil.move( log_file , 'OPT.STOPPED/{}'.format(prev_zmat_msg) )
       except( shutil.Error ):
         pass 
       sp.call( 'rm slurm* err.gms status.csv log*', shell = True )
       
       gms_obj.write_input_file_ZMAT( zmat_dict = last_zmat, header_dict = gms_inp_dict, msg=last_msg )
     
       print( 'submitting new job' )
       sp.call( 'rm status.csv', shell=True )
       sp.call( 'sbatch -J {} submit_gamess.sh {}'.format( inp_file, inp_file ), shell=True )


if __name__ == '__main__':
  main()

