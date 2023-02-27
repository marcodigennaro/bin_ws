#!/home/mdi0316/anaconda3/bin/python

import json, os, sys
import shutil
import subprocess as sp

scripts_dir = '/home/mdi0316/FUNCTIONS'
classes_dir = '/home/mdi0316/CLASSES'
sys.path.insert(0, scripts_dir)
sys.path.insert(0, classes_dir)

import GAMESS
from Functions import create_missing_zmat_labels

work_dir = '/data/mdi0316/WORK'
mono_json = os.path.join( work_dir, 'monomers_mdi0316.json' )
print( 'Reading: {}'.format( mono_json ))
with open(mono_json,'r') as json_file:
   mono_dict = json.load(json_file)

def main():
    cwd_list = os.getcwd().split('/')
    cwd = os.getcwd()
    print( 'checking ', cwd )
    inp_file = [ f for f in os.listdir() if f.endswith('inp') ][0]
    log_file = [ f for f in os.listdir() if f.startswith('log') ][0]

    if cwd_list[6] == 'MONOMERS':

       print( 'monomers optimization' )
       mono_label = cwd_list[8]
       mono_v = mono_dict[mono_label]
       mono_natom  = mono_v['nat']
       mono_charge = mono_v['charge']
       mono_mult   = mono_v['mult']
       gms_basis = cwd_list[9]
       gms_funct = cwd_list[10]
        
       gms_obj = GAMESS.GAMESS( inp_name = inp_file, out_name = log_file, run_dir = cwd,
                                icharge = mono_charge,
                                run_type = 'OPTIMIZE', post_scf = 'DFTTYP', 
                                functional = gms_funct, basis = gms_basis,
                                mult = mono_mult, 
                                natoms = mono_natom ) 

    elif cwd_list[6] == 'DIMERS':
       print( 'dimers optimization' )
       
    elif cwd_list[6] == 'CONVERGENCE': 
       print( 'convergence' )
       gms_basis = cwd_list[9]
       gms_funct = cwd_list[10]

       gms_obj = GAMESS.GAMESS( inp_name = inp_file, out_name = log_file, run_dir = cwd,
                                run_type = 'OPTIMIZE', post_scf = 'DFTTYP', 
                                functional = gms_funct, basis = gms_basis,
                                natoms = 24 ) 


    gms_inp_dict, gms_out_dict, gms_scf, gms_geom = gms_obj.get_job_results()
    gms_inp_dict.pop( 'DATA' )
    gms_inp_dict.pop( 'ZMAT' )
    
    last_iteration = sorted(list(gms_out_dict['OPTIMIZE'].keys()))[-1]
    last_zmat = gms_out_dict['OPTIMIZE'][last_iteration]['ZMAT'] 
   
    os.makedirs( 'OPT.STOPPED' )
    shutil.move( inp_file , 'OPT.STOPPED' )
    shutil.move( log_file , 'OPT.STOPPED' )
    sp.call( 'rm slurm* err.gms status.csv', shell = True )
    
    gms_obj.write_input_file_ZMAT( zmat_dict = last_zmat, header_dict = gms_inp_dict, 
                                   msg='prev.stopped.opt' )

    sp.call( 'sbatch -J {} submit_gamess.sh {}'.format( inp_file, inp_file ), shell=True )


if __name__ == '__main__':
  main()

