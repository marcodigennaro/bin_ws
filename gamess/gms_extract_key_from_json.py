#!/home/mdi0316/anaconda3/bin/python

import json, os, sys

scripts_dir = '/home/mdi0316/FUNCTIONS'
classes_dir = '/home/mdi0316/CLASSES'
zmat_converter_dir = '/home/mdi0316/CLASSES/zmatrix-master'
sys.path.insert(0, scripts_dir)
sys.path.insert(0, classes_dir)
sys.path.insert(0, zmat_converter_dir)

import GAMESS

from gms_write_json import read_input

def main():

  with open( 'gms.json' , 'r' ) as jf:
     json_dict = json.load( jf )

  run_dir = os.getcwd()
  obj_inp = sys.argv[1]
  obj_key = sys.argv[2]
  inp_dict = read_input( obj_inp )
  
  tmp_coords = inp_dict['CONTRL']['COORD'] 
  tmp_runtyp = inp_dict['CONTRL']['RUNTYP']
  if 'DFTTYP' in inp_dict['CONTRL'].keys():
     tmp_postscf = 'DFTTPY'
     tmp_postscf_lab = 'DFT'
  else:
     tmp_postscf = 'NONE'
     tmp_postscf_lab = 'NONE'

  obj_calc = GAMESS.GAMESS_calculation( obj_inp, run_dir = run_dir, natoms = 24, 
                                        runtyp = tmp_runtyp, post_scf = tmp_postscf, coordinates = tmp_coords )
  print( obj_calc )
  out_dict = json_dict[tmp_postscf_lab][tmp_runtyp[:3]]
  if list(out_dict.keys()) == ['NSERCH', 'FINAL']:
     out_dict = dict(out_dict['FINAL'])

 
  for key, value in out_dict.items():
      if key == obj_key:
         print(value)

if __name__ == '__main__':
  main()

