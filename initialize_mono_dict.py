#!/home/mdi0316/anaconda3/bin/python

import os, shutil
import subprocess as sp
import json
import getpass
user = getpass.getuser()
work_dir = '/data/{}/WORK'.format(user)
mono_json = os.path.join( work_dir, 'monomers_{}.json'.format( user ) )

init_mono_dict = {

## ANIONS
   'Cl' : {'nat' :  1, 'charge': 0, 'mult' : 2, 'scftyp' : 'ROHF', 
                  'composition' : { 'Cl': 1 }, 
                  'at.orbitals' : { 'p' : 1 } },
    'F' : {'nat' :  1, 'charge': 0, 'mult' : 2, 'scftyp' : 'ROHF',
                  'composition' : { 'F' : 1 }, 
                  'at.orbitals' : { 'p' : 1 } },
    'S' : {'nat' :  1, 'charge': 0, 'mult' : 3, 'scftyp' : 'ROHF',
                  'composition' : { 'S' : 1 }, 
                  'at.orbitals' : { 'p' : 1 } },
    'P' : {'nat' :  1, 'charge': 0, 'mult' : 4, 'scftyp' : 'ROHF',
                  'composition' : { 'P' : 1 }, 
                  'at.orbitals' : { 'p' : 1 } },
  'PF6' : {'nat' :  7, 'charge':-1, 'mult' : 1, 'scftyp' : 'RHF', 
                 'composition' : { 'F' : 6, 'P' : 1 },
                 'at.orbitals' : { 'p' : 7 }},
  'SCN' : {'nat' :  3, 'charge':-1, 'mult' : 1, 'scftyp' : 'RHF', 
                 'composition' : { 'C' : 1, 'N' : 1, 'S' : 1},
                 'at.orbitals' : { 'p' : 3 }},
  'BF4' : {'nat' :  5, 'charge':-1, 'mult' : 1, 'scftyp' : 'RHF', 
                 'composition' : { 'B' : 1, 'F' : 4 },
                 'at.orbitals' : { 'p' : 3 }},
  'DEP' : {'nat' : 20, 'charge':-1, 'mult' : 1, 'scftyp' : 'RHF',
                 'composition' : { 'H' : 11, 'C' : 4, 'O' : 4, 'P' : 1 },
#                 'at.orbitals' : { 's' : 3 }
                  },
 'Tf2N' : {'nat' : 15, 'charge':-1, 'mult': 1, 'scftyp' : 'RHF',
                 'composition' : { 'C' : 2, 'N' : 1, 'O' : 4, 'S' : 2, 'F' : 6 }
                  },

## CATIONS

 'EMIM' : {'nat' : 19, 'charge': 1, 'mult' : 1, 'scftyp' : 'RHF',   #EMIM
                 'composition' : { 'H' : 11, 'C' : 6 , 'N' : 2 },
                 'at.orbitals' : { 's' : 11, 'p' : 8 }},

'C1MIM' : {'nat' : 16, 'charge': 1, 'mult' : 1, 'scftyp' : 'RHF',   #MMIM 
                 'composition' : { 'H' : 9, 'C' : 5, 'N' : 2 },
                 'at.orbitals' : { 's' : 9, 'p' : 7 }},

'C2MIM' : {'nat' : 19, 'charge': 1, 'mult' : 1, 'scftyp' : 'RHF',   #EMIM
                 'composition' : { 'H' : 11, 'C' : 6 , 'N' : 2 },
                 'at.orbitals' : { 's' : 11, 'p' : 8 }},

'C3MIM' : {'nat' : 22, 'charge': 1, 'mult' : 1, 'scftyp' : 'RHF', 
                 'composition' : { 'H' : 13, 'C' : 7 , 'N' : 2 },
                 'at.orbitals' : { 's' : 13, 'p' : 9 }},

'C4MIM' : {'nat' : 25, 'charge': 1, 'mult' : 1, 'scftyp' : 'RHF',   #BMIM
                 'composition' : { 'H' : 15, 'C' : 8 , 'N' : 2 },
                 'at.orbitals' : { 's' : 15, 'p' : 10 }},

'C5MIM' : {'nat' : 28, 'charge': 1, 'mult' : 1, 'scftyp' : 'RHF', 
                 'composition' : { 'H' : 17, 'C' : 9 , 'N' : 2 },
                 'at.orbitals' : { 's' : 17, 'p' : 11 }},

'C6MIM' : {'nat' : 31, 'charge': 1, 'mult' : 1, 'scftyp' : 'RHF', 
                 'composition' : { 'H' : 19, 'C' : 10 , 'N' : 2 },
                 'at.orbitals' : { 's' : 19, 'p' : 12 }},

'C7MIM' : {'nat' : 34, 'charge': 1, 'mult' : 1, 'scftyp' : 'RHF', 
                 'composition' : { 'H' : 21, 'C' : 11 , 'N' : 2 },
                 'at.orbitals' : { 's' : 21, 'p' : 13 }},

'C8MIM' : {'nat' : 37, 'charge': 1, 'mult' : 1, 'scftyp' : 'RHF', 
                 'composition' : { 'H' : 23, 'C' : 12 , 'N' : 2 },
                 'at.orbitals' : { 's' : 23, 'p' : 14 }}
}


#all_monomers = list(init_mono_dict.keys()) 
#all_monomers.remove('C2MIM')
##   all_monomers = ['EMIM']
#
#if os.path.exists( mono_json ):
#   existing_dict = json.load( open(mono_json, 'r') )
#   for mono in all_monomers:
#     sp.call( './monomers.py {}'.format( mono ), shell=True )
#
#else:
#   existing_dict = init_mono_dict

with open( mono_json, 'w+') as fp:
    json.dump(init_mono_dict, fp)

shutil.copy2( mono_json, '/home/{}/Inputfiles/GAMESS'.format(user) )


