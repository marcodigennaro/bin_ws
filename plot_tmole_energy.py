#!/home/mdi0316/miniconda3/envs/kubas/bin/python

'''
quick print of energies from TURBOMOLE
'''

from turbomoleio.output.data import ScfEnergiesData
from turbomoleio.output.files import ScfOutput
from turbomoleio.output.files import JobexOutput

from monty.json import MSONable
from dataclasses import dataclass
from pathlib import Path
import os

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from pprint import pprint
import subprocess as sp

if len(sys.argv) == 2:
   run_dir = sys.argv[1]
else:
   run_dir = os.getcwd()

print(run_dir)
@dataclass
class TMOLE_class(MSONable):

    job_last = Path( run_dir ).joinpath( 'job.last' )
    ene_file = Path( run_dir ).joinpath( 'energy' )

    def read_turbomole_out_file(self):
        if Path(self.job_last).exists():
            with open(self.job_last, 'r') as f:
                out_lines = f.readlines()
            for o_line in out_lines:
                if 'd s c f - program' in o_line:
                    self.tmol_exe = 'dscf'
                    break
                if 'PROGRAM RIPER' in o_line:
                    self.tmol_exe = 'riper'
                    break
                if 'r i d f t' in o_line:
                    self.tmol_exe = 'ridft'
                    break
                if 'r d g r a d' in o_line:
                    self.tmol_exe = 'rdgrad'
                    break
            print(f'{self.tmol_exe} job found')

    '''this works for 0D/dscf - missing for riper'''
    def read_scf_energies(self):
        if self.tmol_exe in ['dscf', 'ridft', 'rdgrad']:
            #self.iterations = self.scf_output.as_dict()['scf']['iterations']
            self.scf_status = self.iterations['converged']
            ##equivalently
            #self.jbx_output = JobexOutput.from_file(self.job_last)
            #iterations = self.jbx_output.energy.scf.iterations
        if self.tmol_exe == 'riper':
            print('No reading implemented in turbomoleio')
            #return False, False

    def read_energy_file(self):
        self.jobex_energy = []
        if Path(self.ene_file).exists():
            with open(self.ene_file, 'r') as f:
                ene_lines = f.readlines()
            for count, line in enumerate(ene_lines):
                if line.startswith('$'):
                    continue
                iter, scf, scfkin, scfpot = line.split()
                self.jobex_energy.append(float(scf))
        return

    def read_riper_out(self):
        with open(self.job_last, 'r') as f:
            riper_lines = f.readlines()
        self.riper_status = False
        self.band_gap = False

        if any( [ 'PROGRAM RIPER' in line for line in riper_lines ] ):
            pass
        else:
            return 'Not a RIPER'

        for line in riper_lines[-200:]:
            if 'Band gap     ' in line:
                self.band_gap = float(line.split()[3])
                break
            if 'SCF converged within' in line:
                self.riper_status = 'converged'
            elif 'SCF HAS NOT converged within' in line:
                self.riper_status = 'not.converged'
                break

        df = pd.DataFrame()
        scf_iter_lines = [c for (c, l) in enumerate(
            riper_lines[:-30]) if "   SCF iteration " in l]
        last_ene = math.nan
        for count in scf_iter_lines:
            scf_iteration, damp_value, tot_energy, rms_density = [
                False for _ in range(4)]
            scf_iteration = riper_lines[count].split()[3]
            for jj in range(30):
                jj_line = riper_lines[count+jj]
                if 'new damping factor' in jj_line:
                    damp_value = float(jj_line.split()[4])
                elif 'TOTAL ENERGY' in jj_line:
                    tot_energy = float(jj_line.split()[4])
                elif 'RMS of difference density' in jj_line:
                    #rms_density = float(jj_line.split()[5].replace('D','E'))
                    rms_density = jj_line.split()[5].replace('D', 'E')
            tmp_dict = {'SCF': scf_iteration,
                        'damp': damp_value,
                        'tot_ene': tot_energy,
                        'last_ene': last_ene,
                        'rms_den': rms_density}
            last_ene = tot_energy
            tmp_df = pd.DataFrame([tmp_dict])
            df = pd.concat([df, tmp_df],  axis=0, ignore_index=True)

        if not df.empty:
            df['DE'] = df['tot_ene'] - df['last_ene']

        return df


    def avogadro(self):
        os.chdir(run_dir)
        sp.call('t2x coord > coord.xyz', shell=True)
        sp.call('avogadro coord.xyz', shell=True)
        return
        

tmole_obj = TMOLE_class()

if Path(tmole_obj.job_last).exists():
   tmole_obj.read_turbomole_out_file()
   
   tmole_obj.read_energy_file() 
   print('Jobex energy step = ', len(tmole_obj.jobex_energy))
   
   if tmole_obj.tmol_exe == 'riper': 
      riper_df = tmole_obj.read_riper_out()
      
      fig, ax = plt.subplots(2,1,sharex=True)
      ax[0].plot(riper_df['tot_ene'], label='tot_ene')
      ax[1].plot(riper_df['damp'], label='damp')
      plt.legend()
      plt.show()
   
   else: 
      fig, ax = plt.subplots()
      ax.plot(tmole_obj.jobex_energy)
      plt.legend()
      plt.show()


tmole_obj.avogadro()
