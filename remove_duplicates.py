#!/home/mdi0316/miniconda3/envs/gec/bin/python

import sys
import subprocess as sp

cmd = 'fdupes -r ./ > duplicates.txt'
sp.call( cmd, shell=True)

with open('duplicates.txt', 'r') as rf:
  readlines = rf.readlines()

with open('to_remove.txt', 'w+') as wf:
   for cc, ll in enumerate(readlines):
     if readlines[cc+1] == '\n':
        pass
     elif readlines[cc] == '\n':
        pass
     else:
        wf.write(f'rm {ll}')

