#!/usr/bin/env python
'''
Created on 29 Jun 2016

@author: COSMOlogic

Converts cif files to turbomole format
Needs openbabel installation: https://sourceforge.net/projects/openbabel/
and obabel binary (obabel) in path environmental variable
-no support for symmetry currently!


'''
import re,os.path,sys
import subprocess

def cif2tm(filename,basedir=''):
    cif_file = basedir  +filename
    
    #parse cif add H etc, only for pdb!
    subprocess.call(["obabel","-icif",cif_file,"-opdb","-O","tmp.pdb","--fillUC","keepconnect"])
    
    #convert to coord
    subprocess.call(["obabel","-ipdb","tmp.pdb","-otmol","-O","coord"])
    
    #direct
    #subprocess.call(["obabel","-icif",cif_file,"-otmol","-O","coord"])
    
    #read cellparameter
    with open(cif_file, 'r') as f:
        tmp = f.read().replace("\n","")
    a = float(re.findall("_cell_length_a\s+(\d+\.\d+)",tmp,re.IGNORECASE)[0])
    b = float(re.findall("_cell_length_b\s+(\d+\.\d+)",tmp,re.IGNORECASE)[0])
    c = float(re.findall("_cell_length_c\s+(\d+\.\d+)",tmp,re.IGNORECASE)[0])
    alpha = float(re.findall("_cell_angle_alpha\s+(\d+\.?\d+)",tmp,re.IGNORECASE)[0])
    beta = float(re.findall("_cell_angle_beta\s+(\d+\.?\d+)",tmp,re.IGNORECASE)[0])
    gamma = float(re.findall("_cell_angle_gamma\s+(\d+\.?\d+)",tmp,re.IGNORECASE)[0])
    symmetry = re.findall("_symmetry_cell_setting\s+([a-z]+)",tmp,re.IGNORECASE)
    
    #converts cell parameters to a.u.
    r_bohr = 0.529177210
    a = a / r_bohr
    b = b /r_bohr
    c = c/ r_bohr    
    cellp = [str(x) for x in [a,b,c,alpha,beta,gamma]]

    #check if control contains keywords
    if os.path.isfile("control"):
        with open("control",'r') as f:
            tmp = f.read()
        
        if re.findall("\$periodic",tmp):
            sys.stderr.write("WARNING: Keyword $periodic already exists\n")
    
        if re.findall("\$cell",tmp):
            sys.stderr.write("WARNING: Keyword $cell already exists\n")
    
    #open control file
    with open("control",'a') as f:
        f.write("$coord file=coord\n")
        f.write("$periodic    3\n")
        f.write("$cell\n")
        f.write(" ".join(cellp)+"\n")

if __name__=="__main__":
    if len(sys.argv)<2:
        print "USAGE: cif2tm cif_file"
    else:
        cif2tm(sys.argv[1])

    
        
    
    
    
