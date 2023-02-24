def get_opt_dict( self, save_zmat=False, out_file=None, 
                  optimize = True, mp2 = False, ccsdt = False, 
                  dft = True, eda = False, hessian = True ):
    
    if os.path.exists( self.nserch_json_file ):
       print_tab(4, 'all zmatrix already saved to disk')
       save_zmat = False

    out_lines = open( self.out_file, 'r', encoding = "ISO-8859-1" ).readlines()

    # common
    mulliken_lines = []
    distance_lines = []

    out_dict = {}
    opt_dict = {}
    mp2_dict = {}
    hes_dict = {}
    eda_dict = {}
    ccsdt_dict = {} 

    begin_geom_lines  = []
    located_geom_line = False
    aborted_optimization = False

    for o_count, o_line in enumerate( out_lines ):

        # COMMON
        dens_conv_search = re.search( self.density_conv_msg , o_line )
        ener_conv_search = re.search( self.energy_conv_msg  , o_line )
        diis_conv_search = re.search( self.diis_conv_msg    , o_line )
        fail_scf_search  = re.search( self.no_scf_conv_msg  , o_line )
        tot_ener_search  = re.search( self.tot_ener_msg     , o_line )

        if dens_conv_search or ener_conv_search:
           out_dict['SCF'] = 'CONVERGED'
        if diis_conv_search:
           #out_dict['DIIS'] = 'CONVERGED'
           out_dict['SCF'] = 'CONVERGED'
           if self.verbose:
              warnings.warn('DIIS converged not SCF')
        if fail_scf_search:
           out_dict['SCF'] = 'UNCONVERGED'
        if tot_ener_search:
           out_dict['TOT.EN.'] = float(o_line.strip().split('=')[1])*Ha2eV
           out_dict['INT.EN.'] = float(o_line.strip().split('=')[1])*Ha2eV - self.zero_energy

        mulliken_search     = re.search( self.mulliken_msg, o_line )
        distance_search     = re.search( self.distance_msg, o_line )
        if mulliken_search:
           mulliken_lines.append( (o_count+1, o_line) )
        if distance_search:
           distance_lines.append( (o_count+1, o_line) )
       
        # HESSIAN
        if heshian:
           ene_not_conv_search = re.search( ' ENERGY DID NOT CONVERGE...ABORTING HESSIAN', o_line )
           if ene_not_conv_search:
              hes_dict['ENE'] = 'NOT.CONVERGED'

        # CCSD(T)
        if ccsdt:
           stop_match = re.search( 'SUMMARY OF RESULTS', o_line )
           conv_match = re.search( 'THE    CCSD     ITERATIONS HAVE CONVERGED', o_line )
           ccsdt_en_match1 = re.search( 'CCSD\(T\) ENERGY', o_line )
           ccsdt_en_match2 = re.search( 'CCSD\[T\] ENERGY', o_line )
           ccsd_en_match   = re.search( 'CCSD    ENERGY', o_line )
           mbpt2_en_match  = re.search( 'MBPT\(2\) ENERGY', o_line )
           refer_en_match  = re.search( 'REFERENCE ENERGY', o_line )

           if ccsdt_en_match1:
                ccsdt_dict['CCSD(T)'] = { 'EN.' : Ha2eV*float(o_line.split()[2]), 'CORR.E.' : Ha2eV*float(o_line.split()[4]) }
           elif ccsdt_en_match2:
                ccsdt_dict['CCSD[T]'] = { 'EN.' : Ha2eV*float(o_line.split()[2]), 'CORR.E.' : Ha2eV*float(o_line.split()[4]) }
           elif ccsd_en_match  :
                ccsdt_dict['CCSD']    = { 'EN.' : Ha2eV*float(o_line.split()[2]), 'CORR.E.' : Ha2eV*float(o_line.split()[4]) }
           elif mbpt2_en_match :
                ccsdt_dict['MBPT(2)'] = { 'EN.' : Ha2eV*float(o_line.split()[2]), 'CORR.E.' : Ha2eV*float(o_line.split()[4]) }
           elif refer_en_match :
                ccsdt_dict['REF.EN.'] = Ha2eV*float(o_line.split()[2])
           elif conv_match:
                ccsdt_dict['SCF'] = 'CONVERGED'
             
        # MP2
        if mp2:
           en_mp2_search   = re.search( 'E\(MP2\)', o_line )
           if en_mp2_search:
              mp2_dict['MP2.EN.'] = Ha2eV * float(o_line.split('=')[1])

        # EDA
        if eda:
           own_basis_search = re.search( 'OWN BASIS SET', o_line )
           all_basis_search = re.search( 'ALL BASIS SET', o_line )
           if own_basis_search:
              own_basis_line = o_count
           if all_basis_search:
              all_basis_line = o_count
   
        # DFT-OPT
        if dft and optimize:
           located_geom_search = re.search( self.ok_geom_conv_msg, o_line )
           beg_geom_search     = re.search( self.beg_geom_msg, o_line )
           aborted_opt_search  = re.search( self.opt_abort_msg, o_line )
           if beg_geom_search:
              begin_geom_lines.append( (o_count, o_line) )
           if located_geom_search:
              located_geom_line = o_count + 1
           if aborted_opt_search:
              aborted_optimization = True

    ## END OF READING LOOP

    if eda:
       own_chunk = out_lines[ own_basis_line + 2 : own_basis_line + 7 ]
       all_chunk = out_lines[ all_basis_line + 2 : all_basis_line + 7 ]
       
       eda_dict['OWN.BS'] = {}
       eda_dict['OWN.BS']['ES.'] = Ha2eV * float(own_chunk[0].split()[3])
       eda_dict['OWN.BS']['EX.'] = Ha2eV * float(own_chunk[1].split()[3])
       eda_dict['OWN.BS']['REP.'] = Ha2eV * float(own_chunk[2].split()[3])
       eda_dict['OWN.BS']['POL.'] = Ha2eV * float(own_chunk[3].split()[3])
       eda_dict['OWN.BS']['INT.EN.'] = Ha2eV * float(own_chunk[4].split()[7])

       eda_dict['ALL.BS'] = {}
       eda_dict['ALL.BS']['ES.'] = Ha2eV * float(all_chunk[0].split()[3])
       eda_dict['ALL.BS']['EX.'] = Ha2eV * float(all_chunk[1].split()[3])
       eda_dict['ALL.BS']['REP.'] = Ha2eV * float(all_chunk[2].split()[3])
       eda_dict['ALL.BS']['POL.'] = Ha2eV * float(all_chunk[3].split()[3])
       eda_dict['ALL.BS']['INT.EN.'] = Ha2eV * float(all_chunk[4].split()[7])

    if opt:
       nserch_dict = {}
       for count, (start_line, beg_line_ii) in enumerate(begin_geom_lines):
           ii_dict = {'start.line' : start_line }
           try:
             end_line = begin_geom_lines[count+1][0] - 1
             ii_dict['end.line'] = end_line 
             ii_dict['chunk.size'] = end_line - start_line 
           except(IndexError): ## final iteration 
             if located_geom_line :
                end_line = located_geom_line
                ii_dict['chunk.size'] = located_geom_line - start_line 
             else:
                end_line = -1
                ii_dict['chunk.size'] = len(out_lines) - start_line 
           ii_dict['end.line'] = end_line
           ii_chunk = out_lines[ start_line: end_line]
           ii_nserch = int(beg_line_ii.split('=')[1].replace('...','').strip())
           ii_dict['NSERCH'] = ii_nserch
           ii_dict['SCF'] = 'unknown'
           for ii_line in ii_chunk:
               dens_conv_search = re.search( self.density_conv_msg , ii_line )
               ener_conv_search = re.search( self.energy_conv_msg  , ii_line )
               diis_conv_search = re.search( self.diis_conv_msg    , ii_line )
               fail_scf_search  = re.search( self.no_scf_conv_msg  , ii_line )
               end_geom_search  = re.search( self.end_geom_msg     , ii_line )

               if dens_conv_search or ener_conv_search:
                  ii_dict['SCF'] = 'CONVERGED'
               if diis_conv_search:
                  ii_dict['DIIS'] = 'CONVERGED'
               if fail_scf_search:
                  ii_dict['SCF'] = 'UNCONVERGED'
               if end_geom_search:
                  ii_dict['TOT.EN.'] = float(ii_line.strip().replace('NSERCH:','').split()[2])*Ha2eV
                  ii_dict['INT.EN.'] = float(ii_line.strip().replace('NSERCH:','').split()[2])*Ha2eV - self.zero_energy
           if save_zmat:
              ii_dict['ATOMS'], ii_dict['CART.COORDS.'], ii_dict['ZMAT'] = self.read_out_chunk( ii_chunk )
           nserch_dict[ii_nserch] = ii_dict 

       opt_dict['NSERCH'] = nserch_dict 

       if save_zmat:
          print_tab(4, 'saving all zmatrix to disk')
          with open( self.nserch_json_file, 'w+' ) as jf:
               json.dump( nserch_dict, jf )

       # initiate final dict
       if located_geom_line:
          
          final_dict  = ii_dict
          final_dict['GEOM.'] = 'LOCATED'

          final_chunk = ii_chunk
          final_dict['ATOMS'], final_dict['CART.COORDS.'], final_dict['ZMAT'] = self.read_out_chunk( final_chunk )
   
          final_dict['MULL.CHARGES'], final_dict['CHARG.CAT.'], final_dict['CHARG.ANI.'] = self.read_mulliken_charges( mulliken_lines )
          #final_dict['INTERNUCL.DISTANCES'] = self.read_internuclear_distances( distance_lines )
       else:
          # equilibrium not located
          length_dict = len(opt_dict['NSERCH'])
          ## read dict from last converged iteration
          final_dict = opt_dict['NSERCH'][length_dict-2]
          geometry_dict = opt_dict['NSERCH'][length_dict-1]
          ## read geometry error from last chunk
          geometry_chunk = out_lines[ geometry_dict[ 'start.line' ] : -1 ]
          for ii_line in geometry_chunk:
              no_geom_conv_search    = re.search(self.no_geom_conv_msg, ii_line)
              if   no_geom_conv_search:
                   final_dict['GEOM.'] = 'NOT.LOCATED'
          if aborted_optimization:
             final_dict['GEOM.'] = 'ABORTED'

       opt_dict['FINAL'] = final_dict

    if dft and optimize:
       out_dict['OPTIMIZE'] = opt_dict
    if mp2:
       out_dict['MP2'] = mp2_dict
    if hes:
       out_dict['HES'] = hes_dict
    if eda:
       out_dict['EDA'] = eda_dict
    if ccsdt:
       out_dict['CCSDT'] = ccsdt_dict

    out_dict['MULL.CHARGES'], out_dict['CHARG.CAT.'], out_dict['CHARG.ANI.'] = self.read_mulliken_charges( mulliken_lines )
    out_dict['INTERNUCL.DISTANCES'] = self.read_internuclear_distances( distance_lines )
    return out_dict
