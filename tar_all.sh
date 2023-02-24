#!/usr/bin/sh

src='/data/mdi0316/WORK'
dest='/storage/mdi0316/AC2T'

basis_list='N311 PBE0'
basis_list='N311'
funct_list='B3LYP'

for mono in EMIM PF6 BF4; do
    for basis in $basis_list; do 
        for funct in $funct_list; do 

             start_time=`date +%s`
             opt_dir=MONOMERS/$mono/$basis/$funct/OPT/DFT
             dft_dir=MONOMERS/$mono/$basis/$funct/ENE/DFT
             mp2_dir=MONOMERS/$mono/$basis/$funct/ENE/MP2

             mkdir -p $dest/$opt_dir
             mkdir -p $dest/$dft_dir
             mkdir -p $dest/$mp2_dir

             rsync -ra --update --exclude 'slurm*' --exclude 'err.gms' $src/$opt_dir/ $dest/$opt_dir --delete
             rsync -ra --update --exclude 'slurm*' --exclude 'err.gms' $src/$dft_dir/ $dest/$dft_dir --delete
             rsync -ra --update --exclude 'slurm*' --exclude 'err.gms' $src/$mp2_dir/ $dest/$mp2_dir --delete

             echo $basis $funct $mono $(expr `date +%s` - $start_time) secs. 
        done
    done
done

for dime in EMIM_PF6 EMIM_BF4; do

    csv_dir=DIMERS/$dime/CSV/$basis/$funct
    mkdir -p $dest/$csv_dir 
    rsync -ra --update --exclude '*csv' --exclude 'slurm*' --exclude 'err.gms' $src/$csv_dir/  $dest/$csv_dir  --delete    

    for T in T_90 T_5  T_45  T_135  T_175; do 
       for P in P_0 P_90 P_135  P_180  P_225  P_270  P_315  P_45; do
         for basis in $basis_list; do 
           for funct in $funct_list; do 

             start_time=`date +%s`
              
             opt_dir=DIMERS/$dime/RUNS/SCAN/OPT/DFT/$basis/$funct
             dft_dir=DIMERS/$dime/RUNS/SCAN/ENE/OPT_from_DFT_N311_B3LYP/DFT/$basis/$funct
             mp2_dir=DIMERS/$dime/RUNS/SCAN/ENE/OPT_from_DFT_N311_B3LYP/MP2/$basis

             TP_dir=$T/$P

             mkdir -p $dest/$opt_dir/$TP_dir 
             mkdir -p $dest/$dft_dir/$TP_dir 
             mkdir -p $dest/$mp2_dir/$TP_dir 

             for tmp_dir in $opt_dir $dft_dir $mp2_dir; do 
                 rsync -ra --update --exclude '*csv' --exclude 'slurm*' --exclude 'err.gms' \
                       $src/$tmp_dir/$TP_dir/ $dest/$tmp_dir/$TP_dir --delete
#             rsync -ra --update --exclude '*csv' --exclude 'slurm*' --exclude 'err.gms' $src/$dft_dir/$TP_dir/ $dest/$dft_dir/$TP_dir --delete
#             rsync -ra --update --exclude '*csv' --exclude 'slurm*' --exclude 'err.gms' $src/$mp2_dir/$TP_dir/ $dest/$mp2_dir/$TP_dir --delete
                 done
             echo $basis $funct $T $P $(expr `date +%s` - $start_time) secs. 
             sleep 1
             done
          done
        done
    done
    tar -zcvf AC2T/DIMERS/ene_${dime}.tar.gz AC2T/DIMERS/$dime/RUNS/SCAN/ENE/
done


