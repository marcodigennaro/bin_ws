#!/usr/bin/sh

src='/data/mdi0316/WORK/IONIC_LIQUIDS'
dest='/storage/mdi0316/IONIC_LIQUIDS'
mkdir -p $dest/CSV
mkdir -p $dest/CSV/CHARGES
mkdir -p $dest/CSV/CHARGES/GAMESS
mkdir -p $dest/CSV/CHARGES/NWCHEM
mkdir -p $dest/CSV/ENERGIES
mkdir -p $dest/CSV/ENERGIES/GAMESS
mkdir -p $dest/CSV/ENERGIES/NWCHEM
mkdir -p $dest/TAR

read -p 'rsync_csv: ' rsync_csv
read -p 'rsync_mono: ' rsync_mono
read -p 'rsync_dime: ' rsync_dime
read -p 'rsync_convergence: ' rsync_convergence
read -p 'rsync_charges: ' rsync_charges
read -p 'tar_all:' tar_all

#if [ $# -eq 0 ]; then
#   rsync_csv=1
#   rsync_mono=1
#   rsync_dime=0
#   rsync_convergence=1
#   rsync_charges=0
#   tar_all=1
#else
#   rsync_csv=$1
#   rsync_mono=$2
#   rsync_dime=$3
#   rsync_convergence=$4
#   rsync_charges=$5
#   tar_all=$6
#fi

find_all_files()
{
  echo "find $1 -mindepth $2 -maxdepth $2 -name 'log*' -o -name 'gms*inp' -o -name 'gms*dat' -o -name 'gms*wfn' | grep $3 >> $4 "
  find $1 -mindepth $2 -maxdepth $2 -name 'log*' -o -name '*inp' -o -name '*dat' -o -name '*wfn' | grep -v -e FAILED -e STOPPED -e CHECK | grep $3 >> /storage/mdi0316/$4
}

rsync_from_file()
{
echo reading $1
echo $1 is `wc -l $1` lines long
for file_abs_path in `cat $1`; do 
  split_path=`echo $file_abs_path | sed "s/\// /g" `
  dest_dir="/storage/"
  dest_dir=` echo $file_abs_path | sed 's/data/storage/g; s/WORK//g' | sed 's/\// /g' | awk '{$NF="";print $0}' | sed 's/ /\//g; s/storage/\/storage/g'`
  echo $file_abs_path
  echo $dest_dir
  mkdir -p $dest_dir
  echo "rsync -ra --update $file_abs_path $dest_dir"
  rsync -ra --update $file_abs_path $dest_dir
  done
}

## CSV 
if [ $rsync_csv = 1 ]; then
  dest_dir="/storage/mdi0316/IONIC_LIQUIDS/CSV/CONVERGENCE"
  mkdir -p $dest_dir
  for file_name in `ls /data/mdi0316/WORK/IONIC_LIQUIDS/CONVERGENCE/T_90_P_90_R_3.5 | grep csv`; do
    echo $file_name
    rsync -ra --update /data/mdi0316/WORK/IONIC_LIQUIDS/CONVERGENCE/T_90_P_90_R_3.5/$file_name $dest_dir
    done
  
  dest_dir="/storage/mdi0316/IONIC_LIQUIDS/CSV/MONOMERS/"
  mkdir -p $dest_dir
  for file_name in `ls /data/mdi0316/WORK/IONIC_LIQUIDS/MONOMERS/CSV/`; do
    echo $file_name
    rsync -ra --update /data/mdi0316/WORK/IONIC_LIQUIDS/MONOMERS/CSV/$file_name $dest_dir
    done
  
  for file_name in scan_opt.csv scan_err.csv scan_dft_ene.csv scan_dft_eda.csv; do
    for ii in `find /data/mdi0316/WORK/IONIC_LIQUIDS/DIMERS/*/CSV -name $file_name` ; do
        dest_dir=`echo $ii | sed "s/\// /g" | awk '{print "/storage/mdi0316/IONIC_LIQUIDS/CSV/"$5"/"$6}'`
        dest_file=`echo $ii $file_name | sed "s/\// /g" | awk '{print $8"_"$9"_"$NF}'`
        mkdir -p $dest_dir
        echo "rsync -rva --update $ii $dest_dir/$dest_file"
        rsync -ra --update $ii $dest_dir/$dest_file
      done
    done
  
  for file_name in scan_mp2_ene.csv scan_mp2_eda.csv; do
    for ii in `find /data/mdi0316/WORK/IONIC_LIQUIDS/DIMERS/*/CSV -name $file_name` ; do
        dest_dir=`echo $ii | sed "s/\// /g" | awk '{print "/storage/mdi0316/IONIC_LIQUIDS/CSV/"$5"/"$6}'`
        dest_file=`echo $ii $file_name | sed "s/\// /g" | awk '{print $8"_"$NF}'`
        mkdir -p $dest_dir
        echo "rsync -ra --update $ii $dest_dir/$dest_file"
        rsync -ra --update $ii $dest_dir/$dest_file
      done
    done
fi
## CSV 

## MONOMERS inp/log files
if [ $rsync_mono = 1 ]; then
  if [ -f "$dest/all_gamess_mono_data" ]; then
    echo $dest/all_gamess_mono_data exists
    else
    echo printing $dest/all_gamess_mono_data
    find_all_files $src/MONOMERS 7 'OPT/DFT'   $dest/all_gamess_mono_data 
    find_all_files $src/MONOMERS 7 'ENE/DFT'   $dest/all_gamess_mono_data 
    find_all_files $src/MONOMERS 7 'ENE/MP2'   $dest/all_gamess_mono_data 
    find_all_files $src/MONOMERS 7 'ENE/CCSDT' $dest/all_gamess_mono_data 
    fi
  rsync_from_file $dest/all_gamess_mono_data 
  fi

## DIMERS inp/log files
if [ $rsync_dime = 1 ]; then
  if [ -f "$dest/all_gamess_dime_data" ]; then
    echo $dest/all_gamess_dime_data exists
    else
    echo printing $dest/all_gamess_dime_data
    find_all_files $src/DIMERS 11 'OPT/DFT'                          $dest/all_gamess_dime_data 
    find_all_files $src/DIMERS 12 'ENE/OPT_from_DFT_N311_B3LYP/DFT'  $dest/all_gamess_dime_data
    find_all_files $src/DIMERS 11 'ENE/OPT_from_DFT_N311_B3LYP/MP2'  $dest/all_gamess_dime_data
    fi
  rsync_from_file $dest/all_gamess_dime_data 
  fi

## CONVERGENCE inp/log files
if [ $rsync_convergence = 1 ]; then
  if [ -f "$dest/all_gamess_convergence_data" ]; then
    echo $dest/all_gamess_convergence_data 
    else 
    echo printing $dest/all_gamess_convergence_data 
    find_all_files $src/CONVERGENCE/T_90_P_90_R_3.5 4 'CONVERGENCE' $dest/all_gamess_convergence_data 
    find_all_files $src/CONVERGENCE/T_90_P_90_R_3.5 5 'CONVERGENCE' $dest/all_gamess_convergence_data
    find_all_files $src/CONVERGENCE/T_90_P_90_R_3.5 6 'CONVERGENCE' $dest/all_gamess_convergence_data
    find_all_files $src/CONVERGENCE/T_90_P_90_R_3.5 7 'CONVERGENCE' $dest/all_gamess_convergence_data
    fi
  rsync_from_file $dest/all_gamess_convergence_data
  fi

## CHARGES begins ##
if [ $rsync_charges = 1 ]; then
  #GAMESS
  if [ -f "$dest/all_gamess_charges" ]; then
    echo $dest/all_gamess_charges exists
    else
    echo printing $dest/all_gamess_charges
    echo "find ${src}/DIMERS -name 'atomic_charges.csv' | grep "T_90\/P_90" > $dest/all_gamess_charges"
    find ${src}/DIMERS -name 'atomic_charges.csv' | grep "T_90\/P_90" > $dest/all_gamess_charges
    fi
  echo reading $dest/all_gamess_charges
  echo $dest/all_gamess_charges is `wc -l $dest/all_gamess_charges` lines long
  for ii in `cat $dest/all_gamess_charges`; do
      new_ii=`echo $ii | sed 's/\// /g'| awk '{print $5" "$10" "$11" "$12" "$13" "$14" "$15" "$16}' | sed 's/.csv /.csv/g; s/ /_/g'`
      new_charge_file=$dest/CSV/CHARGES/GAMESS/$new_ii
      rsync -ra --update $ii $new_charge_file
      done
  #NWCHEM
  if [ -f "$dest/all_nwchem_charges" ]; then
    echo $dest/all_nwchem_charges exists
    else
    echo printing $dest/all_nwchem_charges
    echo "find ${src}/NWCHEM/DIMERS -name 'atomic_charges.csv' | grep "T_90\/P_90" > $dest/all_nwchem_charges"
    find ${src}/NWCHEM/DIMERS -name 'atomic_charges.csv' | grep "T_90\/P_90" > $dest/all_nwchem_charges
    fi
  echo reading $dest/all_nwchem_charges
  echo $dest/all_nwchem_charges is `wc -l $dest/all_nwchem_charges` lines long
  for ii in `cat $dest/all_nwchem_charges`; do
      new_ii=`echo $ii | sed 's/\// /g' | awk '{print $6" "$8" "$9" "$10" "$11" "$12" "$13" "$14}' | sed 's/.csv /.csv/g; s/ /_/g'`
      new_charge_file=$dest/CSV/CHARGES/NWCHEM/$new_ii
      rsync -ra --update $ii $new_charge_file
      done
fi
## CHARGES ends  ##

## TAR ALL
if [ $tar_all = 1 ]; then 
  echo 'Tarring all dirs'
  if [ -f "$dest/all_dft_P_dirs" ]; then
    echo $dest/all_dft_P_dirs exists
    else
    echo printing $dest/all_dft_P_dirs
    find IONIC_LIQUIDS -name "P_*" | grep ENE | grep -v MP2 > $dest/all_dft_P_dirs
    find IONIC_LIQUIDS -name "P_*" | grep ENE | grep    MP2 > $dest/all_mp2_P_dirs
  fi
  ## DFT
  #for P_dir in `cat $dest/all_dft_P_dirs`; do
  #  dest_name=`echo $P_dir | sed "s/\// /g" | awk '{print $3"_DFT_"$9"_"$10"_"$11"_"$12}'`
  #  dest_obj=$dest/TAR/${dest_name}.tar.gz
  #  echo "tar -zcvf $dest_obj $P_dir"
  #  tar -zcvf $dest_obj $P_dir 
  #  done
  # MP2
  for P_dir in `cat $dest/all_mp2_P_dirs`; do
    dest_name=`echo $P_dir | sed "s/\// /g" | awk '{print $3"_MP2_"$9"_"$10"_"$11}'`
    dest_obj=$dest/TAR/${dest_name}.tar.gz
    echo "tar -zcvf $dest_obj $P_dir"
    tar -zcvf $dest_obj $P_dir 
    done
  fi
## TAR ALL

#rm $dest/all_gamess_mono_data
#rm all_gamess_dime_data
#rm $dest/all_gamess_charges
#rm all_nwchem_charges
#rm all_dft_P_dirs
#rm $dest/all_gamess_convergence_data
