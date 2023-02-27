#!/usr/bin/sh

mono_dir=/data/mdi0316/WORK/MONOMERS
for basis in `ls EMIM`; do
  bas_dir=`echo $mono_dir/EMIM/$basis`
  for func in `ls $bas_dir`; do
    func_dir=`echo $bas_dir/$func`
    for calc in `ls $func_dir`; do
      calc_dir=`echo $func_dir/$calc`
      for typ in `ls $calc_dir`; do
        run_dir=`echo $calc_dir/$typ`
        cd $run_dir
        all_emim_files=`ls | grep emim`
        c2mim_dir=`echo $run_dir | sed "s/EMIM/C2MIM/g"` 
        mkdir -p $c2mim_dir
        for file in $all_emim_files; do
          c2mim_file=`echo $file | sed "s/emim/c2mim/g"` 
          abs_c2mim_path=`echo $c2mim_dir/$c2mim_file`
          cp $file $abs_c2mim_path
          done
        done
      pwd
      done
    done
  done
