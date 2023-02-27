#!/usr/bin/sh

mkdir CHECK

inp_file=`ls ./ | grep inp`

echo "input file = $inp_file"
sed "s/EXETYP=RUN/EXETYP=CHECK/g" $inp_file > CHECK/$inp_file
cp submit_gamess.sh CHECK
cd CHECK

sbatch -J $inp_file submit_gamess.sh $inp_file

cd -
