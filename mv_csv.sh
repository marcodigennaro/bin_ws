#!/usr/bin/sh

now=`date +'%Y_%m_%d_%T'`

mkdir -p old_CSV/$now

mv  DIMERS/EMIM_BF4/CSV/N311/B3LYP/scan_err.csv        old_CSV/$now/EMIM_BF4_N311_B3LYP_scan_err.csv
mv  DIMERS/EMIM_BF4/CSV/N311/PBE0/scan_err.csv         old_CSV/$now/EMIM_BF4_N311_PBE0_scan_err.csv
mv  DIMERS/EMIM_BF4/CSV/N311/M11/scan_err.csv          old_CSV/$now/EMIM_BF4_N311_M11_scan_err.csv
mv  DIMERS/EMIM_BF4/CSV/N311/wB97x-D/scan_err.csv      old_CSV/$now/EMIM_BF4_N311_wB97x-D_scan_err.csv
mv  DIMERS/EMIM_BF4/CSV/APCseg-1/wB97x-D/scan_err.csv  old_CSV/$now/EMIM_BF4_APCseg-1_wB97x-D_scan_err.csv
mv  DIMERS/EMIM_BF4/CSV/APCseg-1/PBE0/scan_err.csv     old_CSV/$now/EMIM_BF4_APCseg-1_PBE0_scan_err.csv
mv  DIMERS/EMIM_BF4/CSV/APCseg-1/M11/scan_err.csv      old_CSV/$now/EMIM_BF4_APCseg-1_M11_scan_err.csv
mv  DIMERS/EMIM_BF4/CSV/APCseg-1/B3LYP/scan_err.csv    old_CSV/$now/EMIM_BF4_APCseg-1_B3LYP_scan_err.csv
mv  DIMERS/EMIM_BF4/CSV/STO/PBE0/scan_err.csv          old_CSV/$now/EMIM_BF4_STO_PBE0_scan_err.csv
mv  DIMERS/EMIM_BF4/CSV/STO/B3LYP/scan_err.csv         old_CSV/$now/EMIM_BF4_STO_B3LYP_scan_err.csv
mv  DIMERS/EMIM_BF4/CSV/STO/wB97x-D/scan_err.csv       old_CSV/$now/EMIM_BF4_STO_wB97x-D_scan_err.csv
mv  DIMERS/EMIM_BF4/CSV/STO/M11/scan_err.csv           old_CSV/$now/EMIM_BF4_STO_M11_scan_err.csv
mv  DIMERS/EMIM_BF4/CSV/CCQ/B3LYP/scan_err.csv         old_CSV/$now/EMIM_BF4_CCQ_B3LYP_scan_err.csv
mv  DIMERS/EMIM_BF4/CSV/CCQ/PBE0/scan_err.csv          old_CSV/$now/EMIM_BF4_CCQ_PBE0_scan_err.csv
mv  DIMERS/EMIM_BF4/CSV/CCQ/M11/scan_err.csv           old_CSV/$now/EMIM_BF4_CCQ_M11_scan_err.csv
mv  DIMERS/EMIM_BF4/CSV/CCQ/wB97x-D/scan_err.csv       old_CSV/$now/EMIM_BF4_CCQ_wB97x-D_scan_err.csv
mv  DIMERS/EMIM_PF6/CSV/N311/B3LYP/scan_err.csv        old_CSV/$now/EMIM_PF6_N311_B3LYP_scan_err.csv
