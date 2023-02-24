#!/usr/bin/sh

echo removing all $1 log/csv files
rm -r /data/mdi0316/WORK/MOFS/QML9/CSV/partial_csv/*${1}*
rm -r /data/mdi0316/WORK/MOFS/QML9/LOGS/*${1}*
scancel `squeue -u mdi0316 | grep dsgdb9ns | awk '{print $1}'`
