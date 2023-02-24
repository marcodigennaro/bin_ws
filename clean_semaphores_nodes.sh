#!/usr/bin/sh

##Developed By Kamal maiti, 
##check if root can run it.
#if [[ $EUID -ne 0 ]]; then
#   exit 1
# else

#collect all semaphore ID
for SEMID in `ipcs -s|egrep -v -e "Semaphore|key"|sed '/^$/d'|awk '{print $2}'|sort -u`
   do
   #GETPID of semaphore
   PID=`ipcs -s -i $SEMID | tail -2 | head -1 | awk '{print $NF}'`
   #GET PROCESS ID
   #Ignore process ID 0 which is main process & test PID greater than 0
   echo $SEMID $PID

   if [ $PID -gt 0 ]; then

   #Test of PID exists in process list, if exits then don't do anything.

   if ps -p $PID > /dev/null; then
     #running process are
     echo "$SEMID   $PID" &>/dev/null; else
     # dead process are, kill corresponding semaphore of related PID is not exisitng.
     echo "$SEMID   $PID" &>/dev/null 

    #cleaning semaphore of dead process:
    echo "ipcrm -s $SEMID"
    ipcrm -s $SEMID
   fi
   fi
   done
#fi

ME=`whoami`

IPCS_S=`ipcs -s | egrep "0x[0-9a-f]+ [0-9]+" | grep $ME | cut -f2 -d" "`
IPCS_M=`ipcs -m | egrep "0x[0-9a-f]+ [0-9]+" | grep $ME | cut -f2 -d" "`
IPCS_Q=`ipcs -q | egrep "0x[0-9a-f]+ [0-9]+" | grep $ME | cut -f2 -d" "`

for id in $IPCS_M; do
  ipcrm -m $id;
done

for id in $IPCS_S; do
  ipcrm -s $id;
done

for id in $IPCS_Q; do
  ipcrm -q $id;
done
