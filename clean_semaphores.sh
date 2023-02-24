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

    #cleaning semaphore of dead process :
    ipcrm -s $SEMID
   fi
   fi
   done
#fi
