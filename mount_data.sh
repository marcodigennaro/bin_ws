#!/bin/bash

#psw TmE2025
#sudo mount -v -t cifs //10.100.192.1/data /data -o user=$(whoami),uid=$UID,domain=TUX-NET
#sudo mount -t cifs -o user=mdi0316,vers=1.0 //10.100.192.1/data /data
sudo mount -v -t cifs //10.100.192.1/data /data -o user=$(whoami),uid=$UID,vers=1.0
