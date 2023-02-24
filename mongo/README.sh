# ###############
# Load modules
###############

0) in ~/.bashrc

export LC_ALL=C
module use /home/mdi0316/local/modulefiles
module load mongo/5.1.1

########################################
# Starting and configuring the database:
########################################

1) Starting mongod without any admin/password:
mongod --dbpath ~/mongodb/db/ --port 27017

2) Connect to the database:
mongo --port 27017

3) Setup admin user and password in the mongo shell:
use admin
db.createUser({ user : 'mdi0316', pwd : 'sfws4438EYRW', roles: [{role:"userAdminAnyDatabase", db:"admin" }]})

4) stop db in 1)

5) Starting mongod as a daemon (fork)
mongod --fork --bind_ip_all --logpath mongodb/log/mongod.log --port 27017 --dbpath mongodb/db --auth

6) Now, admin can connect using:
mongo --port 27017 -u "mdi0316" -p 'sfws4438EYRW' --authenticationDatabase "admin"

7) In the same mongo shell
# create user and password for the database:
use ionic_liquids
db.createUser( { user : 'mdi0316', pwd : 'rhew3621KRHE', roles: [ { role: "dbOwner", db: "ionic_liquids" } ] } )

###########################
# kill previous instances #
###########################

sudo lsof -iTCP -sTCP:LISTEN -n -P | grep mongo
ps -ef | grep mongo
ls -l /proc/12054/exe
sudo kill 12054

###########################
# check firewall on cluster
###########################

systemctl status ufw
systemctl stop us

netstat -nplt | grep 27017
lsof | grep 27017

telnet 10.100.192.1 27017

###########################

mongod --dbpath ~/mongodb/ --bind_ip_all
mongo "mongodb://10.100.192.47:27017"
mongo "mongodb://10.100.192.1:27017"
mongod --dbpath ~/DATA/mongodb --bind_ip 10.100.192.47

###########################
start instance (after mongod deamon is running)
###########################

mongo fireworks_marcodigennaro -u marcodigennaro -p 'M@rc0!D1!G3nn@r0'
mongo admin -u mdg_admin -p 'n1md@!gdM'

#create admin
/home/dwa7192/TEMPORARY/local/setup_database/README.txt

