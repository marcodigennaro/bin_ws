#!/usr/bin/sh

#kill `ps -ef | grep mongo | awk '{print $2}'`
version=$1 #5.1.1
mongod=/home/mdi0316/anaconda3/envs/mongo_${version}/bin/mongod

mongo_dir=/home/mdi0316/mongodb/${version}/
db_path=${mongo_dir}/db
log_path=${mongo_dir}/log

rm -r $mongo_dir

$mongod --fork --bind_ip_all --logpath $log_path --port 27017 --dbpath $db_path --auth 
