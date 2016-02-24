#! /usr/bin/env bash

export APP_HOME=/home/sandip/deepdive-0.7.1/app/bshvl
#export DEEPDIVE_HOME=`cd ${APP_HOME}/../../../deepdive; pwd`

# Database Configuration
export PGUSER="sandip"
export PGPASSWORD=""
export PGDATABASE="bshvl"
export PGHOST="localhost"
export PGPORT=5432

#export PGUSER=${PGUSER:-`whoami`}
#export PGPASSWORD=${PGPASSWORD:-}
#export PGDATABASE=${PGDATABASE}
#export PGHOST=${PGHOST}
#export PGPORT=${PGPORT}
export PARALLELISM="4"
export MEMORY="64g"
export JAVA_OPTS="-Xmx"$MEMORY
export SBT_OPTS="-Xmx"$MEMORY

export LD_LIBRARY_PATH=[DEEPDIVE_HOME]/util/dw_linux/lib/protobuf/lib:[DEEPDIVE_HOME]/util/dw_linux/lib64
#export PYTHONPATH=$DEEPDIVE_HOME/ddlib/:$PYTHONPATH

#cd $DEEPDIVE_HOME
cd $APP_HOME

deepdive initdb
deepdive run

#$DEEPDIVE_HOME/sbt/sbt "run -c $APP_HOME/application.conf" $@
