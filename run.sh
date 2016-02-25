#! /usr/bin/env bash

export APP_HOME=/home/sandip/deepdive-0.7.1/app/bshvl

echo "Running application at $APP_HOME"

# Database Configuration
export PGUSER=`whoami`
export PGPASSWORD=""
export PGDATABASE="bshvl"
export PGHOST="localhost"
export PGPORT=5432

export PARALLELISM="4"
export MEMORY="64g"
export JAVA_OPTS="-Xmx"$MEMORY
export SBT_OPTS="-Xmx"$MEMORY

export LD_LIBRARY_PATH=[DEEPDIVE_HOME]/util/dw_linux/lib/protobuf/lib:[DEEPDIVE_HOME]/util/dw_linux/lib64
#export PYTHONPATH=$DEEPDIVE_HOME/ddlib/:$PYTHONPATH

cd $APP_HOME

deepdive initdb
deepdive run
