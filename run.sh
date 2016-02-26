#! /usr/bin/env bash

export APP_HOME=/home/sandip/deepdive-0.7.1/app/bshvl
export PATH=~/local/bin:$PATH

cd $APP_HOME

# check if util symlink is valid (set to ~/local/util if not)
if [ ! -e "util" ] ; then
    echo "util directory missing; symlinking ~/local/util";
    rm util && ln -s ~/local/util .;
fi

echo "Running DeepDive located at `which deepdive`"
echo "Running DD application at $APP_HOME"

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

deepdive initdb
deepdive run
