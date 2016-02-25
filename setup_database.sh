#!/bin/bash

set -e

PGUSER=`whoami`
PGDATABASE="bshvl"
PGHOST="localhost"
PGPORT=5432

echo "Creating database superuser $PGUSER... will require sudo privileges"
sudo -u postgres createuser -s $PGUSER

echo

echo "Creating database $PGDATABASE"
createdb -U $PGUSER -p $PGPORT -h $PGHOST $PGDATABASE

echo
echo "Done setting up database"