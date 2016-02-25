# bshvl (BioShovel)

DeepDive application for biological text mining

## Setup

*tested on Ubuntu 14.04 with PostgreSQL 9.3*

* Install DeepDive v0.7.1 and PostgreSQL 9.3+ using the installation guide at http://deepdive.stanford.edu/
* Set `$DEEPDIVE_HOME` environment variable to DeepDive application directory in `.bashrc` or elsewhere
* Clone repository to `$DEEPDIVE_HOME/app/bshvl`
* Create DB and DB superuser by running the bash script `setup_database.sh` **(requires sudo privileges)**
* Set up data directory
    * Download PLOS Open Dataset corpus **(CoNLL format)** from http://deepdive.stanford.edu/opendata/#plos-public-library-of-science
    * Extract file and move `plos_full_conll` directory to `data/`
* Set `$APP_HOME` to `app/bshvl`
* Replace `util` with a symlink to `~/local/util` by running `cd $APP_HOME && ln -s ~/local/util .`
* Run DeepDive by running `./run.sh`