#!/bin/bash

export TMPDIR=/localscratch/barbara.hoellbacher
export _JAVA_OPTIONS="${_JAVA_OPTIONS} -Djava.io.tmpdir=$TMPDIR"

njobs=$1

eval "$(conda shell.bash hook)"
conda activate gr_env

profile='--profile src/profiles'
if type sbatch >/dev/null 2>&1
then
  profile="$profile/slurm"
elif type qsub >/dev/null 2>&1
then
  profile="$profile/sge"
else
  if [[ ! -z $njobs ]]
  then
    if [[ $njobs =~ ^[0-9]+$ ]]
    then
      echo -e "\033[0;31mERROR:\nNo workload manager detected. Please, run the workflow in a system that allows job submissions to a cluster.\033[0m"
      exit 1
    else
      echo -e "\033[0;33mWARNING:\nNo workload manager detected. In production mode the workflow must be executed in a system that allows job submissions to a cluster.\033[0m"
    fi
  fi
fi


### Snakemake commands of interest ###
## --forcerun <rule>
## --rerun-incomplete
## --unlock


cmd="snakemake --configfile config/sample_config.yaml -s Snakefile --use-conda"
if [[ -z $njobs ]] || [[ ! $njobs =~ ^[0-9]+$ ]]
then
  echo -e "\033[0;33mWARNING:\nYou are just running a test. If you want to actually execute the workflow run this script passing the number of cores/jobs as the first parameter: \n\tsh $0 njobs [rule1 rule2 ... rulen]\033[0m\n"
  cmd="$cmd -np $@"
else
  cmd="$cmd $profile -j $@"
fi

eval "$cmd"

