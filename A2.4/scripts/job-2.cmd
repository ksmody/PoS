#!/bin/bash
# DO NOT USE environment = COPY_ALL
#@ job_name = A2.4-03
#@ job_type = parallel
#@ class = test
#@ node = 1
#@ total_tasks = 2
#@ island_count = 1
#@ energy_policy_tag = NONE
#@ network.MPI = sn_all,not_shared,us
#@ wall_clock_limit = 00:30:00
#@ initialdir = $(home)/A2.4/code3
#@ output = $(home)/A2.4/code3/jobs/job$(jobid).out
#@ error = $(home)/A2.4/code3/jobs/job$(jobid).err
#@ notification=always
#@ notify_user=kunal.mody@tum.de
#@ queue
. /etc/profile
. /etc/profile.d/modules.sh
#setup of environment
perf_off
module unload mpi.intel
module load mpi.ibm
module load scalasca
module load metis
module load scorep
module unload papi
module load papi/5.3
module load cube

export MP_TASK_AFFINITY=cpu
export SCOREP_ENABLE_TRACING="true"
export SCOREP_ENABLE_PROFILING="true"
export SCOREP_METRIC_PAPI=PAPI_FP_OPS

for input in "pent" "cojack" "drall"
do
	for ptype in "oneread" "allread"
	do
		for partition in "classic" "nodal" "dual"
		do
			for run in "run1" "run2" "run3"
			do
export SCOREP_EXPERIMENT_DIRECTORY="SCOREP/01/${input}/${ptype}/${partition}/02/${run}"
mpiexec -n 2 ./gccg $input $partition $ptype
echo "\n $input $partition $ptype 2"
			done
		done
	done
done