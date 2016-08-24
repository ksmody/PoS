#!/bin/bash
# DO NOT USE environment = COPY_ALL
#@ job_type = MPICH
#@ class = test
#@ node = 1
#@ total_tasks = 16
#@ island_count = 1
#@ energy_policy_tag = NONE
#@ wall_clock_limit = 00:20:00
#@ job_name = test
#@ network.MPI = sn_all,not_shared,us
#@ initialdir = $(home)/code
#@ output = job$(jobid).out
#@ error = job$(jobid).err
#@ notification=always
#@ notify_user=kunal.mody@tum.de
#@ queue
. /etc/profile
. /etc/profile.d/modules.sh
#setup of environment
module unload mpi.ibm
module load mpi.intel
module load papi
module load metis
perf_off

################ 4 PROCS ########################
## 
########### DRALL -- ALL VARIATIONS
mpiexec -n 4 ./gccg drall.geo.bin dual oneread


mpiexec -n 4 ./gccg drall.geo.bin nodal oneread


mpiexec -n 4 ./gccg drall.geo.bin classic oneread


mpiexec -n 4 ./gccg drall.geo.bin dual allread


mpiexec -n 4 ./gccg drall.geo.bin nodal allread


mpiexec -n 4 ./gccg drall.geo.bin classic allread