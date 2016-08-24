#!/bin/bash
# DO NOT USE environment = COPY_ALL
#@ job_type = MPICH
#@ class = test
#@ node = 1
#@ total_tasks = 16
#@ island_count = 1
#@ energy_policy_tag = NONE
#@ wall_clock_limit = 00:30:00
#@ job_name = test
#@ network.MPI = sn_all,not_shared,us
#@ initialdir = $(home)/code
#@ output = $(home)/code/Jobs/job$(jobid).out
#@ error = $(home)/code/Jobs/job$(jobid).err
#@ notification=always
#@ notify_user=srdjan.krivokapic@tum.de
#@ queue
. /etc/profile
. /etc/profile.d/modules.sh
#setup of environment
module unload mpi.ibm
module load mpi.intel
module load papi
module load metis
perf_off

mpiexec -n 1 ./gccg cojack classic oneread
mpiexec -n 2 ./gccg cojack dual oneread
mpiexec -n 4 ./gccg cojack dual oneread
mpiexec -n 8 ./gccg cojack dual oneread
mpiexec -n 16 ./gccg cojack dual oneread