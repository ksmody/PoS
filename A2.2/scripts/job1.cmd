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
#@ initialdir = $(home)/2.2
#@ output = $(home)/2.2/jobs/job$(jobid).out
#@ error = $(home)/2.2/jobs/job$(jobid).err
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

mpiexec -n 9 ./gccg drall dual oneread
mpiexec -n 9 ./gccg pent nodal oneread
mpiexec -n 9 ./gccg pent classic oneread
