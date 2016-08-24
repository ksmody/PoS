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

mpiexec -n 16 ./gccg drall dual oneread
mpiexec -n 16 ./gccg drall nodal allread
mpiexec -n 16 ./gccg drall classic oneread

mpiexec -n 16 ./gccg pent dual allread
mpiexec -n 16 ./gccg pent nodal allread
mpiexec -n 16 ./gccg pent classic oneread

mpiexec -n 16 ./gccg cojack dual oneread
mpiexec -n 16 ./gccg cojack nodal oneread
mpiexec -n 16 ./gccg cojack classic allread