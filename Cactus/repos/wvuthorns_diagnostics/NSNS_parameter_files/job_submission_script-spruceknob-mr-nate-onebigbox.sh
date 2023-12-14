#!/bin/sh

#This is an example script for executing generic jobs with 
# the use of the command 'qsub <name of this script>'


#These commands set up the Grid Environment for your job.  Words surrounding by a backet ('<','>') should be changed
#Any of the PBS directives can be commented out by placing another pound sign in front
#example
##PBS -N name
#The above line will be skipped by qsub because of the two consecutive # signs 

# Specify job name
#PBS -N G2e14GTRDN

# Specify the resources need for the job
# Walltime is specified as hh:mm:ss (hours:minutes:seconds)
####PBS -l nodes=2:ppn=<number_of_processors_per_node,walltime=<time_needed_by_job>
### set vmem equal to the number of nodes times 64GB (medium memory node)
#####PBS -l nodes=2:ppn=16,walltime=1:00:00,pvmem=3gb
###PBS -l nodes=2:ppn=32,walltime=1:00:00,pvmem=32gb
###PBS -l nodes=1:ppn=32,walltime=1:00:00,pvmem=2gb
####PBS -l nodes=8:ppn=32,walltime=1:00:00,pvmem=2gb

##PBS -l nodes=4:ppn=4,walltime=70:00:00,pvmem=12gb
#PBS -l nodes=8:ppn=4,walltime=4:00:00,pvmem=15gb
#PBS -n


# Specify when Moab should send e-mails 'ae' below user will
# receive e-mail for any errors with the job and/or upon completion
# If you don't want e-mails just comment out these next two PBS lines
#PBS -m e
#PBS -M noemail@hpc.wvu.edu 

# Specify the queue to execute task in. Current options can be found by excuting the command qstat -q at the terminal
#PBS -q standby
##PBS -q comm_mmem_week
###PBS -q debug
####PBS -q stmcwilliams
#PBS -W x=nmatchpolicy:exactnode

module load compilers/intel/13.0.1     
module load mpi/intel/4.1.0.024       
module load libraries/hdf5/1.8.13_intel

which mpirun

cd $HOME
rm -f hostfile 
touch hostfile
cat $PBS_NODEFILE | sort -u >> hostfile.$PBS_JOBID

# get number of MPI processes 
#export NPROCS=`wc -l hostfile | awk ' { print $1 } '`h
export NPROCS=16  # This is total MPI processes so if 2 nodes and 4 MPI processes per node this variable will be 8

# Enter your command below with arguments just as if you where going to execute on the command line
# It is generally good practice to issue a 'cd' command into the directory that contains the files
# you want to use or use full path names
cd /scratch/zbetienne/G2_I14vs14_D5R33_60km-biggerbox-gitrdone-eta0.14/
OMP_NUM_THREADS=8 mpirun -envall -n $NPROCS -machinefile $HOME/hostfile.$PBS_JOBID ./cactus_etilgrmhd -reo nsns_test.par
