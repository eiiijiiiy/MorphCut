#!/bin/sh
#SBATCH --job-name=5                       # Job name
#SBATCH --mail-type=END,FAIL                       # Mail events
#SBATCH --mail-user=yijiewu@connect.hku.hk      # Update your email address
#SBATCH --time=0-01:00:00                          # Wall time limit (days-hrs:min:sec)
#SBATCH --qos=normal
#SBATCH --partition=intel
#SBATCH --nodes=1                                  # Number of Nodes
#SBATCH --ntasks-per-node=1                        # Number of tasks (MPI processes)
#SBATCH --cpus-per-task=32                         # Number of threads per task (OMP threads)
#SBATCH --output=%x.out.%j                         # Standard output & error files


module load gmp mpfr

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

/home/yijiewu/code/pp/build/bcd /scr/u/yijiewu/BCD/test_sample_fp_opt/5/B35427175530106011.off