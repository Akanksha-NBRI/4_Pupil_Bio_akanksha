#!/bin/bash -l
#!/bin/bash -l
#SBATCH -D /scratch/Meth/  
#SBATCH -o /scratch/Meth/logs/Meth-Log-%j.txt
#SBATCH -e /scratch/Meth/logs/Meth-Log-%j.err
#SBATCH -t 10:00:00
#SBATCH -J Meth
#SBATCH --partition=smp
#SBATCH --mem 96g
#SBATCH --mail-type=ALL # if you want emails, otherwise remove
#SBATCH --account=ag-stetter
#SBATCH --mail-user=asingh3@uni-koeln.de # receive an email with updates

module unload compiler/GCC
module load lang/R/

cd Test_Summary/
#R --vanilla -f PMP.R
R --vanilla -f PMP_Rerun.R






