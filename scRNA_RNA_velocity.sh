#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=80G
#SBATCH --time=50:00:00
#SBATCH --output=scRNA_RNA_velocity.sh.%A_%a.stdout
#SBATCH -p all,abild
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh


start=`date +%s`

CPU=$SLURM_NTASKS
if [ ! $CPU ]; then
   CPU=2
fi

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
    N=1
fi

echo "CPU: $CPU"
echo "N: $N"

sample=`cat samples.list | head -n $N | tail -n 1`
fq_dir=../Preprocess/VG_CAMA1_D11_count03/
GTF=/home/jichen/Projects/Database/Reference/refdata-cellranger-hg19-3.0.0/genes/genes.gtf
TenXfolder=$fq_dir/$sample/

echo "Analyze spliced and unspliced reads and genrate a loom file for velocityR"
export PATH=$PATH:/home/jichen/software/BETSY/install/envs/scRNA_velocyto/bin/
export PYTHONPATH=/home/jichen/software/BETSY/install/envs/scRNA_velocyto/lib/python3.6/site-packages/
module load samtools

velocyto run10x --samtools-threads $CPU --samtools-memory 500 $TenXfolder $GTF
 
echo ""
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA/lib/R/library
module load hdf5/1.10.3

/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_RNA_velocity_run_velocityR.R


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

