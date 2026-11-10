#$ -q all.q
#$ -j y
#$ -cwd
#$ -S /bin/bash
#$ -terse
#$ -pe smp 16
#$ -N cellranger_mkfastq

module load cellranger_8.0.0

cellranger mkfastq \
--samplesheet=SampleSheet.csv \
--run=/igm/runs/IGM_Seq19/240815_A01088_0816_AHHG2GDRX5 >& bclconv.out
