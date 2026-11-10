
#$ -j y
#$ -cwd
#$ -S /bin/bash
#$ -N cellranger_count
#$ -q all.q
#$ -pe smp 24
#$ -V

cd /igm/projects/250305_GSL-SH-4314/count

module purge

module load cellranger_8.0.0

cellranger multi \
--id=FF033_HS06 \
--csv=FF033_HS06.csv \
--localcores=24 \
--localmem=64 \
