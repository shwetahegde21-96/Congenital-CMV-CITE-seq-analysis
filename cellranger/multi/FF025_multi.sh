
#$ -j y
#$ -cwd
#$ -S /bin/bash
#$ -N cellranger_count
#$ -q all.q
#$ -pe smp 24
#$ -V

cd /igm/projects/240813_GSL-MS-4043/count

module purge

module load cellranger_8.0.0

cellranger multi \
--id=FF025 \
--csv=FF025_multiplex_template.csv \
--localcores=24 \
--localmem=64 \
