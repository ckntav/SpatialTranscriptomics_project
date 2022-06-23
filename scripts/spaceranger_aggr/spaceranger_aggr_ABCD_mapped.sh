#!/bin/sh
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --account=def-stbil30
#SBATCH --mail-user=christophe.tav@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name=spaceranger_aggr_ABCD_mapped
#SBATCH --output=/home/chris11/projects/def-stbil30/chris11/SpatialTranscriptomics_project/output/log/spaceranger_aggr_ABCD_mapped_log.txt

spaceranger aggr --id=aggr_ABCD_mapped \
                 --csv=/home/chris11/projects/def-stbil30/chris11/SpatialTranscriptomics_project/input/ABCD_aggrcsv.csv \
                 --normalize=mapped