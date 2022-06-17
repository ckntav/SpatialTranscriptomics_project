#!/bin/sh
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --account=def-stbil30
#SBATCH --mail-user=christophe.tav@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name=spaceranger_count_A1
#SBATCH --output=/home/chris11/projects/def-stbil30/chris11/SpatialTranscriptomics_project/output/log/spaceranger_count_A1_log.txt


spaceranger count --reorient-images --id=sampleA1 --transcriptome=/home/chris11/projects/def-stbil30/chris11/SpatialTranscriptomics_project/input/refdata-gex-GRCh38-2020-A --fastqs=/home/chris11/projects/def-stbil30/chris11/SpatialTranscriptomics_project/raw/fastq --sample=TEP175_PG1_A1 --image=/home/chris11/projects/def-stbil30/chris11/SpatialTranscriptomics_project/raw/image/ST01-20210127-09-Stitching-01_s2-A1-1079-PG1.tif --slidefile=/home/chris11/projects/def-stbil30/chris11/SpatialTranscriptomics_project/input/V10J27-102.gpr --slide=V10J27-102 --area=A1 --localcores=32 --localmem=128
