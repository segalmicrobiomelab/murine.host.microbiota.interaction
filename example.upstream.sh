#!/bin/bash
#SBATCH --job-name=example.upstream							# Job name
#SBATCH --mail-type=BEGIN,END,FAIL       						# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=benjamin.wu@nyumc.org					        # Where to send mail
#SBATCH --ntasks=16                    							# Run on a single CPU
#SBATCH --mem=64gb                     							# Job memory request
#SBATCH --time=72:05:00              							# Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log   							# Standard output and error log


module add miniconda2/4.5.4
source activate qiime2-2019.7


qiime tools import \
  --type EMPPairedEndSequences \
  --input-path ~/QIIME2_2_import \
  --output-path ~/QIIME2_3_demux/paired-end-demux.qza

qiime demux emp-paired \
--i-seqs ~/QIIME2_3_demux/paired-end-demux.qza \
--o-per-sample-sequences ~/QIIME2_3_demux/demultiplex.qza \
--m-barcodes-file ~/map/MSQ.99.Map.txt \
--m-barcodes-column BarcodeSequence \
--p-rev-comp-barcodes \
--p-rev-comp-mapping-barcodes \
--verbose \
--o-error-correction-details ~/QIIME2_3_demux/error.qza

qiime demux summarize \
  --i-data ~/QIIME2_3_demux/demultiplex.qza \
  --o-visualization ~/QIIME2_3_demux/demultiplex.qzv

### DADA2

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ~/QIIME2_3_demux/demultiplex.qza \
  --p-trim-left-f 13 \
  --p-trim-left-r 13 \
  --p-trunc-len-f 140 \
  --p-trunc-len-r 140 \
  --p-n-threads 0 \
  --o-table ~/QIIME2_4_DADA2/table.qza \
  --o-representative-sequences ~/QIIME2_4_DADA2/rep-seqs.qza \
  --o-denoising-stats ~/QIIME2_4_DADA2/denoising-stats.qza

qiime feature-table summarize \
  --i-table ~/QIIME2_4_DADA2/table.qza \
  --o-visualization ~/QIIME2_4_DADA2/table.qzv \
  --m-sample-metadata-file ~/map/MSQ.99.Map.txt

qiime feature-table tabulate-seqs \
  --i-data ~/QIIME2_4_DADA2/rep-seqs.qza \
  --o-visualization ~/QIIME2_4_DADA2/rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file ~/QIIME2_4_DADA2/denoising-stats.qza \
  --o-visualization ~/QIIME2_4_DADA2/denoising-stats.qzv

### Quality Control Step 

qiime quality-control exclude-seqs \
  --i-query-sequences  ~/QIIME2_4_DADA2/rep-seqs.qza \
  --i-reference-sequences   ~/greengenes/99_otus.qza \
  --p-method vsearch \
  --p-perc-identity 0.99 \
  --p-perc-query-aligned 0.99 \
  --p-threads 16 \
  --o-sequence-hits  ~/QIIME2_4_DADA2/hits_quality.qza \
  --o-sequence-misses  ~/QIIME2_4_DADA2/misses_quality.qza \
  --verbose

qiime feature-table filter-features \
  --i-table ~/QIIME2_4_DADA2/table.qza \
  --m-metadata-file ~/QIIME2_4_DADA2/misses_quality.qza  \
  --o-filtered-table ~/QIIME2_4_DADA2/no-miss-table.qza \
  --p-exclude-ids

qiime feature-table summarize \
  --i-table ~/QIIME2_4_DADA2/no-miss-table.qza  \
  --o-visualization ~/QIIME2_4_DADA2/no-miss-table.qzv \
  --m-sample-metadata-file ~/map/MSQ.99.Map.txt

qiime feature-table tabulate-seqs \
  --i-data ~/QIIME2_4_DADA2/hits_quality.qza \
  --o-visualization ~/QIIME2_4_DADA2/hits_quality.qza
