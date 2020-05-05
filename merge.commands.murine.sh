#!/bin/bash
#SBATCH --job-name=MSQ105.new								# Job name
#SBATCH --mail-type=BEGIN,END,FAIL       						# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=benjamin.wu@nyumc.org					        # Where to send mail
#SBATCH --ntasks=16                    							# Run on a single CPU
#SBATCH --mem=64gb                     							# Job memory request
#SBATCH --time=72:05:00              							# Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log   							# Standard output and error log


module add miniconda2/4.5.4
source activate qiime2-2019.7

mkdir QIIME2_8_taxonomy
mkdir QIIME2_6_tree

qiime feature-table merge \
  --i-tables ~/QIIME2_4_DADA2/no-miss-table.qza \
  --i-tables ~/QIIME2_4_DADA2/no-miss-table.qza \
  --i-tables ~/QIIME2_4_DADA2/no-miss-table.qza \
  --i-tables ~/QIIME2_4_DADA2/no-miss-table-dada2_norev.qza \
  --i-tables ~/QIIME2_4_DADA2/no-miss-table-dada2_norev.qza \
  --o-merged-table ~/QIIME2/02.merge.mice/QIIME2_4_DADA2/merge.table.qza

qiime feature-table merge-seqs \
  --i-data ~/QIIME2_4_DADA2/hits_quality.qza \
  --i-data ~/MSQ100/QIIME2_4_DADA2/hits_quality.qza \
  --i-data ~/QIIME2_4_DADA2/hits_quality.qza \
  --i-data ~/QIIME2_4_DADA2/hits_quality.qza \
  --i-data ~/QIIME2_4_DADA2/hits_quality.qza \
  --o-merged-data ~/02.merge.mice/QIIME2_4_DADA2/merge.rep-seqs.qza

qiime feature-table filter-samples \
  --i-table ~/QIIME2_4_DADA2/merge.table.qza \
  --m-metadata-file ~/QIIME2/01.map/MSQ.mouse.master.map.txt \
  --o-filtered-table ~/02.merge.mice/QIIME2_4_DADA2/filtered.merge.table.qza

qiime feature-table filter-seqs \
  --i-data ~/QIIME2_4_DADA2/merge.rep-seqs.qza \
  --i-table ~/QIIME2_4_DADA2/filtered.merge.table.qza \
  --o-filtered-data ~/02.merge.mice/QIIME2_4_DADA2/filtered.merge.rep-seqs.qza

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ~/02.merge.mice/QIIME2_4_DADA2/filtered.merge.rep-seqs.qza \
  --o-alignment ~/02.merge.mice/QIIME2_6_tree/aligned-rep-seqs.qza \
  --o-masked-alignment ~/02.merge.mice/QIIME2_6_tree/masked-aligned-rep-seqs.qza \
  --o-tree ~/02.merge.mice/QIIME2_6_tree/unrooted-tree.qza \
  --o-rooted-tree ~/02.merge.mice/QIIME2_6_tree/rooted-tree.qza

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny ~/02.merge.mice/QIIME2_6_tree/rooted-tree.qza \
  --i-table ~/02.merge.mice/QIIME2_4_DADA2/filtered.merge.table.qza \
  --p-sampling-depth 1000 \
  --m-metadata-file ~/QIIME2/01.map/MSQ.mouse.master.map.txt \
  --output-dir ~/02.merge.mice/QIIME2_7_core_metrics

qiime feature-classifier classify-sklearn \
  --i-classifier ~/02.merge.mice/greengenes/99_classifier.qza \
  --i-reads ~/02.merge.mice/QIIME2_4_DADA2/filtered.merge.rep-seqs.qza \
  --o-classification ~/02.merge.mice/QIIME2_8_taxonomy/taxonomy.qza

qiime metadata tabulate \
  --m-input-file ~/02.merge.mice/QIIME2_8_taxonomy/taxonomy.qza \
  --o-visualization ~/02.merge.mice/QIIME2_8_taxonomy/taxonomy.qzv

qiime taxa barplot \
  --i-table ~/02.merge.mice/QIIME2_4_DADA2/filtered.merge.table.qza \
  --i-taxonomy ~/02.merge.mice/QIIME2_8_taxonomy/taxonomy.qza \
  --m-metadata-file ~/01.map/MSQ.mouse.master.map.txt \
  --o-visualization ~/02.merge.mice/QIIME2_8_taxonomy/taxa-bar-plots.qzv