## Overview

Long reads sequencing has shown the tremendous potential to resolve genomic assembly challenges e.g., long repeat regions and structural variants. The characteristics of 
longreads, on the one hand, could greatly improve the contiguity of the genome assemblies, on the other hand, assembled contigs of the error-prone long reads have a much 
higher error rates at the nucleotide level compared with the high-accuracy short-read based assemblies. Recently, a few assemblers have been developed to take the advantage
of the long reads, including the hybrid assembler Unicycler, OPERA-MS and long-read assembler, Canu and Flye, and have been documented by various studies. However, each 
assembler has its own limitations, e.g., feasibility and computing resources demand, when used as an assembler for metagenomes of high-complexity environmental samples. 
So, there is an urgent need to develop a new workflow that could resolve the above issues and be used for the high-complexity microbiome.

We firstly developed the hybrid assembly workflow for high-complexity systems and integrated the workflow into the Hierarchical Clustering Based Hybrid Assembly (HCBHA) 
approach to reconstruct more high-quality genomes. Totally, seven steps and various tools were integrated into the hybrid assembly workflow to reconstruct high-quality 
genomes from the AS sample. 

### Hybrid assembly workflow
```
mkdir 01-Pre  02-Flye  03-Initial-binning  04-Initial_bins_polish  05-Re-binning  06-Re-assembly  07-Final-binning
```
#### Step 1, Metagenomic Dataset prepare
\# both Illumina short reads and ONT long reads need prepared. e.g., test_S1_1.fastq and test_S1_2.fastq, test_S2_1.fastq and test_S2_2.fastq for short reads and test_lr.fastq for long reads
\# if your sample was highly sequenced, I strongly suggested you to subsample your data accrodingly using `seqtk` or `seqkit`. Then store the sequences to the 01-Pre folder

```
mkdir 01-Pre H1-SRs && cp test_S1_*.fastq test_S2_*.fastq 01-Pre/H1-SRs && cp test_lr.fastq 01-Pre
cat test_*_1.fastq > test_1.fastq
cat test_*_2.fastq > test_2.fastq
```
#### Step 2, Nanopore long-read assembly using Flye

\# Flye assembly
  ```
  cd 02-Flye
  flye --nano-raw ../01-Pre/test_lr.fastq -t 40 -i 5 -g 5m -o flye-polish --meta
  ```
\# filter flye-assembled contigs length using 1 Kbp
  ```
  cat flye-polish/assembly.fasta | seqkit seq -m 1000 -o flye_assembled_len1K.fasta
  ```
#### Step 3, Initial binning

\# Initial binning to reconstruct raw bins, which might contain many mis-binnings with low-accuracy contigs
  ```
  cd 03-Initial-binnig
  metawrap binning -o INITIAL_BINNING -t 40 -a ../02-Flye/flye_assembled_len1K.fasta --metabat2 --maxbin2 ../01-Pre/H1-SRs/test-*fastq 
  metawrap bin_refinement -o BIN_REFINEMENT -c 50 -x 10 -t 40 -A INITIAL_BINNING/metabat2_bins/ -B INITIAL_BINNING/maxbin2_bins
  ```
\# prepare the refined bins in the folder binsAB
  ```
  mkdir Refined_bins 
  cp BIN_REFINEMENT/work_files/binsAB/Refined_*fa Refined_bins/
  cd Refined_bins
  ls *fa | while read line; do echo "sed -i 's/^>/>${line}_/g' $line "; done > rename-header.cmd
  parallel -j40 < rename-header.cmd
  cat *fa > Refined-bins.fasta && cp Refined-bins.fasta ../../04-Initial_bins_polish
  ls *fa > Refined-binsID && cp Refined-binsID ../../04-Initial_bins_polish
  ```
#### Step 4 Initial bins polish
\# this step was to correct LRs assembled contigs grouped from the initial binning step, including mapping, filtering, ID extraction, Seq extraction and polish using unicycler
  ```
  cd ../../04-Initial_bins_polish
  ```
\# LRs mapping and filtering 70 && 70
  ```
  minimap2 -x map-ont -t 40 Refined-bins.fasta ../01-Pre/test_lr.fastq > Refined-bins-mapping.lr.paf 
  awk -F'[\t]' '($4-$3+1)/$2 >=0.70 && $27<0.30 {print $1"\t"$6}' Refined-bins-mapping.lr.paf > filtered_70_70-lr.paf
  ```
\# SRs mapping and filtering 80 && 80
  ```
  minimap2 -x sr -t 40 Refined-bins.fasta ../01-Pre/test_1.fastq ../01-Pre/test_2.fastq   > Refined-bins-mapping.sr.paf
  awk '($13~"tp:A:P")  {print $0}' Refined-bins-mapping.sr.paf | awk '($4-$3+1)/$2 >=0.80 && $10/$11 >=0.80 {print $1"\t"$6}' > filtered_80_80-sr.paf
  ```
\# Bin clustering
  ```
  mkdir Bin-cluster && cd Bin-cluster &&  mkdir ID LRs SRs
  ```
\# Bin ID extraction
  ```
  cat ../Refined-binsID | while read line; do echo "grep '^>' ../../03-Initial-binning/Refined_bins/$line | sed -e 's/^>//g' > ID/${line}_ID "; done > binID-extra.cmd
  parallel -j40 < binID-extra.cmd
  ```
\# for each refined bin, extract the LRs ID which mapped onto the given bin, and then extract the corresponding sequences using the ID
  ```
  cat ../Refined-binsID | while read line; do echo "grep ${line} ../filtered_70_70-lr.paf | awk '{print \$1}' > LRs/${line}_70_70.lr.ID "; done > LRs-eachBin-ID-extra.cmd
  parallel -j40 < LRs-eachBin-ID-extra.cmd 
  
  ls -lht LRs/*ID |awk '$5!=0' |awk '{print $NF}'| cut -d"/" -f2 | cut -d_ -f1-2 > Non-zero-Refined-binsID

  cat Non-zero-Refined-binsID  | while read line; do echo "seqtk subseq ../../01-Pre/test_lr.fastq LRs/${line}_70_70.lr.ID > LRs/${line}_70_70.lr.fastq "; done >LRs-eachBin-Seq-extra.cmd
  parallel -j40 < LRs-eachBin-Seq-extra.cmd
  ```
\# for each refined bin, extract the SRs ID which mapped onto the given bin, and then extract the corresponding sequences using the ID, for SRs, as it might extract too much 
coverage of certain bins, so only randomly selected 1 million paired SRs would be used to do the polishing
   ```
   cat Non-zero-Refined-binsID |  while read line; do echo "grep ${line} ../filtered_80_80-sr.paf | awk '{print \$1}' > SRs/${line}_80_80.sr.ID" ;done > SRs-eachBin-ID-extra.cmd
   parallel -j40 < SRs-eachBin-ID-extra.cmd
   cd SRs && mkdir tmp
   ls *ID | while read line; do echo "cut -d/ -f1 ${line} | sort | uniq > tmp/${line}"; done > de-duplicate.cmd
   parallel -j40 < de-duplicate.cmd
   cd ../
   ```
\# 1M subsampled Paired SRs extraction
   ```
   cat Non-zero-Refined-binsID | while read line; do echo "cat SRs/tmp/${line}_80_80.sr.ID | shuf | head -1000000 > SRs/${line}_80_80.sr_1M.ID" ; done > SRs-1M-extra.cmd
   parallel -j40 < SRs-1M-extra.cmd

   cat Non-zero-Refined-binsID | while read line; do echo "seqtk subseq ../../01-Pre/test_1.fastq SRs/${line}_80_80.sr_1M.ID  > SRs/${line}_80_80.sr_1.fastq "; done > SRs-eachBin-seq_1-extra.cmd
   cat Non-zero-Refined-binsID | while read line; do echo "seqtk subseq ../../01-Pre/test_2.fastq SRs/${line}_80_80.sr_1M.ID  > SRs/${line}_80_80.sr_2.fastq "; done > SRs-eachBin-seq_2-extra.cmd

   parallel -j40 < SRs-eachBin-seq_1-extra.cmd
   parallel -j40 < SRs-eachBin-seq_2-extra.cmd
   ```

\# polishing using Unicycler
   ```
   mkdir Unicycler

   cat Non-zero-Refined-binsID | while read line; do echo "unicycler-runner.py --no_correct -1 SRs/${line}_80_80.sr_1.fastq -2 SRs/${line}_80_80.sr_2.fastq -l LRs/${line}_70_70.lr.fastq -t 15 --min_fasta_length 1000 -o Unicycler/${line}-unicycler"; done > Initial_bins_polish.cmd   
   parallel -j3 < Initial_bins_polish.cmd 
   ```
\#collect all polished initial bins and rename; pwd: Initial_bins_polish/Bin-cluster
   ```
   mkdir polish-01
   find -name assembly.fasta > tmpID
   cut -d/ -f3-5 tmpID | tr "/" "_" | sed -e 's/^/polish-01\//g' | paste tmpID - | sed -e 's/^/cp /g' > tmp-cp.cmd
   parallel -j40 < tmp-cp.cmd
   cd polish-01
   ls *fasta | cut -d- -f1 | while read line; do echo "sed -i 's/^>/>polished_01_${line}_/g' ${line}-unicycler_assembly.fasta "; done  > rename-header.cmd 
   parallel -j10 < rename-header.cmd
   cat *fasta > ../polished_01.fasta
   
   cd ../
   cp polished_01.fasta ../../05-Re-binning/
  ```
#### Step 5 Re-binning
\# This step using MetaWRAP: the polished raw bins from the above step, might contain mis-grouping, so all the polished initial bins would be concatenate as assembled contigs 
(polished_01.fasta) combined with short reads, to do the binning again by MetaWRAP; pwd: Re-binning
  ```
  metawrap binning -o INITIAL_BINNING -t 40 -a polished_01.fasta --metabat2 --maxbin2 ../01-Pre/H1-SRs/test-*fastq 
  metawrap bin_refinement -o BIN_REFINEMENT -c 50 -x 10 -t 40 -A INITIAL_BINNING/metabat2_bins/ -B INITIAL_BINNING/maxbin2_bins
  ```
  \# concatenate all the bins (-c 50 -x 10) and move to the next step
  ```
  cd BIN_REFINEMENT/metawrap_50_10_bins  
  cat *fa > ../../re-bins.fasta
  ls *fa > ../../re-binsID
  cd ../../
  cp re-bins.fasta re-binsID ../06-Re-assembly
  ```
#### Step 6 Re-assembly
\#Similar with Step 4 Mapping, filtering, ID extraction, Seq extraction and polish using unicycler

\# LRs mapping and filtering 70 && 70
  ```
  minimap2 -x map-ont -t 40 re-bins.fasta ../01-Pre/test_lr.fastq > re-bins-mapping.lr.paf 
  awk -F'[\t]' '($4-$3+1)/$2 >=0.70 && $12>=30 && $27<0.30 {print $1"\t"$6}' re-bins-mapping.lr.paf > filtered_70_70-lr.paf
  ```
\# SRs mapping and filtering 80 && 80
  ```
  minimap2 -x sr -t 40 re-bins.fasta ../01-Pre/test_1.fastq ../01-Pre/test_2.fastq  > re-bins-mapping.sr.paf
  awk '($13~"tp:A:P")  {print $0}' re-bins-mapping.sr.paf | awk '($4-$3+1)/$2 >=0.80 && $10/$11 >=0.80 {print $1"\t"$6}' > filtered_80_80-sr.paf
  ```
\# Bin clustering
  ```
  mkdir Bin-cluster && cd Bin-cluster &&  mkdir ID LRs SRs
  ```
\# Bin ID extraction
  
  ```cat ../re-binsID | while read line; do echo "grep '^>' ../../05-Re-binning/BIN_REFINEMENT/metawrap_50_10_bins/$line | sed -e 's/^>//g' > ID/${line}_ID "; done > binID-extra.cmd
  parallel -j40 < binID-extra.cmd
  ```
\# for each rebined bin, extract the LRs ID which mapped onto the given bin, and then extract the corresponding sequences using the ID
  ```cat ../re-binsID | while read line; do echo "grep -wFf ID/${line}_ID ../filtered_70_70-lr.paf | awk '{print \$1}' > LRs/${line}_70_70.lr.ID "; done > LRs-eachBin-ID-extra.cmd
  parallel -j40 < LRs-eachBin-ID-extra.cmd 
   
  ls -lht LRs/*ID |awk '$5!=0' | awk '{print $NF}'| cut -d"/" -f2 | cut -d_ -f1 > Non-zero-Refined-binsID
   
  cat Non-zero-Refined-binsID | while read line; do echo "seqtk subseq ../../01-Pre/test_lr.fastq LRs/${line}_70_70.lr.ID > LRs/${line}_70_70.lr.fastq "; done >LRs-eachBin-Seq-extra.cmd
  parallel -j40 < LRs-eachBin-Seq-extra.cmd
```
\# each Bin mapped SRs extraction: firstly extract mapped all SRs ID then paired SRs sequences
  ```
  cat Non-zero-Refined-binsID |  while read line; do echo "grep -wFf ID/${line}_ID ../filtered_80_80-sr.paf | awk '{print \$1}' > SRs/${line}_80_80.sr.ID" ;done > SRs-eachBin-ID-extra.cmd
  parallel -j40 < SRs-eachBin-ID-extra.cmd

   cd SRs && mkdir tmp
   ls *ID | while read line; do echo "cat ${line} | sort | uniq > tmp/${line}"; done > de-duplicate.cmd
   parallel -j40 < de-duplicate.cmd
   cd ../
```
\# 1M subsampled Paired SRs extraction
 ```
 cat Non-zero-Refined-binsID | while read line; do echo "cat SRs/tmp/${line}_80_80.sr.ID | shuf | head -1000000 > SRs/${line}_80_80.sr_1M.ID" ; done > SRs-1M-extra.cmd
 parallel -j40 < SRs-1M-extra.cmd
 cat Non-zero-Refined-binsID | while read line; do echo "seqtk subseq ../../01-Pre/test-sr_1.fastq SRs/${line}_80_80.sr_1M.ID  > SRs/${line}_80_80.sr_1.fastq "; done > SRs-eachBin-seq_1-extra.cmd
 cat Non-zero-Refined-binsID | while read line; do echo "seqtk subseq ../../01-Pre/test-sr_2.fastq SRs/${line}_80_80.sr_1M.ID  > SRs/${line}_80_80.sr_2.fastq "; done > SRs-eachBin-seq_2-extra.cmd
  
 parallel -j40 < SRs-eachBin-seq_1-extra.cmd
 parallel -j40 < SRs-eachBin-seq_2-extra.cmd
```
\# polishing
 ```
 mkdir Unicycler

 cat Non-zero-Refined-binsID | while read line; do echo "unicycler-runner.py --no_correct -1 SRs/${line}_80_80.sr_1.fastq -2 SRs/${line}_80_80.sr_2.fastq -l LRs/${line}_70_70.lr.fastq -t 15 --min_fasta_length 1000 -o Unicycler/${line}-unicycler"; done > re-assembly.cmd
 parallel -j3 < re-assembly.cmd
```
\## prepare for the next step
 ```
 mkdir polish_02
 cd Unicycler 
 find -name assembly.fasta > tmp-ID
 cat tmp-ID | cut -d/ -f2-3 | tr "/" "_" | sed -e 's/^/..\/polish_02\//g' | paste tmp-ID - | sed -e 's/^/cp /g' > cp.cmd
 parallel -j30 < cp.cmd
```
\## rename header of the each bin
 ```
 cd ../polish_02
 ls *fasta | cut -d- -f1 | cut -d. -f1-2 |while read line; do echo "sed -i 's/^>/>polished_02_${line}_/g' ${line}.fa-unicycler_assembly.fasta"; done > rename-header.cmd
 parallel -j40 < rename-header.cmd
 cat *fasta > ../polished_02.fasta
 cp ../polished_02.fasta ../../../07-Final-binning
```
#### Step 7 Final binning
```
metawrap binning -o INITIAL_BINNING -t 40 -a polished_02.fasta --metabat2 --maxbin2 ../01-Pre/H1-SRs/test-*fastq
metawrap bin_refinement -o BIN_REFINEMENT -c 50 -x 10 -t 40 -A INITIAL_BINNING/metabat2_bins/ -B INITIAL_BINNING/maxbin2_bins
```
#########################################################################################################

* After the hybrid assembled genomes were retrieved from the system, then select the qualified MAGs, in our study, we choosed the MAGs with completeness >=90%, 
contamination <=10% and contig contig <=30. To facilitate the reconstrcution of the remaining community members, we need to take out the short and long reads
that assigned to the qualified MAGs using minimap2, then repeat the above hybrid assembly process



