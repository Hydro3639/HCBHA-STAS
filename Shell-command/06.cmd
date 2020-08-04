# LRs mapping and filtering 70 && 70
  minimap2 -x map-ont -t 40 re-bins.fasta ../01-Pre/test_lr.fastq > re-bins-mapping.lr.paf 
  awk -F'[\t]' '($4-$3+1)/$2 >=0.70 && $12>=30 && $27<0.30 {print $1"\t"$6}' re-bins-mapping.lr.paf > filtered_70_70-lr.paf

# SRs mapping and filtering 80 && 80
  minimap2 -x sr -t 40 re-bins.fasta ../01-Pre/test_1.fastq ../01-Pre/test_2.fastq  > re-bins-mapping.sr.paf
  awk '($13~"tp:A:P")  {print $0}' re-bins-mapping.sr.paf | awk '($4-$3+1)/$2 >=0.80 && $10/$11 >=0.80 {print $1"\t"$6}' > filtered_80_80-sr.paf

# Bin clustering
  mkdir Bin-cluster && cd Bin-cluster &&  mkdir ID LRs SRs
# Bin ID extraction
  cat ../re-binsID | while read line; do echo "grep '^>' ../../05-Re-binning/BIN_REFINEMENT/metawrap_50_10_bins/$line | sed -e 's/^>//g' > ID/${line}_ID "; done > binID-extra.cmd
  parallel -j40 < binID-extra.cmd
  
# for each rebined bin, extract the LRs ID which mapped onto the given bin, and then extract the corresponding sequences using the ID
  cat ../re-binsID | while read line; do echo "grep -wFf ID/${line}_ID ../filtered_70_70-lr.paf | awk '{print \$1}' > LRs/${line}_70_70.lr.ID "; done > LRs-eachBin-ID-extra.cmd
  parallel -j40 < LRs-eachBin-ID-extra.cmd 
   
  ls -lht LRs/*ID |awk '$5!=0' | awk '{print $NF}'| cut -d"/" -f2 | cut -d_ -f1 > Non-zero-Refined-binsID
   
  cat Non-zero-Refined-binsID | while read line; do echo "seqtk subseq ../../01-Pre/H5-LRs.fastq LRs/${line}_70_70.lr.ID > LRs/${line}_70_70.lr.fastq "; done >LRs-eachBin-Seq-extra.cmd
  parallel -j40 < LRs-eachBin-Seq-extra.cmd

# each Bin mapped SRs extraction: firstly extract mapped all SRs ID then paired SRs sequences
  cat Non-zero-Refined-binsID |  while read line; do echo "grep -wFf ID/${line}_ID ../filtered_80_80-sr.paf | awk '{print \$1}' > SRs/${line}_80_80.sr.ID" ;done > SRs-eachBin-ID-extra.cmd
  parallel -j40 < SRs-eachBin-ID-extra.cmd

   cd SRs && mkdir tmp
   ls *ID | while read line; do echo "cat ${line} | sort | uniq > tmp/${line}"; done > de-duplicate.cmd
   parallel -j40 < de-duplicate.cmd
   cd ../
# 1M subsampled Paired SRs extraction
   cat Non-zero-Refined-binsID | while read line; do echo "cat SRs/tmp/${line}_80_80.sr.ID | shuf | head -1000000 > SRs/${line}_80_80.sr_1M.ID" ; done > SRs-1M-extra.cmd
   parallel -j40 < SRs-1M-extra.cmd

   cat Non-zero-Refined-binsID | while read line; do echo "seqtk subseq ../../01-Pre/H5-UASB_1.fastq SRs/${line}_80_80.sr_1M.ID  > SRs/${line}_80_80.sr_1.fastq "; done > SRs-eachBin-seq_1-extra.cmd
   cat Non-zero-Refined-binsID | while read line; do echo "seqtk subseq ../../01-Pre/H5-UASB_2.fastq SRs/${line}_80_80.sr_1M.ID  > SRs/${line}_80_80.sr_2.fastq "; done > SRs-eachBin-seq_2-extra.cmd

   parallel -j40 < SRs-eachBin-seq_1-extra.cmd
   parallel -j40 < SRs-eachBin-seq_2-extra.cmd

# polishing
   mkdir Unicycler

   cat Non-zero-Refined-binsID | while read line; do echo "unicycler-runner.py --no_correct -1 SRs/${line}_80_80.sr_1.fastq -2 SRs/${line}_80_80.sr_2.fastq -l LRs/${line}_70_70.lr.fastq -t 15 --min_fasta_length 1000 -o Unicycler/${line}-unicycler"; done > re-assembly.cmd
   parallel -j3 < re-assembly.cmd

## prepare for the next step
   mkdir polish_02
   cd Unicycler 
   find -name assembly.fasta > tmp-ID
   cat tmp-ID | cut -d/ -f2-3 | tr "/" "_" | sed -e 's/^/..\/polish_02\//g' | paste tmp-ID - | sed -e 's/^/cp /g' > cp.cmd
   parallel -j30 < cp.cmd
## rename header of the each bin
   cd ../polish_02
   ls *fasta | cut -d- -f1 | cut -d. -f1-2 |while read line; do echo "sed -i 's/^>/>polished_02_${line}_/g' ${line}.fa-unicycler_assembly.fasta"; done > rename-header.cmd
   parallel -j40 < rename-header.cmd
   cat *fasta > ../polished_02.fasta
   cp ../polished_02.fasta ../../../07-Final-binning
