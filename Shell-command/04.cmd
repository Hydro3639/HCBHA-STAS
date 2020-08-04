# this step was to correct LRs assembled contigs grouped from the initial binning step, including mapping, filtering, ID extraction, Seq extraction and polish using unicycler
  cd ../../04-Initial_bins_polish

# LRs mapping and filtering 70 && 70
  minimap2 -x map-ont -t 40 Refined-bins.fasta ../01-Pre/test_lr.fastq > Refined-bins-mapping.lr.paf 
  awk -F'[\t]' '($4-$3+1)/$2 >=0.70 && $12>=30 && $27<0.30 {print $1"\t"$6}' Refined-bins-mapping.lr.paf > filtered_70_70-lr.paf

# SRs mapping and filtering 80 && 80
  minimap2 -x sr -t 40 Refined-bins.fasta ../01-Pre/test_1.fastq ../01-Pre/test_2.fastq   > Refined-bins-mapping.sr.paf
  awk '($13~"tp:A:P")  {print $0}' Refined-bins-mapping.sr.paf | awk '($4-$3+1)/$2 >=0.80 && $10/$11 >=0.80 {print $1"\t"$6}' > filtered_80_80-sr.paf

# Bin clustering
  mkdir Bin-cluster && cd Bin-cluster &&  mkdir ID LRs SRs
# Bin ID extraction
  cat ../Refined-binsID | while read line; do echo "grep '^>' ../../03-Initial-binning/Refined_bins/$line | sed -e 's/^>//g' > ID/${line}_ID "; done > binID-extra.cmd
  parallel -j40 < binID-extra.cmd

# for each refined bin, extract the LRs ID which mapped onto the given bin, and then extract the corresponding sequences using the ID
  cat ../Refined-binsID | while read line; do echo "grep ${line} ../filtered_70_70-lr.paf | awk '{print \$1}' > LRs/${line}_70_70.lr.ID "; done > LRs-eachBin-ID-extra.cmd
  parallel -j40 < LRs-eachBin-ID-extra.cmd 
   
  ls -lht LRs/*ID |awk '$5!=0' |awk '{print $NF}'| cut -d"/" -f2 | cut -d_ -f1-2 > Non-zero-Refined-binsID

  cat Non-zero-Refined-binsID  | while read line; do echo "seqtk subseq ../../01-Pre/H5-LRs.fastq LRs/${line}_70_70.lr.ID > LRs/${line}_70_70.lr.fastq "; done >LRs-eachBin-Seq-extra.cmd
  parallel -j40 < LRs-eachBin-Seq-extra.cmd

# for each refined bin, extract the SRs ID which mapped onto the given bin, and then extract the corresponding sequences using the ID, for SRs, as it might extract too much coverage of certain bins, so only randomly selected 1 million paired SRs would be used to do the polishing
   cat Non-zero-Refined-binsID |  while read line; do echo "grep ${line} ../filtered_80_80-sr.paf | awk '{print \$1}' > SRs/${line}_80_80.sr.ID" ;done > SRs-eachBin-ID-extra.cmd
   parallel -j40 < SRs-eachBin-ID-extra.cmd
   cd SRs && mkdir tmp
   ls *ID | while read line; do echo "cut -d/ -f1 ${line} | sort | uniq > tmp/${line}"; done > de-duplicate.cmd
   parallel -j40 < de-duplicate.cmd
   cd ../
# 1M subsampled Paired SRs extraction
   cat Non-zero-Refined-binsID | while read line; do echo "cat SRs/tmp/${line}_80_80.sr.ID | shuf | head -1000000 > SRs/${line}_80_80.sr_1M.ID" ; done > SRs-1M-extra.cmd
   parallel -j40 < SRs-1M-extra.cmd

   cat Non-zero-Refined-binsID | while read line; do echo "seqtk subseq ../../01-Pre/test_1.fastq SRs/${line}_80_80.sr_1M.ID  > SRs/${line}_80_80.sr_1.fastq "; done > SRs-eachBin-seq_1-extra.cmd
   cat Non-zero-Refined-binsID | while read line; do echo "seqtk subseq ../../01-Pre/test_2.fastq SRs/${line}_80_80.sr_1M.ID  > SRs/${line}_80_80.sr_2.fastq "; done > SRs-eachBin-seq_2-extra.cmd

   parallel -j40 < SRs-eachBin-seq_1-extra.cmd
   parallel -j40 < SRs-eachBin-seq_2-extra.cmd


# polishing using Unicycler
   mkdir Unicycler

   cat Non-zero-Refined-binsID | while read line; do echo "unicycler-runner.py --no_correct -1 SRs/${line}_80_80.sr_1.fastq -2 SRs/${line}_80_80.sr_2.fastq -l LRs/${line}_70_70.lr.fastq -t 15 --min_fasta_length 1000 -o Unicycler/${line}-unicycler"; done > Initial_bins_polish.cmd   
   parallel -j3 < re-assembly.cmd

#collect all polished initial bins and rename; pwd: Initial_bins_polish/Bin-cluster
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
