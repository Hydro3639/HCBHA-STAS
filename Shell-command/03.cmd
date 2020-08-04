# Initial binning to reconstruct raw bins, which might contain many mis-binnings with low-accuracy contigs
  cd 03-Initial-binnig
  metawrap binning -o INITIAL_BINNING -t 40 -a ../02-Flye/flye_assembled_len1K.fasta --metabat2 --maxbin2 ../01-Pre/H1-SRs/test-*fastq 
  metawrap bin_refinement -o BIN_REFINEMENT -c 50 -x 10 -t 40 -A INITIAL_BINNING/metabat2_bins/ -B INITIAL_BINNING/maxbin2_bins

# prepare the refined bins in the folder binsAB
  mkdir Refined_bins 
  cp BIN_REFINEMENT/work_files/binsAB/Refined_*fa Refined_bins/
  cd Refined_bins
  ls *fa | while read line; do echo "sed -i 's/^>/>${line}_/g' $line "; done > rename-header.cmd
  parallel -j40 < rename-header.cmd
  cat *fa > Refined-bins.fasta && cp Refined-bins.fasta ../../04-Initial_bins_polish
  ls *fa > Refined-binsID && cp Refined-binsID ../../04-Initial_bins_polish
