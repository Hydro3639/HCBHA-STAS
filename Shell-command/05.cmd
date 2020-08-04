# This step using MetaWRAP: the polished raw bins from the above step, might contain mis-grouping, so all the polished initial bins would be concatenate as assembled contigs (polished_01.fasta) combined with shor reads from UASB, to do the binning again by MetaWRAP; pwd: Re-binning
  metawrap binning -o INITIAL_BINNING -t 40 -a polished_01.fasta --metabat2 --maxbin2 ../01-Pre/H1-SRs/test-*fastq 
  metawrap bin_refinement -o BIN_REFINEMENT -c 50 -x 10 -t 40 -A INITIAL_BINNING/metabat2_bins/ -B INITIAL_BINNING/maxbin2_bins

  # concatenate all the bins (-c 50 -x 10) and move to the next step
  cd BIN_REFINEMENT/metawrap_50_10_bins  
  cat *fa > ../../re-bins.fasta
  ls *fa > ../../re-binsID
  cd ../../
  cp re-bins.fasta re-binsID ../06-Re-assembly
