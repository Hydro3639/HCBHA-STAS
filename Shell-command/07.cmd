 metawrap binning -o INITIAL_BINNING -t 40 -a polished_02.fasta --metabat2 --maxbin2 ../01-Pre/H1-SRs/test-*fastq
 metawrap bin_refinement -o BIN_REFINEMENT -c 50 -x 10 -t 40 -A INITIAL_BINNING/metabat2_bins/ -B INITIAL_BINNING/maxbin2_bins
