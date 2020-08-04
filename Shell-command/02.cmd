# Flye assembly
  cd 02-Flye
  flye --nano-raw ../01-Pre/test_lr.fastq -t 40 -i 5 -g 5m -o flye-polish --meta
## filter flye-assembled contigs length using 1 Kbp
  cat flye-polish/assembly.fasta | seqkit seq -m 1000 -o flye_assembled_len1K.fasta
