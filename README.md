### Brief introduction
* HCBHA stands for the Hierachical Clustering Based Hybrid Assembly, which could be applied for the `complete` and  `high-quality/high-contiguity` genome reconstruction from  highly complex ecosystems by integrating the iterative strategy (`hybrid assembly, bining and reads mapping`). The genome reconstruction methods require both the short reads and long reads from the same sample. Please note that the iterative strategy is well described in our previous work [`High-Quality Bacterial Genomes of a Partial-Nitritation/Anammox System by an Iterative Hybrid Assembly Method`](https://github.com/Hydro3639/Iterative-Hybrid-Assembly-for-enrichment-system "https://github.com/Hydro3639/Iterative-Hybrid-Assembly-for-enrichment-system")

* The HCBHA workflow is a `haplotype-resolved genome reconstruction approach`. The key idea is phasing short and long reads into candidate bacterial/archaeal haplotypes and then hybrid assembling these genomes individually. We designed a seven-step workflow to achieve this goal: 1) short and long reads preparation; 2) long-read-only *de novo* assembly, aiming to generate high-contiguity error-prone contigs; 3) initial binning, aiming to obtain candidate genome bins; 4) initial bins polish, aiming to improve contigs accuracy; 5) re-binning, aiming to improve bins accuracy; 6) re-assembly, aiming to further improve contigs contiguity and 7) final binning, aiming to obtain the high-accurate MAGs.

* For instructions, usage and explanations of the HCBHA approach, please refer to the [`HCBHA wiki`](https://github.com/Hydro3639/HCBHA-STAS/blob/master/HCBHA%20wiki.md "HCBHA wiki")!
* Aside from the evaluation of the proposed iterative HCBHA approach using the Mock dasets, we also demostrated the performance of the approach using a highly complex metagenomic sample, the Activated Sludge sample from Shatin wastewater treatment plant, Hong Kong.
* If you are using the HCBHA workflow for the genome reconstruction, please cite the corresponding reference papers, particullary the paper described the [metaFlye](https://github.com/fenderglass/Flye), [Unicycler](https://github.com/rrwick/Unicycler)
  * [metaFlye: scalable long-read metagenome assembly using repeat graphs](https://www.nature.com/articles/s41592-020-00971-x)
  * [Unicycler: Resolving bacterial genome assemblies from short and long sequencing reads](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005595) <br>
  * [MetaWRAP—a flexible pipeline for genome-resolved metagenomic data analysis](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0541-1)
  * [Minimap2: pairwise alignment for nucleotide sequences](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778)
  * [Fast gapped-read alignment with Bowtie 2](https://www.nature.com/articles/nmeth.1923)
  * [seqtk](https://github.com/lh3/seqtk) (https://github.com/lh3/seqtk)
  * [SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation](https://github.com/shenwei356/seqkit)
#### Citing our paper
* If you found this iterative haplotype-resolved framework useful in your research, a citation would be appreciated! <br>
* Charting the complexity of an activated sludge microbiome through a hybrid sequencing strategy <br>

* [High-Quality Bacterial Genomes of a Partial-Nitritation/Anammox System by an Iterative Hybrid Assembly Method](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00937-3) <br>


