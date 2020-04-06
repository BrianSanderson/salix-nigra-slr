## Probe candidate polymorphism

A Python script to screen SNP and indel polymorphism in a VCF, used here to design a targeted sequence capture array

### Python 3 dependencies

* pyvcf (imported as vcf)
* tqdm

### Required parameters

* -c, --candidate_file: path to a tab-delimited text file that describes coordinates of interest for probe design, with the following columns:

  0. chromosome or contig name
  1. the start coordinate
  2. the stop coordinate
  3. the length of the region
  4. the name of the locus (gene, etc.)

* -i, --vcf_file: path a VCF that summarizes the variants from a set of alignments (in this case we used samtools mpileup). Can be in compressed (e.g. .vcf.gz)

*  -o, --out_file: path to write a tab-delimited file that describes the number of SNPs and indels in each candidate region

### Optional parameters

* --num_cand: the number of lines in the candidate file, which is used to provide more a informative status bar

* --num_vcf: the number of entries in the VCF, which is used to provide a more informative status bar
