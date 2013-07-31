Running SVfinder
================
###STEP1: Classifying mapped paired reads/tags
    python classify_reads.py paired_samfile target_samfile max_inster_size

###STEP2: identifying inter-chromosomal translocation
	python sv_for_diff_chr.py diffchr_pairs.txt output.txt read_len(optional, default 100) search_len(optional, default 1000) cluster_cutoff(optional, default 1)

###STEP3: identifiying intra-chromosomal translocation
	python sv_within_chr.py longisize.txt output.txt read_len(optional, default 100) search_len(optional, default 1000) cluster_cutoff(optional, default 1)

Getting Annotation
==================
To use addcatann.pl, you need to specify
 
Please provide the following:

	(1) species (human or mouse),	
	(2) upstream and downstream distal distance limit (in kb),
	(3) upstream proxi distance limit (in kb),	
	(4) downstream proxi distance limit (in kb),
	(5) peak file name,
	(6) column that contains peak summit (5 for hpeak output, 0 if no summit provided),
	(7) reference TSS, TES, ATG, GTA location (sorted) and gene names file
	(8) output file name that contains distances and names.
	 
The syntax is
	 
	perl addcatann.pl human 10 2 2 [input file in “chr start end” format] 0 gene_annotation_file [output file name]
	 
An example (in hg19) is
	 
	perl addcatann.pl human 100 2 2 mof3.cd4.out.hpeak.out 5 hg19.ucsc.gene.txt output_file

Gene annotaion file format (seven columns):

	chr_num  strand  txSTART  txEND  cdsSTART  cdsEND  geneName
	
*chr_num column should only contain the number, no 'chr' label. Use 23 instead of X and 24 for Y.*
