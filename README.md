Running SVfinder
================
###STEP1: Classifying mapped paired reads/tags
    python classify_reads.py paired_samfile target_samfile max_inster_size

###STEP2: identifying inter-chromosomal translocation
	python sv_for_diff_chr.py diffchr_pairs.txt output.txt read_len(optional, default 100) search_len(optional, default 1000) cluster_cutoff(optional, default 1)

###STEP3: identifiying intra-chromosomal translocation
	python sv_within_chr.py longisize.txt output.txt read_len(optional, default 100) search_len(optional, default 1000) cluster_cutoff(optional, default 1)

Annotation Script
=================
To use addcatann, you need to specify
 
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
	 
	perl /home/prj/compbio/software/utl/addcatann.pl human 10 2 2 [input file in “chr start end” format] 0 /home/prj/compbio/data/annotations/hg19/hg19.refseq.unique.annot [output file name]
	 
An example (in hg18) is
	 
	perl addcatann.pl human 100 2 2 mof3.cd4.out.hpeak.out 5 /home/prj/compbio/data/annotations/hg18/hg18.refseq.unique.annot mof3.cd4.annot
