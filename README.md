Running SVfinder
================
###STEP 0: Install Python Packages
    Scipy (http://scipy.org/)
    Numpy (http://www.numpy.org/)
    HTSeq (http://www-huber.embl.de/users/anders/HTSeq)
###STEP 1: Classifying mapped paired reads/tags
    python classify_reads.py paired_samfile target_samfile max_inster_size

###STEP 2: identifying inter/intra-chromosomal genomic rearrangement
	python SVfinder.py -i <discordant_reads.sam> -o <output.txt> [opts] 
#### Options:
	-n <int>    :cutoff of number of discordant pais to define a cluster  (default:2)'
	-l <int>    :extention length to join overlaped reads together (default:1000)'
	-r <int>    :read length (default:100)'
	-g <int>    :gene annotation file (default:hg19.ucsc.gene.txt)'
	-h          :produce this menu'

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
	
*chr_num column should only contain the number, no 'chr' label. Use 23 instead of X and 24 for Y. Chrom should noly be in 1-24.*
