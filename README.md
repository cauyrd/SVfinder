Introduction
------------
SVfinder provides genome-wide detection of structural variants from next generation paired-end sequencing reads.

Pre-installtalation
-------------------
Python Packages:
* Scipy (http://scipy.org/)
* Numpy (http://www.numpy.org/)
* HTSeq (http://www-huber.embl.de/users/anders/HTSeq)

Running SVfinder
----------------
### command-line usage
	python SVfinder.py -i <input_mapped_reads.sam> -o <output.txt> [opts] 
#### Options:
	-n <int>    :cutoff of number of discordant pais to define a cluster  (default:2)
	-l <int>    :extention length to join overlaped reads together (default:1000)
	-r <int>    :read length (default:100)
	-g <int>    :gene annotation file (default:hg19.ucsc.gene.txt)
	-h          :produce this menu

Output
-------------
The output files include a output summary file and a output BED file. the discordant reads are listed in the folder with suffix '_svreads'.

A. The summary file consists of follwoing columns:
   1. SV type (SVfinder support Insertion, Deletion, Inversion, Intra- and Inter-chromosome translocation)
   2. # of discordant read pairs
   3. Chromosome 1
   4. Position 1 start
   5. Position 1 end
   6. Orientation 1
   7. Chromosome 2
   8. Position 2 start
   9. Position 2 end
   10. Orientation 2 
   11. Annotation 1
   12. Annotation 2
   13. Putative gene fusion

B. The BED file consists of follwoing columns:
   1. Chromesome
   2. Position start
   3. Position end
   4. partener chr+start+end+strand
   5. Orientation
