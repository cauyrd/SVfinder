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
#### SV type
#### # of discordant read pairs
#### Chromosome 1
#### Position 1 start
#### Position 1 end
#### Orientation 1
#### Chromosome 2
#### Position 2 start
#### Position 2 end
#### Orientation 2 
#### Annotation 1
#### Annotation 2
#### Putative gene fusion

B. The BED file consists of follwoing columns:
#### Chromesome
#### Position start
#### Position end
#### partener chr+start+end+strand
#### Orientation
