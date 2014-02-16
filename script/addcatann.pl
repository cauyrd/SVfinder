#!/usr/bin/perl -w
#use strict;
# Steve Qin 7/01/2011
# modified from douann.pl on 10/08/2011/
# to add annotation categoryies to the peak list.

$argc = @ARGV;
$argc == 8 || die "Please provide the following: 
(1) species (human or mouse), 
(2) upstream and downstream distal distance limit (in kb),
(3) upstream proxi distance limit (in kb),
(4) downstream proxi distance limit (in kb),
(5) peak file name, 
(6) column that contains peak summit (5 for hpeak output, 0 if no summit provided),
(7) reference TSS, TES, ATG, GTA location (sorted) and gene names file 
(8) output file name that contains distances and names.\n";

$spe = $ARGV[0];
$dis = $ARGV[1];
$up = $ARGV[2];
$down0 = $ARGV[3];
$inputname = $ARGV[4];
$summitcol = $ARGV[5];
$refname = $ARGV[6];
$outputname = $ARGV[7];
open(REFFILE,$refname);

$distallimit = 1000 * $dis;
$upproxilimit = 1000 * $up;
$downproxilimit = 1000* $down0;

$count = 0;
if($spe eq "human")
{
	$chrcount = 24;
}
elsif($spe eq "mouse")
{
        $chrcount = 21;
}
else
{
	print "Sorry, only human and mouse please.\n";
       	exit(0);
}

my @counterbychr;
for($k=0;$k<$chrcount;$k++)
{
	$counterbychr[$k] = 0;
}

$fiveUTRlength = 0;
$codinglength = 0;
$threeUTRlength = 0;
while($datline = <REFFILE>)
{
        chop($datline);
#        print OUTPUTFILE "$datline";
	@items = split(" ",$datline);
	$chr = $items[0];
	$chr =~ s/chr//;
	if($chr eq "X")
	{
		if($spe eq "human")
	        {
			$refchr[$count] = 23;
		}
		elsif($chr eq "mouse")
                {        
                        $refchr[$count] = 20;
                }
		else
		{
			print "Sorry, only human and mouse please.\n";
			exit(0);
		}
	}
	elsif($chr eq "Y")
        {
		 if($spe eq "human")
                {        
                        $refchr[$count] = 24;
                }
                elsif($chr eq "mouse")
                {
                        $refchr[$count] = 21;
                }
                else
                {
                        print "Sorry, only human and mouse please.\n";
                        exit(0);
                }
        }
	else
	{
		$refchr[$count] = $chr;
	}
	$temchr = $refchr[$count] - 1;
        $refstrand[$temchr][$counterbychr[$temchr]] = $items[1];
	$reftss[$temchr][$counterbychr[$temchr]] = $items[2];
        $reftes[$temchr][$counterbychr[$temchr]] = $items[3];              
        $refatg[$temchr][$counterbychr[$temchr]] = $items[4];              
        $refgta[$temchr][$counterbychr[$temchr]] = $items[5];              
	$refnames[$temchr][$counterbychr[$temchr]] = $items[6];
	if($items[1] eq "+")
	{
		$fiveUTRlength = $fiveUTRlength + $items[4] - $items[2];
		$threeUTRlength = $threeUTRlength + $items[3] - $items[5];
	}
	else
	{
		$fiveUTRlength = $fiveUTRlength + $items[3] - $items[5];
                $threeUTRlength = $threeUTRlength + $items[4] - $items[2];
	}
	$codinglength = $codinglength + $items[5] - $items[4];
	$counterbychr[$temchr]++;
#	print "strand = $items[1];tss = $items[2];atg = $items[4].\n";
#        print "gta = $items[5];tes = $items[3].\n";
#        print "5' UTR length= $fiveUTRlength; coding length = $codinglength; 3' UTR length = $threeUTRlength.\n";
#	if($count ==1)
#	{
#		exit(0);
#	}
	$count ++;
}
print "There are $chrcount chromosomes in the reference genome.\n";
for($k=0;$k<$chrcount;$k++)                  
{                
	print "For ",$k+1,"th chromosome, there are $counterbychr[$k] genes in the reference.\n";
}  
print "There are $count unique genes in the reference.\n";
close(REFFILE);
print "Total length of 5' distal regions is ",$count * 98," kb.\n";
print "Total length of 5' proximal regions is ",$count * 2," kb.\n";
print "Total length of 5' UTR regions is ",0.001*$fiveUTRlength," kb.\n";
print "Total length of coding regions is ",0.001*$codinglength," kb.\n";
print "Total length of 3' UTR regions is ",0.001*$threeUTRlength," kb.\n";
print "Total length of 3' proximal regions is ",$count * 2," kb.\n";
print "Total length of 3' distal regions is ",$count * 98," kb.\n";
#exit(0);

$names[0]  = "5' distal";
$names[1]  = "5' proximal";
$names[2]  = "5' UTR";
$names[3]  = "coding";
$names[4]  = "3' UTR";
$names[5]  = "3' proximal";
$names[6]  = "3' distal";

open(INPUTFILE,$inputname);
open(OUTPUTFILE,">".$outputname);
open(OUTFILE1,">".$outputname.".5proxi");
open(OUTFILE2,">".$outputname.".coding");
my @anncount;
for($j=0;$j<8;$j++)
{
	$anncount[$j] = 0;
}
$totalcount = 0;
while($datline = <INPUTFILE>)
{
        chop($datline);
#        print "$datline";
	print OUTPUTFILE "$datline\t";
        @items = split(" ",$datline);
	$chr = $items[0];
	$totalcount ++;
	if($summitcol == 0)
	{
		$summit = int(0.5*($items[1] + $items[2]));
	}
	else
	{
		$summit = $items[1] + $items[$summitcol-1];
	}
#	print "$summitcol\t$summit\n";
#	exit(0);
	$current = $chr -1;
	for($j=0;$j<8;$j++)
	{
        	$select[$j] = 0;
	}
	for($j=0;$j<7;$j++)
        {
                $refseqnameofselect[$j] = " ";                 
        } 
	$hit = 0;
	for($j=0;$j < $counterbychr[$current];$j++)
	{
		if(($summit >= ($reftss[$current][$j]-$distallimit )) && ($summit <= ($reftes[$current][$j] + $distallimit )))
        	{
                        $hit = 1;
##### 
##### if "+" strand
#####
			if($refstrand[$current][$j] eq "+")
			{
				if(($summit >= ($reftss[$current][$j] - $distallimit ))&&($summit < ($reftss[$current][$j] - $upproxilimit)))
		                {
					$select[0] = 1;
					$refseqnameofselect[0] = $refnames[$current][$j];
				}
                                elsif(($summit >= ($reftss[$current][$j] - $upproxilimit ))&&($summit < $reftss[$current][$j]))
                                {
					$select[1] = 1;
                                        $refseqnameofselect[1] = $refnames[$current][$j];  
                                }
				elsif(($summit >= $reftss[$current][$j])&&($summit < $refatg[$current][$j]))       
                                {
					$select[2] = 1;  
                                        $refseqnameofselect[2] = $refnames[$current][$j];
                                }
				elsif(($summit >= $refatg[$current][$j])&&($summit < $refgta[$current][$j]))      
                                {
					$select[3] = 1;  
                                        $refseqnameofselect[3] = $refnames[$current][$j];
                                }
                                elsif(($summit >= $refgta[$current][$j])&&($summit < $reftes[$current][$j])) 
                                {
                   			$select[4] = 1;  
                                        $refseqnameofselect[4] = $refnames[$current][$j];
                                }
                                elsif(($summit >= $reftes[$current][$j])&&($summit < ($reftes[$current][$j] + $downproxilimit ))) 
                                {
					$select[5] = 1;  
                                        $refseqnameofselect[5] = $refnames[$current][$j];
                                }
                                elsif(($summit >= ($reftes[$current][$j] + $downproxilimit ))&&($summit <= ($reftes[$current][$j] + $distallimit))) 
                                {
					$select[6] = 1;  
                                        $refseqnameofselect[6] = $refnames[$current][$j];
                                }
			}
##### 
##### if "-" strand
#####
                        elsif($refstrand[$current][$j] eq "-")
                        {
                                if(($summit >= ($reftes[$current][$j]+$upproxilimit ))&&($summit <= ($reftes[$current][$j] + $distallimit )))
                                {
					$select[0] = 1;  
                                        $refseqnameofselect[0] = $refnames[$current][$j];
                                }
                                elsif(($summit >= $reftes[$current][$j])&&($summit < ($reftes[$current][$j] + $upproxilimit )))
                                {
					$select[1] = 1;  
                                        $refseqnameofselect[1] = $refnames[$current][$j];
                                }
                                elsif(($summit >= $refgta[$current][$j])&&($summit < $reftes[$current][$j]))       
                                {
					$select[2] = 1;  
                                        $refseqnameofselect[2] = $refnames[$current][$j];
                                }
                                elsif(($summit >= $refatg[$current][$j])&&($summit < $refgta[$current][$j]))      
                                {
					$select[3] = 1;  
                                        $refseqnameofselect[3] = $refnames[$current][$j];
                                }
                                elsif(($summit >= $reftss[$current][$j])&&($summit < $refatg[$current][$j])) 
                                {
					$select[4] = 1;  
                                        $refseqnameofselect[4] = $refnames[$current][$j];
                                }
                                elsif(($summit >= ($reftss[$current][$j] - $downproxilimit ))&&($summit < $reftss[$current][$j])) 
                                {
					$select[5] = 1;  
                                        $refseqnameofselect[5] = $refnames[$current][$j];
                                }
                                elsif(($summit >= ($reftss[$current][$j] - $distallimit ))&&($summit < ($reftes[$current][$j] - $downproxilimit ))) 
                                {
					$select[6] = 1;  
                                        $refseqnameofselect[6] = $refnames[$current][$j];
                                }
                        }
			else
			{
				print "strand error.\n";
				exit(0);
			}
		}
	}
	if ($hit == 0)
	{
		$anncount[7] ++;
	}
	else
	{
		$selectsum = 0;
                for($j=0;$j<7;$j++)                
                {
			$selectsum = $selectsum + $select[$j];
		}
		if($selectsum == 0)
                {
			print "selectsum error.\n";
		}
		elsif($selectsum ==1)
		{
			for($j=0;$j<7;$j++)
			{
				if ($select[$j] == 1)
		        	{
                			$anncount[$j] ++;
        			}
			}
		}
		else
                {
			if (($select[0] == 1)&&($select[6] == 0))
			{
				for($j=1;$j<6;$j++)            
	                        {        
        	                        if ($select[$j] == 1)
                	                {       
                        	                $anncount[$j] = $anncount[$j] + 1/($selectsum -1);
                                	}       
 	                        }
			}
			elsif (($select[0] == 0)&&($select[6] == 1))
                        {
                                for($j=1;$j<6;$j++)
                                {
                                        if ($select[$j] == 1)
                                        {
                                                $anncount[$j] = $anncount[$j] + 1/($selectsum -1);
                                        }      
                                }       
                        }
			elsif (($select[0] == 1)&&($select[6] == 1)&&($selectsum >2))
                        {
                                for($j=1;$j<6;$j++)
                                {
                                        if ($select[$j] == 1)
                                        {
                                                $anncount[$j] = $anncount[$j] + 1/($selectsum -2);
                                        }
                                }
                        }         
			else
			{
	                        for($j=0;$j<7;$j++)        
        	                {
                	                if ($select[$j] == 1)
                        	        {        
                                	        $anncount[$j] = $anncount[$j] + 1/$selectsum;
                                	}
				}
                        }
                }     
	}
	$occupied = "a";
	for($j=0;$j<7;$j++)
        {                
#      	        print OUTPUTFILE "$select[$j] ";
		if($select[$j] == 1)
		{
			$occupied = "b";
			print OUTPUTFILE "$names[$j]-$refseqnameofselect[$j]\t";
#			print "$names[$j] $refseqnameofselect[$j]\t";
#			exit(0);
		}
        }
	if($occupied eq "a") 
        {
                print OUTPUTFILE "intergenic\n"; 
        }
	else
	{
	        print OUTPUTFILE "\n";
	}
	if($select[1] == 1)
        {
                        print OUTFILE1 "$refseqnameofselect[1]\n";
#                       print "$refseqnameofselect[1]\n";
#                       exit(0);
        }
        if($select[3] == 1)
        {
                        print OUTFILE2 "$refseqnameofselect[3]\n";           
#                       print "$refseqnameofselect[3]\n";           
#                       exit(0);
        }    
}
close(INPUTFILE);
close(OUTPUTFILE);       
close(OUTFILE1);
close(OUTFILE2);
print "There are $totalcount peaks.\n";
print "There are $anncount[0] peaks that are 5' dismal.\n";
print "There are $anncount[1] peaks that are 5' proximal.\n";            
print "There are $anncount[2] peaks that are 5' UTR.\n";            
print "There are $anncount[3] peaks that are coding.\n";            
print "There are $anncount[4] peaks that are 3' UTR.\n";            
print "There are $anncount[5] peaks that are 3' proximal.\n";            
print "There are $anncount[6] peaks that are 3' dismal.\n";            
print "There are $anncount[7] peaks that are intergenic.\n";            








