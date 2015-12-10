#!/usr/bin/perl
# Program to score fasta sequence files using the method of Berg and von Hippel (1987)
# And promoter scoring method of Rhodius & Mutalik (2010) and Rhodius et al (2012)
#
# Specific for sigma 70 promoters
# Based on all s70 promoters -35, spacer, -10, disc, start for all promoters
# Has specific PWMs and spacer penalties
# -35 PWM, spacer penalty, -10 PWM, disc penalty, start PWM
# Note no combined S+D penalty as these spacer lengths are NOT correlated

# Default scores just top strand
# Default outputs highest scoring promoter if there are overlaps that share -10 or -35 motif
# Threshold for -10/-35 motifs greater than z-score cutoff threshold (user input)
# Option to score upstream (UP-element) sequences (user input)
# Option to score both strands (user input)
# Option to just give highest scoring promoter per sequence (user input)
# Option to screen promoters by their total zscore (user input)
# Option to output all promoter hits (even if overlapping) (user input)

# Search format is <PWM1><spacer1><PWM2><spacer2><PWM3>

# THIS PROGRAM REQUIRES INPUT FILES IN THE <DATAFILES> DIRECTORY
# s70_model_key.txt descriptor file of promoter model file names of format
# <Model 1 name><PWM1><PWM2><PWM3><spacer1><spacer2>
#
# options.txt lists options for program
# s70_promoter_stats.txt lists average and standard deviation scores for s70 promoter library
#
# Also requires following data files
# Tab delimited Position Weight Matrix files (*PWM.txt) 1-3
# Tab delimited spacer and discriminator length frequency files 1-2 (length<tab>frequency)
#
# Promoter sequence input file (FASTA format or single line sequence)

# THIS PROGRAM CREATES OUTPUT FILES 
# Predictions of matching sites

use Getopt::Long; #Use Getopt module for commandline options

use Cwd;
$dir=cwd(); 	#Gets current working directory

# Input commandline options
$seqname=(); # Input sequence filename
$outputname=(); # Output filename
$zscore=(); # Z-score cutoff for PWM motif search cutof threshold. Must be a number
$pzscore=(); # Promoter Z-score cutoff threshold for total promoter z-score. Must be a number
$both=(); # Flag for both strand search (Default is top strand only)
$up=(); # Flag to score upstream sequence
$high=(); # Flag to output just highest scoring promoter for each sequence
$all=(); # Flag to output all promoter hits including shared -10/-35 motifs
GetOptions('f=s'=>\$seqname,'o=s'=>\$outputname,'z=f'=>\$zscore,'p=f'=>\$pzscore,
	'b'=>\$both,'u'=>\$up,'h'=>\$high,'a'=>\$all);

unless (($seqname) and ($outputname) and ($zscore)) { #Tests for required inputs
	chdir("Datafiles")|| die "Can't chdir to Datafiles: $!";
	open(OPTIONS,"<options.txt") || die "Can't open file $options.txt: $!";
	print <OPTIONS>,"\n";
	chdir"$dir";
	exit;
} 

# Obtain promoter model stats file
chdir("Datafiles") || die "Can't chdir to Datafiles: $!";
open(STATS, "<s70_promoter_stats.txt") || die "Can't open file s70_promoter_stats.txt: $!";
$i=0;
while ($statsin=<STATS>) {
	chomp $statsin;
	$stats[$i]=[split'\t',$statsin];
	++$i;
}
close (STATS);

# Obtain promoter model key descriptor file
open(MODEL, "<s70_model_key.txt") || die "Can't open file s70_model_key.txt: $!";
$modelin=<MODEL>;
chomp $modelin;
@model=split'\t',$modelin;
close (MODEL);

# Opening and storing Position Weight Matrices
for ($b=0;$b<4;++$b) { # Loop for each PWM
	open(PWM, "<$model[$b+1]") || die "Can't open file $model[$b+1]: $!";
	$i=0;
	while ($pwmin=<PWM>) {
		chomp $pwmin;
		$pwm[$b][$i]=[split'\t',$pwmin]; #Stores PWMs into pos 0-3
		if ($i=~/1/) { # To obtain width of array for pwm length
			@pwmwidth=split'\t',$pwmin; #Stores width of PWM array into pos 0-3
			$pwmlength[$b]=@pwmwidth-1; # length of pwm site
		}
		++$i;
	}
	close (PWM);
}

# Calculating motif cutoff thresholds and pwm lengths
@threshold=();
for ($a=0;$a<4;++$a) { # For each PWM
	$threshold[$a]=$pwm[$a][5][1]+($zscore*$pwm[$a][6][1]); # Threshold for PWM
}

# Opening and storing spacer lengths
for ($b=0;$b<2;++$b) { # Loop for each spacer
open(SPACER, "<$model[$b+5]") || die "Can't open file $model[$b+5]: $!";
	$i=0;
	while ($spacerin=<SPACER>) {
		chomp $spacerin;
		$spacer[$b][$i]=[split'\t',$spacerin]; #Stores spacers into pos 0-1
			++$i;
		}
	$range[$b]=$i; #Stores range of spacer lengths
	close (SPACER);
}

# Find highest spacer freq
for ($a=0;$a<2;++$a) { # For each spacer model
	$i=0;
	$highfreq[$a]=0;
	while ($spacer[$a][$i][0]) {
		if ($spacer[$a][$i][1]>$highfreq[$a]) {
			$highfreq[$a]=$spacer[$a][$i][1];
		}
		++$i;
	}
}	
chdir"$dir"; #returns to parent directory

# Opening and storing Sequence file
open(PROMOTER, "<$seqname") || die "Can't open file $seqname: $!";
$i=0;
@promoter=();
while ($promoterin=<PROMOTER>) {
	chomp $promoterin;
	@promoter[$i]=$promoterin;
	++$i;
}
close (PROMOTER);

# Calculating scores for each sequence
@scores=(); # Initiating scores array
@prediction=(); # Initiating output predictions array
$fasta=(); #Sets sequence input flag as single line input
SEQ: for ($i=0;$i<@promoter;++$i) {
	if ($promoter[$i]=~/\>/) { #Screening Fasta header lines containing ">"
		$fasta=1; # Flags sequence as FASTA format
		next SEQ;
	}
	$promoter[$i]=~tr/gatc/GATC/; # Converting small letters into large
	$sequence=(); # Resetting $sequence
	$sequence=$promoter[$i]; # Assigning sequence

	# Grabbing sequence name
	if ($fasta) { #If FASTA format
		$name=$promoter[$i-1]; # FASTA sequence name
	} else { 
		$name=">".$seqname; # Assigning filename as sequence name
	}
	SCORES(); # Scores calculation subroutine

	if ($both) { #If scoring both strands
		REVERSE(); #Generate reverse complementary sequence
		# Grabbing sequence name and tag with "_r" for reverse
		if ($fasta) { #If FASTA format
			$name=$promoter[$i-1]."_r"; # FASTA sequence name
		} else { 
			$name=">".$seqname."_r"; # Assigning filename as sequence name
		}
		SCORES(); # Scores calculation subroutine
	}
}

RESULT(); # Print out @prediction
print"complete\n";

# Subroutines from here

sub SCORES { # For each sequence, calculating scores for each search window
	$highscore=-99; #Resetting highscore tally
 	$seqpos=@prediction; # Size of @prediction for sequence
	@position=(); #Resetting position
	$sequence=~tr/GATC/1234/; # Converts sequence into positional numbers for @pwm
	$seqlength=length($sequence);

	for ($s=0;$s<$range[0];++$s) { # For each spacer length
		for ($t=0;$t<$range[1];++$t) { # For each discriminator length
			$search=($seqlength-$pwmlength[0]-$pwmlength[1]-$pwmlength[2]-$spacer[0][$s][0]-$spacer[1][$t][0])+1; # Calculating no of search windows
			if ($up) { # If scoring upstream sequence
				$start=18; #Offset start search position for -35 to position 18 to enable upstream search
			}
			WINDOW: for ($a=$start;$a<$search;++$a) { # Looping for each search window
				@score=(); #Resets subsite scores
				# 35PWM (PWM1)
				$window[0]=substr($sequence,$a,$pwmlength[0]); # Extracts from sequence the residue numbers for each search window
				$x=0; #Initiating counter
				while ($window[0]) {
					$residue=chop($window[0]);
					# Scores residue in pwmmatrix from 3' end and sums into @score
					$score[0]=$score[0]+$pwm[0][$residue][$pwmlength[0]-$x];
					++$x; # Increment counter
				}
				unless ($score[0]>$threshold[0]) { # If not >threshold then next search window
					next WINDOW;
				}
				# 10PWM (PWM2)
				$b=$a+$pwmlength[0]+$spacer[0][$s][0]; # 2nd PWM search position
				$window[1]=substr($sequence,$b,$pwmlength[1]); # Extracts from sequence the residue numbers for each search window
				$x=0; #Initiating counter
				while ($window[1]) {
					$residue=chop($window[1]);
					# Scores residue in pwmmatrix from 3' end and sums into @score
					$score[1]=$score[1]+$pwm[1][$residue][$pwmlength[1]-$x];
					++$x; # Increment counter
				}
				unless ($score[1]>$threshold[1]) { # If not >threshold then next search window
					next WINDOW;
				}
				# startPWM (PWM3)
				$c=$b+$pwmlength[1]+$spacer[1][$t][0]; # 3rd PWM search position
				$window[2]=substr($sequence,$c,$pwmlength[2]); # Extracts from sequence the residue numbers for each search window
				$x=0; #Initiating counter
				while ($window[2]) {
					$residue=chop($window[2]);
					# Scores residue in pwmmatrix from 3' end and sums into @score
					$score[2]=$score[2]+$pwm[2][$residue][$pwmlength[2]-$x];
					++$x; # Increment counter
				}
				unless ($score[2]>$threshold[2]) { # If not >threshold then next search window
					next WINDOW;
				}
				# spPWM (PWM4)
				$d=$b-$pwmlength[3]; # 4th PWM search position (no threshold for this PWM)
				$window[3]=substr($sequence,$d,$pwmlength[3]); # Extracts from sequence the residue numbers for each search window
				$x=0; #Initiating counter
				while ($window[3]) {
					$residue=chop($window[3]);
					# Scores residue in pwmmatrix from 3' end and sums into @score
					$score[3]=$score[3]+$pwm[3][$residue][$pwmlength[3]-$x];
					++$x; # Increment counter
				}
				# Spacer penalties
				$score[4]=log((($spacer[0][$s][1])+($highfreq[0]*0.005))/(($highfreq[0])+($highfreq[0]*0.005))); #Spacer penalty
				$score[5]=log((($spacer[1][$t][1])+($highfreq[1]*0.005))/(($highfreq[1])+($highfreq[1]*0.005))); #Discriminator penalty
				# Scoring UP-element (distal and proximal subsites from -57 to -37)
				if ($up) { #If scoring upstream sequence (-57 to -37)
					$e=$a-21; # UP-element scoring position (-57)
					if ($e<0) { #If negative pos then reset to 0
						$e=0;
					}
					if ($a<21) { # If upstream sequence shorter than UP-element
						$uplen=$a; #Set to shorter seq
					} else {
						$uplen=21; # Else grab full seq
					}
					$window[4]=substr($sequence,$e,$uplen); #Extracts UP-element sequence
					$window[4]=~tr/1234/GATC/; # Converting numbers back to sequence
					for ($j=0;$j<22-2;++$j) {
						$tri=substr($window[4],$j,3);
						if ($tri=~"AAA") {
							++$score[6]; #Tallying "AAA"
						} elsif ($tri=~"TTT") {
							++$score[6]; #Tallying "TTT"
						}
					}
				}			

				# Calculating total score and total zscore
				$totalscore=$score[0]+$score[1]+$score[2]+$score[3]+$score[4]+$score[5]+$score[6];
				if ($up) { #Calculating zscore if scoring upstream sequences
					$totalzscore=($totalscore-$stats[2][1])/$stats[2][2];
					} else { #Calculating zscore if just scoring core promoter sequences
					$totalzscore=($totalscore-$stats[1][1])/$stats[1][2];
				}

				# If filtering by promoter zscore
				if ($pzscore) { 
					if ($totalzscore<$pzscore) {
						next WINDOW; #If less than cutoff then discard
					}
				}
 
				$predsize=@prediction; # Size of @prediction
				if ($high) { #If only outputting highest score for each sequence
					if ($totalscore>$highscore) { #Storing prediction if highest score
						$highscore=$totalscore; #Resetting highest score
						$m=$seqpos; #Resets $m to end of @prediction for sequence
						PREDICTION(); # Prediction subroutine
					}
				} elsif ($all) { # Output all predictions
					$m=$predsize; #Resets $m to end of @prediction
					PREDICTION(); # Prediction subroutine
				} else { # Removes lower scoring predictions that share -10 or -35 hits
					OVERLAP(); # Overlap subroutine to remove overlaps
				}
			}
		}
	}
}

sub REVERSE { # Reverse complementary sequence subroutine
	$sequence=$promoter[$i]; #Reassigning sequence
	$sequence=~tr/GATC/CTAG/; # Generate complementary nts
	@revcomp=(); #Initiating @reverse comp array
	$a=0; #Initiating counter
	while ($sequence) { #Removing sequence 1 nt at a time from 3'end
		$revcomp[$a]=chop($sequence);
		++$a; #Increment counter
	}
	$sequence=join"",@revcomp; #Assigning sequence reverse order
}

sub OVERLAP { #Discard current hit if overlaps with previous hits and has lower totalscore
	$predsize=@prediction; # Size of @prediction
	for ($m=0;$m<$predsize;++$m) { #For each prediction
		if ($prediction[$m][0]=~/$name/) { #If seq name matches
			if ($prediction[$m][3] eq ($a+1)) { #If PWM1 pos matches
				if ($totalscore>$prediction[$m][17]) {
					PREDICTION(); #Overwrite with higher scoring overlap
				} else {
					next WINDOW; #Lower score so discard
				}
			}
			if ($prediction[$m][8] eq ($b+1)) { #If PWM2 pos matches
				if ($totalscore>$prediction[$m][17]) {
					PREDICTION(); #Overwrite with higher scoring overlap
				} else {
					next WINDOW; #Lower score so discard
				}
			}
		}
	} 
	#No overlap so new prediction
	$m=$predsize; #Resets $m to end of @prediction
	PREDICTION(); #Adds as new prediction
}

sub PREDICTION { # Adding predicted sites to @prediction
	$prediction[$m][0]=$name; # Assigns sequence name
	$prediction[$m][1]=$window[4]; # UP sequence
	$prediction[$m][2]=$score[6]; # UP-element score
	$prediction[$m][3]=$a+1; # PWM1 position (35motif)
	$site=substr($sequence,$a,$pwmlength[0]);
	$site=~tr/1234/GATC/; # Converting numbers back to sequence
	$prediction[$m][4]=$site; # PWM1 sequence
	$prediction[$m][5]=$score[0]; # PWM1 score
	$site=substr($sequence,$d,$pwmlength[3]);
	$site=~tr/1234/GATC/; # Converting numbers back to sequence
	$prediction[$m][6]=$site; # PWM4 sequence (spacer motif)
	$prediction[$m][7]=$score[3]; # PWM4 score
	$prediction[$m][8]=$b+1; # PWM2 position (10motif)
	$site=substr($sequence,$b,$pwmlength[1]);
	$site=~tr/1234/GATC/; # Converting numbers back to sequence
	$prediction[$m][9]=$site; # PWM2 sequence
	$prediction[$m][10]=$score[1]; # PWM2 score
	$prediction[$m][11]=$spacer[0][$s][0]; # Spacer 1 length
	$prediction[$m][12]=$score[4]; # Spacer 1 penalty
	$site=substr($sequence,$c,$pwmlength[2]);
	$site=~tr/1234/GATC/; # Converting numbers back to sequence
	$prediction[$m][13]=$site; # PWM3 sequence (start motif)
	$prediction[$m][14]=$score[2]; # PWM3 score
	$prediction[$m][15]=$spacer[1][$t][0]; # Spacer 2 length
	$prediction[$m][16]=$score[5]; # Spacer 2 penalty
	$prediction[$m][17]=$totalscore; # Total site score
	$prediction[$m][18]=$totalzscore; # Total site zscore
	next WINDOW;
}

# Writing @prediction files in results directory
sub RESULT {
	open(OUTPUT, ">$outputname") || die "Can't create file $outputname: $!";
	print OUTPUT"# Input is:\t$seqname\n";
	print OUTPUT"# Output file is:\t$outputname\n";
	print OUTPUT"# Motif z-score threshold:\t$zscore\n";
	print OUTPUT"# Promoter z-score threshold:\t$pzscore\n";
	print OUTPUT"# Search both strands:\t$both\n";
	print OUTPUT"# Score UP-element:\t$up\n";
	print OUTPUT"# Highest score only:\t$high\n";
	print OUTPUT"# Outputs all predictions:\t$all\n";
	print OUTPUT"# Sequence\tUP sequence\tUP-score\t-35 pos\t-35 seq\t-35 score\tspacer seq\tspacer score\t-10 pos\t-10 seq\t-10 score\tspacer len\tspacer pen\tstart seq\tstart score\tdisc len\tdisc pen\tTotal score\tTotal zscore\n";
	for ($a=0;$a<@prediction;++$a) {
		for ($b=0;$b<@{@prediction[$a]};++$b) {
			print OUTPUT"$prediction[$a][$b]\t";
		}
		print OUTPUT"\n";
	}
	close (OUTPUT);
}