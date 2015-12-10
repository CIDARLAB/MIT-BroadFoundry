http://allaboutbioinfo.blogspot.com/2012/04/estimating-paired-end-read-insert.html


Usage:

getinsertsize.py [ SAM file | -]

or

samtools view [ BAM file ] | getinsertsize.py - 


Detailed Usage:


usage: getinsertsize.py [-h] [--span-distribution-file SPAN_DISTRIBUTION_FILE]
                        [--read-distribution-file READ_DISTRIBUTION_FILE]
                        SAMFILE

Automatically estimate the insert size of the paired-end reads for a given
SAM/BAM file.

positional arguments:
  SAMFILE               Input SAM file (use - from standard input)

optional arguments:
  -h, --help            show this help message and exit
  --span-distribution-file SPAN_DISTRIBUTION_FILE, -s SPAN_DISTRIBUTION_FILE
                        Write the distribution of the paired-end read span
                        into a text file with name SPAN_DISTRIBUTION_FILE.
                        This text file is tab-delimited, each line containing
                        two numbers: the span and the number of such paired-
                        end reads.
  --read-distribution-file READ_DISTRIBUTION_FILE, -r READ_DISTRIBUTION_FILE
                        Write the distribution of the paired-end read length
                        into a text file with name READ_DISTRIBUTION_FILE.
                        This text file is tab-delimited, each line containing
                        two numbers: the read length and the number of such
                        paired-end reads.


Sample output:

Read length: mean 90.6697303194, STD=15.9446036414
Possible read length and their counts:
{108: 43070882, 76: 50882326}
Read span: mean 165.217445903, STD=32.8914834802


Note: If the SAM/BAM file size is too large, it is accurate enough to estimate based on a few reads (like 1 millioin). In this case, you can run the script as follows:

head -n 1000000 [ SAM file ] |  getinsertsize.py -

or

samtools view [ BAM file ] | head -n 1000000 | getinsertsize.py -

Note: According to the SAM definition, the read span "equals the number of bases from the leftmost mapped base to the rightmost mapped base". This span is the distance between two reads in a paired-end read PLUS 2 times read length. Read span is different from the "mate-inner-distance" in Tophat (-r option), which measures only the distance between two reads in a paired-end read.