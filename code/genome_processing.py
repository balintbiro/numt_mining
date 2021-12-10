#import the required modules
import os
import filecmp
from subprocess import call

call("grep -E '>' genome.fa >> headers.txt", shell=True)#get all the headers from the genome
call("grep -E '>*REF' headers.txt >> ref_headers.txt", shell=True)#get the reference headers
if filecmp.cmp('headers.txt','ref_headers.txt')==False:#that means that there are alternative headers
    os.rename('genome.fa','alt_genome.fa')#reanem the original genome.fa file to alt_genome.fa which contains the alternative sequences
    call("awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < alt_genome.fa > processed_genome.fa", shell=True) #transform alt_genome to one liner form
    call("grep -f ref_headers.txt -A1 processed_genome.fa >> genome.fa", shell=True)#get the reference sequences based on the reference headers in the ref_headers.txt file
    call("sed -i '/^-/d' genome.fa", shell=True)##delete the line starting with '-' sign which is is introduced by grep command