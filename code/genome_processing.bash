grep -E '>' genome.fa >> headers.txt #get all the headers from the genome
grep -E '>*REF' headers.txt >> ref_headers.txt #get the reference headers

#if the size of headers != size of ref_headers, get the genome without alternative contigs!
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < genome.fa > processed_genome.fa #transform genome to one liner form
grep -f ref_headers.txt -A1 processed_genome.fa >> new_genome.fa #get the reference sequences
sed -i '/^-/d' new_genome.fa #delete the line starting with '-' sign which is not supported by LASTAL