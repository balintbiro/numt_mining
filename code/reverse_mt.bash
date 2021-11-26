awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' mt.fa >> mt_oneliner.fa
grep '>' mt_oneliner.fa >> r_mt.fa
sed -n 2p mt_oneliner.fa | rev >> r_mt.fa