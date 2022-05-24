from ftplib import FTP

ftp = FTP('ftp.ncbi.nlm.nih.gov')
ftp.login()
ftp.cwd('/genomes/refseq/vertebrate_mammalian/')
ftp.nlst()