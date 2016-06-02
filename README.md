# G-PAnalysis
Various scripts used to analyze genetic and phenotypic growth information for the Genotype-Phenotype Project at Edwards Bioinformatics Lab

Protocol for Blastp Analysis
All code available at https://github.com/brobinit/G-PAnalysis
Acquiring formatted genome sequences
1.	Log into RAST and access your Jobs Overview page. 
2.	View job details for the desired genome. 
3.	From the drowdown menu, select “Spreadsheet (tab-separated text format)” and download. 
Reformat to fasta
4.	Rename extension from .txt to .xls and move to same folder as create_fasta.py.
5.	In terminal, make the folder containing the genome(s).xls and create_fasta.py the active directory.
6.	Create directory “fastas.”
7.	Run:
python create_fasta.py
This will open each .xls file, reformat it to .faa, and save it to the fastas/ directory. 
8.	Copy fasta sequences, format_db.py, run_blastp.py from local server to a new directory on the Edwards Anthill remote server.
On Edwards Anthill remote server
9.	In the newly created directory, create “fastas” directory and copy the created .faa sequences to it.
10.	Create the following directories: db, OG.
11.	Run
python format_db.py
This will use makeblastdb (on the remote server) to create a database in db/ for each .faa file.
12.	Run
python run_blatp.py
This runs an all-versus-all blastp analysis of all the created databases from .faa files. 
13.	Run
perl og.pl fastas/ OG/
This will match orthologous groups from the all-versus-all blastp analysis. 
On local server
14.	Copy OGs.txt file from remote server to local server in the same directory as parse_og.py.
15.	Run
python parse_og.py
This will reformat the orthologous groups identified via og.pl into a more user-friendly format. Also looks at the name of nucleotide sequence function of same orthologous group across analyzed genomes and see if they all match, identified as a yes or no in “Same function?” column. 

