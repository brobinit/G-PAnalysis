# Written by Genivaldo Silva

import os
allfiles=[i for i in os.listdir(os.curdir) if ".faa" in i]
for i in allfiles:
    for j in allfiles:
        print i,j
        os.system("blastp -db db/"+j+" -query "+i+" -out OG/"+i+"-vs-"+j+" -outfmt 6 -evalue 0.00001 -num_threads 24")
