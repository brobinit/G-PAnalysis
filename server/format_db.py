#Written by Genivaldo Silva

import os

for myfile in [i for i in os.listdir(os.curdir) if ".faa" in i]:
    os.system("makeblastdb -in "+myfile+" -dbtype prot  -out db/"+myfile)
    
