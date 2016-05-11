# Written by Genivaldo Silva

import os

for myfile in [i for i in os.listdir(os.curdir) if ".xls" in i]:
    print myfile

    #myfile="Prevotella.xls"
    f=open(myfile)
    f.readline()
    o=open("fastas/"+myfile.replace(".xls","_RAST")+".faa","w+")
    for i in f:
        i=i.replace("\n","").split("\t")

        sequence=i[-1]
        ID=i[3]
        function=i[7].replace(" ","_")
        if len(sequence)>0:
            o.write(">"+ID+"___"+function+"\n"+sequence+"\n")

    f.close()
    o.close()
