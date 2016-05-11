#Written by Genivaldo Silva

#Change to output of og.pl
OG="14051.OGs.txt"

f=open(OG)




h={};K=["127371_RAST.faa","127372_RAST.faa","186156_RAST.faa"];IDs={}

for line in f:
    line=line.replace("\n","").split("\t")
    myOG=line[2]
    fasta=line[1].split("6666666.")[-1]
    function=line[0].split("___")[-1]
    sequenceID=line[0].split("___")[0]

    if fasta not in IDs:
        IDs[fasta]={}
    IDs[fasta][sequenceID]=0
    
    if myOG not in h:
        h[myOG]={x:[] for x in K}

    h[myOG][fasta]+=[function]
f.close()

o=open(OG+"___.xls","w+")

o.write("OG\t"+"\t".join(K)+"\tSame Function?\n")
for i in h:
    temp="";c=[]
    for k in [h[i][x] for x in K]:
        k=";".join(k)
        temp=temp+"\t"+k
        c.append(k)


    if len(set(c))==1:
        c="Yes"
    else:
        c="No"
    o.write(i+temp+"\t"+c+"\n")

o.close()
       