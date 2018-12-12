#normalization script
import sys
inp= open(sys.argv[1], 'r')
output= open(str(sys.argv[1])+'_norm.txt', 'w')

#normalize values and write
def normalize_values(lista):
    gene = lista[0]
    newlist=[]
    #clust= lista[1]
    data= lista[1:]
    datafl=[]
    for l in data:
        try: 
            float(l)
            datafl.append(l)
        except(ValueError):
            datafl.append('NA')
    mindata= min(datafl)
    maxdata= max(datafl)
    if maxdata == 'NA':
        pass
    elif mindata == 'NA':
        pass
    else:
        dem= float(maxdata)-float(mindata)
        if dem != float(0):
            for i in datafl:
                if i != 'NA':
                    norm= (float(i)-float(mindata))/(dem)
                    newlist.append(str(norm))
                else:
                    newlist.append('NA')
            newstr= '\t'.join(newlist)
            output.write('%s\t%s\n' % (gene, newstr))
        else:
            for i in datafl:
                newlist.append(str(0))
            newstr= '\t'.join(newlist)
            output.write('%s\t%s\n' % (gene, newstr))

#write header
header =inp.readline().strip().split('\t')
headerstr= '\t'.join(header)
output.write(headerstr+'\n')
        
#loop through file
for line in inp:
    L= line.strip().split('\t')
    normalize_values(L)
    
inp.close()
output.close()        