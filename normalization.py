#normalization script
import sys
import pandas as pd
import numpy as np
df = pd.read_csv(sys.argv[1], sep='\t', index_col = 0, header=0)

#normalize values and write
def normalize_values(lista):
    #gene = lista[0]
    newlist=[]
    data= lista
    datafl=[]
    datafl = data[np.logical_not(np.isnan(data))] # remove NA's to get floats and min/max
    if len(datafl) > 0: #check list is not empty
        datafl=np.array(datafl, dtype=np.float32) #convert to float in numpy
        mindata= np.amin(datafl)
        maxdata= np.amax(datafl)
        print(mindata, maxdata)
        dem= float(maxdata)-float(mindata)
        if dem != float(0):
            for i in data: #have to go back to original data because have removed NAs which can be placeholders
                    if i != np.nan:
                            norm= (float(i)-float(mindata))/(dem)
                            newlist.append(str(norm))
                    else:
                        newlist.append(np.nan)
            #print(newlist)
            return(newlist)
        else: #replace with NAs if denominator is = 0
            dem=float(0.00001)
            for i in data: #have to go back to original data because have removed NAs which can be placeholders
                    if i != np.nan:
                            norm= (float(i)-float(mindata))/(dem)
                            newlist.append(str(norm))
                    else:
                        newlist.append(np.nan)
            #print(newlist)
            return(newlist)
    else: #replace with NAs if all NAs
        for i in data:
            newlist.append(np.nan)
        #print(newlist)
        return(newlist)

        

#get NAs in np format 
df = df.replace("?",np.nan) 
df = df.replace("NA",np.nan)
df = df.replace("",np.nan)   

col_list= list(df.columns.values)

#loop through file
rows_to_norm = df.index.values.tolist()
cols_to_norm = list(df)

#df[rows_to_norm] = df[rows_to_norm, ].apply(lambda x: (x - x.min()) / (x.max() - x.min()))
#print(df.loc[rows_to_norm,:])

newdf= df.loc[rows_to_norm,:].apply(normalize_values, axis=1, result_type= 'expand')

print(newdf)
newdf.columns = col_list
print(newdf)

newdf.to_csv(path_or_buf=str(sys.argv[1])+'_norm.txt', sep="\t", header=True)       