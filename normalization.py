#normalization script
import sys
import pandas as pd
import numpy as np
df = pd.read_csv(sys.argv[1], sep='\t', index_col = 0, header=None) #dataframe to normalize
direct= str(sys.argv[2]) #by row or col
#normalize values and write
def normalize_values(lista):
    #gene = lista[0]
    newlist=[]
    data= lista
    #print(data)
    datafl=[]
    datafl = data[np.logical_not(pd.isnull(data))] # remove NA's to get floats and min/max
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

#column names as list
#cols_list= list(df.columns.values)
cols_list= list(df[:1])
#get columns to normalize
cols_to_norm= list(df.columns[1:].values)

#get rows to normalize
rows_to_norm = df.index[1:].values.tolist()

#rows_to_norm = df.loc[1:,].tolist()
#cols_to_norm = list(df)
#print(rows_to_norm)
#df[rows_to_norm] = df[rows_to_norm, ].apply(lambda x: (x - x.min()) / (x.max() - x.min()))
#print(df.loc[rows_to_norm,:])

if direct == "row":
    newdf= df.loc[rows_to_norm,:].apply(normalize_values, axis=1, result_type= 'expand')
elif direct == "col":
    newdf= df.loc[rows_to_norm,:].apply(normalize_values, axis=0, result_type= 'expand')
else:
    print("Need second argument: row = normalize rows, col = normalize by column")

newdf.columns = cols_list
print(newdf)
print(cols_to_norm, cols_list)
newdf.to_csv(path_or_buf=str(sys.argv[1])+'_norm.txt', sep="\t", header=True)       