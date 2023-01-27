import pandas as pd
import numpy as np
from collections import Counter
from itertools import chain
from sklearn.metrics import multilabel_confusion_matrix
from sklearn.metrics import classification_report
from sklearn.metrics import accuracy_score

a=pd.read_csv("/Users/islekbro/Desktop/genomize_makale/3_statistics/statistics_inputs/comparison.txt",sep="\t")

a["seq-pathogenicity"]=a["seq-pathogenicity"].str.replace("-","")
a["seq-pathogenicity"]=a["seq-pathogenicity"].str.replace("+","")

a=a.rename({'seq-pathogenicity' : 'seq'}, axis=1)
a=a.rename({'clingen-pathogenicity' : 'clingen'}, axis=1)
a=a.rename({'franklin_pathogenicity' : 'franklin'}, axis=1)


a["seq"]=a["seq"].str.replace(" ","")
a["clingen"]=a["clingen"].str.replace(" ","")
a["franklin"]=a["franklin"].str.replace(" ","")

#I would like to create a third column (column C) which returns the items that exists in column A but does not exist in column B

a = a.replace(r'^\s*$', np.nan, regex=True)
#a.fillna("blank", inplace = True)


#create new df for SEQ by unique analysis comparing to Clingen
df = a[["clingen","seq","franklin"]].copy()
sorter = ["P","LP","VUS","LB","B"]

C=df["clingen"].str.split(pat = ",", expand=True).apply(lambda x : x.value_counts(), axis = 1).fillna(0).astype(int)
C=C.reindex(sorter, axis='columns')
clingen=C.iloc[:].values

S=df["seq"].str.split(pat = ",", expand=True).apply(lambda x : x.value_counts(), axis = 1).fillna(0).astype(int)
S=S.reindex(sorter, axis='columns')
seq=S.iloc[:].values

F=df["franklin"].str.split(pat = ",", expand=True).apply(lambda x : x.value_counts(), axis = 1).fillna(0).astype(int)
F=F.reindex(sorter, axis='columns')
franklin=F.iloc[:].values


seq_report=classification_report(clingen, seq, output_dict=True, zero_division=0, target_names=sorter)
S_df = pd.DataFrame(seq_report).transpose().reset_index()
seq_df = S_df.iloc[:-4].drop(columns="support")
seq_df["platform"]="seq"


franklin_report=classification_report(clingen, franklin, output_dict=True, zero_division=0, target_names=sorter)
F_df = pd.DataFrame(franklin_report).transpose().reset_index()
franklin_df = F_df.iloc[:-4].drop(columns="support")
franklin_df["platform"]="franklin"

final_df=pd.concat([seq_df,franklin_df])


final_df=final_df[(final_df != 0).all(axis=1)]
#final_df=final_df.rename({'index' : 'Evidence Codes'}, axis=1)
final_df.to_csv("/Users/islekbro/Desktop/genomize_makale/3_statistics/statistics_outputs/statistics_variant_pathogenicity_out.csv",sep="\t", index=False)
