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

a=a.rename({'seq-evidence_codes' : 'seq'}, axis=1)
a=a.rename({'clingen-evidence-codes' : 'clingen'}, axis=1)
a=a.rename({'franklin_evidence_codes' : 'franklin'}, axis=1)
a=a.rename({'franklin_pathogenicity' : 'franklin-pathogenicity'}, axis=1)

a["seq"]=a["seq"].str.replace("-","")
a["seq"]=a["seq"].str.replace("+","")


a["clingen"]=a["clingen"].str.replace("_Very_Strong","")
a["clingen"]=a["clingen"].str.replace("_Strong","")
a["clingen"]=a["clingen"].str.replace("_Moderate","")
a["clingen"]=a["clingen"].str.replace("_Supporting","")
a["clingen"]=a["clingen"].str.replace("_Ignore","")
a["clingen"]=a["clingen"].str.replace("_Stand_Alone","")
a["clingen"]=a["clingen"].str.replace("_Very Strong","")
a["clingen"]=a["clingen"].str.replace("_Stand Alone","")

a["seq"]=a["seq"].str.replace(" ","")
a["clingen"]=a["clingen"].str.replace(" ","")
a["franklin"]=a["franklin"].str.replace(" ","")

#I would like to create a third column (column C) which returns the items that exists in column A but does not exist in column B

a = a.replace(r'^\s*$', np.nan, regex=True)
a.fillna("blank", inplace = True)


sorter=["PVS1","PS1","PS2","PS3","PS4","PM1","PM2","PM3","PM4","PM5","PM6","PP1","PP2","PP3","PP4","PP5","BA1","BS1","BS2","BS3","BS4","BP1","BP2","BP3","BP4","BP5","BP6","BP7"]

#create new df for SEQ by unique analysis comparing to Clingen
df = a[["clingen-pathogenicity","seq-pathogenicity","franklin-pathogenicity","seq","clingen","franklin"]].copy()

df_P=df[df['clingen-pathogenicity'] =="P"]
df_LP=df[df['clingen-pathogenicity'] =="LP"]
df_VUS=df[df['clingen-pathogenicity'] =="VUS"]
df_LB=df[df['clingen-pathogenicity'] =="LB"]
df_B=df[df['clingen-pathogenicity'] =="B"]

dfs=[df_P,df_LP,df_VUS,df_LB,df_B]

seqdfs = []
franklindfs = []
for i in dfs:

##############################################################	Clingen	#####################################################################

	#create columns obtained from the counts of comma separated elements in
	C_evidence=i["clingen"].str.split(pat = ",", expand=True).apply(lambda x : x.value_counts(), axis = 1).fillna(0).astype(int)
	
	#getting missing evidence codes
	C_missing=list(set(sorter)-set(list(C_evidence)))
	
	#assign them into a frame as "0"
	C_evidence[C_missing] = pd.DataFrame([[int(0)]*len(C_missing)], index=C_evidence.index)
	
	#reindexing df like sorter's index
	C_evidence=C_evidence.reindex(sorter, axis='columns')
	
	# Extract all rows as binary arrays
	clingen=C_evidence.iloc[:].values



##############################################################	Seq	#####################################################################

	#create columns obtained from the counts of comma separated elements in CS_seq_unique 
	S_evidence=i["seq"].str.split(pat = ",", expand=True).apply(lambda x : x.value_counts(), axis = 1).fillna(0).astype(int)

	#getting missing evidence codes
	S_missing=list(set(sorter)-set(list(S_evidence)))
	
	#assign them into a frame as "0"
	S_evidence[S_missing] = pd.DataFrame([[int(0)]*len(S_missing)], index=S_evidence.index)
	
	#reindexing df like sorter's index
	S_evidence=S_evidence.reindex(sorter, axis='columns')
	
	# Extract all rows as binary arrays
	seq=S_evidence.iloc[:].values


	seq_report=classification_report(clingen, seq, output_dict=True, zero_division=0, target_names=sorter)
	df = pd.DataFrame(seq_report).transpose().reset_index()
	seq_df = df.iloc[:-4].drop(columns="support")
	seq_df["platform"]="seq"
	seq_df["reference"]=str(i["clingen-pathogenicity"].values[0])
	seqdfs.append(seq_df)


##############################################################	Franklin_P	#####################################################################

	#create columns obtained from the counts of comma separated elements in CS_seq_unique 
	F_evidence=i["franklin"].str.split(pat = ",", expand=True).apply(lambda x : x.value_counts(), axis = 1).fillna(0).astype(int)

	#getting missing evidence codes
	F_missing=list(set(sorter)-set(list(F_evidence)))

	#assign them into a frame as "0"
	F_evidence[F_missing] = pd.DataFrame([[int(0)]*len(F_missing)], index=F_evidence.index)

	#reindexing df like sorter's index
	F_evidence=F_evidence.reindex(sorter, axis='columns')

	# Extract all rows as binary arrays
	franklin=F_evidence.iloc[:].values



	franklin_report=classification_report(clingen, franklin, output_dict=True, zero_division=0, target_names=sorter)
	df = pd.DataFrame(franklin_report).transpose().reset_index()
	franklin_df = df.iloc[:-4].drop(columns="support")
	franklin_df["platform"]="franklin"
	franklin_df["reference"]=str(i["clingen-pathogenicity"].values[0])
	franklindfs.append(franklin_df)




final_seqdf = pd.concat(seqdfs)
final_franklindf = pd.concat(franklindfs)
final_df=pd.concat([final_seqdf,final_franklindf])
final_df=final_df[(final_df != 0).all(axis=1)]
#final_df=final_df.rename({'index' : 'Evidence Codes'}, axis=1)
final_df.to_csv("/Users/islekbro/Desktop/genomize_makale/3_statistics/statistics_outputs/statistics_evidence_code_out.csv",sep="\t", index=False)
