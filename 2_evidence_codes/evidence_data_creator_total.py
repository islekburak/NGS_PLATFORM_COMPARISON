import pandas as pd
import numpy as np
from collections import Counter
from itertools import chain


a=pd.read_csv("/Users/islekbro/Desktop/genomize_makale/2_evidence_codes/evidence_inputs/comparison.txt",sep="\t")

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

a = a.replace(r'^\s*$', np.nan, regex=True)
a.fillna("blank", inplace = True)

a["seq"]=a["seq"].str.replace(" ","")
a["clingen"]=a["clingen"].str.replace(" ","")
a["franklin"]=a["franklin"].str.replace(" ","")

#for difference
inter = a.seq.str.split(',').apply(set) - a.clingen.str.split(',').apply(set).values
a['CS_seq_unique'] = inter.str.join(',')

inter1 = a.clingen.str.split(',').apply(set) - a.seq.str.split(',').apply(set).values
a['CS_clingen_unique'] = inter1.str.join(',')


inter2 = a.franklin.str.split(',').apply(set) - a.clingen.str.split(',').apply(set).values
a['CF_franklin_unique'] = inter2.str.join(',')

inter3 = a.clingen.str.split(',').apply(set) - a.franklin.str.split(',').apply(set).values
a['CF_clingen_unique'] = inter3.str.join(',')

#for intersection
inter4 = a.seq.str.split(',').apply(set) - (a.seq.str.split(',').apply(set) - a.clingen.str.split(',').apply(set).values)
a['CS_common'] = inter4.str.join(',')
inter5 = a.franklin.str.split(',').apply(set) - (a.franklin.str.split(',').apply(set) - a.clingen.str.split(',').apply(set).values)
a['CF_common'] = inter5.str.join(',')


a = a.replace(r'^\s*$', np.nan, regex=True)
a.fillna("none", inplace = True)


sorter=["PVS1","PS1","PS2","PS3","PS4","PM1","PM2","PM3","PM4","PM5","PM6","PP1","PP2","PP3","PP4","PP5","BA1","BS1","BS2","BS3","BS4","BP1","BP2","BP3","BP4","BP5","BP6","BP7"]

######################################################	DIFFERENTIAL SPECIFIC ANALYSIS	#########################################################################

##############################################################	CSS_DF	(Specific)	#####################################################################

#create new df for SEQ by unique analysis comparing to Clingen
CSSdf1 = a[["clingen-pathogenicity","seq-pathogenicity","CS_seq_unique"]].copy()
#CSSdf1.drop(CSSdf1.loc[CSSdf1['clingen-pathogenicity']!="P"].index, inplace=True)
CSSdf1.drop(CSSdf1.loc[CSSdf1['CS_seq_unique']=="none"].index, inplace=True)
CSSdf1.drop(CSSdf1.loc[CSSdf1['CS_seq_unique']=="blank"].index, inplace=True)

#create columns obtained from the counts of comma separated elements in CS_seq_unique 
CSSevidence1=CSSdf1["CS_seq_unique"].str.split(pat = ",", expand=True).apply(lambda x : x.value_counts(), axis = 1).fillna(0).astype(int)
#merge dfs for final df
CSS_df1=pd.concat([CSSdf1,CSSevidence1],axis=1)

#list column names for create dictionary with present evidence codes
CSSdf_codes1=list(CSS_df1.columns.values)
#delete pathogenicity values from the list
del CSSdf_codes1[0:3]
#create dictionary from the values extracted from the list, that is created above. #add string of "sum" for the next step. 
CSSdct1 = {CSSdf_codes1[i]: "sum" for i in range(0, len(CSSdf_codes1))}


CSS1=CSS_df1.agg(CSSdct1)
#CSS1=CSS1.reset_index()
#CSS1["info"]="CSS_"+CSS1["clingen-pathogenicity"]
#del CSS1["clingen-pathogenicity"]


CSS1=CSS1.reindex(sorter, axis='columns').fillna(0)
#create the final df (grouped by pathogenicity codes) with the count of evidence codes
#CSS1=CSS_df1.groupby(["clingen-pathogenicity", "seq-pathogenicity"]).agg(CSSdct1)

"""
#CSS_P=CSS1.loc[CSS1['info']=="CSS_P"]
#CSS_LP=CSS1.loc[CSS1['info']=="CSS_LP"]
#CSS_VUS=CSS1.loc[CSS1['info']=="CSS_VUS"]
#CSS_LB=CSS1.loc[CSS1['info']=="CSS_LB"]
#CSS_B=CSS1.loc[CSS1['info']=="CSS_B"]
"""


##############################################################	CSC_DF	(Specific)	#####################################################################

CSCdf1 = a[["clingen-pathogenicity","seq-pathogenicity","CS_clingen_unique"]].copy()
CSCdf1.drop(CSCdf1.loc[CSCdf1['CS_clingen_unique']=="none"].index, inplace=True)
CSCdf1.drop(CSCdf1.loc[CSCdf1['CS_clingen_unique']=="blank"].index, inplace=True)

#create columns obtained from the counts of comma separated elements in CS_seq_unique 
CSCevidence1=CSCdf1["CS_clingen_unique"].str.split(pat = ",", expand=True).apply(lambda x : x.value_counts(), axis = 1).fillna(0).astype(int)
#merge dfs for final df
CSC_df1=pd.concat([CSCdf1,CSCevidence1],axis=1)

#list column names for create dictionary with present evidence codes
CSCdf_codes1=list(CSC_df1.columns.values)
#delete pathogenicity values from the list
del CSCdf_codes1[0:3]
#create dictionary from the values extracted from the list, that is created above. #add string of "sum" for the next step. 
CSCdct1 = {CSCdf_codes1[i]: "sum" for i in range(0, len(CSCdf_codes1))}

CSC1=CSC_df1.agg(CSCdct1)
#CSC1=CSC1.reset_index()
#CSC1["info"]="CSC_"+CSC1["clingen-pathogenicity"]
#del CSC1["clingen-pathogenicity"]


CSC1=CSC1.reindex(sorter, axis='columns').fillna(0)
#create the final df (grouped by pathogenicity codes) with the count of evidence codes
#CSC1=CSC_df1.groupby(["clingen-pathogenicity", "seq-pathogenicity"]).agg(CSCdct1)

"""
CSC_P=CSC1.loc[CSC1['info']=="CSC_P"]
CSC_LP=CSC1.loc[CSC1['info']=="CSC_LP"]
CSC_VUS=CSC1.loc[CSC1['info']=="CSC_VUS"]
CSC_LB=CSC1.loc[CSC1['info']=="CSC_LB"]
CSC_B=CSC1.loc[CSC1['info']=="CSC_B"]
"""


##############################################################	CFF_DF	(Specific)	#####################################################################

#create new df for SEQ by unique analysis comparing to Clingen
CFFdf1 = a[["clingen-pathogenicity","franklin-pathogenicity","CF_franklin_unique"]].copy()
#CSSdf1.drop(CSSdf1.loc[CSSdf1['clingen-pathogenicity']!="P"].index, inplace=True)
CFFdf1.drop(CFFdf1.loc[CFFdf1['CF_franklin_unique']=="none"].index, inplace=True)
CFFdf1.drop(CFFdf1.loc[CFFdf1['CF_franklin_unique']=="blank"].index, inplace=True)

#create columns obtained from the counts of comma separated elements in CS_seq_unique 
CFFevidence1=CFFdf1["CF_franklin_unique"].str.split(pat = ",", expand=True).apply(lambda x : x.value_counts(), axis = 1).fillna(0).astype(int)
#merge dfs for final df
CFF_df1=pd.concat([CFFdf1,CFFevidence1],axis=1)

#list column names for create dictionary with present evidence codes
CFFdf_codes1=list(CFF_df1.columns.values)
#delete pathogenicity values from the list
del CFFdf_codes1[0:3]
#create dictionary from the values extracted from the list, that is created above. #add string of "sum" for the next step. 
CFFdct1 = {CFFdf_codes1[i]: "sum" for i in range(0, len(CFFdf_codes1))}


CFF1=CFF_df1.agg(CFFdct1)
#CFF1=CFF_df1.groupby(["clingen-pathogenicity"]).agg(CFFdct1)
#CFF1=CFF1.reset_index()
#CFF1["info"]="CFF_"+CFF1["clingen-pathogenicity"]
#del CFF1["clingen-pathogenicity"]

CFF1=CFF1.reindex(sorter, axis='columns').fillna(0)

#create the final df (grouped by pathogenicity codes) with the count of evidence codes
#CSS1=CSS_df1.groupby(["clingen-pathogenicity", "seq-pathogenicity"]).agg(CSSdct1)

"""
CFF_P=CFF1.loc[CFF1['info']=="CFF_P"]
CFF_LP=CFF1.loc[CFF1['info']=="CFF_LP"]
CFF_VUS=CFF1.loc[CFF1['info']=="CFF_VUS"]
CFF_LB=CFF1.loc[CFF1['info']=="CFF_LB"]
CFF_B=CFF1.loc[CFF1['info']=="CFF_B"]
"""


##############################################################	CFC_DF	(Specific)	#####################################################################

CFCdf1 = a[["clingen-pathogenicity","franklin-pathogenicity","CF_clingen_unique"]].copy()
CFCdf1.drop(CFCdf1.loc[CFCdf1['CF_clingen_unique']=="none"].index, inplace=True)
CFCdf1.drop(CFCdf1.loc[CFCdf1['CF_clingen_unique']=="blank"].index, inplace=True)

#create columns obtained from the counts of comma separated elements in CS_seq_unique 
CFCevidence1=CFCdf1["CF_clingen_unique"].str.split(pat = ",", expand=True).apply(lambda x : x.value_counts(), axis = 1).fillna(0).astype(int)
#merge dfs for final df
CFC_df1=pd.concat([CFCdf1,CFCevidence1],axis=1)

#list column names for create dictionary with present evidence codes
CFCdf_codes1=list(CFC_df1.columns.values)
#delete pathogenicity values from the list
del CFCdf_codes1[0:3]
#create dictionary from the values extracted from the list, that is created above. #add string of "sum" for the next step. 
CFCdct1 = {CFCdf_codes1[i]: "sum" for i in range(0, len(CFCdf_codes1))}


CFC1=CFC_df1.agg(CFCdct1)
#CFC1=CFC_df1.groupby(["clingen-pathogenicity"]).agg(CFCdct1)
#CFC1=CFC1.reset_index()
#CFC1["info"]="CFC_"+CFC1["clingen-pathogenicity"]
#del CFC1["clingen-pathogenicity"]


CFC1=CFC1.reindex(sorter, axis='columns').fillna(0)
#create the final df (grouped by pathogenicity codes) with the count of evidence codes
#CSC1=CSC_df1.groupby(["clingen-pathogenicity", "seq-pathogenicity"]).agg(CSCdct1)

"""
CFC_P=CFC1.loc[CFC1['info']=="CFC_P"]
CFC_LP=CFC1.loc[CFC1['info']=="CFC_LP"]
CFC_VUS=CFC1.loc[CFC1['info']=="CFC_VUS"]
CFC_LB=CFC1.loc[CFC1['info']=="CFC_LB"]
CFC_B=CFC1.loc[CFC1['info']=="CFC_B"]
"""



#####################################################################################################################################################
##############################################################	C_DF	(Total)	#####################################################################

#create new df for SEQ by unique analysis comparing to Clingen
Cdf1 = a[["clingen-pathogenicity","seq-pathogenicity","clingen"]].copy()
#CSSdf1.drop(CSSdf1.loc[CSSdf1['clingen-pathogenicity']!="P"].index, inplace=True)
Cdf1.drop(Cdf1.loc[Cdf1['clingen']=="none"].index, inplace=True)
Cdf1.drop(Cdf1.loc[Cdf1['clingen']=="blank"].index, inplace=True)

#create columns obtained from the counts of comma separated elements in CS_seq_unique 
Cevidence1=Cdf1["clingen"].str.split(pat = ",", expand=True).apply(lambda x : x.value_counts(), axis = 1).fillna(0).astype(int)
#merge dfs for final df
C_df1=pd.concat([Cdf1,Cevidence1],axis=1)

#list column names for create dictionary with present evidence codes
Cdf_codes1=list(C_df1.columns.values)
#delete pathogenicity values from the list
del Cdf_codes1[0:3]
#create dictionary from the values extracted from the list, that is created above. #add string of "sum" for the next step. 
Cdct1 = {Cdf_codes1[i]: "sum" for i in range(0, len(Cdf_codes1))}


C1=C_df1.agg(Cdct1)
#C1=C_df1.groupby(["clingen-pathogenicity"]).agg(Cdct1)
#C1=C1.reset_index()
#C1["info"]="C_"+C1["clingen-pathogenicity"]
#del C1["clingen-pathogenicity"]

C1=C1.reindex(sorter, axis='columns').fillna(0)
#create the final df (grouped by pathogenicity codes) with the count of evidence codes
#CSS1=CSS_df1.groupby(["clingen-pathogenicity", "seq-pathogenicity"]).agg(CSSdct1)

"""
C_P=C1.loc[C1['info']=="C_P"]
C_LP=C1.loc[C1['info']=="C_LP"]
C_VUS=C1.loc[C1['info']=="C_VUS"]
C_LB=C1.loc[C1['info']=="C_LB"]
C_B=C1.loc[C1['info']=="C_B"]
"""



##############################################################	S_DF	(Total)	#####################################################################

#create new df for SEQ by unique analysis comparing to Clingen
Sdf1 = a[["clingen-pathogenicity","seq-pathogenicity","seq"]].copy()
#CSSdf1.drop(CSSdf1.loc[CSSdf1['clingen-pathogenicity']!="P"].index, inplace=True)
Sdf1.drop(Sdf1.loc[Sdf1['seq']=="none"].index, inplace=True)
Sdf1.drop(Sdf1.loc[Sdf1['seq']=="blank"].index, inplace=True)

#create columns obtained from the counts of comma separated elements in CS_seq_unique 
Sevidence1=Sdf1["seq"].str.split(pat = ",", expand=True).apply(lambda x : x.value_counts(), axis = 1).fillna(0).astype(int)
#merge dfs for final df
S_df1=pd.concat([Sdf1,Sevidence1],axis=1)

#list column names for create dictionary with present evidence codes
Sdf_codes1=list(S_df1.columns.values)
#delete pathogenicity values from the list
del Sdf_codes1[0:3]
#create dictionary from the values extracted from the list, that is created above. #add string of "sum" for the next step. 
Sdct1 = {Sdf_codes1[i]: "sum" for i in range(0, len(Sdf_codes1))}


S1=S_df1.agg(Sdct1)
#S1=S_df1.groupby(["clingen-pathogenicity"]).agg(Sdct1)
#S1=S1.reset_index()
#S1["info"]="S_"+S1["clingen-pathogenicity"]
#del S1["clingen-pathogenicity"]

S1=S1.reindex(sorter, axis='columns').fillna(0)
#create the final df (grouped by pathogenicity codes) with the count of evidence codes
#CSS1=CSS_df1.groupby(["clingen-pathogenicity", "seq-pathogenicity"]).agg(CSSdct1)


"""
S_P=S1.loc[S1['info']=="S_P"]
S_LP=S1.loc[S1['info']=="S_LP"]
S_VUS=S1.loc[S1['info']=="S_VUS"]
S_LB=S1.loc[S1['info']=="S_LB"]
S_B=S1.loc[S1['info']=="S_B"]
"""


##############################################################	F_DF	(Total)	#####################################################################

#create new df for SEQ by unique analysis comparing to Clingen
Fdf1 = a[["clingen-pathogenicity","seq-pathogenicity","franklin"]].copy()
#CSSdf1.drop(CSSdf1.loc[CSSdf1['clingen-pathogenicity']!="P"].index, inplace=True)
Fdf1.drop(Fdf1.loc[Fdf1['franklin']=="none"].index, inplace=True)
Fdf1.drop(Fdf1.loc[Fdf1['franklin']=="blank"].index, inplace=True)

#create columns obtained from the counts of comma separated elements in CS_seq_unique 
Fevidence1=Fdf1["franklin"].str.split(pat = ",", expand=True).apply(lambda x : x.value_counts(), axis = 1).fillna(0).astype(int)
#merge dfs for final df
F_df1=pd.concat([Fdf1,Fevidence1],axis=1)

#list column names for create dictionary with present evidence codes
Fdf_codes1=list(F_df1.columns.values)
#delete pathogenicity values from the list
del Fdf_codes1[0:3]
#create dictionary from the values extracted from the list, that is created above. #add string of "sum" for the next step. 
Fdct1 = {Fdf_codes1[i]: "sum" for i in range(0, len(Fdf_codes1))}


F1=F_df1.agg(Fdct1)
#F1=F_df1.groupby(["clingen-pathogenicity"]).agg(Fdct1)
#F1=F1.reset_index()
#F1["info"]="F_"+F1["clingen-pathogenicity"]
#del F1["clingen-pathogenicity"]


F1=F1.reindex(sorter, axis='columns').fillna(0)
#create the final df (grouped by pathogenicity codes) with the count of evidence codes
#CSS1=CSS_df1.groupby(["clingen-pathogenicity", "seq-pathogenicity"]).agg(CSSdct1)

"""
F_P=F1.loc[F1['info']=="F_P"]
F_LP=F1.loc[F1['info']=="F_LP"]
F_VUS=F1.loc[F1['info']=="F_VUS"]
F_LB=F1.loc[F1['info']=="F_LB"]
F_B=F1.loc[F1['info']=="F_B"]
"""



##for common

#sorter.remove("info")

fin=pd.concat([CSS1,CSC1,CFF1,CFC1,C1,S1,F1],axis=1).reset_index().rename(
	columns={
	'index': 'codes',
	0:"CSS",
	1:"CSC",
	2:"CFF",
	3:"CFC",
	4:"C",
	5:"S",
	6:"F"
	})

sublist=[
"main","total","main","total",
"main","total","main","total",
"main","total","main","total",
"main","total","main","total",
"main","total","main","total",
"main","total","main","total",
"main","total","main","total",
"main","total","main","total",
"main","total","main","total",
"main","total","main","total",
"main","total","main","total",
"main","total","main","total",
"main","total","main","total",
"main","total","main","total"
]



fin["CandS"]=fin["S"]-fin["CSS"]
fin["CandF"]=fin["F"]-fin["CFF"]

fin = fin.astype({
	'CSS': 'int',
	'CSC': 'int',
	'CFF':'int',
	'CFC':'int',
	'C':'int',
	'S':'int',
	'F':'int',
	'CandS':'int',
	'CandF':'int'
	})


CS_common = fin[["CandS"]].copy().rename(columns={ 'CandS': 'count'})
CS_common["codes"]=sorter
CS_common["platform"]="CS"
CS_common["label"]=CS_common["codes"] + "=(" + CS_common["count"].astype(str) + ")"

CF_common = fin[["CandF"]].copy().rename(columns={ 'CandF': 'count'})
CF_common["codes"]=sorter
CF_common["platform"]="CF"
CF_common["label"]=CF_common["codes"] + "=(" + CF_common["count"].astype(str) + ")"

common=pd.concat([CS_common,CF_common])
common.drop(common.loc[common['count']==0].index, inplace=True)
common.reset_index(drop=True, inplace=True)
common.to_csv("/Users/islekbro/Desktop/genomize_makale/2_evidence_codes/evidence_outputs/common_total.csv",sep="\t",index=False)


k=fin.iloc[np.arange(len(fin)).repeat(2)]
####################################################### FOR clingen_vs_seq #################################################

#for CSS
x=fin["CSS"].tolist()
y=fin["CandS"].tolist()
xy = list(zip(x,y))
new=list(chain.from_iterable(xy))
CSS_df = pd.DataFrame(new,columns =['variant_count'])

freqlist=round((k['CSS']/k['S']*100),1).tolist()
frqlist = [0 if x != x else x for x in freqlist]
frqlist=[int(i) for i in frqlist]
CSS_df["freq"]=frqlist
CSS_df["freq"]="~"+CSS_df["freq"].astype(str) + "%"


mainlist=k["CSS"].tolist()
CSS_df["main"]=mainlist

totallist=k["S"].tolist()
CSS_df["total"]=totallist

CSS_df["total_label"] = "(" + CSS_df["total"].astype(str) + ")"
CSS_df[["platforms", "orientation"]]=["Clingen_vs_SEQ","up"]

CSS_df.insert(0,"codes",k["codes"].tolist())
CSS_df["sub"]=sublist



#for CSC
x=fin["CSC"].tolist()
y=fin["CandS"].tolist()
xy = list(zip(x,y))
new=list(chain.from_iterable(xy))
CSC_df = pd.DataFrame(new,columns =['variant_count'])

freqlist=round((k['CSC']/k['C']*100),1).tolist()
frqlist = [0 if x != x else x for x in freqlist]
frqlist=[int(i) for i in frqlist]
CSC_df["freq"]=frqlist
CSC_df["freq"]="~"+CSC_df["freq"].astype(str) + "%"


mainlist=k["CSC"].tolist()
CSC_df["main"]=mainlist

totallist=k["C"].tolist()
CSC_df["total"]=totallist

CSC_df["total_label"] = "(" + CSC_df["total"].astype(str) + ")"
CSC_df[["platforms", "orientation"]]=["Clingen_vs_SEQ","down"]

CSC_df.insert(0,"codes",k["codes"].tolist())
CSC_df["sub"]=sublist


####################################################### FOR clingen_vs_franklin #################################################

#for CFF
x=fin["CFF"].tolist()
y=fin["CandF"].tolist()
xy = list(zip(x,y))
new=list(chain.from_iterable(xy))
CFF_df = pd.DataFrame(new,columns =['variant_count'])

freqlist=round((k['CFF']/k['F']*100),1).tolist()
frqlist = [0 if x != x else x for x in freqlist]
frqlist=[int(i) for i in frqlist]
CFF_df["freq"]=frqlist
CFF_df["freq"]="~"+CFF_df["freq"].astype(str) + "%"


mainlist=k["CFF"].tolist()
CFF_df["main"]=mainlist

totallist=k["F"].tolist()
CFF_df["total"]=totallist

CFF_df["total_label"] = "(" + CFF_df["total"].astype(str) + ")"
CFF_df[["platforms", "orientation"]]=["Clingen_vs_Franklin","up"]

CFF_df.insert(0,"codes",k["codes"].tolist())
CFF_df["sub"]=sublist



#for CFC

x=fin["CFC"].tolist()
y=fin["CandF"].tolist()
xy = list(zip(x,y))
new=list(chain.from_iterable(xy))
CFC_df = pd.DataFrame(new,columns =['variant_count'])

freqlist=round((k['CFC']/k['C']*100),1).tolist()
frqlist = [0 if x != x else x for x in freqlist]
frqlist=[int(i) for i in frqlist]
CFC_df["freq"]=frqlist
CFC_df["freq"]="~"+CFC_df["freq"].astype(str) + "%"


mainlist=k["CFC"].tolist()
CFC_df["main"]=mainlist

totallist=k["C"].tolist()
CFC_df["total"]=totallist

CFC_df["total_label"] = "(" + CFC_df["total"].astype(str) + ")"
CFC_df[["platforms", "orientation"]]=["Clingen_vs_Franklin","down"]

CFC_df.insert(0,"codes",k["codes"].tolist())
CFC_df["sub"]=sublist



final = pd.concat([CSS_df,CSC_df,CFF_df,CFC_df], axis=0)
final.to_csv("/Users/islekbro/Desktop/genomize_makale/2_evidence_codes/evidence_outputs/evidence_total.csv", sep="\t", index=False)
