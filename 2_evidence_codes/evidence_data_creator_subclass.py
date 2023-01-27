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


sorter=["info","PVS1","PS1","PS2","PS3","PS4","PM1","PM2","PM3","PM4","PM5","PM6","PP1","PP2","PP3","PP4","PP5","BA1","BS1","BS2","BS3","BS4","BP1","BP2","BP3","BP4","BP5","BP6","BP7"]

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


CSS1=CSS_df1.groupby(["clingen-pathogenicity"]).agg(CSSdct1)
CSS1=CSS1.reset_index()
CSS1["info"]="CSS_"+CSS1["clingen-pathogenicity"]
del CSS1["clingen-pathogenicity"]


CSS1=CSS1.reindex(sorter, axis='columns').fillna(0)
#create the final df (grouped by pathogenicity codes) with the count of evidence codes
#CSS1=CSS_df1.groupby(["clingen-pathogenicity", "seq-pathogenicity"]).agg(CSSdct1)

CSS_P=CSS1.loc[CSS1['info']=="CSS_P"]
CSS_LP=CSS1.loc[CSS1['info']=="CSS_LP"]
CSS_VUS=CSS1.loc[CSS1['info']=="CSS_VUS"]
CSS_LB=CSS1.loc[CSS1['info']=="CSS_LB"]
CSS_B=CSS1.loc[CSS1['info']=="CSS_B"]



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

CSC1=CSC_df1.groupby(["clingen-pathogenicity"]).agg(CSCdct1)
CSC1=CSC1.reset_index()
CSC1["info"]="CSC_"+CSC1["clingen-pathogenicity"]
del CSC1["clingen-pathogenicity"]


CSC1=CSC1.reindex(sorter, axis='columns').fillna(0)
#create the final df (grouped by pathogenicity codes) with the count of evidence codes
#CSC1=CSC_df1.groupby(["clingen-pathogenicity", "seq-pathogenicity"]).agg(CSCdct1)


CSC_P=CSC1.loc[CSC1['info']=="CSC_P"]
CSC_LP=CSC1.loc[CSC1['info']=="CSC_LP"]
CSC_VUS=CSC1.loc[CSC1['info']=="CSC_VUS"]
CSC_LB=CSC1.loc[CSC1['info']=="CSC_LB"]
CSC_B=CSC1.loc[CSC1['info']=="CSC_B"]




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


CFF1=CFF_df1.groupby(["clingen-pathogenicity"]).agg(CFFdct1)
CFF1=CFF1.reset_index()
CFF1["info"]="CFF_"+CFF1["clingen-pathogenicity"]
del CFF1["clingen-pathogenicity"]

CFF1=CFF1.reindex(sorter, axis='columns').fillna(0)

#create the final df (grouped by pathogenicity codes) with the count of evidence codes
#CSS1=CSS_df1.groupby(["clingen-pathogenicity", "seq-pathogenicity"]).agg(CSSdct1)

CFF_P=CFF1.loc[CFF1['info']=="CFF_P"]
CFF_LP=CFF1.loc[CFF1['info']=="CFF_LP"]
CFF_VUS=CFF1.loc[CFF1['info']=="CFF_VUS"]
CFF_LB=CFF1.loc[CFF1['info']=="CFF_LB"]
CFF_B=CFF1.loc[CFF1['info']=="CFF_B"]



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

CFC1=CFC_df1.groupby(["clingen-pathogenicity"]).agg(CFCdct1)
CFC1=CFC1.reset_index()
CFC1["info"]="CFC_"+CFC1["clingen-pathogenicity"]
del CFC1["clingen-pathogenicity"]


CFC1=CFC1.reindex(sorter, axis='columns').fillna(0)
#create the final df (grouped by pathogenicity codes) with the count of evidence codes
#CSC1=CSC_df1.groupby(["clingen-pathogenicity", "seq-pathogenicity"]).agg(CSCdct1)


CFC_P=CFC1.loc[CFC1['info']=="CFC_P"]
CFC_LP=CFC1.loc[CFC1['info']=="CFC_LP"]
CFC_VUS=CFC1.loc[CFC1['info']=="CFC_VUS"]
CFC_LB=CFC1.loc[CFC1['info']=="CFC_LB"]
CFC_B=CFC1.loc[CFC1['info']=="CFC_B"]




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


C1=C_df1.groupby(["clingen-pathogenicity"]).agg(Cdct1)
C1=C1.reset_index()
C1["info"]="C_"+C1["clingen-pathogenicity"]
del C1["clingen-pathogenicity"]

C1=C1.reindex(sorter, axis='columns').fillna(0)
#create the final df (grouped by pathogenicity codes) with the count of evidence codes
#CSS1=CSS_df1.groupby(["clingen-pathogenicity", "seq-pathogenicity"]).agg(CSSdct1)

C_P=C1.loc[C1['info']=="C_P"]
C_LP=C1.loc[C1['info']=="C_LP"]
C_VUS=C1.loc[C1['info']=="C_VUS"]
C_LB=C1.loc[C1['info']=="C_LB"]
C_B=C1.loc[C1['info']=="C_B"]




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


S1=S_df1.groupby(["clingen-pathogenicity"]).agg(Sdct1)
S1=S1.reset_index()
S1["info"]="S_"+S1["clingen-pathogenicity"]
del S1["clingen-pathogenicity"]

S1=S1.reindex(sorter, axis='columns').fillna(0)
#create the final df (grouped by pathogenicity codes) with the count of evidence codes
#CSS1=CSS_df1.groupby(["clingen-pathogenicity", "seq-pathogenicity"]).agg(CSSdct1)

S_P=S1.loc[S1['info']=="S_P"]
S_LP=S1.loc[S1['info']=="S_LP"]
S_VUS=S1.loc[S1['info']=="S_VUS"]
S_LB=S1.loc[S1['info']=="S_LB"]
S_B=S1.loc[S1['info']=="S_B"]



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


F1=F_df1.groupby(["clingen-pathogenicity"]).agg(Fdct1)
F1=F1.reset_index()
F1["info"]="F_"+F1["clingen-pathogenicity"]
del F1["clingen-pathogenicity"]


F1=F1.reindex(sorter, axis='columns').fillna(0)
#create the final df (grouped by pathogenicity codes) with the count of evidence codes
#CSS1=CSS_df1.groupby(["clingen-pathogenicity", "seq-pathogenicity"]).agg(CSSdct1)

F_P=F1.loc[F1['info']=="F_P"]
F_LP=F1.loc[F1['info']=="F_LP"]
F_VUS=F1.loc[F1['info']=="F_VUS"]
F_LB=F1.loc[F1['info']=="F_LB"]
F_B=F1.loc[F1['info']=="F_B"]



sorter.remove("info")
df = pd.DataFrame(sorter,columns =['codes'])

P=pd.concat([CSS_P,CSC_P,CFF_P,CFC_P,C_P,S_P,F_P],ignore_index=True)
P=P.set_index("info").T.reset_index(drop=True).astype(int)
P = pd.concat([df, P], axis=1)

LP=pd.concat([CSS_LP,CSC_LP,CFF_LP,CFC_LP,C_LP,S_LP,F_LP],ignore_index=True)
LP=LP.set_index("info").T.reset_index(drop=True).astype(int)
LP = pd.concat([df, LP], axis=1)

VUS=pd.concat([CSS_VUS,CSC_VUS,CFF_VUS,CFC_VUS,C_VUS,S_VUS,F_VUS],ignore_index=True)
VUS=VUS.set_index("info").T.reset_index(drop=True).astype(int)
VUS = pd.concat([df, VUS], axis=1)

LB=pd.concat([CSS_LB,CSC_LB,CFF_LB,CFC_LB,C_LB,S_LB,F_LB],ignore_index=True)
LB=LB.set_index("info").T.reset_index(drop=True).astype(int)
LB = pd.concat([df, LB], axis=1)

B=pd.concat([CSS_B,CSC_B,CFF_B,CFC_B,C_B,S_B,F_B],ignore_index=True)
B=B.set_index("info").T.reset_index(drop=True).astype(int)
B = pd.concat([df, B], axis=1)



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


####################################################### FOR P #################################################
kp=P.iloc[np.arange(len(P)).repeat(2)]
#for P_CSS
P["CSS_P_diff"]=P['S_P']-P['CSS_P']

x=P["CSS_P"].tolist()
y=P["CSS_P_diff"].tolist()
xy = list(zip(x,y))
new=list(chain.from_iterable(xy))
CSS_P_df = pd.DataFrame(new,columns =['variant_count'])

freqlist=round((kp['CSS_P']/kp['S_P']*100),1).tolist()
frqlist = [0 if x != x else x for x in freqlist]
frqlist=[int(i) for i in frqlist]
CSS_P_df["freq"]=frqlist
CSS_P_df["freq"]="~"+CSS_P_df["freq"].astype(str) + "%"


mainlist=kp["CSS_P"].tolist()
CSS_P_df["main"]=mainlist

totallist=kp["S_P"].tolist()
CSS_P_df["total"]=totallist

CSS_P_df["total_label"] = "(" + CSS_P_df["total"].astype(str) + ")"
CSS_P_df[["platforms", "orientation"]]=["Clingen_vs_SEQ","up"]

CSS_P_df.insert(0,"codes",kp["codes"].tolist())
CSS_P_df["sub"]=sublist



#for P_CSC
P["CSC_P_diff"]=P['C_P']-P['CSC_P']

x=P["CSC_P"].tolist()
y=P["CSC_P_diff"].tolist()
xy = list(zip(x,y))
new=list(chain.from_iterable(xy))
CSC_P_df = pd.DataFrame(new,columns =['variant_count'])

freqlist=round((kp['CSC_P']/kp['C_P']*100),1).tolist()
frqlist = [0 if x != x else x for x in freqlist]
frqlist=[int(i) for i in frqlist]
CSC_P_df["freq"]=frqlist
CSC_P_df["freq"]="~"+CSC_P_df["freq"].astype(str) + "%"


mainlist=kp["CSC_P"].tolist()
CSC_P_df["main"]=mainlist

totallist=kp["C_P"].tolist()
CSC_P_df["total"]=totallist

CSC_P_df["total_label"] = "(" + CSC_P_df["total"].astype(str) + ")"
CSC_P_df[["platforms", "orientation"]]=["Clingen_vs_SEQ","down"]

CSC_P_df.insert(0,"codes",kp["codes"].tolist())
CSC_P_df["sub"]=sublist


#for P_CFF
P["CFF_P_diff"]=P['F_P']-P['CFF_P']

x=P["CFF_P"].tolist()
y=P["CFF_P_diff"].tolist()
xy = list(zip(x,y))
new=list(chain.from_iterable(xy))
CFF_P_df = pd.DataFrame(new,columns =['variant_count'])

freqlist=round((kp['CFF_P']/kp['F_P']*100),1).tolist()
frqlist = [0 if x != x else x for x in freqlist]
frqlist=[int(i) for i in frqlist]
CFF_P_df["freq"]=frqlist
CFF_P_df["freq"]="~"+CFF_P_df["freq"].astype(str) + "%"


mainlist=kp["CFF_P"].tolist()
CFF_P_df["main"]=mainlist

totallist=kp["F_P"].tolist()
CFF_P_df["total"]=totallist

CFF_P_df["total_label"] = "(" + CFF_P_df["total"].astype(str) + ")"
CFF_P_df[["platforms", "orientation"]]=["Clingen_vs_Franklin","up"]

CFF_P_df.insert(0,"codes",kp["codes"].tolist())
CFF_P_df["sub"]=sublist



#for P_CFC
P["CFC_P_diff"]=P['C_P']-P['CFC_P']

x=P["CFC_P"].tolist()
y=P["CFC_P_diff"].tolist()
xy = list(zip(x,y))
new=list(chain.from_iterable(xy))
CFC_P_df = pd.DataFrame(new,columns =['variant_count'])

freqlist=round((kp['CFC_P']/kp['C_P']*100),1).tolist()
frqlist = [0 if x != x else x for x in freqlist]
frqlist=[int(i) for i in frqlist]
CFC_P_df["freq"]=frqlist
CFC_P_df["freq"]="~"+CFC_P_df["freq"].astype(str) + "%"


mainlist=kp["CFC_P"].tolist()
CFC_P_df["main"]=mainlist

totallist=kp["C_P"].tolist()
CFC_P_df["total"]=totallist

CFC_P_df["total_label"] = "(" + CFC_P_df["total"].astype(str) + ")"
CFC_P_df[["platforms", "orientation"]]=["Clingen_vs_Franklin","down"]

CFC_P_df.insert(0,"codes",kp["codes"].tolist())
CFC_P_df["sub"]=sublist



final_P = pd.concat([CSS_P_df,CSC_P_df,CFF_P_df,CFC_P_df], axis=0)
final_P.to_csv("/Users/islekbro/Desktop/genomize_makale/2_evidence_codes/evidence_outputs/evidence_P.csv", sep="\t", index=False)

##for common P

CS_P_common = P[["CSS_P_diff"]].copy().rename(columns={ 'CSS_P_diff': 'count'})
CS_P_common["codes"]=sorter
CS_P_common["platform"]="CS_P"
CS_P_common["label"]=CS_P_common["codes"] + "=(" + CS_P_common["count"].astype(str) + ")"

CF_P_common = P[["CFF_P_diff"]].copy().rename(columns={ 'CFF_P_diff': 'count'})
CF_P_common["codes"]=sorter
CF_P_common["platform"]="CF_P"
CF_P_common["label"]=CF_P_common["codes"] + "=(" + CF_P_common["count"].astype(str) + ")"

common_P=pd.concat([CS_P_common,CF_P_common])
common_P.drop(common_P.loc[common_P['count']==0].index, inplace=True)
common_P.reset_index(drop=True, inplace=True)
common_P.to_csv("/Users/islekbro/Desktop/genomize_makale/2_evidence_codes/evidence_outputs/common_P.csv",sep="\t",index=False)



####################################################### FOR LP #################################################
kLP=LP.iloc[np.arange(len(LP)).repeat(2)]
#for LP_CSS
LP["CSS_LP_diff"]=LP['S_LP']-LP['CSS_LP']

x=LP["CSS_LP"].tolist()
y=LP["CSS_LP_diff"].tolist()
xy = list(zip(x,y))
new=list(chain.from_iterable(xy))
CSS_LP_df = pd.DataFrame(new,columns =['variant_count'])

freqlist=round((kLP['CSS_LP']/kLP['S_LP']*100),1).tolist()
frqlist = [0 if x != x else x for x in freqlist]
frqlist=[int(i) for i in frqlist]
CSS_LP_df["freq"]=frqlist
CSS_LP_df["freq"]="~"+CSS_LP_df["freq"].astype(str) + "%"


mainlist=kLP["CSS_LP"].tolist()
CSS_LP_df["main"]=mainlist

totallist=kLP["S_LP"].tolist()
CSS_LP_df["total"]=totallist

CSS_LP_df["total_label"] = "(" + CSS_LP_df["total"].astype(str) + ")"
CSS_LP_df[["platforms", "orientation"]]=["Clingen_vs_SEQ","up"]

CSS_LP_df.insert(0,"codes",kLP["codes"].tolist())
CSS_LP_df["sub"]=sublist


#for LP_CSC
LP["CSC_LP_diff"]=LP['C_LP']-LP['CSC_LP']

x=LP["CSC_LP"].tolist()
y=LP["CSC_LP_diff"].tolist()
xy = list(zip(x,y))
new=list(chain.from_iterable(xy))
CSC_LP_df = pd.DataFrame(new,columns =['variant_count'])

freqlist=round((kLP['CSC_LP']/kLP['C_LP']*100),1).tolist()
frqlist = [0 if x != x else x for x in freqlist]
frqlist=[int(i) for i in frqlist]
CSC_LP_df["freq"]=frqlist
CSC_LP_df["freq"]="~"+CSC_LP_df["freq"].astype(str) + "%"


mainlist=kLP["CSC_LP"].tolist()
CSC_LP_df["main"]=mainlist

totallist=kLP["C_LP"].tolist()
CSC_LP_df["total"]=totallist

CSC_LP_df["total_label"] = "(" + CSC_LP_df["total"].astype(str) + ")"
CSC_LP_df[["platforms", "orientation"]]=["Clingen_vs_SEQ","down"]

CSC_LP_df.insert(0,"codes",kLP["codes"].tolist())
CSC_LP_df["sub"]=sublist


#for LP_CFF
LP["CFF_LP_diff"]=LP['F_LP']-LP['CFF_LP']

x=LP["CFF_LP"].tolist()
y=LP["CFF_LP_diff"].tolist()
xy = list(zip(x,y))
new=list(chain.from_iterable(xy))
CFF_LP_df = pd.DataFrame(new,columns =['variant_count'])

freqlist=round((kLP['CFF_LP']/kLP['F_LP']*100),1).tolist()
frqlist = [0 if x != x else x for x in freqlist]
frqlist=[int(i) for i in frqlist]
CFF_LP_df["freq"]=frqlist
CFF_LP_df["freq"]="~"+CFF_LP_df["freq"].astype(str) + "%"


mainlist=kLP["CFF_LP"].tolist()
CFF_LP_df["main"]=mainlist

totallist=kLP["F_LP"].tolist()
CFF_LP_df["total"]=totallist

CFF_LP_df["total_label"] = "(" + CFF_LP_df["total"].astype(str) + ")"
CFF_LP_df[["platforms", "orientation"]]=["Clingen_vs_Franklin","up"]

CFF_LP_df.insert(0,"codes",kLP["codes"].tolist())
CFF_LP_df["sub"]=sublist



#for LP_CFC
LP["CFC_LP_diff"]=LP['C_LP']-LP['CFC_LP']

x=LP["CFC_LP"].tolist()
y=LP["CFC_LP_diff"].tolist()
xy = list(zip(x,y))
new=list(chain.from_iterable(xy))
CFC_LP_df = pd.DataFrame(new,columns =['variant_count'])

freqlist=round((kLP['CFC_LP']/kLP['C_LP']*100),1).tolist()
frqlist = [0 if x != x else x for x in freqlist]
frqlist=[int(i) for i in frqlist]
CFC_LP_df["freq"]=frqlist
CFC_LP_df["freq"]="~"+CFC_LP_df["freq"].astype(str) + "%"


mainlist=kLP["CFC_LP"].tolist()
CFC_LP_df["main"]=mainlist

totallist=kLP["C_LP"].tolist()
CFC_LP_df["total"]=totallist

CFC_LP_df["total_label"] = "(" + CFC_LP_df["total"].astype(str) + ")"
CFC_LP_df[["platforms", "orientation"]]=["Clingen_vs_Franklin","down"]

CFC_LP_df.insert(0,"codes",kLP["codes"].tolist())
CFC_LP_df["sub"]=sublist



final_LP = pd.concat([CSS_LP_df,CSC_LP_df,CFF_LP_df,CFC_LP_df], axis=0)
final_LP.to_csv("/Users/islekbro/Desktop/genomize_makale/2_evidence_codes/evidence_outputs/evidence_LP.csv", sep="\t", index=False)



##for common LP

CS_LP_common = LP[["CSS_LP_diff"]].copy().rename(columns={ 'CSS_LP_diff': 'count'})
CS_LP_common["codes"]=sorter
CS_LP_common["platform"]="CS_LP"
CS_LP_common["label"]=CS_LP_common["codes"] + "=(" + CS_LP_common["count"].astype(str) + ")"

CF_LP_common = LP[["CFF_LP_diff"]].copy().rename(columns={ 'CFF_LP_diff': 'count'})
CF_LP_common["codes"]=sorter
CF_LP_common["platform"]="CF_LP"
CF_LP_common["label"]=CF_LP_common["codes"] + "=(" + CF_LP_common["count"].astype(str) + ")"

common_LP=pd.concat([CS_LP_common,CF_LP_common])
common_LP.drop(common_LP.loc[common_LP['count']==0].index, inplace=True)
common_LP.reset_index(drop=True, inplace=True)
common_LP.to_csv("/Users/islekbro/Desktop/genomize_makale/2_evidence_codes/evidence_outputs/common_LP.csv",sep="\t",index=False)

####################################################### FOR VUS #################################################
kVUS=VUS.iloc[np.arange(len(VUS)).repeat(2)]
#for VUS_CSS
VUS["CSS_VUS_diff"]=VUS['S_VUS']-VUS['CSS_VUS']

x=VUS["CSS_VUS"].tolist()
y=VUS["CSS_VUS_diff"].tolist()
xy = list(zip(x,y))
new=list(chain.from_iterable(xy))
CSS_VUS_df = pd.DataFrame(new,columns =['variant_count'])

freqlist=round((kVUS['CSS_VUS']/kVUS['S_VUS']*100),1).tolist()
frqlist = [0 if x != x else x for x in freqlist]
frqlist=[int(i) for i in frqlist]
CSS_VUS_df["freq"]=frqlist
CSS_VUS_df["freq"]="~"+CSS_VUS_df["freq"].astype(str) + "%"


mainlist=kVUS["CSS_VUS"].tolist()
CSS_VUS_df["main"]=mainlist

totallist=kVUS["S_VUS"].tolist()
CSS_VUS_df["total"]=totallist

CSS_VUS_df["total_label"] = "(" + CSS_VUS_df["total"].astype(str) + ")"
CSS_VUS_df[["platforms", "orientation"]]=["Clingen_vs_SEQ","up"]

CSS_VUS_df.insert(0,"codes",kVUS["codes"].tolist())
CSS_VUS_df["sub"]=sublist


#for VUS_CSC
VUS["CSC_VUS_diff"]=VUS['C_VUS']-VUS['CSC_VUS']

x=VUS["CSC_VUS"].tolist()
y=VUS["CSC_VUS_diff"].tolist()
xy = list(zip(x,y))
new=list(chain.from_iterable(xy))
CSC_VUS_df = pd.DataFrame(new,columns =['variant_count'])

freqlist=round((kVUS['CSC_VUS']/kVUS['C_VUS']*100),1).tolist()
frqlist = [0 if x != x else x for x in freqlist]
frqlist=[int(i) for i in frqlist]
CSC_VUS_df["freq"]=frqlist
CSC_VUS_df["freq"]="~"+CSC_VUS_df["freq"].astype(str) + "%"


mainlist=kVUS["CSC_VUS"].tolist()
CSC_VUS_df["main"]=mainlist

totallist=kVUS["C_VUS"].tolist()
CSC_VUS_df["total"]=totallist

CSC_VUS_df["total_label"] = "(" + CSC_VUS_df["total"].astype(str) + ")"
CSC_VUS_df[["platforms", "orientation"]]=["Clingen_vs_SEQ","down"]

CSC_VUS_df.insert(0,"codes",kVUS["codes"].tolist())
CSC_VUS_df["sub"]=sublist


#for VUS_CFF
VUS["CFF_VUS_diff"]=VUS['F_VUS']-VUS['CFF_VUS']

x=VUS["CFF_VUS"].tolist()
y=VUS["CFF_VUS_diff"].tolist()
xy = list(zip(x,y))
new=list(chain.from_iterable(xy))
CFF_VUS_df = pd.DataFrame(new,columns =['variant_count'])

freqlist=round((kVUS['CFF_VUS']/kVUS['F_VUS']*100),1).tolist()
frqlist = [0 if x != x else x for x in freqlist]
frqlist=[int(i) for i in frqlist]
CFF_VUS_df["freq"]=frqlist
CFF_VUS_df["freq"]="~"+CFF_VUS_df["freq"].astype(str) + "%"


mainlist=kVUS["CFF_VUS"].tolist()
CFF_VUS_df["main"]=mainlist

totallist=kVUS["F_VUS"].tolist()
CFF_VUS_df["total"]=totallist

CFF_VUS_df["total_label"] = "(" + CFF_VUS_df["total"].astype(str) + ")"
CFF_VUS_df[["platforms", "orientation"]]=["Clingen_vs_Franklin","up"]

CFF_VUS_df.insert(0,"codes",kVUS["codes"].tolist())
CFF_VUS_df["sub"]=sublist



#for VUS_CFC
VUS["CFC_VUS_diff"]=VUS['C_VUS']-VUS['CFC_VUS']

x=VUS["CFC_VUS"].tolist()
y=VUS["CFC_VUS_diff"].tolist()
xy = list(zip(x,y))
new=list(chain.from_iterable(xy))
CFC_VUS_df = pd.DataFrame(new,columns =['variant_count'])

freqlist=round((kVUS['CFC_VUS']/kVUS['C_VUS']*100),1).tolist()
frqlist = [0 if x != x else x for x in freqlist]
frqlist=[int(i) for i in frqlist]
CFC_VUS_df["freq"]=frqlist
CFC_VUS_df["freq"]="~"+CFC_VUS_df["freq"].astype(str) + "%"


mainlist=kVUS["CFC_VUS"].tolist()
CFC_VUS_df["main"]=mainlist

totallist=kVUS["C_VUS"].tolist()
CFC_VUS_df["total"]=totallist

CFC_VUS_df["total_label"] = "(" + CFC_VUS_df["total"].astype(str) + ")"
CFC_VUS_df[["platforms", "orientation"]]=["Clingen_vs_Franklin","down"]

CFC_VUS_df.insert(0,"codes",kVUS["codes"].tolist())
CFC_VUS_df["sub"]=sublist



final_VUS = pd.concat([CSS_VUS_df,CSC_VUS_df,CFF_VUS_df,CFC_VUS_df], axis=0)
final_VUS.to_csv("/Users/islekbro/Desktop/genomize_makale/2_evidence_codes/evidence_outputs/evidence_VUS.csv", sep="\t", index=False)


##for common VUS

CS_VUS_common = VUS[["CSS_VUS_diff"]].copy().rename(columns={ 'CSS_VUS_diff': 'count'})
CS_VUS_common["codes"]=sorter
CS_VUS_common["platform"]="CS_VUS"
CS_VUS_common["label"]=CS_VUS_common["codes"] + "=(" + CS_VUS_common["count"].astype(str) + ")"

CF_VUS_common = VUS[["CFF_VUS_diff"]].copy().rename(columns={ 'CFF_VUS_diff': 'count'})
CF_VUS_common["codes"]=sorter
CF_VUS_common["platform"]="CF_VUS"
CF_VUS_common["label"]=CF_VUS_common["codes"] + "=(" + CF_VUS_common["count"].astype(str) + ")"

common_VUS=pd.concat([CS_VUS_common,CF_VUS_common])
common_VUS.drop(common_VUS.loc[common_VUS['count']==0].index, inplace=True)
common_VUS.reset_index(drop=True, inplace=True)
common_VUS.to_csv("/Users/islekbro/Desktop/genomize_makale/2_evidence_codes/evidence_outputs/common_VUS.csv",sep="\t",index=False)


####################################################### FOR LB #################################################
kLB=LB.iloc[np.arange(len(LB)).repeat(2)]
#for LB_CSS
LB["CSS_LB_diff"]=LB['S_LB']-LB['CSS_LB']

x=LB["CSS_LB"].tolist()
y=LB["CSS_LB_diff"].tolist()
xy = list(zip(x,y))
new=list(chain.from_iterable(xy))
CSS_LB_df = pd.DataFrame(new,columns =['variant_count'])

freqlist=round((kLB['CSS_LB']/kLB['S_LB']*100),1).tolist()
frqlist = [0 if x != x else x for x in freqlist]
frqlist=[int(i) for i in frqlist]
CSS_LB_df["freq"]=frqlist
CSS_LB_df["freq"]="~"+CSS_LB_df["freq"].astype(str) + "%"


mainlist=kLB["CSS_LB"].tolist()
CSS_LB_df["main"]=mainlist

totallist=kLB["S_LB"].tolist()
CSS_LB_df["total"]=totallist

CSS_LB_df["total_label"] = "(" + CSS_LB_df["total"].astype(str) + ")"
CSS_LB_df[["platforms", "orientation"]]=["Clingen_vs_SEQ","up"]

CSS_LB_df.insert(0,"codes",kLB["codes"].tolist())
CSS_LB_df["sub"]=sublist


#for LB_CSC
LB["CSC_LB_diff"]=LB['C_LB']-LB['CSC_LB']

x=LB["CSC_LB"].tolist()
y=LB["CSC_LB_diff"].tolist()
xy = list(zip(x,y))
new=list(chain.from_iterable(xy))
CSC_LB_df = pd.DataFrame(new,columns =['variant_count'])

freqlist=round((kLB['CSC_LB']/kLB['C_LB']*100),1).tolist()
frqlist = [0 if x != x else x for x in freqlist]
frqlist=[int(i) for i in frqlist]
CSC_LB_df["freq"]=frqlist
CSC_LB_df["freq"]="~"+CSC_LB_df["freq"].astype(str) + "%"


mainlist=kLB["CSC_LB"].tolist()
CSC_LB_df["main"]=mainlist

totallist=kLB["C_LB"].tolist()
CSC_LB_df["total"]=totallist

CSC_LB_df["total_label"] = "(" + CSC_LB_df["total"].astype(str) + ")"
CSC_LB_df[["platforms", "orientation"]]=["Clingen_vs_SEQ","down"]

CSC_LB_df.insert(0,"codes",kLB["codes"].tolist())
CSC_LB_df["sub"]=sublist


#for LB_CFF
LB["CFF_LB_diff"]=LB['F_LB']-LB['CFF_LB']

x=LB["CFF_LB"].tolist()
y=LB["CFF_LB_diff"].tolist()
xy = list(zip(x,y))
new=list(chain.from_iterable(xy))
CFF_LB_df = pd.DataFrame(new,columns =['variant_count'])

freqlist=round((kLB['CFF_LB']/kLB['F_LB']*100),1).tolist()
frqlist = [0 if x != x else x for x in freqlist]
frqlist=[int(i) for i in frqlist]
CFF_LB_df["freq"]=frqlist
CFF_LB_df["freq"]="~"+CFF_LB_df["freq"].astype(str) + "%"


mainlist=kLB["CFF_LB"].tolist()
CFF_LB_df["main"]=mainlist

totallist=kLB["F_LB"].tolist()
CFF_LB_df["total"]=totallist

CFF_LB_df["total_label"] = "(" + CFF_LB_df["total"].astype(str) + ")"
CFF_LB_df[["platforms", "orientation"]]=["Clingen_vs_Franklin","up"]

CFF_LB_df.insert(0,"codes",kLB["codes"].tolist())
CFF_LB_df["sub"]=sublist



#for LB_CFC
LB["CFC_LB_diff"]=LB['C_LB']-LB['CFC_LB']

x=LB["CFC_LB"].tolist()
y=LB["CFC_LB_diff"].tolist()
xy = list(zip(x,y))
new=list(chain.from_iterable(xy))
CFC_LB_df = pd.DataFrame(new,columns =['variant_count'])

freqlist=round((kLB['CFC_LB']/kLB['C_LB']*100),1).tolist()
frqlist = [0 if x != x else x for x in freqlist]
frqlist=[int(i) for i in frqlist]
CFC_LB_df["freq"]=frqlist
CFC_LB_df["freq"]="~"+CFC_LB_df["freq"].astype(str) + "%"


mainlist=kLB["CFC_LB"].tolist()
CFC_LB_df["main"]=mainlist

totallist=kLB["C_LB"].tolist()
CFC_LB_df["total"]=totallist

CFC_LB_df["total_label"] = "(" + CFC_LB_df["total"].astype(str) + ")"
CFC_LB_df[["platforms", "orientation"]]=["Clingen_vs_Franklin","down"]

CFC_LB_df.insert(0,"codes",kLB["codes"].tolist())
CFC_LB_df["sub"]=sublist



final_LB = pd.concat([CSS_LB_df,CSC_LB_df,CFF_LB_df,CFC_LB_df], axis=0)
final_LB.to_csv("/Users/islekbro/Desktop/genomize_makale/2_evidence_codes/evidence_outputs/evidence_LB.csv", sep="\t", index=False)



##for common LB

CS_LB_common = LB[["CSS_LB_diff"]].copy().rename(columns={ 'CSS_LB_diff': 'count'})
CS_LB_common["codes"]=sorter
CS_LB_common["platform"]="CS_LB"
CS_LB_common["label"]=CS_LB_common["codes"] + "=(" + CS_LB_common["count"].astype(str) + ")"

CF_LB_common = LB[["CFF_LB_diff"]].copy().rename(columns={ 'CFF_LB_diff': 'count'})
CF_LB_common["codes"]=sorter
CF_LB_common["platform"]="CF_LB"
CF_LB_common["label"]=CF_LB_common["codes"] + "=(" + CF_LB_common["count"].astype(str) + ")"

common_LB=pd.concat([CS_LB_common,CF_LB_common])
common_LB.drop(common_LB.loc[common_LB['count']==0].index, inplace=True)
common_LB.reset_index(drop=True, inplace=True)
common_LB.to_csv("/Users/islekbro/Desktop/genomize_makale/2_evidence_codes/evidence_outputs/common_LB.csv",sep="\t",index=False)


####################################################### FOR B #################################################
kB=B.iloc[np.arange(len(B)).repeat(2)]
#for B_CSS
B["CSS_B_diff"]=B['S_B']-B['CSS_B']

x=B["CSS_B"].tolist()
y=B["CSS_B_diff"].tolist()
xy = list(zip(x,y))
new=list(chain.from_iterable(xy))
CSS_B_df = pd.DataFrame(new,columns =['variant_count'])

freqlist=round((kB['CSS_B']/kB['S_B']*100),1).tolist()
frqlist = [0 if x != x else x for x in freqlist]
frqlist=[int(i) for i in frqlist]
CSS_B_df["freq"]=frqlist
CSS_B_df["freq"]="~"+CSS_B_df["freq"].astype(str) + "%"


mainlist=kB["CSS_B"].tolist()
CSS_B_df["main"]=mainlist

totallist=kB["S_B"].tolist()
CSS_B_df["total"]=totallist

CSS_B_df["total_label"] = "(" + CSS_B_df["total"].astype(str) + ")"
CSS_B_df[["platforms", "orientation"]]=["Clingen_vs_SEQ","up"]

CSS_B_df.insert(0,"codes",kB["codes"].tolist())
CSS_B_df["sub"]=sublist


#for B_CSC
B["CSC_B_diff"]=B['C_B']-B['CSC_B']

x=B["CSC_B"].tolist()
y=B["CSC_B_diff"].tolist()
xy = list(zip(x,y))
new=list(chain.from_iterable(xy))
CSC_B_df = pd.DataFrame(new,columns =['variant_count'])

freqlist=round((kB['CSC_B']/kB['C_B']*100),1).tolist()
frqlist = [0 if x != x else x for x in freqlist]
frqlist=[int(i) for i in frqlist]
CSC_B_df["freq"]=frqlist
CSC_B_df["freq"]="~"+CSC_B_df["freq"].astype(str) + "%"


mainlist=kB["CSC_B"].tolist()
CSC_B_df["main"]=mainlist

totallist=kB["C_B"].tolist()
CSC_B_df["total"]=totallist

CSC_B_df["total_label"] = "(" + CSC_B_df["total"].astype(str) + ")"
CSC_B_df[["platforms", "orientation"]]=["Clingen_vs_SEQ","down"]

CSC_B_df.insert(0,"codes",kB["codes"].tolist())
CSC_B_df["sub"]=sublist


#for B_CFF
B["CFF_B_diff"]=B['F_B']-B['CFF_B']

x=B["CFF_B"].tolist()
y=B["CFF_B_diff"].tolist()
xy = list(zip(x,y))
new=list(chain.from_iterable(xy))
CFF_B_df = pd.DataFrame(new,columns =['variant_count'])

freqlist=round((kB['CFF_B']/kB['F_B']*100),1).tolist()
frqlist = [0 if x != x else x for x in freqlist]
frqlist=[int(i) for i in frqlist]
CFF_B_df["freq"]=frqlist
CFF_B_df["freq"]="~"+CFF_B_df["freq"].astype(str) + "%"


mainlist=kB["CFF_B"].tolist()
CFF_B_df["main"]=mainlist

totallist=kB["F_B"].tolist()
CFF_B_df["total"]=totallist

CFF_B_df["total_label"] = "(" + CFF_B_df["total"].astype(str) + ")"
CFF_B_df[["platforms", "orientation"]]=["Clingen_vs_Franklin","up"]

CFF_B_df.insert(0,"codes",kB["codes"].tolist())
CFF_B_df["sub"]=sublist



#for B_CFC
B["CFC_B_diff"]=B['C_B']-B['CFC_B']

x=B["CFC_B"].tolist()
y=B["CFC_B_diff"].tolist()
xy = list(zip(x,y))
new=list(chain.from_iterable(xy))
CFC_B_df = pd.DataFrame(new,columns =['variant_count'])

freqlist=round((kB['CFC_B']/kB['C_B']*100),1).tolist()
frqlist = [0 if x != x else x for x in freqlist]
frqlist=[int(i) for i in frqlist]
CFC_B_df["freq"]=frqlist
CFC_B_df["freq"]="~"+CFC_B_df["freq"].astype(str) + "%"


mainlist=kB["CFC_B"].tolist()
CFC_B_df["main"]=mainlist

totallist=kB["C_B"].tolist()
CFC_B_df["total"]=totallist

CFC_B_df["total_label"] = "(" + CFC_B_df["total"].astype(str) + ")"
CFC_B_df[["platforms", "orientation"]]=["Clingen_vs_Franklin","down"]

CFC_B_df.insert(0,"codes",kB["codes"].tolist())
CFC_B_df["sub"]=sublist



final_B = pd.concat([CSS_B_df,CSC_B_df,CFF_B_df,CFC_B_df], axis=0)
final_B.to_csv("/Users/islekbro/Desktop/genomize_makale/2_evidence_codes/evidence_outputs/evidence_B.csv", sep="\t", index=False)


##for common B

CS_B_common = B[["CSS_B_diff"]].copy().rename(columns={ 'CSS_B_diff': 'count'})
CS_B_common["codes"]=sorter
CS_B_common["platform"]="CS_B"
CS_B_common["label"]=CS_B_common["codes"] + "=(" + CS_B_common["count"].astype(str) + ")"

CF_B_common = B[["CFF_B_diff"]].copy().rename(columns={ 'CFF_B_diff': 'count'})
CF_B_common["codes"]=sorter
CF_B_common["platform"]="CF_B"
CF_B_common["label"]=CF_B_common["codes"] + "=(" + CF_B_common["count"].astype(str) + ")"

common_B=pd.concat([CS_B_common,CF_B_common])
common_B.drop(common_B.loc[common_B['count']==0].index, inplace=True)
common_B.reset_index(drop=True, inplace=True)
common_B.to_csv("/Users/islekbro/Desktop/genomize_makale/2_evidence_codes/evidence_outputs/common_B.csv",sep="\t",index=False)
