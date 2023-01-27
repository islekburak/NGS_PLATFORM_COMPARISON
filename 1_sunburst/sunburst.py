######################################

import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

a=pd.read_csv("/Users/islekbro/Desktop/genomize_makale/1_sunburst/sunburst_inputs/comparison.txt",sep="\t")

a["seq-pathogenicity"]=a["seq-pathogenicity"].str.replace("-","")
a["seq-pathogenicity"]=a["seq-pathogenicity"].str.replace("+","")


a=a.rename({'clingen-pathogenicity' : 'clingen'}, axis=1)
a=a.rename({'seq-pathogenicity' : 'seq'}, axis=1)
a=a.rename({'franklin_pathogenicity' : 'franklin'}, axis=1)

a["seq"]=a["seq"].str.replace(" ","")
a["clingen"]=a["clingen"].str.replace(" ","")
a["franklin"]=a["franklin"].str.replace(" ","")


df = a[['clingen', 'seq', 'franklin']].copy()
df = df.replace(r'^\s*$', np.nan, regex=True)
df.fillna("unclassified", inplace = True)



Sdf=df.groupby(["clingen","seq"]).size().reset_index().rename(columns={0:"count"})
Sdf["platform"]="GENOMIZE"


Fdf=df.groupby(["clingen","franklin"]).size().reset_index().rename(columns={0:"count"})
Fdf["platform"]="FRANKLIN"


fig1=px.sunburst(
	Sdf,
	path=["platform","clingen","seq"],
	values="count",
	color="seq",
	color_discrete_sequence=px.colors.qualitative.Pastel,
	template='seaborn',
	title="Pathogenicity Classification of Seq (Compared to Clingen)"
	)


#figure 2
fig2=px.sunburst(
	Fdf,
	path=["platform","clingen","franklin"],
	values="count",
	color="franklin",
	color_discrete_sequence=px.colors.qualitative.Pastel,
	template='ggplot2',
	title="Specific Evidence Codes to Franklin (Compared to Clingen)"
	)


labels1=list(fig1["data"][0]["labels"])
labels2=list(fig2["data"][0]["labels"])
colordict={

"P":"#ff876a",
"LP":"#ffbbac",
"VUS":"#7eaadb",
"LB":"#aae4bd",
"B":"#7fcd83",
"unclassified":"#feffff",

"GENOMIZE":"#bae1ff",
"FRANKLIN":"#ffb3ba"
}

fontcolordict={

"P":"#36454f",
"LP":"#36454f",
"VUS":"#36454f",
"LB":"#36454f",
"B":"#36454f",
"unclassified":"#feffff",

"GENOMIZE":"#36454f",
"FRANKLIN":"#36454f"
}

color_first=[colordict.get(item,item)  for item in labels1]
color_second=[colordict.get(item,item)  for item in labels2]

"""
fontcolor_first=[fontcolordict.get(item,item)  for item in labels1]
fontcolor_second=[fontcolordict.get(item,item)  for item in labels2]
"""

fontcolor_first=len(labels1[:-6])*["blue"]+len(labels1[-6:-1])*["darkgreen"]+["black"]
indices1  = [index for (index, item) in enumerate(labels1) if item == "unclassified"]
for i in indices1:
    fontcolor_first[i] = "#feffff"

fontcolor_second=len(labels2[:-6])*["red"]+len(labels2[-6:-1])*["green"]+["black"]
indices2  = [index for (index, item) in enumerate(labels2) if item == "unclassified"]
for i in indices2:
    fontcolor_second[i] = "#feffff"


fonts1=(len(labels1)-6)*["Arial"]+6*["Arial Black"]
sizes1=len(labels1[:-6])*[15]+len(labels1[-6:-1])*[20]+[25]


fonts2=(len(labels2)-6)*["Arial"]+6*["Arial Black"]
sizes2=len(labels2[:-6])*[15]+len(labels2[-6:-1])*[20]+[25]


fig1.update_traces(
	go.Sunburst(
		hovertemplate='<b>%{label}</b> <br><b>Code Count:</b> %{value}<br><b>Flow:</b> <i>%{id}</i>'
		),
	insidetextorientation='horizontal',
	textfont=dict(family=fonts1, size=sizes1, color=fontcolor_first),
	marker=dict(colors=color_first)
	
	)

fig2.update_traces(
	go.Sunburst(
		hovertemplate='<b>%{label}</b> <br><b>Code Count:</b> %{value}<br><b>Flow:</b> <i>%{id}</i>'
		),
	insidetextorientation='horizontal',
	textfont=dict(family=fonts2, size=sizes2, color=fontcolor_second),
	marker=dict(colors=color_second)
	
	)


#update main figure layout
fig1.update_layout(
	title=dict(
		font_size=25,
		text="GENOMIZE-SEQ",
		font_family="Arial Black",
		font_color="black",
		pad=dict(t=100)
		),
	#grid= dict(columns=2, rows=1),
	hoverlabel=dict(
        bgcolor="whitesmoke",
        font_size=15
        #font_family="Rockwell"
        ),
	hoverlabel_bordercolor="blue",
	hoverlabel_font_color="black",
	hoverlabel_font_family="Arial",
	spikedistance=10,
	uniformtext=dict(minsize=10, mode=False),
	#margin=dict(t=50, b=10, r=100, l=100)
	#paper_bgcolor='whitesmoke'
	)


fig2.update_layout(
	title=dict(
		font_size=25,
		text="FRANKLIN",
		font_family="Arial Black",
		font_color="black",
		pad=dict(t=100)
		),
	#grid= dict(columns=2, rows=1),
	hoverlabel=dict(
        bgcolor="whitesmoke",
        font_size=15
        #font_family="Rockwell"
        ),
	hoverlabel_bordercolor="blue",
	hoverlabel_font_color="black",
	hoverlabel_font_family="Arial",
	spikedistance=10,
	uniformtext=dict(minsize=10, mode=False),
	#margin=dict(t=50, b=10, r=100, l=100),
	legend=dict(bordercolor="black")

	#paper_bgcolor='whitesmoke'
	)

fig1.data[0].marker.line.width = 1.5
fig1.data[0].marker.line.color = "black"
fig2.data[0].marker.line.width = 1.5
fig2.data[0].marker.line.color = "black"

"""
fig1.add_annotation(
    text=f'Azura',
    xanchor='right',
    x=0.41, y=0.21,
    ax=-100, ay=90)


#fig.show()
"""
fig1.write_html("/Users/islekbro/Desktop/genomize_makale/1_sunburst/sunburst_outputs/Genomize_variants.html")
fig2.write_html("/Users/islekbro/Desktop/genomize_makale/1_sunburst/sunburst_outputs/Franklin_variants.html")
