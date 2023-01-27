import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
data = pd.read_csv("/Users/islekbro/Desktop/genomize_makale/3_statistics/statistics_outputs/statistics_evidence_code_total_out.csv", sep="\t")
sns.set_style("darkgrid", {"grid.color": ".6", "grid.linestyle": ":"})



#for F1
f1=sns.catplot(
    data=data, palette=["#bae1ff", "#ffb3ba"],
    x="index", y="f1-score",
    hue="platform",
    kind="bar",
    aspect=1.5,
    height=4,
    edgecolor="black",
).set_axis_labels("Evidence Codes", "F1 Scores").set(xlim=(-1,19),ylim=(0, 1.1), yticks=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]).tick_params(axis='x', rotation=90, bottom=True)


f1.set(title="Total Evidence Codes F1 Score Comparison \nfor Different Platforms (Based on Clingen)")
f1.savefig('/Users/islekbro/Desktop/genomize_makale/3_statistics/statistics_outputs/statistics_plot_outs/total_evidence_codes/f1_evidence_total_code_plot.pdf',dpi=300)   # save the figure to file
plt.close()



#for Precision
p=sns.catplot(
    data=data, palette=["#bae1ff", "#ffb3ba"],
    x="index", y="f1-score",
    hue="platform",
    kind="bar",
    aspect=1.5,
    height=4,
    edgecolor="black",
).set_axis_labels("Evidence Codes", "Precision Scores").set(xlim=(-1,19),ylim=(0, 1.1), yticks=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]).tick_params(axis='x', rotation=90, bottom=True)


p.set(title="Total Evidence Codes Precision Score Comparison \nfor Different Platforms (Based on Clingen)")
p.savefig('/Users/islekbro/Desktop/genomize_makale/3_statistics/statistics_outputs/statistics_plot_outs/total_evidence_codes/precision_evidence_total_code_plot.pdf',dpi=300)   # save the figure to file
plt.close()



#for Recall
r=sns.catplot(
    data=data, palette=["#bae1ff", "#ffb3ba"],
    x="index", y="f1-score",
    hue="platform",
    kind="bar",
    aspect=1.5,
    height=4,
    edgecolor="black",
).set_axis_labels("Evidence Codes", "Recall Scores").set(xlim=(-1,19),ylim=(0, 1.1), yticks=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]).tick_params(axis='x', rotation=90, bottom=True)


r.set(title="Total Evidence Codes Recall Score Comparison \nfor Different Platforms (Based on Clingen)")
r.savefig('/Users/islekbro/Desktop/genomize_makale/3_statistics/statistics_outputs/statistics_plot_outs/total_evidence_codes/recall_evidence_total_code_plot.pdf',dpi=300)   # save the figure to file
plt.close()
