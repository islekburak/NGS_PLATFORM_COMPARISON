import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
data = pd.read_csv("/Users/islekbro/Desktop/genomize_makale/3_statistics/statistics_outputs/statistics_variant_pathogenicity_out.csv", sep="\t")
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
).set_axis_labels("Evidence Codes", "F1 Scores").set(xlim=(-1,5),ylim=(0, 1.1), yticks=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]).tick_params(axis='x', bottom=True)

f1.set(title="F1 Score Comparison for Different Platforms \n(Based on Clingen)")
f1.savefig('/Users/islekbro/Desktop/genomize_makale/3_statistics/statistics_outputs/statistics_plot_outs/variant_pathogenicity/f1_variant_pathogenicity_plot.pdf',dpi=300)   # save the figure to file
plt.close()




p=sns.catplot(
    data=data, palette=["#bae1ff", "#ffb3ba"],
    x="index", y="precision",
    hue="platform",
    kind="bar",
    aspect=1.5,
    height=4,
    edgecolor="black",
).set_axis_labels("Evidence Codes", "Precision Scores").set(xlim=(-1,5),ylim=(0, 1.1), yticks=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]).tick_params(axis='x', bottom=True)


p.set(title="Precision Score Comparison for Different Platforms \n(Based on Clingen)")
p.savefig('/Users/islekbro/Desktop/genomize_makale/3_statistics/statistics_outputs/statistics_plot_outs/variant_pathogenicity/precision_variant_pathogenicity_plot.pdf',dpi=300)   # save the figure to file
plt.close()




r=sns.catplot(
    data=data, palette=["#bae1ff", "#ffb3ba"],
    x="index", y="recall",
    hue="platform",
    kind="bar",
    aspect=1.5,
    height=4,
    edgecolor="black",
).set_axis_labels("Evidence Codes", "Recall Scores").set(xlim=(-1,5),ylim=(0, 1.1), yticks=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]).tick_params(axis='x', bottom=True)


r.set(title="Recall Score Comparison for Different Platforms \n(Based on Clingen)")
r.savefig('/Users/islekbro/Desktop/genomize_makale/3_statistics/statistics_outputs/statistics_plot_outs/variant_pathogenicity/recall_variant_pathogenicity_plot.pdf',dpi=300)   # save the figure to file
plt.close()