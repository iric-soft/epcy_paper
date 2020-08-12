import os
import h5py as h5
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#import umap

hdf = pd.HDFStore('./data/STAG2/SEQC/hdf_counts.h5', mode='r')
counts = hdf.get('counts')
hdf.close()

path_design = os.path.join("./data/design/STAG2_granulo/all")
if not os.path.exists(path_design):
    os.makedirs(path_design)


path_matrix = os.path.join("./data/STAG2_granulo/SEQC/")
if not os.path.exists(path_matrix):
    os.makedirs(path_matrix)

counts.index.unique('sample_ID')


# granulo ko vs wt version
items = [2, 6, 7, 16, 19, 24, 25]
df_selected = counts.query("community in @items")


df_tsne = df_selected.index.to_frame()
df_tsne = df_tsne.reset_index(drop=True)
sns_plot = sns.scatterplot(
    x="tsne_x", y="tsne_y", data=df_tsne,
    hue="STAG2_del", size=1, linewidth=0,
    alpha=0.5)

fig_file = os.path.join(path_matrix, "tsne.pdf")
sns_plot.figure.savefig(fig_file)
plt.close()

items_drop = [
    'sample_ID', 'community', 'tsne_x',
    'tsne_y', 'tsnedm_x', 'tsnedm_y'
]
df_selected.index = df_selected.index.droplevel(items_drop)

df_design = df_selected.index.to_frame()
df_design = df_design.reset_index(drop=True)
df_design.columns = ["sample", "subgroup"]
df_design.to_csv(
    os.path.join(path_design, "design.tsv"),
    sep="\t",
    index=False
)


df_selected.index = df_selected.index.droplevel(['STAG2_del'])
df_matrix = df_selected.transpose()
df_matrix = df_matrix.reset_index()
df_matrix = df_matrix.rename(columns={'index': 'ID'})

df_matrix.to_csv(
    os.path.join(path_matrix, "readcounts.xls"),
    sep="\t",
    index=False
)

#df_tsne_gene = df_selected.reset_index()

#fig = plt.figure(figsize=(14, 7), dpi=80, facecolor='w', edgecolor='k')

#area1 = np.ma.masked_where(r < r0, area)
#area2 = np.ma.masked_where(r >= r0, area)
#plt.scatter(x, y, s=area1, marker='^', c=c)
#plt.scatter(x, y, s=area2, marker='o', c=c)
#plt.scatter(
#    "tsne_x", "tsne_y", data=df_tsne_gene,
#    c="IFITM1", cmap='Reds', marker="STAG2_del")
#plt.colorbar()
#plt.show()


#sns_plot = sns.scatterplot(
#    x="tsne_x", y="tsne_y", data=df_tsne_gene,
#    style="STAG2_del", hue="IFITM1",
#)

#fig_file = os.path.join(path_matrix, "tsne_IFITM1.pdf")
#sns_plot.figure.savefig(fig_file)
#plt.close()


# ko vs wt version
path_design = os.path.join("./data/design/STAG2_ko/all")
if not os.path.exists(path_design):
    os.makedirs(path_design)


path_matrix = os.path.join("./data/STAG2_ko/SEQC/")
if not os.path.exists(path_matrix):
    os.makedirs(path_matrix)


items = ['XA1', 'XG3']
df_selected = counts.query("sample_ID in @items")

items_drop = [
     'sample_ID', 'community', 'tsne_x',
     'tsne_y', 'tsnedm_x', 'tsnedm_y'
]
df_selected.index = df_selected.index.droplevel(items_drop)

df_design = df_selected.index.to_frame()
df_design = df_design.reset_index(drop=True)
df_design.columns = ["sample", "subgroup"]
df_design.to_csv(
    os.path.join(path_design, "design.tsv"),
    sep="\t",
    index=False
)


df_selected.index = df_selected.index.droplevel(['STAG2_del'])
df_matrix = df_selected.transpose()
df_matrix = df_matrix.reset_index()
df_matrix = df_matrix.rename(columns={'index': 'ID'})

df_matrix.to_csv(
    os.path.join(path_matrix, "readcounts.xls"),
    sep="\t",
    index=False
)


# all ko vs wt version
path_design = os.path.join("./data/design/STAG2_ko_all/all")
if not os.path.exists(path_design):
    os.makedirs(path_design)


path_matrix = os.path.join("./data/STAG2_ko_all/SEQC/")
if not os.path.exists(path_matrix):
    os.makedirs(path_matrix)

df_selected = counts

items_drop = [
     'sample_ID', 'community', 'tsne_x',
     'tsne_y', 'tsnedm_x', 'tsnedm_y'
]
df_selected.index = df_selected.index.droplevel(items_drop)


df_design = df_selected.index.to_frame()
df_design = df_design.reset_index(drop=True)
df_design.columns = ["sample", "subgroup"]
df_design.to_csv(
    os.path.join(path_design, "design.tsv"),
    sep="\t",
    index=False
)

df_selected.index = df_selected.index.droplevel(['STAG2_del'])
df_matrix = df_selected.transpose()
df_matrix = df_matrix.reset_index()
df_matrix = df_matrix.rename(columns={'index': 'ID'})

df_matrix.to_csv(
    os.path.join(path_matrix, "readcounts.xls"),
    sep="\t",
    index=False
)
