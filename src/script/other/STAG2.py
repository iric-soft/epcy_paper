import os
import h5py as h5
import pandas as pd

hdf = pd.HDFStore('./data/STAG2/SEQC/hdf_counts.h5', mode='r')
counts = hdf.get('counts')
hdf.close()

path_design = os.path.join("./data/design/STAG2/all")
if not os.path.exists(path_design):
    os.makedirs(path_design)


path_matrix = os.path.join("./data/STAG2/SEQC/")
if not os.path.exists(path_matrix):
    os.makedirs(path_matrix)

counts.index.unique('sample_ID')

items = ['XA1', 'XG3']
df_selected = counts.query("sample_ID in @items")

items_drop = [
    'sample_ID', 'community', 'tsne_x',
    'tsne_y', 'tsnedm_x', 'tsnedm_y'
]
df_selected.index = df_selected.index.droplevel(items_drop)

df_design = df_selected.index.to_frame()
df_design = df_design.reset_index(drop=True)
df_design.columns = ["ID", "subgroup"]
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
