#!/usr/bin/env python

import scanpy as sc
import tangram as tg
import pandas as pd
import sys


print(sys.argv[1])
adata_st = sc.read_visium(path=sys.argv[1])
adata_sc= sc.read(filename=sys.argv[2])

adata_sc.X=adata_sc.raw.X.copy()


tg.pp_adatas(adata_sc,adata_st,genes=None)

ad_map = tg.map_cells_to_space(
                   adata_sc, 
                   adata_st,         
                   mode='clusters',
                   cluster_label='seurat_clusters_tangram')


tg.project_cell_annotations(ad_map, adata_st, annotation="seurat_clusters_tangram")


adata_st.obsm['tangram_ct_pred'].to_csv(sys.argv[3])