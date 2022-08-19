#!/usr/bin/env python

import scanpy as sc
import tangram as tg
import pandas as pd
import sys
import matplotlib.pyplot as plt
import torch



print(sys.argv[1])
adata_st = sc.read_visium(path=sys.argv[1])
adata_sc= sc.read(filename=sys.argv[2])

adata_sc.X=adata_sc.raw.X.copy()


tg.pp_adatas(adata_sc,adata_st,genes=None)


if torch.cuda.is_available():
       ad_map = tg.map_cells_to_space(
                   adata_sc, 
                   adata_st,device="cuda:0")
else:
       ad_map = tg.map_cells_to_space(
                   adata_sc, 
                   adata_st)


ad_ge = tg.project_genes(ad_map, adata_sc)

genes=["LYVE1","F13A1","FOLR2","SELENOP","APOE","SLC40A1","C1QB","DAB2","PDK4","SPP1","ACP5","CD9","FCER1A","CD1C","CLEC10A","HSPA6","DNAJB1",
       "HSPA1B","S100A8","S100A9","S100A12","EREG","G0S2","FCN1","CCL20","IL1B","IL23A","CXCL10","CXCL9","GBP1","CDC1C","CCL3L1",
       "CCL3","CCL4L2","MT1X","MT1E","CTSL","RGS1","FOS","DNASE1L3","MMP9","LYZ","AREG","VCAN","HSPA1A","MARCO","COLEC12"]

genes=list(set(list(ad_ge.var["features"].to_numpy())) & set(genes))
genes=list(map(lambda x: x.lower(),genes))


expid=list(adata_st.uns['spatial'])[0]

scale_fac=adata_st.uns['spatial'][expid]['scalefactors']['tissue_hires_scalef']

pl=tg.plot_genes_sc(genes, adata_measured=adata_st, adata_predicted=ad_ge, perc=0.02,scale_factor=scale_fac,spot_size=40,return_figure=True)

plt.rcParams['figure.figsize'] = [15, 15]
plt.rcParams['figure.dpi'] = 150


pl.savefig(sys.argv[3], format="pdf", bbox_inches="tight")
