
from collections import defaultdict
from yaml import load

files, = glob_wildcards("{sample}/outs/filtered_feature_bc_matrix.h5")


rule all:
    input:
        expand("pdf/{sample}.umap.pdf",sample=files)



rule umap:
    input:
        "{sample}/outs/filtered_feature_bc_matrix.h5"
    output:
        "pdf/{sample}.umap.pdf",
        "rds/{sample}.rds"
    threads: 2
    resources:
        mem_mb=5000
    shell:
        "scripts/spatial-pipeline-stage1.R {wildcards.sample}"
