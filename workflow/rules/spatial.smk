
from collections import defaultdict
from yaml import load

files, = glob_wildcards("data/{sample}/outs/filtered_feature_bc_matrix.h5")



rule rds:
    input:
        "data/{sample}/outs/filtered_feature_bc_matrix.h5"
    output:
        protected("rds/{sample}.rds")
    threads: 2
    resources:
        mem_mb=2500
    shell:
        "workflow/scripts/spatial-rds.R {wildcards.sample}"


rule clustree:
    input:
        "rds/{sample}.rds"
    output:
        "{sample}/clusteringTree/clusteringTree-{sample}.pdf"
    shell:
        "workflow/scripts/clusteringTree.R {wildcards.sample}"

rule imagefix:
    input:
        "data/{sample}/outs/spatial/tissue_lowres_image.png"
    output:
        "{sample}/TissueImage/{sample}.png",
        "data/{sample}/outs/spatial/grayscale.png"

    shell:
        """
        convert {input} -colorspace HSI -channel blue -separate +channel {output[0]}
        convert {input} -colorspace HSI -channel blue -separate +channel {output[1]}
        """

rule qc:
    input:
        "rds/{sample}.rds",
        "data/{sample}/outs/spatial/grayscale.png"
    output:
        "{sample}/technicals/{sample}.n_counts.pdf"
    
    shell:
        """
        workflow/scripts/spatial-qc.R {wildcards.sample}
        """

rule umap:
    input:
        "rds/{sample}.rds",
        "data/{sample}/outs/spatial/grayscale.png"
    output:
        "{sample}/resolution-{res}/{sample}.umap.spatial.pdf"
    run:
        shell("workflow/scripts/spatial-umap.R {wildcards.sample} {wildcards.res}")


rule spatialmetrics:
    input:
        "rds/{sample}.rds"
    output:
        "{sample}/resolution-{res}/{sample}.number-of-cells-per-cluster.xlsx"
    shell:
        """
        workflow/scripts/spatial-metrics.R {wildcards.sample} {wildcards.res}
        """

rule spatialfeatures:
    input:
        "rds/{sample}.rds"
    output:
        "{sample}/{sample}.spatial_markers.xlsx",
        directory("{sample}/spatial-markers")
    shell:
        """
        mkdir -p {output[1]}
        workflow/scripts/spatial-markers.R {wildcards.sample}
        """
    



