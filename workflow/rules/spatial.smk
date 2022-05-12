
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
        "data/{sample}/outs/spatial/tissue_fixed.png"

    shell:
        """
        convert {input} -colorspace HCL -channel R -evaluate set 67% +channel -colorspace sRGB {output}
        """

rule imagetissue:
    input:
        "data/{sample}/outs/spatial/tissue_lowres_image.png"
    output:
        "{sample}/TissueImage/{sample}.png"
    shell:
        """
        convert {input} -colorspace HCL -channel R -evaluate set 67% +channel -colorspace sRGB {output}

        """

rule qc:
    input:
        "rds/{sample}.rds",
        "data/{sample}/outs/spatial/tissue_fixed.png"
    output:
        "{sample}/technicals/{sample}.n_counts.pdf",
        "{sample}/technicals/{sample}.normalization.pdf"
    
    shell:
        """
        workflow/scripts/spatial-qc.R {wildcards.sample}
        """

rule umap:
    input:
        "rds/{sample}.rds",
        "data/{sample}/outs/spatial/tissue_fixed.png"
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

rule clustermarkers:
    input:
        "rds/{sample}.rds"
    output:
        "{sample}/resolution-{res}/{sample}.positive-markers-forAllClusters.xlsx",
        "{sample}/resolution-{res}/{sample}.all-markers-forAllClusters.xlsx"
    shell:
        """
        workflow/scripts/spatial-cluster-markers.R {wildcards.sample} {wildcards.res}
        """

rule clustermarkerplots:
    input:
        "rds/{sample}.rds",
        "data/{sample}/outs/spatial/tissue_fixed.png",
        "{sample}/resolution-{res}/{sample}.positive-markers-forAllClusters.xlsx"
    output:
        directory("{sample}/resolution-{res}/markers")
    shell:
        """
        workflow/scripts/spatial-cluster-markerplots.R {wildcards.sample} {wildcards.res}
        """

rule spatialfeatures:
    input:
        "rds/{sample}.rds",
        "data/{sample}/outs/spatial/tissue_fixed.png"
    output:
        "{sample}/spatial-markers/{sample}.spatial_markers.xlsx",
        directory("{sample}/spatial-markers/plots")
    shell:
        """
        mkdir -p {output[1]}
        workflow/scripts/spatial-spatial-markers.R {wildcards.sample}
        """

rule sc_cluster:
    input:
        "scrna/{datafile}.rds"
    output:
        "scrna/{datafile}.cluster_markers.xlsx"
    shell:
        """
        workflow/scripts/spatial-sc-cluster.R {wildcards.datafile}
        """



rule spotlight:
    input:
        "scrna/{datafile}.rds",
        "rds/{sample}.rds",
        "scrna/{datafile}.cluster_markers.xlsx"

    output:
        "rds_decon/{datafile}/{sample}.rds"
    shell:
        """
        workflow/scripts/spatial-spotlight.R {wildcards.sample} {wildcards.datafile}
        """


rule spotlightpdf:
    input:
        "scrna/{datafile}.rds",
        "rds_decon/{datafile}/{sample}.rds"
    output:
        "{sample}/deconvolution/spotlight/{sample}-{datafile}-spotlight.pdf"
    shell:
        """
        workflow/scripts/spatial-spotlight-pdf.R {wildcards.sample} {wildcards.datafile}
        """

rule seuratdecon:
    input:
        "scrna/{datafile}.rds",
        "rds/{sample}.rds"
    output:
        "{sample}/deconvolution/seurat/{sample}-{datafile}-seurat.pdf"
    shell:
        """
        workflow/scripts/spatial-seurat-decon.R {wildcards.sample} {wildcards.datafile}
        """

rule gbm:
    input:
        "rds/{sample}.rds",
        "models/{modelfile}.rds"

    output:
        "{sample}/deconvolution/gbm/{sample}-{modelfile}-gbm.pdf"
    shell:
        """
        workflow/scripts/spatial-gbmtest.R {wildcards.sample} {wildcards.modelfile}
        """