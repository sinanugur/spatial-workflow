
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
        "results/{sample}/clusteringTree/clusteringTree-{sample}.pdf"
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
        "results/{sample}/TissueImage/{sample}.png"
    shell:
        """
        convert {input} -colorspace HCL -channel R -evaluate set 67% +channel -colorspace sRGB {output}

        """

rule qc:
    input:
        "rds/{sample}.rds",
        "data/{sample}/outs/spatial/tissue_fixed.png"
    output:
        "results/{sample}/technicals/{sample}.n_counts.pdf",
        "results/{sample}/technicals/{sample}.normalization.pdf"
    
    shell:
        """
        workflow/scripts/spatial-qc.R {wildcards.sample}
        """

rule umap:
    input:
        "rds/{sample}.rds",
        "data/{sample}/outs/spatial/tissue_fixed.png"
    output:
        "results/{sample}/resolution-{res}/{sample}.umap.spatial.pdf"
    shell:
        """
        workflow/scripts/spatial-umap.R {wildcards.sample} {wildcards.res}
        """


rule spatialmetrics:
    input:
        "rds/{sample}.rds"
    output:
        "results/{sample}/resolution-{res}/{sample}.number-of-cells-per-cluster.xlsx"
    shell:
        """
        workflow/scripts/spatial-metrics.R {wildcards.sample} {wildcards.res}
        """

rule clustermarkers:
    input:
        "rds/{sample}.rds"
    output:
        "results/{sample}/resolution-{res}/{sample}.positive-markers-forAllClusters.xlsx",
        "results/{sample}/resolution-{res}/{sample}.all-markers-forAllClusters.xlsx"
    shell:
        """
        workflow/scripts/spatial-cluster-markers.R {wildcards.sample} {wildcards.res}
        """

rule clustermarkerplots:
    input:
        "rds/{sample}.rds",
        "data/{sample}/outs/spatial/tissue_fixed.png",
        "results/{sample}/resolution-{res}/{sample}.positive-markers-forAllClusters.xlsx"
    output:
        directory("results/{sample}/resolution-{res}/markers")
    shell:
        """
        workflow/scripts/spatial-cluster-markerplots.R {wildcards.sample} {wildcards.res}
        """

rule spatialfeatures:
    input:
        "rds/{sample}.rds",
        "data/{sample}/outs/spatial/tissue_fixed.png"
    output:
        "results/{sample}/spatial-markers/{sample}.spatial_markers.xlsx",
        directory("results/{sample}/spatial-markers/plots")
    shell:
        """
        mkdir -p {output[1]}
        workflow/scripts/spatial-spatial-markers.R {wildcards.sample}
        """

rule selected_markers:
    input:
        "rds/{sample}.rds",
        "data/{sample}/outs/spatial/tissue_fixed.png"
    output:
        directory("results/{sample}/selected-markers/plots")
    shell:
        """
        mkdir -p {output}
        workflow/scripts/spatial-selected-markers.R {wildcards.sample}
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

rule rctd:
    input:
        "scrna/{datafile}.rds",
        "rds/{sample}.rds"

    output:
        "rds_rctd/{datafile}/{sample}.rds"
    shell:
        """
        workflow/scripts/spatial-rctd.R {wildcards.sample} {wildcards.datafile}
        """

rule rctdpdf:
    input:
        "rds_rctd/{datafile}/{sample}.rds"

    output:
        "results/{sample}/deconvolution/rctd/{sample}-{datafile}-rctd.pdf"
    shell:
        """
        workflow/scripts/spatial-rctd-pdf.R {wildcards.sample} {wildcards.datafile}
        """


rule spotlightpdf:
    input:
        "scrna/{datafile}.rds",
        "rds_decon/{datafile}/{sample}.rds",
        "data/{sample}/outs/spatial/tissue_fixed.png"
    output:
        "results/{sample}/deconvolution/spotlight/{sample}-{datafile}-spotlight.pdf"
    shell:
        """
        workflow/scripts/spatial-spotlight-pdf.R {wildcards.sample} {wildcards.datafile}
        """

rule seuratdecon:
    input:
        "scrna/{datafile}.rds",
        "rds/{sample}.rds",
        "data/{sample}/outs/spatial/tissue_fixed.png"
    output:
        "results/{sample}/deconvolution/seurat/{sample}-{datafile}-seurat.pdf"
    shell:
        """
        workflow/scripts/spatial-seurat-decon.R {wildcards.sample} {wildcards.datafile}
        """

rule gbm:
    input:
        "rds/{sample}.rds",
        "models/{modelfile}.rds",
        "data/{sample}/outs/spatial/tissue_fixed.png"

    output:
        "results/{sample}/deconvolution/gbm/{sample}-{modelfile}-gbm.pdf"
    shell:
        """
        workflow/scripts/spatial-gbmtest.R {wildcards.sample} {wildcards.modelfile}
        """