    
workflow scanpymm {
    Map[String, String] data_map

    scatter (pair in data_map) {
        call run_scanpy {
            input: infile=pair.right,
            dname=pair.left
        }
    }
}


task run_scanpy {
    File infile
    String dname
    command {
        ln -s ${infile} ${dname}.h5
        python3 /usr/local/bin/scanpy_processing.py ${dname}.h5
        mv figures/* ./
    }

    runtime {
        docker: "mmingay2/scanpy@sha256:71b8143c80f0e7d8599bdf819fe08be2c3f0f6117d58a1b3b40796d07e808047"
    }

    output {
        File count_matrix_h5ad = "${dname}_results.h5ad"
        File gene_rank_png = "rank_genes_groups_leiden_${dname}_gene_rank.png"
        File umap_png = "umap_${dname}_umap.png"
    }
}