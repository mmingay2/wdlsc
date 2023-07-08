    
workflow scanpymm {
    Map[String, String] data_map

    scatter (pair in data_map) {
        call get_files {
                input: dlink=pair.right,
                dname=pair.left
            }
        call run_scanpy {
            input: infile=get_files.h5file,
            dname=pair.left
        }
    }
}



task get_files {
    # download link
    String dlink
    # file basename
    String dname

    command {
        # download file
        wget ${dlink} -O ${dname}.h5
    }

    runtime {
        docker: "mmingay2/scanpy@sha256:71b8143c80f0e7d8599bdf819fe08be2c3f0f6117d58a1b3b40796d07e808047"
    }

    output {
        File h5file = "${dname}.h5"
    }
}


task run_scanpy {
    # input file name
    File infile
    # 
    String dname
    command {
        python3 /usr/local/bin/scanpy_processing.py ${infile}
        mv figures/* ./
    }

    runtime {
        docker: "mmingay2/scanpy@sha256:71b8143c80f0e7d8599bdf819fe08be2c3f0f6117d58a1b3b40796d07e808047"
    }

    output {
        File count_matrix_h5ad  = "${dname}_results.h5ad"
        File umap_png = "rank_genes_groups_leiden_${dname}_gene_rank.png"
        File gene_rank_png = "umap_${dname}_umap.png"
    }
}