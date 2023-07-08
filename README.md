# Analyze gene-count matrices of human bone marrow single-cell data using SCANPY

Steps from instructions / the Scanpy tutorial are grouped into functions and run sequentially.

See the plots subdirectory for the UMAP and gene rank plots.

# Steps to Run the wdl pipeline

## Clone this repo

```
git clone https://github.com/mmingay2/scwdl.git
cd scwdl
```

## Set up Cromwell / Java / Docker

1. Check if Java is install by running `java --version`. If this does not return a version, install Java environment (e.g. OpenJDK).

2. If you don't already have [Docker](https://docs.docker.com/get-docker/), install it on your computer.

3. Download [Cromwell](https://github.com/broadinstitute/cromwell).

```
wget https://github.com/broadinstitute/cromwell/releases/download/85/cromwell-85.jar
```

## Run the workflow with Cromwell:

### Option 1: Download input data with `wget` in WDL workflow.

Run `scanpy_dl.wdl` to download the necessary data before analyzing it.

```
java -jar cromwell-85.jar run scanpy_dl.wdl -i inputs/scanpy_inputs_dl.json
```

### Option 2: Download inputs before running WDL workflow

#### Download Data

```
mkdir data
for i in {1..8}; do wget "https://storage.googleapis.com/terra-featured-workspaces/Cumulus/cellranger_output/MantonBM"${i}"_HiSeq_1/raw_feature_bc_matrix.h5" -O "data/raw_feature_bc_matrix_${i}.h5"; done
```

#### Run Workflow

Run `scanpy_img.wdl` and include`-Dconfig.file=custom.conf` to use local data after downloading the data.

```
java -Dconfig.file=custom.conf -jar cromwell-85.jar run scanpy_img.wdl -i inputs/scanpy_inputs_img.json
```

> Note the scanpy python script `scanpy_processing.py` is located and run from inside the docker image in both of the above WDL workflows.

