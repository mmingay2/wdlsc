import numpy as np
import pandas as pd
import scanpy as sc
import sys
import os


# get arguments
infile = str(sys.argv[1])

# get output file strings
bfile = os.path.basename(infile)
bname = bfile.split(".")[0]
results_file, umap_plot, generank_plot = bname + "_results.h5ad", "_" + bname + "_umap.png", "_" + bname + "_gene_rank.png"


def load_data(infile):
	# load file
	adata = sc.read_10x_h5(infile)
	# make gene names unique
	adata.var_names_make_unique()
	return(adata)
	
def filter_data(adata, xmin_genes=500, xmax_genes=6000, cell_pct=0.05, mito_thresh=5):
	# get cells with more than 500 genes expressed
	sc.pp.filter_cells(adata, min_genes=xmin_genes)
	# get cells with less than 6000 genes expressed
	sc.pp.filter_cells(adata, max_genes=xmax_genes)
	ncells = adata.shape[0]
	sc.pp.filter_genes(adata, min_cells=ncells*cell_pct)
	get_qc_mt(adata)
	# remove cells in which more than 5% of expressed genes are mitochondrial genes
	adata = adata[adata.obs.pct_counts_mt < mito_thresh, :]
	return(adata)
		
def extract_variable_genes(adata, xn_top_genes=2000):
	sc.pp.highly_variable_genes(adata, n_top_genes=2000)
	adata = adata[:, adata.var.highly_variable]#; print(adata.shape)
	return(adata)
	
def normalize_data(adata, xtarget_sum=1e5):
	adata = sc.pp.normalize_total(adata, target_sum=xtarget_sum, copy=True)
	sc.pp.log1p(adata)
	return(adata)

def gene_rank(adata, xn_genes=25):
	# rank genes with Wilcoxin (Mann Whitney) test
	sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
	return(None)

def cluster_data(adata, xn_comps=50, scalemax=10, xnn=100, xn_pcs=50, leiden_res=1.3):
	# regress out total_counts and pct_counts_mt columns
	adata = sc.pp.regress_out(adata.copy(), ['total_counts', 'pct_counts_mt'], copy=True)
	# scale adata
	sc.pp.scale(adata, max_value=scalemax)
	sc.tl.pca(adata, n_comps=xn_comps, svd_solver='arpack')
	# Get Nearest Neighborhood graph of 100 neighbors from the 50 PCs
	sc.pp.neighbors(adata, n_neighbors=xnn, n_pcs=xn_pcs)
	# Leiden clustering on PCA embedding with resolution 1.3.
	######################## on PCA embedding 
	sc.tl.leiden(adata, resolution=leiden_res)
	# Calculate UMAP embedding 
	sc.tl.umap(adata)
	sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
	######################## on PCA embedding 
	return(adata)

def get_qc_mt(adata):
	# annotate mitochondrial genes (containing 'MT-')
	adata.var['mt'] = adata.var_names.str.startswith('MT-')
	# calculate the percentage expressed genes that are mitochondrial genes
	sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
	return(None)
	

#########
if __name__ == "__main__":
	adata = load_data(infile)
	fdata = filter_data(adata, xmin_genes=500, xmax_genes=6000, cell_pct=0.05, mito_thresh=5)
	ndata = normalize_data(fdata, xtarget_sum=1e5)
	xdata = extract_variable_genes(ndata, xn_top_genes=2000)
	cdata = cluster_data(xdata, xn_comps=50, scalemax=10, xnn=100, xn_pcs=50, leiden_res=1.3)
	# save umap plot
	sc.pl.umap(cdata, color=['leiden'], save=umap_plot)
	# save gene rank plot
	sc.pl.rank_genes_groups(cdata, n_genes=25, sharey=False, save=generank_plot)
	# save to specified output file
	cdata.write(results_file)

