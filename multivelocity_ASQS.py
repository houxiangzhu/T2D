#1 Process the RNA data
import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd

pd.set_option('display.max_columns', None)

# get current directory
os. getcwd()

# load sparse matrix:
X = io.mmread("counts.mtx")

# create anndata object
adata = anndata.AnnData(
    X=X.transpose().tocsr()
)

# load cell metadata:
cell_meta = pd.read_csv("metadata.csv")

# load gene names:
with open("gene_names.csv", 'r') as f:
    gene_names = f.read().splitlines()

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names

# load dimensional reduction:
pca = pd.read_csv("pca.csv")
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

# plot a UMAP colored by sampleID to test:
sc.pl.umap(adata, color=['celltype_final'], frameon=False, save=True)

# save dataset as anndata format
adata.write('RNA_data.h5ad')

# Load spliced and unspliced counts matrices
import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad
import os

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2

# reload dataset if needed
adata = sc.read_h5ad('RNA_data.h5ad')

loom_files = os.listdir("looms")

# create an empty anndata object
ldatas = ad.AnnData()

flag = 1

for loom in loom_files:
    # load loom files for spliced/unspliced matrices
    ldata_path = "looms/"+ loom
    ldata = scv.read(ldata_path)

    # rename barcodes in order to merge:
    barcodes = [bc.replace(":", "_") for bc in ldata.obs.index.tolist()]
    barcodes = [bc.replace("x", "-1") for bc in barcodes]
    ldata.obs.index = barcodes

    # make variable names unique
    ldata.var_names_make_unique()

    # concatenate the looms
    if flag == 1:
        ldatas = ldata
    else:
        ldatas = ldatas.concatenate(ldata, index_unique=None)
    flag = flag + 1

# merge matrices into the original adata object
adata = scv.utils.merge(adata, ldatas)

# plot umap to check
sc.pl.umap(adata, color='Majority_RNA10', frameon=False, title='', save='_celltypes.pdf')

# save dataset as anndata format
adata.write('RNA_data_with_looms.h5ad')

#2 Process the ATAC data
import os
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import matplotlib.pyplot as plt
import anndata as ad

scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo')
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_rows', 200)
np.set_printoptions(suppress=True)

samples = os.listdir("cellranger_output")

# create an empty anndata object
adata_atac_merged = ad.AnnData()

flag = 1

for sample in samples:
    mat_path = 'cellranger_output/'+sample+'/filtered_feature_bc_matrix/'
    adata_atac = sc.read_10x_mtx(mat_path, var_names='gene_symbols', cache=False, gex_only=False)
    adata_atac = adata_atac[:,adata_atac.var['feature_types'] == "Peaks"]

    # We aggregate peaks around each gene as well as those that have high correlations with promoter peak or gene expression.
    # Peak annotation contains the metadata for all peaks.
    # Feature linkage contains pairs of correlated genomic features.
    anno_path = 'cellranger_output/'+sample+'/atac_peak_annotation.tsv'
    link_path = 'cellranger_output/'+sample+'/feature_linkage.bedpe'
    adata_atac = mv.aggregate_peaks_10x(adata_atac,
                                        anno_path,
                                        link_path,
                                        verbose=True)

    # We normalize aggregated peaks with TF-IDF.
    mv.tfidf_norm(adata_atac)

    # modify cell barcodes
    adata_atac.obs.index = sample + "_" + adata_atac.obs.index

    # Merge multiple adata
    if flag == 1:
        adata_atac_merged = adata_atac
    else:
        adata_atac_merged = adata_atac_merged.concatenate(adata_atac, index_unique=None)
    flag = flag + 1

# save dataset as anndata format
adata_atac_merged.write('adata_atac_merged.h5ad')

#################################################################
#3 Multivelocity analysis
#################################################################
import os
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import matplotlib.pyplot as plt
import anndata as ad
import subprocess
import seaborn as sns
from statannotations.Annotator import Annotator

scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo')
scv.settings.plot_prefix = ""
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_rows', 200)
np.set_printoptions(suppress=True)

#celltypes = ['Beta1', 'Beta2', 'Beta3', 'Beta4', 'Beta5', 'Beta6']
celltypes = ['Activated stellate', 'Quiescent stellate']

phenos = ['LN', 'OW', 'OB', 'T2D']

#xlims = [-6, 5.5]
#ylims = [-4, 7]
xlims = [1.5, 5]
ylims = [11, 13]

for pheno in phenos:
    # 3.1 reload RNA data
    adata_rna = sc.read_h5ad('RNA_data_with_looms.h5ad')

    cur_pheno = [pheno]

    adata_rna = adata_rna[adata_rna.obs['sample_group'].isin(cur_pheno)]

    adata_rna = adata_rna[adata_rna.obs['celltype_final'].isin(celltypes)]

    # Top 1000 variable genes are used for downstream analyses.
    scv.pp.filter_and_normalize(adata_rna, min_shared_counts=10, n_top_genes=1000)

    # 3.2 reload ATAC data
    adata_atac = sc.read_h5ad('adata_atac_merged.h5ad')

    # 3.3 Finding shared barcodes and features between RNA and ATAC
    shared_cells = pd.Index(np.intersect1d(adata_rna.obs_names, adata_atac.obs_names))
    shared_genes = pd.Index(np.intersect1d(adata_rna.var_names, adata_atac.var_names))

    adata_rna = sc.read_h5ad('RNA_data_with_looms.h5ad')
    adata_rna = adata_rna[shared_cells, shared_genes]
    adata_atac = adata_atac[shared_cells, shared_genes]

    scv.pp.normalize_per_cell(adata_rna)
    scv.pp.log1p(adata_rna)
    scv.pp.moments(adata_rna, n_pcs=30, n_neighbors=50)

    # 3.4 Smoothing gene aggregagted peaks by neighbors
    # Write out filtered cells and prepare to run Seurat WNN --> R script can be found on Github.
    adata_rna.obs_names.to_frame().to_csv('filtered_cells.txt', header=False, index=False)

    # call R to generate the following .txt files
    subprocess.call(
        "/usr/local/bin/Rscript /Users/biocore/Documents/Projects/scATACseq/CRI-BIO-842/RNA_Velocity/MultiVelo_Seurat_WNN.R",
        shell=True)

    # Read in Seurat WNN neighbors.
    nn_idx = np.loadtxt("nn_idx.txt", delimiter=',')
    nn_dist = np.loadtxt("nn_dist.txt", delimiter=',')
    nn_cells = pd.Index(pd.read_csv("nn_cells.txt", header=None)[0])

    # Make sure cell names match.
    np.all(nn_cells == adata_atac.obs_names)

    mv.knn_smooth_chrom(adata_atac, nn_idx, nn_dist)

    ### Run sevelo and multivelo models
    #3.5 Running multi-omic dynamical model
    # This will take a while. Parallelization is high recommended.

    # RNA velo
    adata_result = mv.recover_dynamics_chrom(adata_rna,
                                             adata_atac,
                                             max_iter=5,
                                             init_mode="invert",
                                             verbose=False,
                                             parallel=True,
                                             save_plot=False,
                                             rna_only=True,
                                             fit=True,
                                             n_anchors=500
                                             )

    # 3.6 Computing velocity stream and latent time
    mv.velocity_graph(adata_result)
    mv.latent_time(adata_result)

    stream_name = 'RNAVelo_' + '_'.join(celltypes) + '_' + '_'.join(cur_pheno) + '_embedding_stream.svg'
    mv.velocity_embedding_stream(adata_result,
                                 basis='umap',
                                 legend_loc='none',
                                 size=100,
                                 linewidth=2,
                                 arrow_size=3,
                                 alpha=1,
                                 color='celltype_final',
                                 title='',
                                 save=stream_name,
                                 xlim=xlims,
                                 ylim=ylims,
                                 palette={
                                     "Activated stellate": "#d62728",
                                     "Quiescent stellate": "#2ca02c"
                                 }
                                 # palette={
                                 #    "Beta1": "#1f77b4",
                                 #    "Beta2": "#ff7f0e",
                                 #    "Beta3": "#2ca02c",
                                 #    "Beta4": "#d62728",
                                 #    "Beta5": "#9467bd",
                                 #    "Beta6": "#8c564b"
                                 # }
                                 )

    latent_name = 'RNAVelo_' + '_'.join(celltypes) + '_' + '_'.join(cur_pheno) + '_latent_time.svg'
    scv.pl.scatter(adata_result,
                   basis='umap',
                   color='latent_time',
                   color_map='gnuplot',
                   size=100,
                   title='',
                   save=latent_name,
                   xlim=xlims,
                   ylim=ylims
                   )



    # 3.7 latent time swarmplot plot
    # sns.histplot(data=adata_result.obs[adata_result.obs["celltype_final"]=="Activated stellate"], x="latent_time")
    # sns.histplot(data=adata_result.obs[adata_result.obs["celltype_final"] == "Quiescent stellate"], x="latent_time")

    adata_result = adata_result[(adata_result.obs['UMAP_1'] >= 1.5) & (adata_result.obs['UMAP_1'] <= 5) & (adata_result.obs['UMAP_2'] >= 11) & (adata_result.obs['UMAP_2'] <= 13)]

    plt.figure(figsize=(10, 7), dpi=100)

    order = ['Quiescent stellate', 'Activated stellate']

    ax = sns.swarmplot(data=adata_result.obs,
                       x='latent_time',
                       y='celltype_final',
                       order=order,
                       # cut=0,
                       orient='h',
                       zorder=1,
                       palette={
                           "Activated stellate": "#d62728",
                           "Quiescent stellate": "#2ca02c"
                       }
                       )

    # ax.set_title("Life Expectancy By Country")
    # ax.set_ylabel("Cell Type", fontsize=15)
    # ax.set_xlabel("Latent Time", fontsize=15)
    # ax.tick_params(axis='x', labelsize=20)
    # ax.tick_params(axis='y', labelsize=15)
    # ax.set_yticklabels(ax.get_yticks(), size=20)

    pairs = [("Activated stellate", "Quiescent stellate")]

    annotator = Annotator(ax, pairs, data=adata_result.obs, x='latent_time', y='celltype_final',
                          order=order,
                          # cut=0,
                          orient='h')
    annotator.configure(test='Mann-Whitney', text_format='full', loc='inside')
    annotator.apply_and_annotate()  # ADD here
    # annotator.apply_test().annotate(line_offset_to_group=0.2, line_offset=0.1)

    ax = sns.pointplot(ax=ax,
                       data=adata_result.obs.groupby('celltype_final', as_index=False).median(),
                       x='latent_time',
                       y='celltype_final',
                       order=order,
                       orient='h',
                       join=False,
                       zorder=100
                       )

    ax.set_ylabel("Cell Type", fontsize=15)
    ax.set_xlabel("Latent Time", fontsize=15)

    plt.savefig('figures/RNAVelo_' + '_'.join(celltypes) + '_' + '_'.join(cur_pheno) + '_latent_time_stat.svg')
    plt.clf()

    adata_result.obs.groupby('celltype_final', as_index=False).median().to_csv('RNAVelo_' + '_'.join(celltypes) + '_' + '_'.join(cur_pheno) + '.csv')

    #multivelo
    adata_result = mv.recover_dynamics_chrom(adata_rna,
                                             adata_atac,
                                             max_iter=5,
                                             init_mode="invert",
                                             verbose=False,
                                             parallel=True,
                                             save_plot=False,
                                             rna_only=False,
                                             fit=True,
                                             n_anchors=500
                                            )

    #3.6 Computing velocity stream and latent time
    mv.velocity_graph(adata_result)
    mv.latent_time(adata_result)

    stream_name = 'MultiVelo_' + '_'.join(celltypes) + '_' + '_'.join(cur_pheno) + '_embedding_stream.svg'
    mv.velocity_embedding_stream(adata_result,
                                 basis='umap',
                                 legend_loc='none',
                                 size=100,
                                 linewidth=2,
                                 arrow_size=3,
                                 alpha=1,
                                 color='celltype_final',
                                 title='',
                                 save=stream_name,
                                 xlim=xlims,
                                 ylim=ylims,
                                 palette={
                                     "Activated stellate": "#d62728",
                                     "Quiescent stellate": "#2ca02c"
                                 }
                                 #palette={
                                 #    "Beta1": "#1f77b4",
                                 #    "Beta2": "#ff7f0e",
                                 #    "Beta3": "#2ca02c",
                                 #    "Beta4": "#d62728",
                                 #    "Beta5": "#9467bd",
                                 #    "Beta6": "#8c564b"
                                 #}
                                 )

    latent_name = 'MultiVelo_' + '_'.join(celltypes) + '_' + '_'.join(cur_pheno) + '_latent_time.svg'
    scv.pl.scatter(adata_result,
                   basis='umap',
                   color='latent_time',
                   color_map='gnuplot',
                   size=100,
                   title='',
                   save=latent_name,
                   xlim=xlims,
                   ylim=ylims
                   )

    # 3.7 latent time swarmplot plot
    # sns.histplot(data=adata_result.obs[adata_result.obs["celltype_final"]=="Activated stellate"], x="latent_time")
    # sns.histplot(data=adata_result.obs[adata_result.obs["celltype_final"] == "Quiescent stellate"], x="latent_time")

    adata_result = adata_result[(adata_result.obs['UMAP_1'] >= 1.5) & (adata_result.obs['UMAP_1'] <= 5) & (adata_result.obs['UMAP_2'] >= 11) & (adata_result.obs['UMAP_2'] <= 13)]

    plt.figure(figsize=(10, 7), dpi=100)

    order = ['Quiescent stellate', 'Activated stellate']

    ax = sns.swarmplot(data=adata_result.obs,
                       x='latent_time',
                       y='celltype_final',
                       order=order,
                       # cut=0,
                       orient='h',
                       zorder=1,
                       palette={
                           "Activated stellate": "#d62728",
                           "Quiescent stellate": "#2ca02c"
                       }
                       )

    # ax.set_title("Life Expectancy By Country")
    #ax.set_ylabel("Cell Type", fontsize=15)
    #ax.set_xlabel("Latent Time", fontsize=15)
    # ax.tick_params(axis='x', labelsize=20)
    # ax.tick_params(axis='y', labelsize=15)
    # ax.set_yticklabels(ax.get_yticks(), size=20)

    pairs = [("Activated stellate", "Quiescent stellate")]

    annotator = Annotator(ax, pairs, data=adata_result.obs, x='latent_time', y='celltype_final',
                          order=order,
                          # cut=0,
                          orient='h')
    annotator.configure(test='Mann-Whitney', text_format='full', loc='inside')
    annotator.apply_and_annotate()  # ADD here
    # annotator.apply_test().annotate(line_offset_to_group=0.2, line_offset=0.1)

    ax = sns.pointplot(ax=ax,
                       data=adata_result.obs.groupby('celltype_final', as_index=False).median(),
                       x='latent_time',
                       y='celltype_final',
                       order=order,
                       orient='h',
                       join=False,
                       zorder=100
                       )

    ax.set_ylabel("Cell Type", fontsize=15)
    ax.set_xlabel("Latent Time", fontsize=15)

    plt.savefig('figures/MultiVelo_' + '_'.join(celltypes) + '_' + '_'.join(cur_pheno) + '_latent_time_stat.svg')
    plt.clf()

    adata_result.obs.groupby('celltype_final', as_index=False).median().to_csv('MultiVelo_' + '_'.join(celltypes) + '_' + '_'.join(cur_pheno) + '.csv')