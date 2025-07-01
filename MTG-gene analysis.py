import numpy as np
import pandas as pd
import scanpy as sc
import os
from scipy import io
from scipy.sparse import csr_matrix
import shutil
import gzip

input_fold = "C:/78"
h5ad_files = [file for file in os.listdir(input_fold) if file.endswith(".h5ad")]
for file_name in h5ad_files:
    input_file = os.path.join(input_fold,file_name)
    adata = sc.read_h5ad(input_file)

    if adata.raw is not None:
      adata = adata.raw.to_adata()
    cell_type_col = "Subclass"
    selected_cell = adata.obs[cell_type_col] == "L2/3 IT"
    adata = adata[selected_cell].copy()

    if not isinstance(adata.X, csr_matrix):
        adata.X = csr_matrix(adata.X)

    file_name = os.path.basename(input_file)
    file_index = file_name.split("_")[0]
    output_dir = os.path.join("C:/Seurat_fold(MTG)", file_index + "_L2_3_IT")
    os.makedirs(output_dir, exist_ok=True)

    metadata_file = os.path.join(output_dir, "metadata.tsv.gz")
    adata.obs.to_csv(metadata_file, sep="\t", compression="gzip")

    barcodes_file = os.path.join(output_dir, "barcodes.tsv.gz")
    with gzip.open(barcodes_file, "wt") as f:
        for item in adata.obs_names:
            f.write(item + "\n")

    features_file = os.path.join(output_dir, "features.tsv.gz")
    with gzip.open(features_file, "wt") as f:
        for item in adata.var_names:
            f.write(f"{item}\t{item}\n")

    matrix_file = os.path.join(output_dir, "matrix.mtx")
    io.mmwrite(matrix_file, adata.X.T)

    with open(matrix_file, "rb") as f_in:
        with gzip.open(matrix_file + ".gz", "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    os.remove(matrix_file)

    print(f"数据已成功导出到 {output_dir}")



