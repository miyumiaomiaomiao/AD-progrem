import os
import scanpy as sc
import gc
import pandas as pd

data_folder = r"C:\rnaseq folder"
result_data={}

for file_name in os.listdir(data_folder):
    file_path = os.path.join(data_folder, file_name)
    adata = sc.read_h5ad(file_path,backed='r')
    
    if 'Subclass' in adata.obs.columns:
        subclass_counts = adata.obs['Subclass'].value_counts()
        total_cells_in_file = adata.shape[0]   
    percentage_dict = (subclass_counts / total_cells_in_file * 100).to_dict()
    result_data[file_name] = percentage_dict
    
    del adata
    gc.collect()
    
result_df = pd.DataFrame.from_dict(result_data, orient='index')
    
print(result_df)
desktop_path = os.path.join(os.path.expanduser("~"), "Desktop")
output_path = os.path.join(desktop_path, 'subclass_statistics_wide.xlsx')
result_df.to_excel(output_path, index=False)


