import pandas as pd
import umap.umap_ as umap
patho = pd.read_csv("patho_per_sample.txt", sep = "\t", index_col = 0)
patho = patho.loc[[not("#" in i and "K" in i) for i  in patho.index], ]
# load sample meta data
meta = pd.read_csv("~/aw27/all_paerug/v2_all_studies/collect_meta/samples_qc_filtered_and_annotated_per_patient.txt", sep = "\t", encoding = 'iso-8859-1', index_col = 11)
meta  = meta.loc[:, []]
patho = meta.join(patho, how = "inner")
fit = umap.UMAP(metric = "jaccard")
umap_transform = fit.fit_transform(patho)
umap_transform_df = pd.DataFrame(umap_transform, index = patho.index, columns = ["umap.1", "umap.2"])
umap_transform_df.to_csv("umap_coords_jaccard_per_patient.txt")
