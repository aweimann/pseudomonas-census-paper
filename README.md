# Pseudomonas Census Paper

## Instructions
The scripts in this repository rely on the data in the associated Zenodo repository (https://doi.org/10.5281/zenodo.10600286) and are described roughly in the intended order of execution.

## Genomics
hc_clust.R - use hierarchical clustering of SNP distances between isolate genomes to define genomic clusters.
subsample_alignment_horizontally.py - subset FASTA alignment to samples of interest.  
transmission.R - infer potential transmission links based on genomic relatedness.    
cf_per_clone.R - get clone CF/non-CF proportions (CF vs non-CF patients) and compare with surveillance data.

## Dating
plotDCCBayesianSkylinePlot.R - based on C. Ruis 20. Aggregates Skyline model population size data across clones and outputs plots of inferred clone historical population sizes.  
modify_beast_xml.py - modifies BEAST XML template generated with BEAUti to prepare input XMLs for BEAST runs for multiple clones.  
root_age.py - from Treeannotator annotated trees, alignment and sample dates work out root age for input tree.    

## Ancestral genome analysis
root_tree.py - use outgroup to root global Pseudomonas tree.  
parsimony.py - implementation of parsimony ancestral character-state reconstruction of binary traits for gene presence/absence or indels.  
reduce_graph.py - subset Panaroo final_graph.gml to samples of interest.  
subset_gpa.py/subset_gpa_csv.py - subset Panaroo gene presence/absence (csv/Rtab) according to samples of interest.  
get_clone_representative.py - get ancestral genome representatives based on parsimony ancestral character state reconstruction.  
extract_events_fasta.py - extract FASTA file with nucleotide/protein sequences for events reconstructed.   
extract_fasta.py - extract multi-FASTA file for Panaroo gene ids of interest from gene_data.csv.  
emapper_to_clones.py - parse emapper output for functional enrichment analysis.  
emapper_stats.R - test functional enrichment of COG categories in epidemic/sporadic clones and along CF proportions.

## Gene expression analysis
cf_association_analysis.R - Perform gene expression CF association analysis.

## Isolate phenotyping analysis
analyse_assay_measurements.R - Plot virulence factor measurements and test association of virulence factor expression with CF proportions.    
plot_abundance_individua.R - Plot macrophage survival across clones with variable CF proportions and test differential survival.

## Pathoadaptation
get_node_transmission_type.py - annotate nodes/branches in clone trees based on (ancestral) infection type and transmissibility (needs correct columns from samples_qc_filtered_and_annotated_v2.txt as inputs).    
mutational_distribution.R - aggregate variant effect annotation across clones and perform mutational burden test.    
manhattan_plot.R - make Manhattan and QQ-plot (run mutational_distribution.R).
process_pathoadaptive_mutations.py - infer sequence of mutations in evolutionary time (only use mutations in mutational burden test hits as input if used for subsequent analyses.).    
plot_gene2patho_position.R - plot gene specific burden over evolutionary time.    
plot_patho_per_sample.py - get UMAP of pathoadaptive space.   
plot_patho_per_sample.R - UMAP scatter plot (run plot_patho_per_sample.py first).    
test_transmissibility_and_cfness.R - test mutational burden test hits for over/underrepresentation of CF vs non-CF, and transmitted vs untransmitted mutations (run mutational_distribution.R and get_node_transmission_type.py first). 
 
