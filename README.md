# pseudomonas-census-paper

## Genomics
hc_clust.R - use hierarchical clustering of SNP distances between isolate genomes to define genomic clusters.  
subsample_alignment_horizontally.py - subset FASTA alignment to samples of interest.  

## Dating

## Pathoadaptation
get_node_transmission_type.py - annotate nodes/branches in clone trees based on (ancestral) infection type and transmissibility  
mutational_distribution.R - aggregate variant effect annotation across clones and perform mutation burden test.
 
## Ancestral genome analysis
get_clone_representative.py - get ancestral genome representatives based on parsimony ancestral character state reconstruction.  
emapper_to_clones.py - parse emapper output for functional enrichment analysis.  
extract_events_fasta.py - extract FASTA file with nucleotide/protein sequences for events reconstructed.   
extract_fasta.py - extract multi-FASTA file for Panaroo gene ids of interest from gene_data.csv.  
root_tree.py - use outgroup to root global Pseudomonas tree.  
parsimony.py - implementation of parsimony ancestral character-state reconstruction of binary traits for gene presence/absence or indels.  
reduce_graph.py - subset Panaroo final_graph.gml to samples of interest.  
subset_gpa.py/subset_gpa_csv.py - subset Panaroo gene presence/absence (csv/Rtab) according to samples of interest.  
