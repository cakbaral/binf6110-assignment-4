# binf6110-assignment-4

## Introduction

## Methods
All clusters were then annotated and 

## Results
### Dimensionality Reduction and Initial Clustering of Mouse Cells
Dimensionality reduction of mouse cell data from the Seurat object was carried out using PC analysis, to identify relevant PCs to retain for clustering.  UMAP embedding was then used to visualize all of the distinct clusters and their spatial separation with the applied PCs.  When looking at the elbow plot in Figure 1 after applying log-normalization and scaling, a total 20 PCs were identified in total from the data.  However, the variance in cells appeared to decrease at around 15 PCs, as evident by the elbow point levelling off at around PC15.  As such, up to 15 PCs were applied to further downstream analysis for UMAP embedding.

<img width="1482" height="1007" alt="Mouse_Elbow_Plot" src="https://github.com/user-attachments/assets/871c79d6-5f04-4681-ae07-af5578cfbdcb" />

**Figure 1**: Elbow plot of principle components (PCs) for 156,572 cells collected from normal mice and mice infected with nasal influenza.  A total of 20 PCs were identified after scaling and normalizing of data.  The elbow point appeared to occur at around 15 PCs.

Looking at the initial UMAP for the mouse cells (Figure 2), a total of 37 distinct clusters were identified.  These clusters appeared mostly heterogenous from each other, with certain clusters appearing to form larger groups than others, while still remaining distinct between each other.  Overlap between clusters did appear to be present within the clusters upon closer inspection, but these instances were usually minor and, for the most part, did not have a major impact on heterogeneity.  Cluster sizes were also varied, though certain groups appeared to have similar sizes between each other, only differing by a few tens or hundreds of cells.  Cluster 0, which represented the largest cluster, reported a population of 24,976 cells.  The smallest cluster, which was represented by cluster 36, reported a population of 245 cells. 

<img width="1465" height="962" alt="UMAP_Mouse_Influenza" src="https://github.com/user-attachments/assets/f6466b54-2320-4ca8-a0ce-2970287d5ecb" />

**Figure 2**: UMAP embedding of 156,572 mouse cells across 37 mostly heterogenous clusters.  Cells were included from both normal mice and mice infected with nasal influenza across three different regions of the nasal mucosa, including the respiratory mucosa (RM), olfactory mucosa (OM), and lateral nasal gland (LM).  They were collected across a time span of 2 weeks, collected at five different endpoints (Naive, D02, D05, D08, D14).

### Annotation of Mouse Cell Clusters
To identify cell type/identity for each cluster, all cells and clusters were annotated using SingleR (v2.8.0).  Looking at the annotated UMAP embedding at the cell level (Figure 3), a total of 18 different cell types were identified from the entire cell community of 156,572 cells.

<img width="822" height="491" alt="Mouse_Annotated_UMAP_Cell" src="https://github.com/user-attachments/assets/1ec84166-c8e5-4b87-a59e-87383c30f511" />

**Figure 3**: Annotated UMAP embedding of all 156,572 mouse cells using the SingleR (v2.8.0) program at the cell level.  A total of 18 different cell types were identified from the entire cell community.

Then, looking at the annotated UMAP embedding at the cluster level (Figure 4), which looked more at the more dominant cell types across clusters, this list was consequently reduced to 10 different cell types overall.  Neurons were the most prominent type of cell, representing about 32% of the cell community, with a total frequency of 49,658 cells identified across 5 clusters.  Meanwhile, fibroblasts were the second most prominent cell type, representing about 28% of the cell community, with a total of 27,829 cells spread across 8 clusters, the most number of clusters out of all the categories.  Lastly, the least prominent cell type amongst the 10 dominant ones were T cells, only representing about 0.2% of the community, only being represented by cluster 36, which, again, had a total frequency of 245 cells.

<img width="822" height="491" alt="Mouse_Annotated_UMAP_Cluster" src="https://github.com/user-attachments/assets/44e879d5-16d6-49a9-a4eb-c699ebe6eb28" />

**Figure 4**: Annotated UMAP embedding of all 37 mouse cell clusters using the SingleR (v2.8.0) program at the cluster level.  There were a total of 10 dominant cell types identified, with each category represented by 1-8 clusters.

To check the accuracy of the cluster level automatic annotations, cluster 4, which identified its dominant cell type as B cells, was selected for verification using manual annotation with the scientific literature.  The top three significant gene markers identified for this cluster, based on average log2Fold changes and adjusted p values, were Iglc2 (avg_log2FC = 8.88), Fcer2a (avg_log2FC = 8.60), and Iglc1 (avg_log2FC = 8.56).  The true adjusted p-values for each gene were much smaller and drastically approached 0.

### Genes of Interest

<img width="1920" height="1112" alt="Feature_Plot_Gene_Immune_Cells" src="https://github.com/user-attachments/assets/b5433449-ce46-49c4-8af3-44fff80b3251" />

**Figure 5**: UMAP-based feature plots of the top gene markers for cell types belonging to immune cells.  In order, these genes are the top markers for: (a) *B cells*; (b) *Granulocytes*; (c) *Macrophages*; (d) *Monocytes*; (e) *Natural Killer (NK) cells*; and (f) *T cells*.

<img width="1920" height="1112" alt="Feature_Plot_Gene_Structural_Specialized_Cells" src="https://github.com/user-attachments/assets/c76d5ba8-321e-4152-8e0e-0f2b66332e13" />

**Figure 6**: UMAP-based feature plots of the top gene markers for cell types belonging to other cells in the immune response, such as more structural and specialized cells.  In order, these genes are the top markers for: (a) *Endothelial cells*; (b) *Epithelial cells*; (c) *Fibroblasts*; and (d) *Neurons*.

## Discussion

## References
Kazer, S. W., Match, C. M., Langan, E. M., Messou, M. A., LaSalle, T. J., O'Leary, E., Marbourg, J., Naughton, K., von Andrian, U. H., & Ordovas-Montanes, J. (2024). Primary nasal influenza infection rewires tissue-scale memory response dynamics. *Immunity*, *57*(8), 1955-1974.e8. https://doi.org/10.1016/j.immuni.2024.06.005
