# binf6110-assignment-4

## Introduction

## Methods

## Results
### Dimensionality Reduction and Initial Clustering of Mouse Cells
Dimensionality reduction of mouse cell data from the Seurat object was carried out using PC analysis to identify relevant PCs to retain for clustering by .  UMAP embedding was then used to visualize all of the distinct clusters in space with the applied PCs.  When looking at the elbow plot in Figure 1 after applying log-normalization and scaling, 20 PCs were identified in total from the data.  However, the variance in cells appeared to decrease at around 15 PCs, as evident by the elbow point levelling off at around PC15.  As such, up to 15 PCs were applied to further downstream analysis for UMAP embedding.

<img width="1482" height="1007" alt="Mouse_Elbow_Plot" src="https://github.com/user-attachments/assets/871c79d6-5f04-4681-ae07-af5578cfbdcb" />

**Figure 1**: Elbow plot of principle components (PCs) for 156,572 cells collected from normal mice and mice infected with nasal influenza.  A total of 20 PCs were identified after scaling and normalizing of data.  The elbow point appeared to occur at around 15 PCs

Looking at the initial UMAP for the mouse cells (Figure 2), a total of 37 distinct clusters were identified.  These clusters appeared mostly heterogenous from each other, with certain clusters appearing to form larger groups than others, while still remaining distinct between each other.  Overlap between clusters did appear to be present within the clusters upon closer inspection, but these instances were usually minor and, for the most part, did not have a major impact on heterogeneity.  Cluster sizes were also varied, though certain groups appeared to have similar sizes between each other, only differing by a few tens or hundreds of cells.  Cluster 0, which represented the largest cluster, reported a cluster size of 24,976 cells.  The smallest cluster, which was represented by cluster 36, reported a size of 245 clusters. 

<img width="1465" height="962" alt="UMAP_Mouse_Influenza" src="https://github.com/user-attachments/assets/f6466b54-2320-4ca8-a0ce-2970287d5ecb" />

**Figure 2**: UMAP embedding of 156,572 mouse cells across 37 mostly heterogenous clusters.  Cells were included from both normal mice and mice infected with nasal influenza across three different regions of the nasal mucosa, including the respiratory mucosa (RM), olfactory mucosa (OM), and lateral nasal gland (LM).  They were collected across a time span of 2 weeks, collected at five different endpoints (Naive, D02, D05, D08, D14).

## Discussion

## References
