# binf6110-assignment-4

## Introduction

One such study, for example, which was conducted by Kazer et al. in 2024, has even attempted to .

Single cell RNA-sequencing (scRNA-seq) has often been at the heart of understanding many respiratory infections, especially influenza (Medaglia et al., 2022).  It has often been used to investigate elements of host pathogenicity, treatment, and immune response at the cellular level, such as combination drug therapy, innate immunity or transcriptional markers in various cell types indicative of disease progression (Medaglia et al., 2022, Zhang et al., 2023, Kelly et al., 2022).  One of the most widely used programs that has been used to conduct many different types of scRNA-seq analyses is the *Seurat* package in R (Satija et al., 2015).  When it comes to the visualization of these clusters, the most prominent method that has been greatly accepted today is the Uniform Manifold Approximation and Projection (UMAP) (Becht et al., 2018).  UMAP embedding has typically been most useful in illustrating high-dimensional data for cell clusters in a low-dimensional space for many scRNA-seq analyses, being more robust at creating highly reproducible visualizations that take less computing time (Becht et al., 2018, Xiang et al., 2021).  However, the most noteworthy aspect about this method that has made it more valuable than prior visualization methods, most notably, t-distributed stochastic neighbor embedding (t-SNE), is its capability of determining cell clusters that best reflect and maintain the global cell populations from the samples and/or organisms they were obtained from (Becht et al., 2018, Xiang et al., 2021).  Unfortunately, one key limitation that remains in this environment of low-dimensional analysis of high-dimensional data, both in Seurat and UMAP, is that embedding these clusters requires cell data to be log-normalized and scaled, which can often underestimate the nature of the cell populations that raw unfiltered counts provide (Verma & Engelhardt, 2020).

The goal of this analysis was to explore the different types of and using single-cell RNA sequencing 


## Methods
Raw single-cell RNA sequencing data was obtained and adopted from the study conducted by Kazer et al. in 2024 and transformed into a Seurat object, created by Dr. Elias Taylor at University of Guelph, to import into R (v4.4.3) (2026).  In this Seurat object, mice had been infected with the Influenza A virus, and a total of 156,572 cells were collected in the nasal mucosa found in the upper respiratory tract (Kazer et al., 2024).  Over a post-infection period of about 14 days, cells were collected during five specific timepoints, at Day 0 (Naive), Day 2 (D02), Day 5 (D05), Day 8 (D08), and Day 14 (D14) (Kazer et al., 2024).  In addition, during these timepoints, the cells were also collected from across three different types of tissue in the nasal mucosa, which were the respiratory mucosa (RM), olfactory mucosa (OM), and lateral nasal gland (LM) (Kazer et al., 2024).  Once collected, the data was pre-processed using Cell Ranger (v), and a Seurat object was created using the *Seurat* package (v) in R (v.4.4.3) (Satija et al., 2015, Hao et al., 2023, Taylor, 2026).

All clusters were then automatically annotated for their most dominant cell types by using the *SingleR* package (v2.8.0), both at the cell level and the cluster level (Aran et al., 2019).  As a reference dataset, the MouseRNAseqData was downloaded and used from the *celldex* package (v1.16.0), adopted from an epigenomic and transcriptomic mice study conducted by Benayoun et al. (2019, Aran et al., 2019).  First, to explore all the different cell types across the entire cell community, the normalized RNA assay was extracted from the UMAP Seurat object and ran through SingleR.  Then, for cluster level annotation, the normalized RNA assays of the UMAP Seurat object were pseudobulked and used in SingleR.  Both levels of annotation were then visualized using separate annotated UMAPs of the cells and clusters.  Finally, to verify the accuracy of these annotations by SingleR using manual annotation, one cluster was selected, and its top three gene markers were used to conduct a literature search.

To determine the genes of interest, which were defined as the top significant gene markers for each annotated cell type, the *FindAllMarkers* function in the *Seurat* package was first used on the clustered cells to obtain all gene markers for each cluster as a starting directory (Satija et al., 2015, Hao et al., 2023).  This data was then filtered to include only the top gene marker per cluster, with the following parameters: an average log2Fold change greater than 0.5 (*avg_log2FC* > 0.5), an adjusted p-value less than 0.05 (*p_val_adj* < 0.05), at least 25% coverage of cells in the cluster (*pct.1* > 0.25), and no more than 10% coverage of cells in other clusters (*pct.2* < 0.10).  These parameters were selected to facilitate more discovery of important and distinct gene markers, while still finding those that were more unique to specific cell types, as multiple stricter parameters, such as pct.1, often limited the number of top genes per cluster or resulted in a cluster not reporting any genes.  Once these top gene markers were identified per cluster, this list was then further filtered to identify the topmost gene for each cell type, given that each annotation was still represented by multiple clusters and thus multiple top genes, whilst still retaining the same parameters.  The expression patterns of these top genes were then visualized and examined using UMAP embedding as feature plots to justify that these markers were unique to each of their respective cell types.

Finally, to understand how gene expression differed between tissue types in the nasal mucosa, as well as how it differed across different timepoints, both before infection and post-infection, differential expression analysis was carried out on one selected cluster using the DESeq2 package (Love et al., 2014).

## Results
### Dimensionality Reduction and Initial Clustering of Mouse Cells
Dimensionality reduction of mouse cell data from the Seurat object was carried out using PC analysis, to identify relevant PCs to retain for clustering.  UMAP embedding was then used to visualize all of the distinct clusters and their spatial separation with the applied PCs.  When looking at the elbow plot in Figure 1 after applying log-normalization and scaling, a total 20 PCs were identified in total from the data.  However, the variance in cells appeared to decrease at around 15 PCs, as evident by the elbow point levelling off at around PC15.  As such, up to 15 PCs were applied to further downstream analysis for UMAP embedding.

<img width="1482" height="1007" alt="Mouse_Elbow_Plot" src="https://github.com/user-attachments/assets/871c79d6-5f04-4681-ae07-af5578cfbdcb" />

**Figure 1**: Elbow plot of principle components (PCs) for 156,572 cells collected from normal mice and mice infected with nasal influenza.  A total of 20 PCs were identified after scaling and normalizing of data.  The elbow point appeared to occur at around 15 PCs.

Looking at the initial UMAP for the mouse cells (Figure 2), a total of 37 distinct clusters were identified.  These clusters appeared mostly heterogenous from each other, with certain clusters appearing to form larger groups than others, while still remaining distinct between each other.  Overlap between clusters did appear to be present within the clusters upon closer inspection, but these instances were usually minor and, for the most part, did not have a major impact on heterogeneity.  Cluster sizes were also varied, though certain groups appeared to have similar sizes between each other, only differing by a few tens or hundreds of cells.  Cluster 0, which represented the largest cluster, reported a population of 24,976 cells.  The smallest cluster, which was represented by cluster 36, reported a population of 245 cells. 

<img width="1465" height="962" alt="UMAP_Mouse_Influenza" src="https://github.com/user-attachments/assets/f6466b54-2320-4ca8-a0ce-2970287d5ecb" />

**Figure 2**: UMAP embedding of 156,572 mouse cells across 37 mostly heterogenous clusters.  Cells were included from both normal mice and mice infected with nasal influenza across three different regions of the nasal mucosa, including the respiratory mucosa (RM), olfactory mucosa (OM), and lateral nasal gland (LM).  They were collected across a time span of 14 days, collected at five different timepoints (Naive, D02, D05, D08, D14).

### Annotation of Mouse Cell Clusters
To identify cell type/identity for each cluster, all cells and clusters were annotated using SingleR (v2.8.0).  Looking at the annotated UMAP embedding at the cell level (Figure 3), a total of 18 different cell types were identified from the entire cell community of 156,572 cells.  These cell types matched accordingly with the cell types determined in the reference dataset (Benayoun et al., 2019).

<img width="822" height="491" alt="Mouse_Annotated_UMAP_Cell" src="https://github.com/user-attachments/assets/1ec84166-c8e5-4b87-a59e-87383c30f511" />

**Figure 3**: Annotated UMAP embedding of all 156,572 mouse cells using the SingleR (v2.8.0) program at the cell level.  A total of 18 different cell types were identified from the entire cell community.

Then, looking at the annotated UMAP embedding at the cluster level (Figure 4), which looked more at the more dominant cell types across clusters, this list was consequently reduced to 10 different cell types overall.  Neurons were the most prominent type of cell, representing about 32% of the cell community, with a total frequency of 49,658 cells identified across 5 clusters.  Meanwhile, fibroblasts were the second most prominent cell type, representing about 28% of the cell community, with a total of 27,829 cells spread across 8 clusters, the most number of clusters out of all the categories.  Lastly, the least prominent cell type amongst the 10 dominant ones were T cells, only representing about 0.2% of the community, only being represented by cluster 36, which, again, had a total frequency of 245 cells.

<img width="822" height="491" alt="Mouse_Annotated_UMAP_Cluster" src="https://github.com/user-attachments/assets/44e879d5-16d6-49a9-a4eb-c699ebe6eb28" />

**Figure 4**: Annotated UMAP embedding of all 37 mouse cell clusters using the SingleR (v2.8.0) program at the cluster level.  There were a total of 10 dominant cell types identified, with each category represented by 1-8 clusters.

To check the accuracy of the cluster level automatic annotations, cluster 4, which identified its dominant cell type as B cells, was selected for verification using manual annotation with the scientific literature.  The top three significant gene markers identified for this cluster, based on average log2Fold changes and adjusted p values, were Iglc2 (avg_log2FC = 8.88), Fcer2a (avg_log2FC = 8.60), and Iglc1 (avg_log2FC = 8.56).  The true adjusted p-values for each gene were much smaller and drastically approached 0.

### Genes of Interest

After annotating each cluster with SingleR, the top significant gene marker was determined for each of the dominant cell types based on the average log2Fold change and adjusted p-values.  Using UMAP embedding, these genes of interest were then visualized using feature plots, separating them into two groups based on genes belonging to immune cell types (Figure 4) and other cell types belonging to more structural/specialized functions (Figure 5).  The average log2Fold changes ranged from about 4.65, which was exhibited by the Nqo1 gene in neurons, to 11.54, which was exhibited by the Gpihbp1 gene in endothelial cells.  Meanwhile, the adjusted p-values for all top gene markers were much smaller that they all greatly approached 0.  Looking at the feature plots, each gene appeared to show a distinct expression pattern, being highly expressed in individual clusters that differed from each other, without major overlap between the dominant cell types.  This helped to further validate the annotations of cell type identities, with each category being supported by at least one uniquely expressed gene marker.

<img width="1920" height="1112" alt="Feature_Plot_Gene_Immune_Cells" src="https://github.com/user-attachments/assets/b5433449-ce46-49c4-8af3-44fff80b3251" />

**Figure 5**: UMAP-based feature plots of the top gene markers for cell types belonging to immune cells (*avg_log2FC* > 0.5, *p_val_adj* < 0.05).  In order, these genes are the top markers for: (a) *B cells*; (b) *Granulocytes*; (c) *Macrophages*; (d) *Monocytes*; (e) *Natural Killer (NK) cells*; and (f) *T cells*.

<img width="1920" height="1112" alt="Feature_Plot_Gene_Structural_Specialized_Cells" src="https://github.com/user-attachments/assets/c76d5ba8-321e-4152-8e0e-0f2b66332e13" />

**Figure 6**: UMAP-based feature plots of the top gene markers for cell types belonging to other cells in the immune response, such as more structural and specialized cells (*avg_log2FC* > 0.5, *p_val_adj* < 0.05).  In order, these genes are the top markers for: (a) *Endothelial cells*; (b) *Epithelial cells*; (c) *Fibroblasts*; and (d) *Neurons*.

## Discussion

## References
Aran, D., Looney, A. P., Liu, L., Wu, E., Fong, V., Hsu, A., Chak, S., Naikawadi, R. P., Wolters, P. J., Abate, A. R., Butte, A. J., & Bhattacharya, M. (2019). Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. *Nature Immunology*, *20*, 163-172. https://doi.org/10.1038/s41590-018-0276-y

Benayoun, B., Pollina, E. A., Singh, P. P., Mahmoudi, S., Harel, I., Casey, K. M., Dulken, B. W., Kundaje, A., & Brunet, A. (2019). Remodeling of epigenome and transcriptome landscapes with aging in mice reveals widespread induction of inflammatory responses. *Genome Research*, *29*, 697-709. http://genome.cshlp.org/lookup/doi/10.1101/gr.240093.118

Becht, E., McInnes, L., Healy, J., Dutertre, C. A., Kwok, I. W. H., Ng, L. G., Ginhoux, F., & Newell, E. W. (2018). Dimensionality reduction for visualizing single-cell data using UMAP. *Nature Biotechnology*, *37*(1), 38-44. https://doi.org/10.1038/nbt.4314

Hao, Y., Stuart, T., Kowalski, M. H., Choudhary, S., Hoffman, P., Hartman, A., Srivastava, A., Molla, G., Madad, S., Fernandez-Granda, C., & Satija, R. (2023). Dictionary learning for integrative, multimodal and scalable single-cell analysis. *Nature Biotechnology*, *42*, 293-304. https://doi.org/10.1038/s41587-023-01767-y

Kazer, S. W., Match, C. M., Langan, E. M., Messou, M. A., LaSalle, T. J., O'Leary, E., Marbourg, J., Naughton, K., von Andrian, U. H., & Ordovas-Montanes, J. (2024). Primary nasal influenza infection rewires tissue-scale memory response dynamics. *Immunity*, *57*(8), 1955-1974.e8. https://doi.org/10.1016/j.immuni.2024.06.005

Kelly, J. N., Laloli, L., V'kovski, P., Holwerda, M., Portmann, J., Thiel, V., & Dijkman, R. (2022). Comprehensive single cell analysis of pandemic influenza A virus infection in the human airways uncovers cell-type specific host transcriptional signatures relevant for disease progression and pathogenesis. *Frontiers in Immunology*, *13*:978824. https://doi.org/10.3389/fimmu.2022.978824

Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, *15*(12), 550. https://doi.org/10.1186/s13059-014-0550-8

Megdalia, C., Kolpakov, I., Zwygart, A. C. A., Zhu, Y., Constant, S., Huang, S., Cagno, V., Dermitzakis, E. T., Stellaci, F., Xenarios, I., & Tapparel, C. (2022). An anti-influenza combined therapy assessed by single cell RNA-sequencing. *Communications Biology*, *5*, 1075. https://doi.org/10.1038/s42003-022-04013-4

Satija, R., Farrell, J. A., Gennert, D., Schier, A. F., & Regev, A. (2015). Spatial reconstruction of single-cell gene expression data. *Nature Biotechnology*, *33*, 495-502. https://doi.org/10.1038/nbt.3192

Verma, A. & Engelhardt, B. E. (2020). A robust nonlinear low-dimensional manifold for single cell RNA-seq data. *BMC Bioinformatics*, *21*, 324. https://doi.org/10.1186/s12859-020-03625-z

Xiang, R., Wang, W., Yang, L., Wang, S., Xu, C., & Chen, X. (2021). A Comparison for Dimensionality Reduction Methods of Single-Cell RNA-seq Data. *Frontiers in Genetics*, *12*:646936. https://doi.org/10.3389/fgene.2021.646936

Zhang, Y., Zong, L., Zheng, Y., Zhang, Y., Li, N., Li, Y., Jin, Y., Chen, L., Ouyang, J., Bibi, A., Huang, Y., & Xu, Y. (2023). A single-cell atlas of the peripheral immune response in patients with influenza A virus infection. *iScience*, *26*(12), 108507. https://doi.org/10.1016/j.isci.2023.108507
