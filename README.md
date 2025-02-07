# Pleiotropy does not facilitate local adaptation in the Silverleaf sunflower (<i>Helianthus argophyllus</i>)

In this study, we explored connections between pleiotropy and local adaptation in the Texas endemic Silverleaf sunflower, or <i>Helianthus argophyllus</i>. Populations of <i>H. argophyllus</i> exhibit a bimodal life history strategy consisting of tall, late-flowering forms and short, early-flowering forms occurring in close geographical proximity. We hypothesize that the differential expression of life history traits within <i>H. argophyllus</i> populations might be linked to local adaptation and controlled by highly pleiotropic genes. This is related to a study by Hämälä <i>et al</i>, 2020, which observed that adaptive genes in common ragweed (Asteraceae family) that were under selection and highly pleiotropic.

Considering <i>H. argophyllus</i> is also a member of the Asteraceae family, we conducted a similar study to examine the extent of pleiotropy in <i>H. argophyllus</i> adaptive genes. We identified spatial selection, selective sweeps, and tested for local adaptation signatures. We assessed pleiotropy by examining whether genes bearing adaptive mutations were more likely than expected to occur in gene coexpression networks and exhibit higher network connectivity. Our results show that candidate locally adapted genes exhibited significantly lower connectivity than non-adapted genes.

<div style="display: flex; justify-content: center; margin-top: 10px;">
  <img src="Figure_Scripts/PLOS_genetics_figures/Moyers_Reisberg.png" alt="PCAdapt" width="800" />
</div>

<p style="margin-bottom: 5px;">
  <strong>A</strong>. South Texas <em>Helianthus argophyllus</em> populations and  
  <strong>B</strong>: Bimodal distribution of life history traits (from Moyers and Rieseberg, 2016).  
</p>

---------------------------------------------------------------------

### Genetic diversity and population structure

PCA (Fig. 1A) and FST (mean FST = 0.21) showed some genetic differentiation between the sampled populations. PC1 (18.4 PVE) and PC2 (13 PVE) showed a clear differentiation between north and coast samples. Admixture analysis using NGSAdmix suggests a population stratification into two distinct genetic groups at K = 2 representing the two sampled populations (Fig. 1B), with differences between the populations remaining at higher K values.

<div style="display: flex; align-items: center; justify-content: space-between;">
  <img src="Figure_Scripts/PLOS_genetics_figures/pca.tiff" alt="PCA plot" width="400" />
  <img src="Figure_Scripts/PLOS_genetics_figures/Admixture_Plot.tiff" alt="K 2" width="400" />
</div>

<p style="margin-bottom: 5px;">
  <strong>Fig 1</strong>. Genetic differentiation and population structure in <em>Helianthus argophyllus</em>.  
  <strong>A</strong>: Genetic differentiation along the first two PCs. Percentages in brackets show the proportion of variance explained by each PC axis.  
  <strong>B</strong>: Approximate admixture proportions for each sampled individual at three different levels of K.  
</p>

---------------------------------------------------------------------

### Selection scans

Two methods, LFMM 2 and PCAdapt, were used to identify genetic regions driving local adaptation in <i>H. argophyllus</i>. LFMM 2 detected selection based on environment-genotype associations, while PCAdapt used population structure. LFMM analysis, controlling for K = 4 latent factors and using the first five BIOCLIM PCs as predictors, identified 1,438 SNPs across 560 genes as selection candidates. PCAdapt identified 4,968 SNPs across 1,122 genes. Only 12 SNPs in 10 genes overlapped between methods, fewer than expected by chance (p < 2.2e-16). Genes from both methods were considered candidates for local adaptation, while non-outlier genes formed the control group.

---------------------------------------------------------------------

### Proportion of eQTLs and eGenes in selection outliers

To assess if selection outlier loci were more pleiotropic than control genes, associations between SNPs and gene expression were tested for trans-regulatory variants. A total of 42,638 genes were used in the eQTL analysis, controlling for population structure with admixture proportions (K = 4). The analysis identified 1,824 eQTL genes (4.3%) associated with the expression of 3,129 eGenes (7.3%) at an FDR < 0.01. Of the 560 LFMM outlier genes, 114 (20.4%) were eQTLs, and 208 (18.5%) of PCAdapt outliers were eQTLs. Only 11 (2%) LFMM outliers and 24 (2.1%) PCAdapt outliers were eGenes. Both outlier gene sets had significantly more eQTLs and fewer eGenes compared to the rest of the transcriptome (p < 2.2e-16, Fisher’s exact test).

<div style="display: flex; justify-content: center; margin-top: 10px;">
  <img src="Figure_Scripts/PLOS_genetics_figures/eQTL_eGene.tiff" alt="K 2" width="600" />
</div>

<p style="margin-bottom: 5px;">
  <strong>Fig 2</strong>. Proportion of eQTLs and eGenes in candidate adaptive genes identified by LFMM and PCAdapt methods compared to the control. Error bars were calculated at 95% bootstrapped CI.  
</p>

---------------------------------------------------------------------

### Evidence of selection at outlier loci

Genetic differentiation, measured by FST, was significantly higher for LFMM outlier genes compared to control genes, while <i>Dxy</i> values were not significantly different. In contrast, PCAdapt outlier genes showed lower <i>FST</i> and higher <i>Dxy</i> than control genes, indicating distinct signatures between the two methods. Population mutation rate measures (<i>θπ</i> and <i>θW</i>) revealed that PCAdapt outlier genes in the north population had lower <i>θW</i> but no significant difference in θπ. In the coast population, PCAdapt outliers had significantly higher θπ. For LFMM outliers, <i>θπ</i> was higher in the north population but not significantly different in the coast population. LFMM outliers in the coast population had lower <i>θW</i> compared to controls. Further validation using Fay and Wu's <i>H</i> neutrality measure and the Composite Likelihood Ratio (CLR) test showed that outlier genes identified by both methods were significantly enriched for selective sweeps, with more negative <i>H</i> estimates compared to control genes. These results indicate that both LFMM and PCAdapt outlier genes are enriched for selective sweeps in both populations.

<div style="display: flex; align-items: center; justify-content: space-between;">
  <img src="Figure_Scripts/PLOS_genetics_figures/FST_Dxy_plot.tiff" alt="FST" width="800" height="400" />
</div>

<img src="Figure_Scripts/PLOS_genetics_figures/North_Coast_CLR.tiff" alt="CLR" width="800" height="400" />
<img src="Figure_Scripts/PLOS_genetics_figures/North_Coast_FayH.tiff" alt="FayH" width="800" height="400" />

<p style="margin-bottom: 5px;">
  <strong>Fig 3</strong>. Genetic differentiation and selective sweep at adaptive genes.  
  <strong>A</strong>: Weighted FST estimates at genes containing selection outliers, eQTLs, and eGenes compared to the transcriptome-wide background.  
  <strong>B</strong>: Distribution of CLR and  
  <strong>C</strong>: *Fay and Wu’s H* values across adaptive genes, eQTLs, eGenes, and the rest of the transcriptome. The horizontal line represents the median value of the transcriptome-wide background for all estimated values. 
</p>

---------------------------------------------------------------------

### Connectivity of selection outliers

Gene co-expression analysis identified 44 modules, with a median of 41 genes per module. Three modules were significantly enriched for LFMM genes, while no modules were enriched for PCAdapt genes or genes identified by both methods. Connectivity analysis revealed that PCAdapt outlier genes and eQTLs had higher connectivity than control genes (p < 2.2e-16), whereas LFMM outlier genes and eGenes were less connected than control genes (p < 2.2e-16).

<div style="display: flex; justify-content: center; margin-top: 10px;">
<img src="Figure_Scripts/PLOS_genetics_figures/connectivity_plot.tiff" alt="connectivity" width="400" height="400" />
</div>

<p style="margin-bottom: 5px;">
  <strong>Fig 4</strong>. Gene network connectivity at candidate adaptive genes.  Estimated connectivity of candidate adaptive genes, eGenes, and eQTLs compared to the transcriptome-wide control.  
  
---------------------------------------------------------------------

### Annotation of adaptation loci and gene co-expression networks

Our analysis revealed that adaptive genes were significantly enriched in pollen recognition and protein serine/threonine kinase activity across both LFMM and PCAdapt genes. LFMM genes were associated with sulfotransferase, catalase, and oxidoreductase activities, while PCAdapt genes were enriched in signal transduction, defense response, and dehydrogenase activity. Genes identified by both methods were linked to defense response, ADP binding, signal transduction, nuclear pore transport, and spliceosome activity. Additionally, WGCNA modules with an overrepresentation of LFMM genes were primarily involved in core metabolic pathways, including histidine biosynthesis, de novo purine biosynthesis, pyruvate oxidation, triacylglycerol biosynthesis, and leucine biosynthesis.

<div style="display: flex; align-items: center; justify-content: center; gap: 20px;">
  <img src="Figure_Scripts/PLOS_genetics_figures/LFMM_gene_over_representation_barplot.tiff" alt="LFMM" width="400" />
  <img src="Figure_Scripts/PLOS_genetics_figures/PCAdapt_gene_over_representation_barplot.tiff" alt="PCAdapt" width="400" />
</div>

<div style="display: flex; justify-content: center; margin-top: 10px;">
  <img src="Figure_Scripts/PLOS_genetics_figures/LFMM_PCAdapt_gene_over_representation_barplot.tiff" alt="PCAdapt" width="400" />
</div>

<p style="margin-bottom: 5px;">
  <strong>Fig 5</strong>. Gene ontology enrichment for biological processes of adaptive genes identified by LFMM (left), PCAdapt (right), and both LFMM and PCAdapt (bottom).
</p> 

<p style="margin-bottom: 20px; margin-top: 20px;">
  <strong>Table 1</strong>. KEGG Module Annotation of the three network modules containing more than expected LFMM adaptive genes. 
</p> 

<div style="clear: both; margin-top: 20px;"> <!-- Ensures table appears at the correct location -->

<table style="width:100%; border-collapse: collapse;">
  <tr>
    <th style="background-color: darkgreen; color: white; padding: 10px; border: 1px solid black; text-align: left;">WGCNA Module</th>
    <th style="background-color: darkgreen; color: white; padding: 10px; border: 1px solid black; text-align: left;">KEGG ID</th>
    <th style="background-color: darkgreen; color: white; padding: 10px; border: 1px solid black; text-align: left;">Description</th>
    <th style="background-color: darkgreen; color: white; padding: 10px; border: 1px solid black; text-align: left;">Adjusted p value</th>
  </tr>
  <tr>
    <td rowspan="2" style="border: 1px solid black; padding: 10px;">Darkgrey</td>
    <td style="border: 1px solid black; padding: 10px;">M00026</td>
    <td style="border: 1px solid black; padding: 10px;">Histidine biosynthesis, PRPP => histidine</td>
    <td style="border: 1px solid black; padding: 10px;">0.00063984</td>
  </tr>
  <tr>
    <td style="border: 1px solid black; padding: 10px;">M00129</td>
    <td style="border: 1px solid black; padding: 10px;">Ascorbate biosynthesis, glucose-1P => ascorbate</td>
    <td style="border: 1px solid black; padding: 10px;">0.0014631</td>
  </tr>
  <tr>
    <td rowspan="4" style="border: 1px solid black; padding: 10px;">Darkturquoise</td>
    <td style="border: 1px solid black; padding: 10px;">M00002</td>
    <td style="border: 1px solid black; padding: 10px;">Glycolysis, core module involving three-carbon compounds</td>
    <td style="border: 1px solid black; padding: 10px;">0.0006257</td>
  </tr>
  <tr>
    <td style="border: 1px solid black; padding: 10px;">M00048</td>
    <td style="border: 1px solid black; padding: 10px;">De novo purine biosynthesis, PRPP + glutamine => IMP</td>
    <td style="border: 1px solid black; padding: 10px;">0.00117279</td>
  </tr>
  <tr>
    <td style="border: 1px solid black; padding: 10px;">M00001</td>
    <td style="border: 1px solid black; padding: 10px;">Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate</td>
    <td style="border: 1px solid black; padding: 10px;">0.00228412</td>
  </tr>
  <tr>
    <td style="border: 1px solid black; padding: 10px;">M00130</td>
    <td style="border: 1px solid black; padding: 10px;">Inositol phosphate metabolism, PI=> PIP2 => Ins(1,4,5)P3 => Ins(1,3,4,5)P4</td>
    <td style="border: 1px solid black; padding: 10px;">0.00288918</td>
  </tr>
  <tr>
    <td rowspan="3" style="border: 1px solid black; padding: 10px;">Green</td>
    <td style="border: 1px solid black; padding: 10px;">M00307</td>
    <td style="border: 1px solid black; padding: 10px;">Pyruvate oxidation, pyruvate => acetyl-CoA</td>
    <td style="border: 1px solid black; padding: 10px;">0.00107016</td>
  </tr>
  <tr>
    <td style="border: 1px solid black; padding: 10px;">M00089</td>
    <td style="border: 1px solid black; padding: 10px;">Triacylglycerol biosynthesis</td>
    <td style="border: 1px solid black; padding: 10px;">0.00954748</td>
  </tr>
  <tr>
    <td style="border: 1px solid black; padding: 10px;">M00432</td>
    <td style="border: 1px solid black; padding: 10px;">Leucine biosynthesis, 2-oxoisovalerate => 2-oxoisocaproate</td>
    <td style="border: 1px solid black; padding: 10px;">0.03465363</td>
  </tr>
</table>

</div> <!-- Ensure this closes properly -->

---------------------------------------------------------------------

### Conclusion

Our results support Fisher's theory, showing that the adaptive genes (LFMM) we found are less likely to affect multiple traits and show signs of positive selection. The reduced connectivity (pleiotropy) in these genes suggests that the two different phenotypes in <i>H. argophyllus</i> have likely reached their best possible adaptation to their environments (optimum). 

---------------------------------------------------------------------

### Methods

- #### Sample Collection and RNA Sequencing
The study used 19 populations of <i>H. argophyllus</i> representing two subpopulations. Seeds were grown in a growth chamber for 3 weeks, and above-ground tissue was frozen for RNA extraction. RNA was sequenced using Illumina's mRNASeq approach on a GAII platform with paired-end sequencing.

- #### Sequence Processing
RNA-seq data were processed with FastQC, trimmed with Trimmomatic, and aligned to the H. annuus reference genome using STAR. Variants were called with Freebayes, filtered with vcftools, and quality control was applied using criteria for read quality, coverage, and missing data.

- #### Genetic Diversity and Population Structure
Genetic diversity metrics (FST, θW, θπ, Tajima’s D) were calculated using ANGSD and pixy. Population structure was assessed using NGSAdmix and various estimators like Dxy, with genetic differentiation measured by FST in a sliding window of 5000 kb.

- #### Selection Outliers
Selection outliers were identified using the LFMM 2 and PCAdapt R packages. LFMM 2 used genotype-environment associations, while PCAdapt focused on population structure to identify locally adapted loci.

- #### eQTL and eGenes
Gene expression was analyzed using FeatureCounts and normalized with Deseq2. eQTLs were identified with MatrixEqtl, controlling for population structure through PCA of expression data and including the first two PCs as covariates.

- #### Selective Sweeps
Selective sweeps were identified using Fay and Wu’s H and CLR methods with SweepFinder2, focusing on selection outlier SNPs and genes. The analysis considered SNPs in selected genomic regions.

- #### Statistical Comparisons
Statistical comparisons of genetic estimators across selection outliers, eQTLs, eGenes, and control genes were performed using the Wilcoxon signed-rank test.

- #### Gene Coexpression Networks
Gene coexpression networks were constructed with the WGCNA R package, using a power threshold of six for adjacency matrices. Connectivity differences between selection outliers, eQTLs, eGenes, and the transcriptome-wide background were analyzed.

- #### Annotation of Genes and Network Modules
Gene and pathway overrepresentation analysis of selection outliers and coexpression modules were performed using clusterProfiler. Gene IDs were extracted using NCBI-blast for alignment to the HanXRQr2.0-SUNRISE genome, and KEGG pathway analysis was conducted using the ‘enrichKEGG’ function.

<div style="text-align: center;">
  <img src="Figure_Scripts/PLOS_genetics_figures/Methods_flowchart.png" alt="LFMM" width="400" />
</div>

### References

1. Tenaillon O. The Utility of Fisher’s Geometric Model in Evolutionary Genetics. *Annu Rev Ecol Evol Syst.* 2014;45:179. doi:[10.1146/ANNUREV-ECOLSYS-120213-091846](https://doi.org/10.1146/ANNUREV-ECOLSYS-120213-091846)

2. Hämälä T, Gorton AJ, Moeller DA, Tiffin P. Pleiotropy facilitates local adaptation to distant optima in common ragweed (*Ambrosia artemisiifolia*). *PLoS Genet.* 2020;16:e1008707. doi:[10.1371/JOURNAL.PGEN.1008707](https://doi.org/10.1371/JOURNAL.PGEN.1008707)

3. Wagner GP, Zhang J. The pleiotropic structure of the genotype–phenotype map: the evolvability of complex organisms. *Nature Reviews Genetics.* 2011;12:3:204–213. doi:[10.1038/nrg2949](https://doi.org/10.1038/nrg2949)

4. Rennison DJ, Peichel CL. Pleiotropy facilitates parallel adaptation in sticklebacks. *Mol Ecol.* 2022;31:1476–1486. doi:[10.1111/MEC.16335](https://doi.org/10.1111/MEC.16335)

5. Edwards AWF. The Genetical Theory of Natural Selection. *Genetics.* 2000;154:1419–1426. doi:[10.1093/GENETICS/154.4.1419](https://doi.org/10.1093/GENETICS/154.4.1419)



