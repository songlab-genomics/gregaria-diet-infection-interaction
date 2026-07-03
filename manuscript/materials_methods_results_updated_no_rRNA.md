# Diet infection transcriptomics manuscript draft

Authors: Emily Baker and Maeva A. Techer

Draft updated from the current no-rRNA workflowR analysis in
`gregaria-diet-infection-interaction`. Values are taken from the latest rendered
outputs generated on 2026-07-03.

## Items to confirm before final submission

- Exact diet formulation represented by diet codes 33, 50, and 83.
- Fungal strain, inoculation dose, inoculation method, and control treatment.
- Confirm the exact fat body dissection description and RNA integrity assay
  used after Qubit quantification.
- Exact fastp version and whether any project-specific changes were made to the
  STAR/featureCounts Snakemake workflow relative to the AGY analysis.
- Whether Metarhizium-mapped reads are available for a formal fungal-load
  covariate or validation figure.

## Materials and Methods

### Experimental design

We tested whether diet modifies the host transcriptional response to fungal
challenge in female *Schistocerca gregaria*. Individuals were assigned to one of
three diet treatments, encoded as 33, 50, and 83, and to one of two infection
treatments, encoded as Control or Infected. All samples in the current metadata
were collected 96 h post-inoculation. Age was controlled by the laboratory
design, reducing the likelihood that age-dependent expression differences drove
the diet or infection contrasts. The current analysis focuses on fat body
RNA-seq because this tissue is central to insect immune, metabolic, storage, and
nutritional physiology.

The full metadata table contained 44 libraries distributed across the six
diet-by-infection groups: diet 33 had seven control and nine infected samples,
diet 50 had six control and nine infected samples, and diet 83 had seven
control and six infected samples. Four control libraries, samples 1036 and 1039
from diet 33 and samples 1007 and 1037 from diet 50, were flagged during
quality-control exploration and evaluated in a separate outlier-sensitivity
analysis.

### Tissue sampling and RNA extraction

Eight to ten females were dissected per condition and treatment in sterile
conditions by a single experimenter. Live individuals were dissected within 20
minutes using cold locust saline and pre-chilled dissection plates to limit RNA
degradation. Fat body tissue was dissected with minimal disturbance to
surrounding tissues, snap-frozen in liquid nitrogen, and stored at -80 degrees C
until molecular processing.

Total RNA was extracted using an automated Promega Maxwell RSC instrument with
the simplyRNA Tissue Kit, following the manufacturer's protocol with minor
modifications. RNA yield was quantified using a Qubit 4 fluorometer and the HS
RNA Assay. RNA quality and integrity were assessed before library selection, and
samples were retained for sequencing based on dissection quality, RNA yield and
integrity, and observations of individual vitality.

### Library preparation and sequencing

Total RNA libraries were prepared in-house using the Illumina Stranded Total
RNA kit with RiboZero rRNA depletion. Libraries were pooled and sequenced by
Novogene on an Illumina NovaSeq X Plus 25B platform, targeting approximately 45
million paired-end reads per library. Residual rRNA signal was quantified from
the pre-filtered gene-count matrix as the percentage of assigned gene-level
counts mapping to annotated rRNA loci. Because this residual signal was
substantial, annotated rRNA loci were removed bioinformatically before
downstream expression analyses.

### Read quality control, mapping, and quantification

RNA-seq preprocessing was performed upstream on the HPC with a custom Snakemake
workflow. Raw paired-end reads were adapter- and quality-trimmed with fastp
using automatic paired-end adapter detection, trimming of the first two bases
from each read, and a minimum retained read length of 50 bp. fastp JSON and HTML
reports were retained for read-quality assessment.

Trimmed reads were aligned to the *S. gregaria* reference genome
GCF_023897955.1 with STAR v2.7.10b. The STAR genome index was built with the
corresponding reference annotation, `sjdbOverhang` set to 149, and
`alignIntronMax` set to 2,500,000. The maximum intron size was increased beyond
the STAR default after assessing species-specific intron lengths in
*Schistocerca* genomes. Alignments were run in two-pass mode with sorted
coordinate BAM output, STAR GeneCounts output, transcriptome BAM output, and BAM
CSI indexing with samtools v1.21. Gene-level counts were generated with
featureCounts from Subread using paired-end fragment counting
(`--countReadPairs`), transcript and exon features, and `gene_id` as the
grouping attribute.

This workflowR repository begins from the resulting gene-level count files
rather than from FASTQ or BAM files. The upstream HPC workflow is summarized in
the workflow figure, while the R Markdown analyses document all downstream
count filtering, quality control, differential expression, and enrichment steps.

Per-sample count files were merged into a master gene-count matrix. Sample names
in the count matrix were converted back to biological sample labels by removing
the `mehreen_` prefix and `_MERGE` suffix, then matched exactly to the metadata
before any plotting or modeling.

### rRNA removal and expression filtering

Because the libraries were generated with RiboZero and rRNA-derived counts are
not reliable evidence for host gene expression, rRNA loci were removed before
all current analyses. rRNA filtering used the project file
`data/excluded_loci/gregaria_rrna_list.txt`. The original matrix contained
35,681 gene rows. The rRNA list contained 12,402 IDs, all of which were present
in the count matrix, leaving 23,279 genes in the active no-rRNA matrix
`data/raw_read_counts/master_counts_no_rRNA.csv`.

Across the original assigned count matrix, rRNA loci accounted for 20.08% of
all assigned gene-level counts. At the library level, the median rRNA fraction
was 20.17%, with values ranging from 7.54% to 63.65%. Thus, rRNA removal was
not only a precaution based on the library protocol, but a measurable filtering
step required to prevent residual rRNA features from contributing to QC,
non-coding RNA summaries, and differential expression tests.

After rRNA removal, library sizes ranged from 777,155 to 8,748,179 assigned
counts, with a median of 3,214,846 assigned counts. Genes with fewer than 10
total counts across the experiment were removed before DESeq2 modeling and VST
exploration, retaining 16,892 genes for the main modeled analyses.

### Quality control and outlier sensitivity

Sample-level quality control used library-size summaries, variance-stabilized
principal component analysis (PCA), sample-distance heatmaps, and heatmaps of
the most variable VST-transformed genes. These plots were used to assess whether
libraries clustered by the experimental design and whether any individual
library was unusually distant from its expected diet-by-treatment group.

The main analyses retain all 44 samples. Samples 1036 and 1039 from diet 33
controls and samples 1007 and 1037 from diet 50 controls were also removed in a
sensitivity analysis because they had atypical global expression patterns. This
sensitivity analysis was not treated as an automatic replacement for the full
dataset because the removed samples were all controls; removing them changes the
group balance as well as the expression structure.

### Differential expression analysis

Differential expression was performed with DESeq2 using raw no-rRNA gene counts.
Genes were considered differentially expressed when the Benjamini-Hochberg
adjusted P-value was below 0.05 and the absolute log2 fold-change was at least
1. Variance-stabilized counts were used only for exploratory PCA and heatmaps,
not as input to DESeq2.

Three complementary DESeq2 models were used. First, a full diet-by-infection
model was fitted:

```r
~ Diet + Treatment + Diet:Treatment
```

In this model, diet 50 and Control were the reference levels. The main infection
coefficient therefore estimates the infected-versus-control response in diet 50,
and the interaction terms test whether the infection response in diets 33 or 83
differs from the diet-50 response.

Second, pairwise diet contrasts were fitted with infection status as a
covariate:

```r
~ Treatment + Diet
```

This model asks whether expression differs among diets after accounting for
infection treatment.

Third, infection effects were tested within each diet separately:

```r
~ Treatment
```

These within-diet models provide the clearest estimates of the infection
response under each nutritional condition. Positive log2 fold-changes indicate
higher expression in infected samples, whereas negative log2 fold-changes
indicate higher expression in controls.

### Gene biotype and candidate-gene annotation

Gene classes were assigned from the reference annotation and local annotation
tables to distinguish coding genes, long non-coding RNAs, and other non-coding
features. Because rRNA loci were removed before modeling, the remaining
"other non-coding" category represents non-rRNA non-coding features retained in
the reference rather than rRNA signal.

Functional annotation used the local eggNOG annotation file for *S. gregaria*
and the corresponding reference GFF. CDS records in the GFF were used to bridge
protein-based eggNOG identifiers to `LOC...` gene IDs in the count matrix. This
lifted eggNOG annotation onto 13,526 genes in the no-rRNA count matrix,
including 6,533 genes with GO terms, 7,412 genes with KEGG orthology IDs, and
4,337 genes with KEGG pathway annotations.

### GO and KEGG enrichment

GO and KEGG pathway enrichment were performed with custom TERM2GENE maps derived
from the local eggNOG annotations. Enrichment was tested for non-directional DEG
sets: all infection DEGs, infection DEGs within each diet, genes shared by at
least two diet-specific infection responses, and genes shared by all three diet
responses. Up- and down-regulated genes were not separated for the current GO
and KEGG enrichment pages, so each enriched term represents the functional
composition of the DEG set as a whole. The enrichment background was restricted
to genes present in the no-rRNA count matrix and represented in the
corresponding annotation map.

### Curated heatmap and publication figures

A curated VST heatmap was generated to summarize interpretable candidate genes.
Genes were prioritized if they were differentially expressed in one or more
within-diet infection contrasts and/or had eggNOG descriptions consistent with
immune function, nutrient storage, lipid transport, or metabolism. Samples were
ordered by infection status and then diet rather than by unsupervised clustering,
so the heatmap directly compares control and infected states across the diet
gradient. Genes were grouped first by DEG status and then by functional theme.

## Results

### rRNA removal refines the active expression universe

The current analysis uses a no-rRNA matrix containing 23,279 genes, compared
with 35,681 genes in the original count matrix. Low-count filtering retained
16,892 genes for DESeq2 modeling. This filtering step is important because
RiboZero libraries can still contain residual rRNA-derived reads, and those rows
could otherwise inflate apparent non-coding signal or distort library-level
summaries.

### Global expression structure follows diet and infection treatment

PCA of the top 500 most variable VST-transformed genes showed strong global
structure across samples. PC1 explained 35.1% of the variance and PC2 explained
28.0%. Samples were not randomly distributed in this reduced expression space;
instead, they showed structure associated with diet and infection status. This
supported the use of both full interaction models and diet-specific infection
contrasts.

### Infection responses differ among diets

Diet-specific infected-versus-control contrasts identified infection-associated
DEGs in all three diets. Diet 50 showed the largest response, with 1,312 DEGs
including 794 genes higher in infected samples and 518 genes higher in controls.
Diet 33 showed 898 DEGs, including 457 higher in infected samples and 441 higher
in controls. Diet 83 showed the smallest response, with 519 DEGs, including 342
higher in infected samples and 177 higher in controls.

Across diets, most infection-associated DEGs were coding genes. The combined
within-diet infection DEG burden included 2,579 coding-gene calls, 144 long
non-coding RNA calls, and six other non-coding calls after rRNA removal. Thus,
the current non-coding signal is not driven by retained rRNA rows.

### The full interaction model supports diet-dependent infection effects

The full diet-by-infection DESeq2 model identified 1,064 genes associated with
the main infection coefficient in the diet-50 reference background, including
734 positive and 330 negative infection coefficients. The diet-33-by-infection
interaction term identified 65 genes, all with negative interaction
coefficients, indicating that a subset of infection responses differed between
diet 33 and the diet-50 reference. The diet-83-by-infection interaction term
identified two genes, both with positive interaction coefficients.

Main diet effects in the full model were much smaller than the infection effect:
three genes were significant for diet 33 versus 50, and 20 genes were
significant for diet 83 versus 50. This pattern suggests that infection is the
dominant driver of differential expression in the current model, while diet
modifies the magnitude and composition of the infection response.

### Pairwise diet effects are detectable but smaller than infection effects

Pairwise diet contrasts, adjusted for infection status, identified fewer DEGs
than the within-diet infection contrasts. Diet 33 versus 50 identified 58 DEGs,
diet 33 versus 83 identified 78 DEGs, and diet 50 versus 83 identified 30 DEGs.
These results indicate that baseline diet-associated transcriptional
differences are present, but they are smaller than the transcriptional response
to fungal challenge in this dataset.

### Outlier removal changes diet-33 and diet-50 infection responses

Removing samples 1036 and 1039 from diet 33 controls and samples 1007 and 1037
from diet 50 controls changed the infection DEG burden in both affected diets.
Diet 83 was unchanged because no diet-83 samples were removed. In diet 33, the
full dataset identified 898 infection DEGs, whereas the no-outlier sensitivity
analysis identified 803 DEGs, including 417 genes higher in infected samples and
386 genes higher in controls. Of the diet-33 DEGs, 677 were retained in both
analyses, 221 were lost after removing the outliers, and 126 were gained.

The effect was larger in diet 50. The full dataset identified 1,312 infection
DEGs in diet 50, whereas the no-outlier sensitivity analysis identified 2,885
DEGs, including 2,453 genes higher in infected samples and 432 genes higher in
controls. Of the diet-50 DEGs, 1,022 were retained in both analyses, 290 were
lost after removing the outliers, and 1,863 were gained.

These results show that infection contrasts in diets 33 and 50 are sensitive to
the flagged control libraries, especially diet 50. Because removal reduces the
diet-33 control group from seven to five samples and the diet-50 control group
from six to four samples, the sensitivity analysis should be presented as
evidence about robustness and uncertainty rather than as a simple replacement
for the full dataset.

### Infection DEGs are enriched for extracellular, microbial-response, and cuticle-associated GO terms

GO enrichment of all infection DEGs highlighted extracellular and immune-related
terms. The strongest terms included extracellular region, extracellular space,
aminoglycan metabolic process, response to bacterium, carbohydrate derivative
binding, peptidoglycan binding, glycosaminoglycan binding, peptidoglycan
metabolic and catabolic processes, response to biotic stimulus, defense response
to Gram-positive bacterium, defense response to bacterium, cuticle development,
and chitin-based cuticle development.

The infection DEG set contained 1,731 unique genes with evidence from at least
one diet-specific infection contrast. Of these, 721 were shared by at least two
diets and 277 were shared by all three diets. This overlap indicates that part
of the host response to fungal challenge is shared across nutritional
conditions, while a substantial fraction remains diet-specific.

### KEGG enrichment highlights immune signaling, metabolism, detoxification, and transport

KEGG enrichment of all infection DEGs identified pathways related to metabolism,
immune signaling, transport, and detoxification. Enriched pathways included
biosynthesis of secondary metabolites, Toll and Imd signaling pathway, arginine
and proline metabolism, glycerolipid metabolism, glycine, serine and threonine
metabolism, carbon metabolism, glutathione metabolism, insect hormone
biosynthesis, xenobiotic and cytochrome P450 metabolism, ABC transporters, and
protein digestion and absorption. The Toll and Imd signaling pathway was also
enriched in diet-specific and shared DEG sets, consistent with activation of
canonical insect immune pathways.

The genes shared by all three diet-specific infection responses were enriched
for bile secretion and ABC transporters. Although these pathway labels are broad
KEGG categories, they point to transport, lipid handling, and detoxification
processes that are biologically plausible for fat body responses to infection
and diet.

### Curated immune, storage, and nutrition candidates support the fat body interpretation

The curated candidate table highlighted immune and nutritional genes that were
both statistically responsive and biologically interpretable. Recurrent
infected-higher immune candidates included peptidoglycan-recognition proteins,
serpin-family genes, lysozyme-like genes, STAT-related genes, C-type lectin or
carbohydrate-recognition-domain genes, and clip-domain serine protease-related
genes. Storage and lipid-transport candidates included hemocyanin-domain genes,
lipid transporter genes, and a nutrient-reservoir gene. Nutrition and metabolism
candidates included amino acid transporters, sodium-ion transporters, lipid
binding genes, trypsin-like genes, and broader transporter genes.

Together, the DEG burden, enrichment results, and curated VST heatmap support a
model in which fungal challenge induces a broad fat body transcriptional
response involving immune recognition, extracellular defense, cuticle-associated
features, transport, detoxification, storage, and metabolic remodeling. Diet 50
currently shows the largest infection-associated transcriptional response, but
the outlier sensitivity analysis indicates that this contrast should be
interpreted with care.

## Figure-ready result statements

**Figure 1. Study design and workflow.** Fat body RNA-seq libraries from three
diet treatments and two infection treatments were processed upstream on the HPC
using a STAR-based RNA-seq workflow, then analyzed in workflowR from merged
gene-count files. rRNA loci were removed before QC, DESeq2, enrichment, and
publication figure generation.

**Figure 2. Count QC and global sample structure.** After rRNA removal, the
active matrix contained 23,279 genes, with 16,892 genes retained after low-count
filtering. VST PCA showed strong sample-level structure, with PC1 and PC2
explaining 35.1% and 28.0% of the variance, respectively.

**Figure 3. Infection DEGs within each diet.** Infection-associated DEG burden
was largest in diet 50, intermediate in diet 33, and smallest in diet 83. Most
DEGs were coding genes, with a smaller but visible long non-coding RNA
component.

**Figure 4. DEG overlap and functional enrichment.** Infection DEG sets shared
721 genes across at least two diets and 277 genes across all three diets.
Enrichment highlighted extracellular, bacterial-response, Toll/Imd, metabolism,
transport, and detoxification terms.

**Figure 5. Curated fat body candidate genes.** Curated immune, storage,
lipid-transport, and nutrition/metabolism genes showed coordinated expression
differences across infection and diet groups, including peptidoglycan-recognition
proteins, serpins, lysozyme-like genes, hemocyanin-domain genes, lipid
transporters, and amino acid transporters.
