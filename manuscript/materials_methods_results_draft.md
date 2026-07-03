# Diet infection transcriptomics manuscript draft

Authors: Emily Baker and Maeva A. Techer

This draft follows the structure and tone of the AGY manuscript model, but it is
based on the current files in `gregaria-diet-infection-interaction`.

## Author notes to confirm before submission

- Confirm the exact artificial diet formulations represented by diets 33, 50,
  and 83, and whether these labels correspond to protein:carbohydrate ratio,
  percentage, or another diet code.
- Confirm the infection protocol for the `Infected` group: fungal strain,
  inoculation method, dose, solvent/control treatment, and whether controls were
  sham-treated.
- Confirm the dissected tissue name. The website currently describes this as
  fat body RNA-seq, but the metadata and count matrix do not encode tissue.
- Confirm RNA extraction, library preparation, sequencing platform, read length,
  mapping, and gene-counting commands. The current repo starts from a
  precomputed gene-count matrix.
- Confirm whether the final manuscript should report the full dataset or the
  no-outlier sensitivity dataset as the main analysis. Current page 06 results
  use the full metadata table; page 09 defines a no-outlier table excluding
  samples 1007 and 1037.
- Confirm whether Metarhizium RNA-seq read mapping/load estimates can be added.
  The current repository does not contain those fungal-load outputs.

## Materials and Methods

### Experimental design

To test whether nutritional environment modifies the transcriptional response
to fungal challenge, female *Schistocerca gregaria* were assigned to one of
three diet treatments and one of two infection treatments. The three diet groups
are encoded in the current metadata as diets 33, 50, and 83. Infection treatment
was encoded as `Control` or `Infected`, and all samples in the current metadata
were collected 96 h post-inoculation. Age was controlled by the laboratory
design, reducing the likelihood that age-dependent expression differences drove
the diet or infection contrasts.

The current expression matrix contains 44 female samples distributed across the
six diet-by-infection groups. In the full metadata table, diet 33 contained
seven control and nine infected samples, diet 50 contained six control and nine
infected samples, and diet 83 contained seven control and six infected samples.
Two diet-50 control samples, 1007 and 1037, were also recorded in a no-outlier
metadata file for sensitivity analyses because they showed atypical global
expression structure during quality-control inspection.

### RNA-seq count matrix and sample matching

Downstream analyses began from a merged gene-level count matrix
(`data/raw_read_counts/master_counts.csv`) and sample metadata
(`data/metadata/mehreen_metadata_fixed.txt`). Count-matrix column names were
converted back to biological sample labels by removing the `mehreen_` prefix and
`_MERGE` suffix, then matched to metadata labels before any modeling. This step
ensured that diet, infection status, and count columns were aligned exactly.

The full count matrix contained 35,681 gene rows and 44 libraries. Library sizes
varied from 1.90 million to 9.49 million assigned counts, with a median of 4.20
million assigned counts. Genes with fewer than 10 total counts across all
samples were excluded before variance-stabilized exploratory plots and DESeq2
modeling, retaining 19,644 genes for the current analyses.

### Quality control and outlier handling

Quality control used library-size summaries, variance-stabilized principal
component analysis (PCA), sample-distance heatmaps, and heatmaps of the most
variable VST-transformed genes. These analyses were used to identify whether any
individual libraries separated strongly from their expected diet-by-treatment
group before interpreting differential expression.

Outlier handling was explicit because the experiment has a modest number of
replicates per group. Two diet-50 control libraries, samples 1007 and 1037, were
recorded as QC outliers in the no-outlier metadata file. Removing these samples
does not imply that age was uncontrolled or that these samples are biologically
uninteresting. Rather, the exclusion provides a sensitivity dataset to assess
whether differential expression patterns in diet 50 depend on two controls with
atypical global expression profiles. The current page 06 infected-versus-control
results were generated from the full metadata table; the no-outlier table is
available for sensitivity reporting.

### Differential expression analyses

Differential expression was performed with DESeq2 using raw gene counts. Genes
were called differentially expressed when the Benjamini-Hochberg adjusted
P-value was less than 0.05 and the absolute log2 fold-change was at least 1.
Variance-stabilized transformed counts were used for PCA and heatmaps, but not
as input to the DESeq2 models.

Three complementary DEG analyses were used. First, a full diet-by-infection
model was fitted using:

```r
~ Diet + Treatment + Diet:Treatment
```

In this model, diet 50 and control treatment were used as references. The main
infection coefficient therefore represents the infected-versus-control response
within diet 50, while the interaction coefficients test whether the infection
response differs in diets 33 or 83 relative to diet 50.

Second, pairwise diet contrasts were tested using a model that included
treatment as a covariate:

```r
~ Treatment + Diet
```

This analysis asked whether genes differed among diets after accounting for
infection status.

Third, infection effects were tested separately within each diet using:

```r
~ Treatment
```

These diet-specific models avoid averaging infection responses across
nutritional conditions and provide the clearest contrasts for asking how each
diet modifies the host transcriptional response to infection. In each contrast,
a positive log2 fold-change indicates higher expression in infected samples,
whereas a negative log2 fold-change indicates higher expression in controls.

### Functional annotation and enrichment

Functional annotation used the local eggNOG annotation file for *S. gregaria*
(`GCF_023897955.1_iqSchGreg1.2_Arthopoda_one2one.emapper.annotations`) together
with the corresponding reference GFF
(`GCF_023897955.1_iqSchGreg1.2_genomic.gff`). Because the count matrix uses
`LOC...` gene IDs and the eggNOG file is protein-query based, CDS records in the
GFF were used to bridge protein IDs to gene IDs. This lifted eggNOG annotation
onto 13,526 genes in the count matrix, including 6,533 genes with GO terms,
7,412 genes with KEGG orthology identifiers, and 4,337 genes with KEGG pathway
identifiers.

GO and KEGG enrichment were performed with `clusterProfiler::enricher()` using
custom TERM2GENE maps built from the local eggNOG annotations. DEG sets were
split by analysis source and expression direction, so genes higher in infected
samples and genes higher in controls were tested separately. The enrichment
background was restricted to genes present in the count matrix and represented
in the corresponding annotation map.

### Curated candidate heatmap

A curated VST heatmap was generated to summarize genes most relevant to the
diet-infection question. The candidate set included the strongest
infected-versus-control DEGs from each diet and direction, plus significant DEGs
whose eggNOG descriptions matched immune-related or nutrition/metabolism-related
keywords. Samples were ordered by infection status first and diet second rather
than clustered, so the figure directly shows infected samples followed by
controls, with diet substructure within each treatment. Genes were ordered by
DEG breadth across diets, then by functional theme.

## Results

### Diet and infection structure host gene expression

Variance-stabilized PCA showed that global expression variation was structured
by the diet-by-treatment design. The first two principal components explained
28.6% and 18.6% of the variance, respectively. When centroids were calculated
for each diet-by-treatment group, infected samples shifted away from their
corresponding controls in all three diets, indicating that fungal challenge was
associated with a broad host transcriptional response. The centroid shifts were
not identical among diets, supporting the use of diet-specific infection
contrasts rather than a single pooled infection model.

### The strongest infection response occurred in diet 50

Diet-specific infected-versus-control DESeq2 models recovered substantial
infection-associated transcriptional changes in all three diets, but the
magnitude of the response differed among diets. Diet 50 showed the largest
number of infection DEGs, with 1,735 significant genes, including 1,177 genes
higher in infected samples and 558 genes higher in controls. Diet 33 showed an
intermediate response, with 1,124 infection DEGs, including 646 higher in
infected samples and 478 higher in controls. Diet 83 showed the smallest
infection response, with 657 infection DEGs, including 442 higher in infected
samples and 215 higher in controls.

Across all three diets, more DEGs were higher in infected samples than in
controls. This pattern suggests that fungal challenge induced a broad activation
of host transcriptional programs, but that the amplitude of the response was
conditioned by diet. The larger DEG burden in diet 50 may reflect a stronger
immune or stress response, whereas the smaller DEG burden in diet 83 may reflect
a dampened response, reduced infection burden, or higher background variability
in that diet group. Because post-mortem fungal loads were not available in the
current repository, these alternatives should be interpreted cautiously until
Metarhizium-mapped read counts can be incorporated.

### The full interaction model supports diet-specific infection responses

The full DESeq2 model tested diet, infection, and diet-by-infection terms in a
single framework. With diet 50 as the reference, the main infection coefficient
identified 1,640 significant genes, of which 1,224 were higher in infected
samples and 416 were higher in controls. This is consistent with the
diet-specific model showing that diet 50 had the largest infected-versus-control
response.

The diet-33 interaction term identified 83 significant genes, mostly with
negative interaction coefficients (80 negative and three positive). This pattern
indicates that a subset of the infection response differed between diet 33 and
the diet-50 reference, largely in the direction of a reduced or reversed
infection effect relative to diet 50. In contrast, the diet-83 interaction term
identified only five significant genes, suggesting that the diet-83 infection
response differed less strongly from the diet-50 reference in this model, even
though the diet-specific DEG burden was lower in diet 83.

### Baseline diet effects were present but smaller than infection effects

Pairwise diet contrasts, adjusted for infection status, recovered fewer DEGs
than the infected-versus-control contrasts. The largest diet contrast was diet
33 versus diet 83, with 313 significant genes, most of which were higher in diet
83 relative to diet 33. Diet 33 versus diet 50 identified 62 DEGs, and diet 50
versus diet 83 identified 43 DEGs. These results suggest that diet alone affects
host transcription, but in the current dataset the infection-associated response
is the dominant source of DEG burden.

### Functional enrichment highlights extracellular and microbial-response terms

Functional enrichment of direction-aware DEG sets highlighted terms consistent
with immune and extracellular responses. Among genes higher in infected samples
across the infection-within-diet analysis, enriched GO terms included
extracellular region, response to bacterium, carbohydrate derivative binding,
glycosaminoglycan binding, aminoglycan metabolic process, and peptidoglycan
binding. Several of these terms were also strongly enriched among genes higher
in infected samples in diet 50, matching the large DEG burden observed in that
diet.

KEGG pathway enrichment was strongest among genes lower in infected samples or
higher in controls in several contrasts, with recurrent enrichment of broad KEGG
pathway categories such as `ko01110`, `ko01130`, and `ko01120`. Because the
current KEGG term map retains pathway IDs rather than descriptive pathway names,
these terms should be interpreted as broad metabolic annotations until the KEGG
labels are expanded in the manuscript tables.

### Curated immune and nutrition candidates concentrate in diet 50

The curated heatmap combined statistical DEG evidence with eggNOG-derived
functional descriptions. Among curated immune and nutrition/metabolism DEGs,
diet 50 again showed the strongest signal, including 14 immune-related and 16
nutrition/metabolism genes higher in infected samples, plus one
nutrition/metabolism gene higher in controls. Diet 33 contributed one
immune-related and one nutrition/metabolism DEG higher in infected samples, and
diet 83 contributed one immune-related and two nutrition/metabolism DEGs higher
in infected samples.

This curated view supports the broader DEG result: the most pronounced
infection-associated transcriptional remodeling is currently observed in diet
50, and many of the annotated candidate genes fall into immune or
metabolic/nutritional categories relevant to the diet-infection interaction.
However, because the curated set is intentionally selected for interpretability,
it should be treated as a candidate-gene summary rather than an unbiased count
of all immune or metabolic genes in the dataset.

## Basic figure text

### Figure 1. Experimental design and global transcriptomic response

Female *S. gregaria* were sampled from three diet treatments and two infection
treatments 96 h post-inoculation. PCA of variance-stabilized counts shows
global expression structure across the six diet-by-treatment groups. Points
represent individual libraries, polygons summarize the breadth of each
diet-by-treatment group, and centroids summarize group means. Arrows connect
control and infected centroids within each diet.

### Figure 2. Diet-specific infection DEGs

Infected-versus-control contrasts were tested separately within each diet using
DESeq2. Bars show the number of significant genes at adjusted P < 0.05 and
absolute log2 fold-change >= 1. Positive bars indicate genes higher in infected
samples; negative bars indicate genes higher in controls. Diet 50 showed the
largest DEG burden, followed by diet 33 and diet 83.

### Figure 3. Curated immune and nutrition DEG heatmap

The heatmap shows row-scaled VST expression for curated infection-responsive
genes. Samples are ordered by infection status and then diet, rather than by
hierarchical clustering. Genes are grouped by DEG breadth and functional theme,
highlighting immune-related and nutrition/metabolism-related candidate genes
identified from eggNOG descriptions.

### Figure 4. Functional enrichment of diet and infection DEG sets

GO and KEGG enrichment were performed on direction-aware DEG sets using local
eggNOG annotations lifted to `LOC...` gene IDs. Dot size indicates the number of
DEGs assigned to each term, and color indicates enrichment strength. Enriched GO
terms among infection-responsive genes included extracellular region, response
to bacterium, glycosaminoglycan binding, and peptidoglycan binding.
