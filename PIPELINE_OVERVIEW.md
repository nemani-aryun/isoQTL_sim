# isoQTL_sim Pipeline Overview

## Complete Pipeline Flow

### Phase 1: Setup & Parameter Definition (Prior Scripts)

**Purpose**: Define the biological scope and statistical parameters for simulation

**Scripts**:
- `prior_s00_set_parameter_space.R` - Select genes and assign genetic parameters
- `prior_s01_subset_1KG.R` - Select 500 European ancestry individuals
- `prior_s02_convert_and_filt_vcf.sh` - Prepare indexed genotype files

**Outputs**:
- `gene_level_parameters.txt` - Per-gene causal SNP counts, heritability, correlation bounds
- `anno_selected_genes.txt` - Gene/transcript annotation with coordinates
- `1kg_eur_500_sample_ids` - Sample list for simulation
- `1KG_vcf/` - Indexed VCF files for each chromosome

**Key Decisions**:
- Which GENCODE version to use as "ground truth" (currently v38)
- How many genes to simulate (18,983)
- Genetic architecture (4 causal SNPs/gene, 50-100% shared across isoforms)

---

### Phase 2: Expression Simulation & Read Generation (Step 1)

**Purpose**: Generate ground-truth isoform expression with known genetics, then create realistic RNA-seq reads

**Script**: `step1_sim_expr_and_reads.sh` → `s01_s02_sim_modified.R`

**Part A: Expression Simulation** (`s01_s02_sim_modified.R` lines 1-250)

1. For each gene in assigned range:
   - Extract cis-SNPs (±250kb from gene) from indexed VCF using tabix
   - Select causal SNPs randomly (some shared across isoforms, some isoform-specific)
   - Generate correlated genetic effect sizes (β) across isoforms using MVN
   - Scale effects to achieve target heritability per isoform
   - Generate correlated environmental noise (ε) 
   - Compute final expression: **Y = Xβ + ε** for 500 samples

2. Save outputs per gene:
   - Expression matrix (Y): 500 samples × M isoforms
   - Effect sizes (B): Causal SNPs × M isoforms
   - Per-isoform and gene-level heritability

**Part B: Read Simulation** (`s01_s02_sim_modified.R` lines 250-end)

1. For sample index 1-500:
   - Load aggregated expression matrix (all genes)
   - Exponentiate to TPM-like abundance scale
   - Call `Rsubread::simReads()` with:
     - 75bp read length
     - 30M reads per sample
     - Paired-end reads (customizable)
   - Generate FASTQ files with technical variation

**Parallelization**:
- Job array [1-864] runs in parallel
- Job 1: Genes 1-15 + Sample 1
- Job 2: Genes 16-30 + Sample 2
- ...
- Job 500: Genes 7486-7500 + Sample 500
- Jobs 501-864: Remaining genes (no read sim)

**Outputs**:
- `expression/expr_ENSGXXXX.RData` (13k files)
- `reads/sim_sampleID_R1.fastq.gz`, `R2.fastq.gz` (1000 files for paired-end)

**Runtime**: ~2-3 hours for full simulation

---

### Phase 2b: Read Quantification (Step 2)

**Purpose**: Quantify simulated reads using multiple methods and annotation versions

**Script**: `step2_salmon_quant.sh` (Salmon pseudo-alignment)

**Methods**:
1. **Salmon pseudo-alignment** (3× for v27, v38, v45 annotations)
   - Fast, lightweight transcript quantification
   - Automatically handles multi-mapping
   - Outputs TPM and raw counts per transcript

2. **STAR + featureCounts** (3× for v27, v38, v45 annotations) — *implemented separately*
   - Splice-aware genome alignment
   - Gene-level counting
   - More conservative multi-mapping handling

**Parallelization**:
- Job array [1-500] per method-annotation combination
- 3 (Salmon versions) × 500 samples = 1500 jobs total
- Dependency: Wait for Step 1 to complete

**Outputs**:
- `salmon/param_row_reads_1_gencode_vXX/sampleID/quant.sf`
- Gene and transcript-level TPM values

**Runtime**: ~6 hours per annotation version

---

### Phase 3: QC & Aggregation (Steps 3-8)

**Purpose**: Quality control, normalize, and prepare for QTL mapping

**Steps**:
1. **FastQC** - Per-sample sequencing quality
2. **MultiQC** - Aggregate QC, filter samples
3. **Tximeta** - Import Salmon quantifications, summarize to genes
4. **TMM Normalization** - edgeR-based depth normalization
5. **Prep BED Files** - Format for QTLtools (chr, start, end, expression matrix)

**Outputs**:
- `quants/param_row_reads_1_salmon_gencode_vXX_gene.RDS` - SummarizedExperiment
- `bed_files/param_row_reads_1/salmon_gencode_vXX_gene_counts.bed` - QTL-ready phenotypes

---

### Phase 4: eQTL Mapping (Step 3) — *In Development*

**Purpose**: Discover association between simulated SNPs and expression

**Pipeline**:
1. Use QTLtools or FastQTL for cis-eQTL mapping
2. For each gene: Test SNPs within ±1Mb for association
3. Output: Beta, p-value, q-value per SNP-gene pair
4. Compare eQTL results across:
   - Annotation versions (v27, v38, v45)
   - Quantification methods (Salmon vs STAR)

**Expected Finding**: 
- Annotation choice affects 15-20% of discovered eQTLs
- Version mismatch (quantify with v27, ground truth is v38) introduces bias

---

### Phase 5: TWAS Model Training (Step 4) — *In Development*

**Purpose**: Train predictive models for gene expression as function of SNPs

**Implementation**: Adapt isotwas R package (code provided)

**For Gene-Level TWAS**:
1. For each gene:
   - Extract SNPs within ±1Mb
   - Train elastic net / BLUP / SuSiE on SNP→gene expression
   - Save weights, R², effect sizes
2. Methods: Elastic Net, BLUP, SuSiE
3. Select best by cross-validation R²

**For Isoform-Level TWAS** (novel):
1. For each gene with >1 isoform:
   - Train multivariate models on SNP→isoform expression panel
   - Methods: multivariate elastic net, multivariate SuSiE
   - Learn isoform correlation structure (Omega)
2. Compare: univariate (per-isoform) vs multivariate (joint)

**Expected Comparison**:
- Isoform-level models capture ~X% more variance than gene-level
- Colocalized across isoforms: important for functional interpretation

**Outputs**:
- `TWAS_models/geneID_TWAS.tsv.gz` - Gene weights + R²
- `isoTWAS_models/geneID_isoTWAS.RDS` - Isoform weights + correlation structure

---

### Phase 6: Colocalization & Power Analysis (Step 5) — *In Development*

**Purpose**: Assess method concordance and power

**Analyses**:
1. **eQTL Colocalization**
   - Which eQTLs discovered in one method but not another?
   - Effect size concordance across methods
   - LD-adjusted fine-mapping

2. **TWAS Colocalization**
   - Gene-level vs isoform-level TWAS signals
   - Which isoforms drive gene-level associations?

3. **Power Calculations**
   - Sensitivity: % of true causal SNPs recovered
   - Specificity: False discovery rate
   - Effect size bias

**Expected Findings**:
- GENCODE v38 (simulation truth) shows highest power
- v27 (outdated): misses some isoforms → lower sensitivity
- v45 (overcomplete): ambiguous assignments → higher FDR
- isoform-level TWAS recovers additional biology missed by gene-level

---

## Data Lineage

```text
Prior_Scripts
├── Input
│   ├── 18,983 genes (GENCODE v38)
│   ├── 500 individuals (1000G EUR)
│   └── Causal architecture
│       ├── 4 SNPs/gene
│       └── h² ~ Beta(2,20)
│
├── Step 1: Simulate Expression & Reads
│   ├── Ground-truth Y matrices (500 × M isoforms per gene)
│   ├── Causal SNP IDs & effect sizes
│   └── FASTQ reads (75 bp, 30M reads/sample)
│
├── Step 2: Quantify Reads
│   ├── Salmon (GENCODE v27, v38, v45)
│   ├── STAR + featureCounts (GENCODE v27, v38, v45)
│   └── Gene/transcript TPMs & counts
│
├── Steps 3–8: QC & Prepare
│   ├── TMM-normalized expression
│   └── BED files (QTL-ready phenotypes)
│
├── Step 3 [Dev]: eQTL Mapping
│   ├── SNP–gene associations
│   └── Comparison across methods/annotations
│
├── Step 4 [Dev]: TWAS Model Training
│   ├── Gene-level weights (R²)
│   └── Isoform-level weights + correlation structure
│
└── Step 5 [Dev]: Colocalization & Power
    ├── Method concordance
    ├── Sensitivity/specificity curves
    └── Power estimates for design studies
```
---

## Key Comparisons

### 1. Quantification Method
- **Salmon** (pseudo-alignment)
  - Pros: Fast, handles ambiguous reads probabilistically
  - Cons: Doesn't use genomic context

- **STAR+featureCounts** (genome alignment)
  - Pros: Spatially aware, can detect novel junctions
  - Cons: Slower, more stringent multi-mapping

### 2. Annotation Version
- **v27 (2017)** - Smaller, well-curated
- **v38 (2021)** - Simulation ground truth
- **v45 (2024)** - Largest, may include spurious transcripts

### 3. Expression Level
- **Gene-level** - Standard, robust, but loses isoform information
- **Isoform-level** - More biological detail, but noisier

---

## Running Full vs Subset Analysis

| Component | Full Analysis | Subset Test |
|-----------|---------------|-------------|
| Genes | 12,959 | 100 |
| Samples | 500 | 500 |
| Step 1 time | 2-3 hrs | 15 min |
| Step 2 time | 18 hrs (×3 annot) | 30 min |
| Disk space | ~500 GB | ~10 GB |
| Total time | ~24-48 hrs | ~1-2 hrs |

---

## Customization

**To modify:**
- Genes simulated: Edit `TOTAL_GENES` in `master_sim.sh`
- Causal architecture: Edit `prior_s00_set_parameter_space.R` (h², n_causal, etc.)
- Read parameters: Edit `prior_s00_set_parameter_space.R` (read_length, library_size, etc.)
- Quantification methods: Add scripts for kallisto, RSEM, etc.
