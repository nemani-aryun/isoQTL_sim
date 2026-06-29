# isoQTL_sim: A Simulation Framework for eQTL/TWAS Power Analysis


## Overview

**isoQTL_sim** is a comprehensive simulation framework for benchmarking RNA-seq quantification methods and powering expression QTL (eQTL) and transcriptome-wide association studies (TWAS). Unlike existing simulators that use synthetic genotypes and gene-level expression, isoQTL_sim integrates **real genotype data from 1000 Genomes** with **biologically-motivated isoform-level expression modeling** to generate realistic RNA-seq data with known genetic architecture.

### Key Features

- **Real Genotype Integration**: Uses actual 1000 Genomes Phase 3 data with realistic linkage disequilibrium structure
- **Isoform-Level Modeling**: Simulates genetic regulation at the transcript level, capturing multi-isoform complexity
- **End-to-End Pipeline**: From expression simulation through quantification to QTL-ready outputs
- **Multi-Method Benchmarking**: Compare Salmon (pseudo-alignment) vs STAR+featureCounts across GENCODE versions
- **TWAS/isoTWAS Support**: Train and evaluate gene-level and isoform-level TWAS models
- **Containerized Deployment**: Reproducible analysis via Docker
- **Extensible Framework**: Easy to add new quantification methods or association tools

### Why isoQTL_sim?

**Existing eQTL simulators have critical limitations:**
1. Generate synthetic SNP dosages from multivariate normal distributions (unrealistic LD structure)
2. Simulate expression at the gene level (doesn't reflect biological regulation at isoform level)
3. Limited ability to power isoform-level TWAS studies

**IsoQTL_sim addresses these gaps by:**
- Using **real genotypes** with authentic correlation structure
- Modeling genetic effects **per isoform** with user-defined causal architecture
- Supporting both **gene-level and isoform-level** association testing and colocalization

---

## Quick Start

### Prerequisites
- R 4.3+
- Modules: tabix, samtools, STAR, salmon
- ~500GB disk space for full simulation (customizable)

### Installation

```bash
git clone https://github.com/nemani-aryun/isoQTL_sim.git
cd isoQTL_sim

# Set up directory structure
mkdir -p logs/{step1_sim,step2_salmon,step3_twas,step4_isotwas}
mkdir -p data/output
```

### Basic Usage

1. **Configure parameters** (edit `master_sim.sh`):
```bash
DIR_PROJECT="/path/to/project"
DIR_GENOTYPES="/path/to/genotypes"
GENES_PER_JOB=15
```

2. **Run simulation pipeline**:
```bash
bash scripts/master_sim.sh
```

3. **Monitor progress**:
```bash
bjobs -A  # View all jobs
tail -f logs/step1_sim/1.log
```

---

## Pipeline Overview

The pipeline consists of three main phases:

### Phase 1: Setup (Prior Scripts)
- Define genetic and RNA-seq parameters
- Select 500 European ancestry individuals from 1000 Genomes
- Prepare indexed VCF genotype files

### Phase 2: Simulation & Quantification
- **Step 1**: Simulate ground-truth isoform expression with real genotypes
- **Step 1**: Generate realistic RNA-seq reads (FASTQ)
- **Step 2**: Quantify reads using Salmon/STAR with multiple GENCODE versions

### Phase 3: Analysis (In Development)
- Train TWAS/isoTWAS prediction models
- Perform eQTL mapping (QTLtools)
- Evaluate colocalization between methods
- Calculate power across annotation versions

---

## Documentation

- **[Getting Started Guide](GETTING_STARTED.md)** - Step-by-step walkthrough
- **[Pipeline Architecture](ARCHITECTURE.md)** - Detailed flowchart and component descriptions
- **[TWAS/isoTWAS Model Training](docs/TWAS_ISOTWAS_TRAINING.md)** - How to train prediction models
- **[Output Interpretation](vignettes/03_output_interpretation.md)** - Understanding results

---

## Citation

If you use isoQTL_sim in your research, please cite:

```bibtex
@article{nemani2026,
  title={A realistic simulation tool for powering eQTL studies and TWAS},
  author={Nemani, Aryun and Head, Taylor and Bhattacharya, Arjun},
  journal={bioRxiv},
  year={2026}
}
```
