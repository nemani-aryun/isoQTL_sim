# Getting Started with isoQTL_sim

## Installation & Setup

### 1. Clone the Repository
```bash
git clone https://github.com/bhattacharya-lab/isoQTL_sim.git
cd isoQTL_sim
```

### 2. Create Required Directories
```bash
mkdir -p {logs/{step1_sim,step2_salmon,step3_twas,step4_isotwas},data/output,temp}
```

### 3. Download Reference Data

**1000 Genomes Phase 3 (PLINK2 Format)**:
```bash
# Download from:
# http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/...
```

**GENCODE Reference Files**:
```bash
# v27, v38, v45 annotations
# Download GTF and FASTA files from: https://www.gencodegenes.org/
```

**Pre-built Salmon Indices**:
```bash
# For each GENCODE version, create Salmon index:
salmon index -t gencode.v38.transcripts.fa -i gencode_v38_index
```

---

## Running Your First Simulation

### Step 1: Prepare Prior Files

First, generate the parameter space and sample list:

```bash
# Set parameter space
Rscript scripts/prior_s00_set_parameter_space.R

# Select 1000 Genomes samples
Rscript scripts/prior_s01_subset_1KG.R

# Convert and filter VCF
bash scripts/prior_s02_convert_and_filt_vcf.sh
```

### Step 2: Configure Master Script

Edit `scripts/master_sim.sh`:

```bash
# Set these paths to your system
DIR_PROJECT="/scratch/project1"
DIR_GENOTYPES="/scratch/project1"  # Can be same as DIR_PROJECT
DIR_SCRIPTS="${DIR_PROJECT}/scripts"

# Customize simulation parameters
GENES_PER_JOB=15        # Number of genes per job (adjust for your HPC)
TOTAL_GENES=12959       # Can subset to test (e.g., 100)
PARAM_ROW_READS=1       # Which read parameter set to use
THREADS=3               # For Salmon quantification
```

### Step 3: Submit to HPC Cluster

```bash
cd $DIR_PROJECT
bash scripts/master_sim.sh
```

This will submit three job arrays:
- **Step 1** (864 jobs): Expression simulation + read generation
- **Step 2** (1500 jobs): Salmon quantification (3 annotations × 500 samples)

Monitor with:
```bash
bjobs -A                          # All your jobs
bjobs -A | grep step              # Just isoQTL_sim jobs
tail -f logs/step1_sim/1.log      # Watch job 1
```

### Step 4: Check Output

After completion, verify outputs:

```bash
# Expression files (~13k files)
ls -lh files_for_analysis/expression/expr_*.RData | head

# FASTQ reads (~1000 files)
ls -lh files_for_analysis/reads/sim_*_R1.fastq.gz | wc -l

# Salmon quantifications
ls -R files_for_analysis/salmon/*/*/quant.sf | head

# QTL-ready BED files
ls -lh files_for_analysis/bed_files/*/
```

---

## Testing on Subset Data

To test the pipeline without running full analysis:

```bash
# Edit scripts/master_sim.sh:
TOTAL_GENES=100         # Test with 100 genes instead of 13k
# Rerun only gene range 1-7 (1 job)
```

Estimated times for 100 genes, 500 samples:
- Step 1: ~2 hours
- Step 2 (Salmon ×3): ~18 hours
- **Total: ~20 hours**

---


## Next Steps

1. ✅ Complete Step 1 & 2 (simulation + quantification)
2. ⏳ Step 3: Train TWAS/isoTWAS models
3. ⏳ Step 4: Run eQTL mapping & colocalization
4. ⏳ Step 5: Generate power curves and comparisons
