#!/usr/bin/env python3
"""Write a small synthetic reference and genotype set for the end-to-end test.

Everything is random. It exercises the pipeline's plumbing and says nothing about
biology. Usage: make_toy_data.py OUT_DIR
"""
import gzip
import random
import subprocess
import sys
from pathlib import Path

N_GENES, N_SAMPLES, SNPS_PER_CHR = 8, 80, 1500
CHROMS = [1, 2]
SPACING = 400_000
GENE_LEN = 3000

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
rng = random.Random(42)
samples = [f"TOY{i:03d}" for i in range(1, N_SAMPLES + 1)]


def seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


gtf, fasta = [], []
for g in range(1, N_GENES + 1):
    chrom = CHROMS[g % len(CHROMS)]
    start = SPACING * ((g + 1) // 2)
    end = start + GENE_LEN
    gene = f"ENSG{g:011d}.1"
    gtf.append(f'chr{chrom}\ttoy\tgene\t{start}\t{end}\t.\t+\t.\t'
               f'gene_id "{gene}"; gene_type "protein_coding"; gene_name "G{g}";')
    for t in range(1, rng.randint(1, 3) + 1):
        tx = f"ENST{g:08d}{t:03d}.1"
        s = seq(rng.randint(1200, 2000))
        wrapped = "\n".join(s[i:i + 60] for i in range(0, len(s), 60))
        fasta.append(f">{tx}|{gene}|-|-|G{g}-{t}|G{g}|{len(s)}|protein_coding|\n{wrapped}")

(out / "toy.gtf").write_text("\n".join(gtf) + "\n")
(out / "toy.fa").write_text("\n".join(fasta) + "\n")

for chrom in CHROMS:
    path = out / f"toy_chr{chrom}.vcf"
    with open(path, "w") as handle:
        handle.write("##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        handle.write(f"##contig=<ID={chrom}>\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
        positions = sorted(rng.sample(range(1000, SPACING * (N_GENES + 2)), SNPS_PER_CHR))
        haps = [[rng.random() < 0.3 for _ in range(2)] for _ in samples]
        for pos in positions:
            # Copy the previous SNP's alleles with occasional switches to create LD.
            maf = rng.uniform(0.1, 0.5)
            haps = [[a if rng.random() > 0.15 else rng.random() < maf for a in h] for h in haps]
            ref, alt = rng.sample("ACGT", 2)
            calls = "\t".join(f"{int(h[0])}|{int(h[1])}" for h in haps)
            handle.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{calls}\n")
    subprocess.run(["bgzip", "-f", str(path)], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", f"{path}.gz"], check=True)
print(f"Toy data written to {out}")
