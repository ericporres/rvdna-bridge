# rvDNA Bridge: 23andMe → RuVector Genomic Analysis

A Python bridge that runs [Reuven Cohen's rvDNA](https://github.com/ruvnet/ruvector/tree/main/examples/dna) genomic analysis pipeline on your 23andMe raw data export. Parses 600K+ genetic markers and produces pharmacogenomic, health variant, and drug metabolism analysis — all locally, in about one second.

## What This Does

rvDNA is a Rust-based genomic analysis toolkit that processes human gene sequences in milliseconds on a CPU. It expects continuous DNA sequences (FASTA format). 23andMe exports sparse SNP genotyping data (rsid + chromosome + position + genotype). Different formats.

This bridge closes the gap. It reads your 23andMe file directly and runs rvDNA's analysis pipeline adapted for genotyping data:

| Stage | What It Does |
|-------|-------------|
| **1. Parse** | Reads the full 23andMe genotype file (v4/v5 format) |
| **2. K-mer similarity** | Genotype-based vector encoding across gene regions (HBB, TP53, BRCA1, CYP2D6, INS) |
| **3. Variant classification** | Homozygous/heterozygous/indel breakdown with het ratio |
| **4. Pharmacogenomics** | CYP2D6 and CYP2C19 star allele calling with CPIC drug recommendations |
| **5. Health variants** | APOE, BRCA1/2, TP53, MTHFR, COMT, OPRM1, and 12+ clinically significant markers |
| **6. Compound analysis** | MTHFR compound status, pain sensitivity profiling |
| **7. Report** | Comprehensive text output saved locally |

## Quick Start

```bash
# Clone the repo
git clone https://github.com/ericporres/rvdna-bridge.git
cd rvdna-bridge

# Run on your 23andMe data
python3 rvdna_bridge.py /path/to/your/genome_file.txt
```

That's it. No dependencies beyond Python 3. No pip install. No API keys. No cloud upload.

## What You Get

The output includes:

**Pharmacogenomics (mirrors rvDNA's pharma.rs):**
- CYP2D6 diplotype and metabolizer phenotype (Normal / Intermediate / Poor / Ultra-Rapid)
- CYP2C19 diplotype and metabolizer phenotype
- CPIC-level drug recommendations with dose adjustment factors

**Health variant analysis:**
- APOE genotype (Alzheimer's risk — ε2/ε3/ε4 determination)
- BRCA1/BRCA2 cancer risk markers
- TP53 p53 Pro72Arg polymorphism
- MTHFR compound status (C677T + A1298C)
- COMT Val158Met (pain/dopamine)
- OPRM1 opioid sensitivity
- SLCO1B1 statin metabolism
- CYP1A2 caffeine metabolism
- Lactase persistence
- And more (OXTR, HTR2A, ANKK1/DRD2, NQO1)

**Compound analysis:**
- MTHFR compound heterozygote assessment with supplementation guidance
- Pain sensitivity profile (COMT + OPRM1 combined)

## Save the Report

```bash
python3 rvdna_bridge.py /path/to/genome.txt > my_analysis.txt
```

## How It Works

The bridge maps 23andMe rsid numbers to the same genomic positions that rvDNA's Rust modules use for analysis:

- **Star allele calling** uses the exact variant definitions from rvDNA's `pharma.rs` — rs3892097 for CYP2D6\*4, rs4244285 for CYP2C19\*2, etc.
- **Phenotype prediction** follows the same activity-score model: sum both allele scores, classify by threshold (>2.0 = Ultra-Rapid, ≥1.0 = Normal, ≥0.5 = Intermediate, <0.5 = Poor)
- **Drug recommendations** mirror rvDNA's CPIC-based lookup table with dose adjustment factors
- **K-mer encoding** adapts rvDNA's FNV-1a rolling hash to work with sparse genotype data in gene regions

## Privacy

Your genome never leaves your machine. The script reads a local file, does math, and prints text. No network calls. No telemetry. No cloud. This is by design — your DNA is the most permanently identifying data you have.

## Requirements

- Python 3.6+
- A 23andMe raw data export (v4 or v5 format — the tab-separated `.txt` file)

To download your 23andMe data: go to [23andMe](https://www.23andme.com/) → Settings → 23andMe Data → Download Raw Data.

## Background

This bridge was built in a single [Claude Cowork](https://claude.ai) session. Claude read through rvDNA's 4,679-line Rust codebase, identified the data format gap between continuous DNA sequences and sparse SNP genotyping, and wrote the Python adapter. The full story is on [Beyond Reason](https://promptedbyeric.substack.com/).

rvDNA is part of [RuVector](https://github.com/ruvnet/ruvector) — Reuven Cohen's unified vector and graph substrate for treating intelligence as structured state.

## Disclaimer

This analysis is for **research and educational purposes only**. It is not a medical diagnosis. Do not make medical decisions based on these results without consulting a healthcare provider or genetic counselor. Pharmacogenomic results should be confirmed with clinical-grade testing before adjusting any medication.

## License

MIT
