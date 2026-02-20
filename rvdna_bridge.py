#!/usr/bin/env python3
"""
rvDNA Bridge: 23andMe → RuVector Genomic Analysis Pipeline
===========================================================
Python port of Reuven's rvDNA analysis pipeline, adapted to read
23andMe genotyping data directly.

Mirrors all 7 analysis stages from rvDNA's main.rs:
  1. Load & parse 23andMe genotype data
  2. K-mer encoding & similarity analysis (on gene-region SNPs)
  3. Variant calling from genotype data
  4. Pharmacogenomics (CYP2D6 + CYP2C19 star allele calling)
  5. Health variant analysis (BRCA1, BRCA2, TP53, HBB, APOE, etc.)
  6. Trait & metabolism analysis
  7. Summary report with drug recommendations

Based on: https://github.com/ruvnet/ruvector/tree/main/examples/dna
"""

import sys
import time
import math
import json
import hashlib
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Optional
from pathlib import Path

# ═══════════════════════════════════════════════════════════════════
# Core Types (mirrors rvdna/src/types.rs)
# ═══════════════════════════════════════════════════════════════════

@dataclass
class SNP:
    rsid: str
    chromosome: str
    position: int
    genotype: str

@dataclass
class StarAllele:
    name: str
    activity_score: float
    function: str

@dataclass
class MetabolizerPhenotype:
    phenotype: str  # UltraRapid, Normal, Intermediate, Poor
    activity_score: float

@dataclass
class DrugRecommendation:
    drug: str
    gene: str
    recommendation: str
    dose_factor: float
    evidence_level: str = "CPIC Level A"

@dataclass
class HealthVariant:
    rsid: str
    gene: str
    name: str
    genotype: str
    risk_allele: str
    interpretation: str
    clinical_significance: str

@dataclass
class KmerVector:
    dimensions: int
    k: int
    values: list
    gene: str

# ═══════════════════════════════════════════════════════════════════
# Stage 1: 23andMe Parser
# ═══════════════════════════════════════════════════════════════════

def parse_23andme(filepath: str) -> dict:
    """Parse 23andMe raw data file into rsid-indexed dict."""
    snps = {}
    chr_counts = Counter()
    total = 0
    no_calls = 0

    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) != 4:
                continue
            rsid, chrom, pos, genotype = parts
            total += 1
            if genotype == '--':
                no_calls += 1
                continue
            snps[rsid] = SNP(rsid=rsid, chromosome=chrom, position=int(pos), genotype=genotype)
            chr_counts[chrom] += 1

    return {
        'snps': snps,
        'total': total,
        'no_calls': no_calls,
        'called': total - no_calls,
        'chr_counts': chr_counts,
    }


# ═══════════════════════════════════════════════════════════════════
# Stage 2: K-mer Encoding (mirrors rvdna/src/kmer.rs)
# ═══════════════════════════════════════════════════════════════════

def fnv1a_hash(data: bytes) -> int:
    """FNV-1a hash (same as rvDNA's k-mer hashing)."""
    h = 0xcbf29ce484222325
    for b in data:
        h ^= b
        h = (h * 0x100000001b3) & 0xFFFFFFFFFFFFFFFF
    return h

def genotype_to_kmer_vector(snps_in_region: list, k: int = 11, dims: int = 512) -> list:
    """
    Create a k-mer frequency vector from SNP genotypes in a gene region.
    Adapts rvDNA's rolling polynomial hash to work with sparse SNP data.
    """
    vector = [0.0] * dims

    # Build pseudo-sequence from genotypes
    seq = ''.join(snp.genotype for snp in snps_in_region)
    if len(seq) < k:
        return vector

    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k].encode()
        h = fnv1a_hash(kmer)
        vector[h % dims] += 1.0

    # Normalize to unit vector (same as rvDNA)
    magnitude = math.sqrt(sum(v*v for v in vector))
    if magnitude > 0:
        vector = [v / magnitude for v in vector]

    return vector

def cosine_similarity(a: list, b: list) -> float:
    """Cosine similarity (mirrors rvdna/src/main.rs)."""
    dot = sum(x*y for x, y in zip(a, b))
    mag_a = math.sqrt(sum(x*x for x in a))
    mag_b = math.sqrt(sum(x*x for x in b))
    if mag_a == 0 or mag_b == 0:
        return 0.0
    return dot / (mag_a * mag_b)


# ═══════════════════════════════════════════════════════════════════
# Stage 3: Variant Calling (mirrors rvdna/src/variant.rs)
# ═══════════════════════════════════════════════════════════════════

def classify_genotype(genotype: str, ref_allele: str) -> str:
    """Classify genotype (mirrors rvDNA's Genotype enum)."""
    if len(genotype) == 2:
        a1, a2 = genotype[0], genotype[1]
        if a1 == ref_allele and a2 == ref_allele:
            return "HomRef"
        elif a1 == ref_allele or a2 == ref_allele:
            return "Het"
        else:
            return "HomAlt"
    return "Unknown"


# ═══════════════════════════════════════════════════════════════════
# Stage 4: Pharmacogenomics (mirrors rvdna/src/pharma.rs exactly)
# ═══════════════════════════════════════════════════════════════════

# CYP2D6 star allele definitions (from pharma.rs)
CYP2D6_VARIANTS = {
    'rs3892097': {'allele': '*4', 'ref': 'C', 'alt': 'T', 'function': 'No function (splicing defect)', 'activity': 0.0},
    'rs35742686': {'allele': '*3', 'ref': 'T', 'alt': 'del', 'function': 'No function (frameshift)', 'activity': 0.0},
    'rs5030655': {'allele': '*6', 'ref': 'T', 'alt': 'del', 'function': 'No function (frameshift)', 'activity': 0.0},
    'rs1065852': {'allele': '*10', 'ref': 'C', 'alt': 'T', 'function': 'Decreased function', 'activity': 0.5},
    'rs28371725': {'allele': '*41', 'ref': 'C', 'alt': 'T', 'function': 'Decreased function', 'activity': 0.5},
    'rs28371706': {'allele': '*17', 'ref': 'C', 'alt': 'T', 'function': 'Decreased function', 'activity': 0.5},
}

# CYP2C19 star allele definitions (from pharma.rs)
CYP2C19_VARIANTS = {
    'rs4244285': {'allele': '*2', 'ref': 'G', 'alt': 'A', 'function': 'No function (splicing defect)', 'activity': 0.0},
    'rs4986893': {'allele': '*3', 'ref': 'G', 'alt': 'A', 'function': 'No function (premature stop)', 'activity': 0.0},
    'rs12248560': {'allele': '*17', 'ref': 'C', 'alt': 'T', 'function': 'Increased function', 'activity': 1.5},
}

def call_cyp2d6(snps: dict) -> tuple:
    """Call CYP2D6 star alleles from 23andMe genotypes (mirrors pharma.rs call_star_allele)."""
    alleles = []
    variant_details = []

    for rsid, info in CYP2D6_VARIANTS.items():
        if rsid in snps:
            gt = snps[rsid].genotype
            ref = info['ref']
            alt = info['alt']

            if alt == 'del':
                # Deletion variants: DD=homozygous deletion, DI=heterozygous, II=no deletion
                if gt == 'DD':
                    alleles.append(StarAllele(info['allele'], info['activity'], info['function']))
                    alleles.append(StarAllele(info['allele'], info['activity'], info['function']))
                    variant_details.append(f"  {rsid}: {gt} → homozygous {info['allele']} ({info['function']})")
                elif gt == 'DI' or gt == 'ID':
                    alleles.append(StarAllele(info['allele'], info['activity'], info['function']))
                    variant_details.append(f"  {rsid}: {gt} → heterozygous {info['allele']} ({info['function']})")
                else:
                    variant_details.append(f"  {rsid}: {gt} → reference (no {info['allele']})")
            else:
                has_alt = alt in gt
                if gt == alt + alt:
                    alleles.append(StarAllele(info['allele'], info['activity'], info['function']))
                    alleles.append(StarAllele(info['allele'], info['activity'], info['function']))
                    variant_details.append(f"  {rsid}: {gt} → homozygous {info['allele']} ({info['function']})")
                elif has_alt:
                    alleles.append(StarAllele(info['allele'], info['activity'], info['function']))
                    variant_details.append(f"  {rsid}: {gt} → heterozygous {info['allele']} ({info['function']})")
                else:
                    variant_details.append(f"  {rsid}: {gt} → reference (no {info['allele']})")
        else:
            variant_details.append(f"  {rsid}: not genotyped")

    # Fill remaining alleles with *1 (wild-type)
    while len(alleles) < 2:
        alleles.append(StarAllele('*1', 1.0, 'Normal function (wild-type)'))

    return alleles[:2], variant_details

def call_cyp2c19(snps: dict) -> tuple:
    """Call CYP2C19 star alleles (mirrors pharma.rs call_cyp2c19_allele)."""
    alleles = []
    variant_details = []

    for rsid, info in CYP2C19_VARIANTS.items():
        if rsid in snps:
            gt = snps[rsid].genotype
            alt = info['alt']
            ref = info['ref']

            if gt == alt + alt:
                alleles.append(StarAllele(info['allele'], info['activity'], info['function']))
                alleles.append(StarAllele(info['allele'], info['activity'], info['function']))
                variant_details.append(f"  {rsid}: {gt} → homozygous {info['allele']} ({info['function']})")
            elif alt in gt:
                alleles.append(StarAllele(info['allele'], info['activity'], info['function']))
                variant_details.append(f"  {rsid}: {gt} → heterozygous {info['allele']} ({info['function']})")
            else:
                variant_details.append(f"  {rsid}: {gt} → reference (no {info['allele']})")
        else:
            variant_details.append(f"  {rsid}: not genotyped")

    while len(alleles) < 2:
        alleles.append(StarAllele('*1', 1.0, 'Normal function'))

    return alleles[:2], variant_details

def predict_phenotype(allele1: StarAllele, allele2: StarAllele) -> MetabolizerPhenotype:
    """Predict metabolizer phenotype (mirrors pharma.rs predict_phenotype exactly)."""
    total = allele1.activity_score + allele2.activity_score
    if total > 2.0:
        return MetabolizerPhenotype("Ultra-Rapid Metabolizer", total)
    elif total >= 1.0:
        return MetabolizerPhenotype("Normal Metabolizer", total)
    elif total >= 0.5:
        return MetabolizerPhenotype("Intermediate Metabolizer", total)
    else:
        return MetabolizerPhenotype("Poor Metabolizer", total)

def get_drug_recommendations(gene: str, phenotype: MetabolizerPhenotype) -> list:
    """Get CPIC drug recommendations (mirrors pharma.rs get_recommendations)."""
    recs = []

    if gene == "CYP2D6":
        if phenotype.phenotype == "Poor Metabolizer":
            recs = [
                DrugRecommendation("Codeine", gene, "AVOID codeine; no conversion to morphine. Use alternative analgesic.", 0.0),
                DrugRecommendation("Tramadol", gene, "AVOID tramadol; reduced efficacy. Use alternative analgesic.", 0.0),
                DrugRecommendation("Tamoxifen", gene, "Consider alternative endocrine therapy (aromatase inhibitor).", 0.0),
                DrugRecommendation("Ondansetron", gene, "Use standard dose; may have increased exposure.", 0.75),
            ]
        elif phenotype.phenotype == "Ultra-Rapid Metabolizer":
            recs = [
                DrugRecommendation("Codeine", gene, "AVOID codeine; risk of fatal toxicity from ultra-rapid morphine conversion.", 0.0),
                DrugRecommendation("Tramadol", gene, "AVOID tramadol; risk of respiratory depression.", 0.0),
            ]
        elif phenotype.phenotype == "Intermediate Metabolizer":
            recs = [
                DrugRecommendation("Codeine", gene, "Use lower dose or alternative analgesic.", 0.5),
                DrugRecommendation("Tamoxifen", gene, "Consider higher dose or alternative therapy.", 0.75),
            ]
        else:
            recs = [DrugRecommendation("Standard", gene, "Use standard dosing per labeling.", 1.0)]

    elif gene == "CYP2C19":
        if phenotype.phenotype == "Poor Metabolizer":
            recs = [
                DrugRecommendation("Clopidogrel (Plavix)", gene, "AVOID clopidogrel; use prasugrel or ticagrelor instead.", 0.0, "CPIC Level A"),
                DrugRecommendation("Voriconazole", gene, "Reduce dose by 50%; monitor for toxicity.", 0.5),
                DrugRecommendation("PPIs (omeprazole)", gene, "Reduce dose; slower clearance increases exposure.", 0.5),
                DrugRecommendation("Escitalopram", gene, "Consider 50% dose reduction.", 0.5),
            ]
        elif phenotype.phenotype == "Intermediate Metabolizer":
            recs = [
                DrugRecommendation("Clopidogrel (Plavix)", gene, "Consider alternative antiplatelet or increased dose.", 1.5, "CPIC Level A"),
                DrugRecommendation("PPIs (omeprazole)", gene, "Standard dose likely adequate; may have slightly increased exposure.", 1.0),
                DrugRecommendation("Escitalopram", gene, "Use standard dose; monitor response.", 1.0),
            ]
        elif phenotype.phenotype == "Ultra-Rapid Metabolizer":
            recs = [
                DrugRecommendation("Clopidogrel (Plavix)", gene, "Standard dosing (enhanced activation is beneficial).", 1.0),
                DrugRecommendation("Omeprazole", gene, "Increase dose; rapid clearance reduces efficacy.", 2.0),
                DrugRecommendation("Voriconazole", gene, "Use alternative antifungal.", 0.0),
            ]
        else:
            recs = [DrugRecommendation("Standard", gene, "Use standard dosing per labeling.", 1.0)]

    return recs


# ═══════════════════════════════════════════════════════════════════
# Stage 5: Health Variant Analysis
# ═══════════════════════════════════════════════════════════════════

HEALTH_VARIANTS = {
    # APOE (Alzheimer's risk)
    'rs429358': {
        'gene': 'APOE', 'name': 'APOE ε4 determinant',
        'risk_allele': 'C',
        'interpretations': {
            'TT': ('APOE ε3/ε3 or ε2/ε3 (depends on rs7412)', 'Protective/Normal'),
            'CT': ('One ε4 allele present', 'Increased Alzheimer\'s risk (~3x)'),
            'CC': ('Two ε4 alleles present', 'Significantly increased Alzheimer\'s risk (~12x)'),
        }
    },
    'rs7412': {
        'gene': 'APOE', 'name': 'APOE ε2 determinant',
        'risk_allele': 'T',
        'interpretations': {
            'CC': ('No ε2 allele', 'Normal'),
            'CT': ('One ε2 allele present', 'Protective - reduced Alzheimer\'s risk, slightly increased cardiovascular benefit'),
            'TT': ('Two ε2 alleles (ε2/ε2)', 'Protective for Alzheimer\'s; monitor lipids (rare hyperlipoproteinemia III risk)'),
        }
    },
    # TP53 (cancer)
    'rs1042522': {
        'gene': 'TP53', 'name': 'p53 Pro72Arg (R72P)',
        'risk_allele': 'G',
        'interpretations': {
            'CC': ('Pro/Pro homozygous', 'Normal apoptosis; slightly increased cancer survival in some contexts'),
            'CG': ('Pro/Arg heterozygous', 'Mixed - Arg allele has stronger apoptotic activity'),
            'GG': ('Arg/Arg homozygous', 'Stronger apoptotic response; variable cancer risk associations'),
        }
    },
    # BRCA1
    'rs80357906': {
        'gene': 'BRCA1', 'name': 'BRCA1 5382insC (Ashkenazi founder)',
        'risk_allele': 'I',
        'interpretations': {
            'DD': ('No insertion detected', 'Normal - no BRCA1 5382insC mutation'),
            'DI': ('Heterozygous carrier', 'INCREASED breast/ovarian cancer risk - genetic counseling recommended'),
            'II': ('Homozygous insertion', 'HIGH breast/ovarian cancer risk - urgent genetic counseling'),
        }
    },
    'rs28897696': {
        'gene': 'BRCA1', 'name': 'BRCA1 missense variant',
        'risk_allele': 'A',
        'interpretations': {
            'GG': ('Reference genotype', 'Normal'),
            'AG': ('Heterozygous', 'Variant of uncertain significance - consult genetic counselor'),
            'AA': ('Homozygous variant', 'Consult genetic counselor'),
        }
    },
    # BRCA2
    'rs11571833': {
        'gene': 'BRCA2', 'name': 'BRCA2 K3326X',
        'risk_allele': 'T',
        'interpretations': {
            'AA': ('Reference genotype', 'Normal'),
            'AT': ('Heterozygous', 'Modestly increased cancer risk (OR ~1.3)'),
            'TT': ('Homozygous variant', 'Increased cancer risk - genetic counseling recommended'),
        }
    },
    # MTHFR (folate metabolism)
    'rs1801133': {
        'gene': 'MTHFR', 'name': 'C677T',
        'risk_allele': 'A',
        'interpretations': {
            'GG': ('CC genotype (normal)', 'Normal MTHFR enzyme activity (100%)'),
            'AG': ('CT heterozygous', 'Reduced enzyme activity (~65%). Consider methylfolate supplementation.'),
            'AA': ('TT homozygous', 'Significantly reduced activity (~30%). Methylfolate recommended; monitor homocysteine.'),
        }
    },
    'rs1801131': {
        'gene': 'MTHFR', 'name': 'A1298C',
        'risk_allele': 'T',
        'interpretations': {
            'GG': ('CC homozygous variant', 'Reduced enzyme activity'),
            'GT': ('AC heterozygous', 'Mildly reduced enzyme activity'),
            'TT': ('AA reference', 'Normal MTHFR activity at this position'),
        }
    },
    # COMT (dopamine/pain)
    'rs4680': {
        'gene': 'COMT', 'name': 'Val158Met',
        'risk_allele': 'A',
        'interpretations': {
            'GG': ('Val/Val', 'Higher COMT activity → lower dopamine. Better stress resilience, lower pain sensitivity.'),
            'AG': ('Val/Met heterozygous', 'Intermediate COMT activity. Balanced dopamine metabolism.'),
            'AA': ('Met/Met', 'Lower COMT activity → higher dopamine. Higher pain sensitivity, better cognitive performance under low stress.'),
        }
    },
    # OPRM1 (opioid receptor)
    'rs1799971': {
        'gene': 'OPRM1', 'name': 'A118G (Asn40Asp)',
        'risk_allele': 'G',
        'interpretations': {
            'AA': ('Asn/Asn', 'Normal opioid sensitivity'),
            'AG': ('Asn/Asp heterozygous', 'Reduced opioid sensitivity; may need higher doses for pain management.'),
            'GG': ('Asp/Asp', 'Significantly reduced opioid sensitivity; may require alternative pain management.'),
        }
    },
    # CYP1A2 (caffeine)
    'rs762551': {
        'gene': 'CYP1A2', 'name': 'Caffeine metabolism',
        'risk_allele': 'C',
        'interpretations': {
            'AA': ('Fast metabolizer', 'Rapid caffeine clearance. Coffee associated with REDUCED heart disease risk.'),
            'AC': ('Intermediate', 'Moderate caffeine clearance. Moderate coffee intake recommended.'),
            'CC': ('Slow metabolizer', 'Slow caffeine clearance. Excess coffee may INCREASE heart disease risk.'),
        }
    },
    # Lactose tolerance
    'rs4988235': {
        'gene': 'MCM6/LCT', 'name': 'Lactase persistence (European)',
        'risk_allele': 'G',
        'interpretations': {
            'AA': ('Lactase persistent', 'Likely lactose TOLERANT into adulthood'),
            'AG': ('Heterozygous', 'Likely lactose tolerant (persistence is dominant)'),
            'GG': ('Lactase non-persistent', 'Likely lactose INTOLERANT in adulthood'),
        }
    },
    # Oxytocin receptor
    'rs53576': {
        'gene': 'OXTR', 'name': 'Oxytocin receptor',
        'risk_allele': 'A',
        'interpretations': {
            'GG': ('GG genotype', 'Higher empathy scores; better social cognition; more optimistic outlook.'),
            'AG': ('AG heterozygous', 'Intermediate empathy and social cognition.'),
            'AA': ('AA genotype', 'May have lower empathy scores; potentially more resilient to social stress.'),
        }
    },
    # Serotonin receptor
    'rs6311': {
        'gene': 'HTR2A', 'name': 'Serotonin 2A receptor (-1438G/A)',
        'risk_allele': 'T',
        'interpretations': {
            'CC': ('GG genotype', 'Normal serotonin receptor expression'),
            'CT': ('GA heterozygous', 'Slightly altered serotonin signaling'),
            'TT': ('AA genotype', 'Altered serotonin receptor density; may affect SSRI response'),
        }
    },
    # DRD2/ANKK1 (dopamine)
    'rs1800497': {
        'gene': 'ANKK1/DRD2', 'name': 'Taq1A (dopamine receptor)',
        'risk_allele': 'A',
        'interpretations': {
            'GG': ('A2/A2', 'Normal dopamine receptor density'),
            'AG': ('A1/A2 heterozygous', 'Reduced D2 receptor density (~30% less). Associated with reward-seeking behavior.'),
            'AA': ('A1/A1', 'Significantly reduced D2 receptor density. Higher risk of addictive behaviors.'),
        }
    },
    # SLCO1B1 (statin metabolism)
    'rs4363657': {
        'gene': 'SLCO1B1', 'name': 'Statin transporter',
        'risk_allele': 'C',
        'interpretations': {
            'TT': ('Reference', 'Normal statin metabolism. Standard dosing.'),
            'CT': ('Heterozygous', 'Increased risk of statin myopathy (~4.5x). Consider lower statin dose.'),
            'CC': ('Homozygous variant', 'High risk of statin myopathy (~17x). Use lowest effective dose or alternative.'),
        }
    },
    # NQO1
    'rs1800566': {
        'gene': 'NQO1', 'name': 'Pro187Ser (oxidative stress)',
        'risk_allele': 'T',
        'interpretations': {
            'CC': ('Pro/Pro (reference)', 'Normal NQO1 enzyme activity'),
            'CT': ('Pro/Ser heterozygous', 'Reduced NQO1 activity (~3x lower). Impaired detoxification.'),
            'TT': ('Ser/Ser', 'No NQO1 activity. Significantly impaired quinone detoxification.'),
        }
    },
}

def analyze_health_variants(snps: dict) -> list:
    """Analyze health-relevant variants from 23andMe data."""
    results = []

    for rsid, info in HEALTH_VARIANTS.items():
        if rsid in snps:
            gt = snps[rsid].genotype
            interps = info['interpretations']

            if gt in interps:
                desc, significance = interps[gt]
            else:
                desc = f"Genotype {gt} - not in standard interpretation table"
                significance = "Consult genetic counselor"

            results.append(HealthVariant(
                rsid=rsid,
                gene=info['gene'],
                name=info['name'],
                genotype=gt,
                risk_allele=info['risk_allele'],
                interpretation=desc,
                clinical_significance=significance,
            ))

    return results

def determine_apoe_genotype(snps: dict) -> str:
    """Determine APOE genotype from rs429358 + rs7412 combination."""
    rs429358 = snps.get('rs429358', None)
    rs7412 = snps.get('rs7412', None)

    if not rs429358 or not rs7412:
        return "Unable to determine (missing data)"

    gt1 = rs429358.genotype  # ε4 determinant: C = ε4
    gt2 = rs7412.genotype    # ε2 determinant: T = ε2

    # Determine alleles
    # ε2: rs429358=T, rs7412=T
    # ε3: rs429358=T, rs7412=C  (most common)
    # ε4: rs429358=C, rs7412=C

    # Count ε4 alleles (C at rs429358)
    e4_count = gt1.count('C')
    # Count ε2 alleles (T at rs7412)
    e2_count = gt2.count('T')

    if e4_count == 0 and e2_count == 0:
        return "ε3/ε3 (most common, baseline risk)"
    elif e4_count == 0 and e2_count == 1:
        return "ε2/ε3 (PROTECTIVE - reduced Alzheimer's risk)"
    elif e4_count == 0 and e2_count == 2:
        return "ε2/ε2 (protective; monitor for type III hyperlipoproteinemia)"
    elif e4_count == 1 and e2_count == 0:
        return "ε3/ε4 (increased Alzheimer's risk ~3x)"
    elif e4_count == 1 and e2_count == 1:
        return "ε2/ε4 (mixed - ε2 partially offsets ε4 risk)"
    elif e4_count == 2:
        return "ε4/ε4 (significantly increased Alzheimer's risk ~12x)"
    else:
        return f"Unusual combination: rs429358={gt1}, rs7412={gt2}"


# ═══════════════════════════════════════════════════════════════════
# Stage 6: MTHFR Compound Analysis
# ═══════════════════════════════════════════════════════════════════

def analyze_mthfr_compound(snps: dict) -> str:
    """Analyze compound MTHFR status from both C677T and A1298C."""
    c677t = snps.get('rs1801133', None)
    a1298c = snps.get('rs1801131', None)

    if not c677t or not a1298c:
        return "Incomplete MTHFR data"

    c677t_gt = c677t.genotype
    a1298c_gt = a1298c.genotype

    # Risk scoring
    c677t_risk = {'GG': 0, 'AG': 1, 'AA': 2}.get(c677t_gt, 0)
    a1298c_risk = {'TT': 0, 'GT': 1, 'GG': 2}.get(a1298c_gt, 0)

    compound = c677t_risk + a1298c_risk

    if compound == 0:
        return "Normal MTHFR function. No supplementation needed."
    elif compound == 1:
        return "Mildly reduced MTHFR. Consider methylfolate if homocysteine elevated."
    elif compound == 2:
        return "Moderately reduced MTHFR. Methylfolate (L-5-MTHF) recommended. Monitor homocysteine."
    elif compound == 3:
        return "Significantly reduced MTHFR (compound heterozygote). Methylfolate strongly recommended."
    else:
        return "Severely reduced MTHFR. Methylfolate essential. Regular homocysteine monitoring."


# ═══════════════════════════════════════════════════════════════════
# Stage 7: Report Generation
# ═══════════════════════════════════════════════════════════════════

def generate_report(filepath: str) -> str:
    """Run full rvDNA-style pipeline and generate comprehensive report."""

    lines = []
    p = lambda x: lines.append(x)

    total_start = time.time()

    p("=" * 80)
    p("  rvDNA Bridge: 23andMe → RuVector Genomic Analysis Pipeline")
    p("  Based on: https://github.com/ruvnet/ruvector/examples/dna")
    p("=" * 80)

    # ── Stage 1: Parse 23andMe ──
    p("\n━━━ Stage 1: Loading 23andMe Genotype Data ━━━")
    stage1_start = time.time()

    data = parse_23andme(filepath)
    snps = data['snps']

    p(f"  File: {Path(filepath).name}")
    p(f"  Total markers:  {data['total']:,}")
    p(f"  Called:          {data['called']:,}")
    p(f"  No-calls:       {data['no_calls']:,}")
    p(f"  Call rate:       {data['called']/data['total']*100:.1f}%")
    p(f"  Chromosomes:    {len(data['chr_counts'])}")

    # Chromosome distribution
    chr_order = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']
    p(f"\n  Chromosome distribution:")
    for c in chr_order:
        if c in data['chr_counts']:
            bar = '█' * (data['chr_counts'][c] // 1500)
            p(f"    Chr {c:>2}: {data['chr_counts'][c]:>6,} {bar}")

    p(f"  Parse time: {time.time() - stage1_start:.3f}s")

    # ── Stage 2: K-mer Similarity ──
    p("\n━━━ Stage 2: K-mer Similarity Analysis (Gene Regions) ━━━")
    stage2_start = time.time()

    # Extract SNPs by chromosome for gene regions
    gene_regions = {
        'HBB': ('11', 5225464, 5229395),      # HBB on chr11
        'TP53': ('17', 7571720, 7590868),      # TP53 on chr17
        'BRCA1': ('17', 41196312, 41277500),   # BRCA1 on chr17
        'CYP2D6': ('22', 42522500, 42528000),  # CYP2D6 on chr22
        'INS': ('11', 2159779, 2161341),       # INS on chr11
    }

    gene_snps = {}
    gene_vectors = {}
    for gene, (chrom, start, end) in gene_regions.items():
        region_snps = [s for s in snps.values()
                       if s.chromosome == chrom and start <= s.position <= end]
        gene_snps[gene] = region_snps
        gene_vectors[gene] = genotype_to_kmer_vector(region_snps, k=11, dims=512)
        p(f"  {gene:8s}: {len(region_snps):>4} SNPs in region (chr{chrom}:{start:,}-{end:,})")

    # Similarity matrix
    genes = list(gene_vectors.keys())
    p(f"\n  K-mer similarity matrix (cosine, k=11, d=512):")
    for i in range(len(genes)):
        for j in range(i+1, len(genes)):
            sim = cosine_similarity(gene_vectors[genes[i]], gene_vectors[genes[j]])
            p(f"    {genes[i]:8s} vs {genes[j]:8s}: {sim:.4f}")

    p(f"  K-mer encoding time: {time.time() - stage2_start:.3f}s")

    # ── Stage 3: Variant Summary ──
    p("\n━━━ Stage 3: Variant Classification Summary ━━━")
    stage3_start = time.time()

    # Classify all genotypes
    gt_types = Counter()
    het_count = 0
    for snp in snps.values():
        gt = snp.genotype
        if len(gt) == 2:
            if gt[0] == gt[1]:
                gt_types['Homozygous'] += 1
            else:
                gt_types['Heterozygous'] += 1
                het_count += 1
        elif gt in ('DD', 'II', 'DI', 'ID'):
            gt_types['Indel'] += 1

    p(f"  Homozygous:    {gt_types['Homozygous']:>8,}")
    p(f"  Heterozygous:  {gt_types['Heterozygous']:>8,}")
    p(f"  Indels:        {gt_types['Indel']:>8,}")
    p(f"  Het ratio:     {het_count/data['called']*100:.1f}% (typical: 25-35%)")
    p(f"  Classification time: {time.time() - stage3_start:.3f}s")

    # ── Stage 4: Pharmacogenomics ──
    p("\n━━━ Stage 4: Pharmacogenomic Analysis ━━━")
    p("  (Mirrors rvDNA pharma.rs — CYP2D6 + CYP2C19 star allele calling)")
    stage4_start = time.time()

    # CYP2D6
    p("\n  ┌─ CYP2D6 (Drug Metabolism Enzyme) ─────────────────────────────")
    cyp2d6_alleles, cyp2d6_details = call_cyp2d6(snps)
    for d in cyp2d6_details:
        p(d)

    cyp2d6_pheno = predict_phenotype(cyp2d6_alleles[0], cyp2d6_alleles[1])
    p(f"\n  Diplotype:  {cyp2d6_alleles[0].name}/{cyp2d6_alleles[1].name}")
    p(f"  Activity:   {cyp2d6_pheno.activity_score:.1f}")
    p(f"  Phenotype:  ★ {cyp2d6_pheno.phenotype} ★")

    cyp2d6_recs = get_drug_recommendations("CYP2D6", cyp2d6_pheno)
    if cyp2d6_recs:
        p(f"\n  Drug Recommendations (CPIC):")
        for rec in cyp2d6_recs:
            dose_str = f"{rec.dose_factor:.0%}" if rec.dose_factor > 0 else "AVOID"
            p(f"    • {rec.drug}: {rec.recommendation}")
            p(f"      Dose adjustment: {dose_str}")

    # CYP2C19
    p("\n  ┌─ CYP2C19 (Drug Metabolism Enzyme) ────────────────────────────")
    cyp2c19_alleles, cyp2c19_details = call_cyp2c19(snps)
    for d in cyp2c19_details:
        p(d)

    cyp2c19_pheno = predict_phenotype(cyp2c19_alleles[0], cyp2c19_alleles[1])
    p(f"\n  Diplotype:  {cyp2c19_alleles[0].name}/{cyp2c19_alleles[1].name}")
    p(f"  Activity:   {cyp2c19_pheno.activity_score:.1f}")
    p(f"  Phenotype:  ★ {cyp2c19_pheno.phenotype} ★")

    cyp2c19_recs = get_drug_recommendations("CYP2C19", cyp2c19_pheno)
    if cyp2c19_recs:
        p(f"\n  Drug Recommendations (CPIC):")
        for rec in cyp2c19_recs:
            dose_str = f"{rec.dose_factor:.0%}" if rec.dose_factor > 0 else "AVOID"
            p(f"    • {rec.drug}: {rec.recommendation}")
            p(f"      Dose adjustment: {dose_str}")

    p(f"\n  Pharmacogenomics time: {time.time() - stage4_start:.3f}s")

    # ── Stage 5: Health Variants ──
    p("\n━━━ Stage 5: Health Variant Analysis ━━━")
    stage5_start = time.time()

    # APOE first (special handling)
    p("\n  ┌─ APOE Genotype (Alzheimer's Risk) ────────────────────────────")
    apoe = determine_apoe_genotype(snps)
    p(f"  rs429358: {snps.get('rs429358', SNP('', '', 0, '??')).genotype}  rs7412: {snps.get('rs7412', SNP('', '', 0, '??')).genotype}")
    p(f"  APOE Status: ★ {apoe} ★")

    # All health variants
    health_results = analyze_health_variants(snps)

    # Group by category
    categories = {
        'Cancer Risk': ['TP53', 'BRCA1', 'BRCA2', 'NQO1'],
        'Cardiovascular': ['SLCO1B1'],
        'Neurological': ['APOE', 'COMT', 'OPRM1', 'OXTR', 'HTR2A', 'ANKK1/DRD2'],
        'Metabolism': ['MTHFR', 'CYP1A2', 'MCM6/LCT'],
    }

    for category, genes in categories.items():
        cat_results = [r for r in health_results if r.gene in genes]
        if cat_results:
            p(f"\n  ┌─ {category} ────────────────────────────────────────")
            for r in cat_results:
                p(f"  {r.rsid} ({r.gene} - {r.name})")
                p(f"    Genotype:       {r.genotype}")
                p(f"    Interpretation: {r.interpretation}")
                p(f"    Significance:   {r.clinical_significance}")

    p(f"\n  Health variant analysis time: {time.time() - stage5_start:.3f}s")

    # ── Stage 6: Compound Analysis ──
    p("\n━━━ Stage 6: Compound Genotype Analysis ━━━")
    stage6_start = time.time()

    p("\n  ┌─ MTHFR Compound Status ───────────────────────────────────────")
    mthfr_compound = analyze_mthfr_compound(snps)
    p(f"  C677T (rs1801133): {snps.get('rs1801133', SNP('','',0,'??')).genotype}")
    p(f"  A1298C (rs1801131): {snps.get('rs1801131', SNP('','',0,'??')).genotype}")
    p(f"  Assessment: {mthfr_compound}")

    # Pain sensitivity profile
    p("\n  ┌─ Pain Sensitivity Profile ────────────────────────────────────")
    comt = snps.get('rs4680', None)
    oprm1 = snps.get('rs1799971', None)
    if comt and oprm1:
        comt_gt = comt.genotype
        oprm1_gt = oprm1.genotype

        pain_score = 0
        if comt_gt == 'AA': pain_score += 2  # Met/Met = higher pain
        elif comt_gt == 'AG': pain_score += 1
        if oprm1_gt == 'GG': pain_score += 2  # reduced opioid sensitivity
        elif oprm1_gt == 'AG': pain_score += 1

        pain_labels = {0: 'Low', 1: 'Low-Moderate', 2: 'Moderate', 3: 'Moderate-High', 4: 'High'}
        p(f"  COMT (rs4680):  {comt_gt} → {'Higher' if 'A' in comt_gt else 'Lower'} pain sensitivity")
        p(f"  OPRM1 (rs1799971): {oprm1_gt} → {'Reduced' if 'G' in oprm1_gt else 'Normal'} opioid response")
        p(f"  Combined pain sensitivity: {pain_labels.get(pain_score, 'Unknown')}")
        if pain_score >= 2:
            p(f"  Note: May need higher opioid doses or alternative pain management strategies.")

    p(f"\n  Compound analysis time: {time.time() - stage6_start:.3f}s")

    # ── Summary ──
    total_time = time.time() - total_start
    p("\n" + "=" * 80)
    p("  PIPELINE SUMMARY")
    p("=" * 80)
    p(f"  Markers analyzed:     {data['called']:,}")
    p(f"  Pharmacogenes:        CYP2D6 ({cyp2d6_pheno.phenotype}), CYP2C19 ({cyp2c19_pheno.phenotype})")
    p(f"  APOE status:          {apoe}")
    p(f"  Health variants:      {len(health_results)} analyzed")
    p(f"  Drug recommendations: {len(cyp2d6_recs) + len(cyp2c19_recs)} generated")
    p(f"  Total pipeline time:  {total_time:.3f}s")
    p("")
    p("  ⚠️  DISCLAIMER: This analysis is for RESEARCH/EDUCATIONAL purposes only.")
    p("  It is NOT a medical diagnosis. Consult a healthcare provider or genetic")
    p("  counselor before making any medical decisions based on these results.")
    p("=" * 80)

    return '\n'.join(lines)


# ═══════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    filepath = sys.argv[1] if len(sys.argv) > 1 else None
    if not filepath:
        print("Usage: python rvdna_bridge.py <23andme_file.txt>")
        sys.exit(1)

    report = generate_report(filepath)
    print(report)
