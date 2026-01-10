# Genetics of growth rate in induced pluripotent stem cells

Phenotype data, summary statistics, and code from [Genetics of growth rate in induced pluripotent stem cells](https://www.biorxiv.org/content/10.1101/2025.07.02.662844v1)


## Phenotype data

gAUC phenotype values for FUSION, iPSCORE, GENESiPS, and HipSci lines can be found under `data/gAUC_values.tsv.gz`.


## Summary statistics

Summary data from the analyses can be found in `summary_statistics`.

### File Descriptions

#### 1. Common variant GWAS: `cvar_gwas.tsv.gz`

A tab-delimited file of common variant GWAS results.

| Column | Description |
| :--- | :--- |
| `variant_id` | SNP identifier |
| `chromosome` | SNP chromosome |
| `base_pair_location` | SNP position |
| `effect_allele` | Effect allele |
| `other_allele` | Non-effect allele |
| `effect_allele_frequency` | Effect allele frequency |
| `beta` | Beta |
| `standard_error` | Standard error |
| `p_value` | Nominal P-value |
| `n` | Sample size |

#### 2. Rare variant GWAS

##### 2a. GMMAT: `rvar_gwas__gmmat.tsv.gz`

A tab-delimited file of rare variant GWAS results using GMMAT.

| Column | Description |
| :--- | :--- |
| `gene` | Gene Ensembl identifier |
| `chromosome` | Gene chromosome |
| `base_pair_start` | Gene starting base pair position |
| `base_pair_end` | Gene ending base pair position |
| `p_value` | P-value from hybrid test |
| `n_variants` | Number of variants in mask |

##### 2b. DeepRVAT: `rvar_gwas__deeprvat.tsv.gz`

A tab-delimited file of rare variant GWAS results using DeepRVAT.

| Column | Description |
| :--- | :--- |
| `gene` | Gene Ensembl identifier |
| `chromosome` | Gene chromosome |
| `p_value` | P-value from hybrid test |

#### 3. Differential gene expression: `dge.tsv.gz`

A tab-delimited file of differential gene expression results.

| Column | Description |
| :--- | :--- |
| `gene` | Gene Ensembl identifier |
| `log2fc` | log_[2] fold-change |
| `standard_error` | Standard error |
| `p_value` | P-value |
| `n` | Sample size |


## Code availability

The following pipelines were used to perform analyses:

1. Common variant GWAS: `scripts/cvas`
2. Rare variant GWAS
    - GMMAT: `scripts/rvas`
    - [DeepRVAT](https://github.com/PMBio/deeprvat)
3. [Differential gene expression analysis](https://github.com/CollinsLabBioComp/nextflow-sc_dge)