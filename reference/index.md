# Package index

## Core workflow

Inspect a new dataset, run guided exploration, and review its state.

- [`explore_genomes()`](https://thierrygosselin.github.io/radr/reference/explore_genomes.md)
  : Explore and filter genomic data
- [`detect_ibm()`](https://thierrygosselin.github.io/radr/reference/detect_ibm.md)
  : Detect identity-by-missingness structure
- [`summarise_genomic_data()`](https://thierrygosselin.github.io/radr/reference/summarise_genomic_data.md)
  : Summarise genomic data
- [`radr_dependencies()`](https://thierrygosselin.github.io/radr/reference/radr_dependencies.md)
  : Check radr dependencies

## Diagnose samples and markers

Investigate missingness, duplicates, mixed samples, marker behaviour,
and data origin before filtering.

- [`detect_all_missing()`](https://thierrygosselin.github.io/radr/reference/detect_all_missing.md)
  : Detect markers with all missing genotypes
- [`detect_allele_problems()`](https://thierrygosselin.github.io/radr/reference/detect_allele_problems.md)
  : Detect alternate allele problems
- [`detect_biallelic_problems()`](https://thierrygosselin.github.io/radr/reference/detect_biallelic_problems.md)
  : Detect biallelic problems
- [`detect_duplicate_genomes()`](https://thierrygosselin.github.io/radr/reference/detect_duplicate_genomes.md)
  : Compute pairwise genome similarity or distance between individuals
  to highligh potential duplicate individuals
- [`detect_het_outliers()`](https://thierrygosselin.github.io/radr/reference/detect_het_outliers.md)
  : Detect heterozygotes outliers and estimate miscall rate
- [`detect_inversions()`](https://thierrygosselin.github.io/radr/reference/detect_inversions.md)
  : Detect candidate inversion-associated genomic regions
- [`genome_scan_context()`](https://thierrygosselin.github.io/radr/reference/genome_scan_context.md)
  : Summarise genomic context around genome-scan signals
- [`detect_mixed_genomes()`](https://thierrygosselin.github.io/radr/reference/detect_mixed_genomes.md)
  : Detect mixed genomes
- [`detect_paralogs()`](https://thierrygosselin.github.io/radr/reference/detect_paralogs.md)
  : Detect paralogs
- [`detect_ref_genome()`](https://thierrygosselin.github.io/radr/reference/detect_ref_genome.md)
  : Detect whether a dataset is reference-guided or de novo assembled

## Filter samples and markers

Apply explicit quality-control and marker-selection decisions to GDS
data.

- [`filter_individuals()`](https://thierrygosselin.github.io/radr/reference/filter_individuals.md)
  : Filter individuals based on genotyping/missingness rate,
  heterozygosity and total coverage
- [`filter_genotyping()`](https://thierrygosselin.github.io/radr/reference/filter_genotyping.md)
  : Filter markers based on genotyping / missing rate
- [`filter_coverage()`](https://thierrygosselin.github.io/radr/reference/filter_coverage.md)
  : Filter markers mean coverage
- [`filter_ma()`](https://thierrygosselin.github.io/radr/reference/filter_ma.md)
  : MAC, MAF and MAD filter
- [`filter_het()`](https://thierrygosselin.github.io/radr/reference/filter_het.md)
  : Heterozygosity filter
- [`filter_fis()`](https://thierrygosselin.github.io/radr/reference/filter_fis.md)
  : Fis filter
- [`filter_hwe()`](https://thierrygosselin.github.io/radr/reference/filter_hwe.md)
  : Filter markers based on Hardy-Weinberg Equilibrium
- [`filter_ld()`](https://thierrygosselin.github.io/radr/reference/filter_ld.md)
  : GBS/RADseq short and long distance linkage disequilibrium pruning
- [`filter_monomorphic()`](https://thierrygosselin.github.io/radr/reference/filter_monomorphic.md)
  : Filter monomorphic markers
- [`filter_common_markers()`](https://thierrygosselin.github.io/radr/reference/filter_common_markers.md)
  : Filter common markers between strata
- [`filter_dart_reproducibility()`](https://thierrygosselin.github.io/radr/reference/filter_dart_reproducibility.md)
  : Filter data based on DArT reproducibility statistics
- [`filter_snp_number()`](https://thierrygosselin.github.io/radr/reference/filter_snp_number.md)
  : Filter SNP number per locus/read
- [`filter_snp_position_read()`](https://thierrygosselin.github.io/radr/reference/filter_snp_position_read.md)
  : Filter markers/SNP based on their position on the read

## Whitelists and blacklists

Retain or exclude selected samples, markers, and individual genotypes.

- [`filter_whitelist()`](https://thierrygosselin.github.io/radr/reference/filter_whitelist.md)
  : Filter dataset with whitelist of markers
- [`filter_blacklist_genotypes()`](https://thierrygosselin.github.io/radr/reference/filter_blacklist_genotypes.md)
  : Filter dataset with blacklist of genotypes
- [`read_blacklist_genotypes()`](https://thierrygosselin.github.io/radr/reference/read_blacklist_genotypes.md)
  : read blacklist of genotypes

## Filter VCF files with bcftools

Apply selected filters directly to VCF files before or outside the GDS
workflow.

- [`filter_genotyping_vcf()`](https://thierrygosselin.github.io/radr/reference/filter_genotyping_vcf.md)
  : Filter SNPs in a VCF based on genotyping / missing rate (bcftools)
- [`filter_mac_vcf()`](https://thierrygosselin.github.io/radr/reference/filter_mac_vcf.md)
  : Filter low-MAC variants in a VCF using bcftools (AC-based)
- [`filter_monomorphic_vcf()`](https://thierrygosselin.github.io/radr/reference/filter_monomorphic_vcf.md)
  : Filter monomorphic SNPs in a VCF using bcftools (AC/AN-based)

## Population-genetic summaries

Calculate diversity, differentiation, relatedness, and private
variation.

- [`allele_frequencies()`](https://thierrygosselin.github.io/radr/reference/allele_frequencies.md)
  : Compute allele frequencies per markers and populations
- [`beta_estimator()`](https://thierrygosselin.github.io/radr/reference/beta_estimator.md)
  : Estimate population-specific allelic FST
- [`ibdg_fh()`](https://thierrygosselin.github.io/radr/reference/ibdg_fh.md)
  : FH measure of IBDg
- [`pi()`](https://thierrygosselin.github.io/radr/reference/pi.md) :
  Nucleotide diversity
- [`private_alleles()`](https://thierrygosselin.github.io/radr/reference/private_alleles.md)
  : Find private alleles
- [`private_haplotypes()`](https://thierrygosselin.github.io/radr/reference/private_haplotypes.md)
  : private haplotypes

## Specialized analyses

Investigate specialized marker classes and run external analytical
workflows.

- [`detect_microsatellites()`](https://thierrygosselin.github.io/radr/reference/detect_microsatellites.md)
  : Detect microsatellites
- [`sexy_markers()`](https://thierrygosselin.github.io/radr/reference/sexy_markers.md)
  : Identify sex-linked markers and reassign genetic sex
- [`run_bayescan()`](https://thierrygosselin.github.io/radr/reference/run_bayescan.md)
  : Run a BayeScan genome scan
- [`check_bayescan()`](https://thierrygosselin.github.io/radr/reference/check_bayescan.md)
  : Locate and validate BayeScan
- [`install_bayescan()`](https://thierrygosselin.github.io/radr/reference/install_bayescan.md)
  : Install BayeScan with Conda or Mamba

## Superseded interfaces

Compatibility names retained for older scripts; use the recommended
replacement documented on each page.

- [`radr_pkg_install()`](https://thierrygosselin.github.io/radr/reference/radr_pkg_install.md)
  : Legacy radr dependency helper
