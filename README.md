<img align="left" src="https://github.com/steinrue/diplo_locus/blob/main/diplolocus_logo.png?raw=true" alt="Logo" width=35%/>

<br/><br/><br/><br/><br/>

# `diplo-locus`: A lightweight toolkit for inference and simulation of time-series genetic data under general diploid selection

This repository hosts a Python CLI tool (`DiploLocus`) and python package (`diplo_locus`) that can be used to compute log-likelihoods for time series genetic data based on the diploid Wright-Fisher diffusion, as well as simulate such data.

To cite our application, use the reference:

>   Cheng, X. and Steinruecken, M. (2023) diplo-locus: A lightweight toolkit for inference and simulation of time-series genetic data under general diploid selection. https://doi.org/10.1101%2F2023.10.12.562101

--------------------

# Table of Contents

<a id="toc"> </a>

 * [Getting Started](#setup)
 * [`DiploLocus` Command-line interface](#CLI)
   * [`likelihood` Mode](#likelihood_CLI)

      <details><summary>Computing likelihood with CLI</summary>

        * [Quick Guide](#quickLL)
        * [Input](#LL_input)
          * [Sample allele counts](#input_samps)
          * [Sampling time points](#input_times)
        * [Specify parameter space for selection coefficients](#LL_grid)
        * [Output](#LL_output)
          * [Log-likelihood surfaces](#LL_likelihood)
          * [Interpolated local max log-likelihoods](#interp)
    
      </details>
    
   * [`simulate` Mode](#sim_CLI)   
      <details><summary>Simulating temporal samples/trajectories with CLI</summary>

        * [Quick Guide](#quickSim)
        * [Input and Output](#sim_IO)
      </details>
 
 * [`diplo_locus` package API](src/README.md)
 * [Examples](examples/README.md)

------------------------------------------------------------------------
# Getting Started
<a id="setup"> </a>

The CLI tools are designed to work in a unix shell-based command line working environment (such as Linux, MacOS, or Windows Sub-Linux system). 

### Install with `pip` from `github.com`

Both the API and CLI can be installed from the github repository using:

```shell
pip install "git+https://github.com/steinrue/diplo_locus@v1.2.0"
```

------------------------------------

# `DiploLocus` Command-line tool
<a id="CLI"> </a>

Once set up is complete, running the main CLI tool without arguments or with `-h` argument will print out
```shell
DiploLocus
#usage: DiploLocus [-h]  ...
#
#DiploLocus CLI
#
#optional arguments:
#  -h, --help  show this help message and exit
#
#DiploLocus Mode:
#  Subcommand to choose whether to compute likelihoods or simulate replicates of time-series
#  samples.
#
#
#    likelihood
#              Computes loglikelihood for time-series samples of independent loci.
#    simulate  Simulate and plot time-series samples and frequency trajectories of independent
#              replicates under diploid selection.
```

## `likelihood` Mode
<a id="likelihood_CLI"> </a>
 
### Quick Guide 
  <a id="quickLL">  </a>
Indicating `likelihood` to the main CLI tool or calling `DiploLocus-likelihood` directly can be used to compute likelihood of temporal data. In other words, the commands below are equivalent.

 ```shell
DiploLocus likelihood
DiploLocus-likelihood
 ````
 Running one of the above commands will list the possible comman dline arguments:
 ```shell
DiploLocus likelihood
#usage: DiploLocus likelihood [-h] --u01 U01 [--u10 U10] --Ne NE [--gen_time GEN_TIME]
#                             (-i INFILE | --vcf VCFFILE | --read_LL_from ONGRID_LL_FILE)
#                             [--info INFOFILE] [--ID_col IDCOLNAME]
#                             [--time_col TIMECOLNAME | --time_ago_col TIMEAGOCOLNAME]
#                             [--inds IND_SUBSET] [--force_hap FORCED_HAPS] [--force_dip FORCED_DIPS]
#                             [--sample_times SAMPLETIMES | --sample_times_ago SAMPLETIMES_AGO]
#                             [--snps SNP_SUBSET] [--minMAF MAF] --init {uniform,initFreq,statBeta}
#                             [--initFreq INITFREQ] [--t0 T0] [--force_t0] [--init_u01 INIT_U01]
#                             [--init_u10 INIT_U10] [--piecewise] [--specify_piece SPECIFY_PIECE]
#                             (--fix_s2 S2 | --linear_s2_range LIN_S2_RANGE | --geom_s2_range GEOM_S2_RANGE)
#                             ([--fix_s1 S1 | --linear_s1_range LIN_S1_RANGE | --geom_s1_range
#                             GEOM_S1_RANGE |]
#                             [--fix_h H | --linear_h_range LIN_H_RANGE | --geom_h_range GEOM_H_RANGE)]
#                             -o OUTPREFIX [--long_LL_table] [--gzip_surface] [--get_on_grid_max]
#                             [--get_off_grid_max] [--get_MLR] [--get_chi2_pval]
#                             [--numStates NUMSTATES] [-v]
#DiploLocus likelihood: error: the following arguments are required: --u01, --Ne, --init, -o/--out_prefix
```

The following arguments are the "must have"s for a complete `DiploLocus likelihood` run using command-line interface:

|                                                        | <b>Required Args</b>                                                                                                          | <b>Optional Args</b>                                                                                                                                                                                                                                                         |
|--------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| <b>Population Parameters</b>                           | `--Ne [pop size]` <br> `--u01 [mut rate 0->1 (per gen/nt)]`                                                                    | `--u10 [mut rate 1->0 (gen/nt)]` <br/> `--gen_time [time]`                                                                                                                                                                                                               |
| <b>Samples </b>  <br/>(VCF)                            | `--vcf <some.vcf(.gz)>` +<br/> `--info <info.table>` <br/>                                                                    | `--[ID/time/time_ago]_col [col name in <info.table>]` <br/> `--force_[hap/dip] {"all"/"none"/"ID1,ID2,..."}` <br/> `--snps {<snp.list>/"SNP1,SNP2,..."}` <br/> `--inds {<sampleID.list>/"ID1,ID2,..."}` <br/> `--snps {<snp.list>/"SNP1,SNP2,..."}` <br/>  `--minMAF [maf]`  |
| <b>Samples </b>  <br/> (allele counts)                 | ` -i <allele.counts>` +<br/> `--sample_times[_ago] <t1,t2,...>`                                                               | `--snps {<snp.list>/"SNP1,SNP2,..."}`<br/> `--minMAF [freq]`                                                                                                                                                                                                                 |
| <b>Initial Condition</b>                               | `--init {"uniform"/"initFreq"/"statBeta"}`                                                                                    | `--initFreq [freq]`<br/> `--init_[u01/u10] [mut rate 0->1/1->0 (gen/nt)]`                                                                                                                                                                                             |
| <b>Selection Parameters</b>                            | Specify `s2` + either `s1` or `h`; <br>`--[linear/geom]_[x]_range="start,end,steps"` <br/> or `--fix_[x] [value/"V1,V2,..."]` | `--piecewise` + `--specify_piece <piece.file>` <br/> `--t0 [init time]` `--force_t0`                                                                                                                                                                                    |
| <b>Output Options</b>                                  | `-o <prefix>` or `--out_prefix <prefix>`                                                                                      |                                                                                                                                                                                                                                                                              |
| <b>Pre-computed On-grid Likelihoods</b><br/>(optional) | `--read_LL_from <LL.table>`                                                                                                   | `--gzip_surface`<br/>`--long_LL_table`<br/> `--get_[on/off]_grid_max`<br/> `--get_MLR` `--get_chi2_pval`                                                                                                                                                                     |


With `-h` or `--help`, the user can retrieve the detailed usage for more information.

<details><summary>Here is the full usage:</summary>


```shell
DiploLocus likelihood -h
#usage: DiploLocus likelihood [-h] --u01 U01 [--u10 U10] --Ne NE [--gen_time GEN_TIME]
#                             (-i INFILE | --vcf VCFFILE | --read_LL_from ONGRID_LL_FILE)
#                             [--info INFOFILE] [--ID_col IDCOLNAME]
#                             [--time_col TIMECOLNAME | --time_ago_col TIMEAGOCOLNAME]
#                             [--inds IND_SUBSET] [--force_hap FORCED_HAPS] [--force_dip FORCED_DIPS]
#                             [--sample_times SAMPLETIMES | --sample_times_ago SAMPLETIMES_AGO]
#                             [--snps SNP_SUBSET] [--minMAF MAF] --init {uniform,initFreq,statBeta}
#                             [--initFreq INITFREQ] [--t0 T0] [--force_t0] [--init_u01 INIT_U01]
#                             [--init_u10 INIT_U10] [--piecewise] [--specify_piece SPECIFY_PIECE]
#                             (--fix_s2 S2 | --linear_s2_range LIN_S2_RANGE | --geom_s2_range GEOM_S2_RANGE)
#                             ([--fix_s1 S1 | --linear_s1_range LIN_S1_RANGE | --geom_s1_range
#                             GEOM_S1_RANGE |]
#                             [--fix_h H | --linear_h_range LIN_H_RANGE | --geom_h_range GEOM_H_RANGE)]
#                             -o OUTPREFIX [--long_LL_table] [--gzip_surface] [--get_on_grid_max]
#                             [--get_off_grid_max] [--get_MLR] [--get_chi2_pval]
#                             [--numStates NUMSTATES] [-v]
#
#optional arguments:
#  -h, --help            show this help message and exit
#  -v, --verbose         Option to write out selection coefficients when computing likelihoods.
#
#Population Parameters:
#  --u01 U01             Mutation rate allele 0 to allele 1 (/site/generation).
#  --u10 U10             Mutation rate allele 1 to allele 0 (/site/generation). Default to be equal to u01.
#  --Ne NE               Effective diploid population size (i.e., total #(allele)=2Ne ).
#  --gen_time GEN_TIME   Average generation time in the same unit as provided to '--sample_times' or
#                        '--sample_times_ago' argument. Default is 1.
#
#Genetic Data Input (choose one):
#  -i INFILE, --infile INFILE
#                        Path and name of the parsed input file.
#  --vcf VCFFILE         Path and name of the vcf file.
#  --read_LL_from ONGRID_LL_FILE
#                        Path and name to pre-computed log-likelihoods. File should be formatted the
#                        same as the on-grid log-likelihood output files generated by DiploLocus.
#
#For VCF input:
#  --info INFOFILE       Path and name of the annotation file. Must be tab-delimited and contain at
#                        least columns ("ID", "Gen_Ago") or ("ID", "Time_Ago") or ("ID", "Gen")
#  --ID_col IDCOLNAME    Name of the column in info table for ID names of the sample (as shown in
#                        VCF). Default is "ID".
#  --time_col TIMECOLNAME
#                        Name of the column in info table for each sample's sampling times (forward;
#                        ascending from past to present). Default name is "Time", unit is the number
#                        of generation unless specified with `--gen_time`.
#  --time_ago_col TIMEAGOCOLNAME
#                        Name of the column in info table for each sample's sampling times (backward;
#                        ascending from present to past). Default name is "Time_Ago", unit is the
#                        number of generation unless specified with `--gen_time`.
#  --inds IND_SUBSET     Comma-separated string or a plain-text file that lists a subset of sample
#                        IDs (among those included in the vcf) to consider.
#  --force_hap FORCED_HAPS
#                        Comma-separated string or a plain-text file that lists IDs (matching VCF
#                        column names) to be considered as haploids even though some may present
#                        diploid genotype calls. If "all", all samples will be considered as
#                        haploids, whose genotypes will be counted as matching haplotype (i.e. half
#                        the alleles). Force quit if any specified individual has heterozygote
#                        variant or any unmentioned diploid-formatted samples lack heterozygote
#                        variant calls. Default is "none".
#  --force_dip FORCED_DIPS
#                        Comma-separated string or a plain-text file that lists samples IDs
#                        (separated by comma) to be considered as diploids even though some may have
#                        loci with only haplotype calls. If "all", all samples will be considered as
#                        diploids, whose haploid GTs will be counted as matching homozygote diploid
#                        GTs (i.e. double-count the allele). Default is "none".
#
#For parsed allele counts (choose one):
#  --sample_times SAMPLETIMES
#                        Specify the time (ascending, from ancient to recent) of the same unit (e.g.,
#                        generations, years, weeks, etc.) when each pool of samples were taken. Must
#                        be positive numbers. Format: "<t1>,<t2>,...,<tk>". Number of time points
#                        must match the number of sample pools.
#  --sample_times_ago SAMPLETIMES_AGO
#                        Specify the time before present (descending, from ancient to recent) of the
#                        same unit (e.g., generations, years, weeks, etc.) when each pool of samples
#                        were taken. Format: "<t1>,<t2>,...,<tk>". Number of time points must match
#                        the number of sample pools.
#
#Input Options:
#  --snps SNP_SUBSET     Comma-separated string or a plain-text file that lists a subset of variant
#                        IDs (among those included in the input file) to consider.
#  --minMAF MAF          Minimum threshold (non-inclusive) for minor allele frequencies in the pooled
#                        samples. SNPs with sub-threshold frequencies will not be considered for
#                        analyses. Default is 0.
#
#Initial Condition:
#  --init {uniform,initFreq,statBeta}
#                        Specify initial condition for computing likelihoods.
#  --initFreq INITFREQ   Required when --init="initFreq". Specify the allele frequency when selection
#                        started.
#  --t0 T0               The time point (in generations) when initial distribution is assumed and selection 
#                        starts. Default is the first sampling time. IMPORTANT: This is interpreted the
#                        same way as the specified sampling times. Thus, if --sample_times or --time_col
#                        is used, this time is regular (forward), whereas with --sample_times_ago or
#                        --time_ago_col, this is time before present (backward).
#  --force_t0            Option to force implement a t0 later than the earliest samples. Data prior
#                        to t0 will be ignored.
#  --init_u01 INIT_U01   Specify mutation rate allele 0 to allele 1 for the initial distribution "statBeta".
#                        Default to be identical to u01.
#  --init_u10 INIT_U10   Specify mutation rate allele 1 to allele 0 for the initial distribution.
#                        Default to be identical to u10.
#
#Optional Mode(s):
#  --piecewise           Assume different selection coefficients in different periods. If selected, user
#                        must provide additional text file, using "--specify_piece", to specify in which
#                        period to vary [s] (and compute likelihooods on the given grid).
#  --specify_piece SPECIFY_PIECE
#                        Name of the file specifying selection periods when "--piecewise" is selected.
#                        It should be tab-delimited with 4 columns, using "x" instead of s1 and s2 to
#                        indicate the coefficients to vary in the respective time period:
#                        <from (forward gen#> <to (forward gen#> > <s1> <s2>
#Parameter Grids (specify s2 and [s1 | h]):
#  --fix_s2 S2           Provide one or more fixed values for s_AA (fitness gain in homozygotes),
#                        separated by comma (without space).
#  --linear_s2_range LIN_S2_RANGE
#                        Specify a linear grid as the parameter space for s_AA value. Format
#                        "(min,max,grid_num)".
#  --geom_s2_range GEOM_S2_RANGE
#                        Specify a geometrically-scaled grid as the parameter space for s_AA value.
#                        Format "(min,max,grid_num)". A closest number to <grid_num> will be
#                        determined such that 10^(-n)-fold values of <min> or <max> are included and
#                        that the smallest absolute value is greater than 1/2Ne. When <min> and <max>
#                        are of opposite signs, s=0 will be included, and two grids with an equal
#                        number of grid points will be formed on either side of zero.
#  --fix_s1 S1           Provide one or more fixed values for s_Aa (fitness gain in heterozygotes),
#                        separated by comma (without space).
#  --linear_s1_range LIN_S1_RANGE
#                        Specify a linear grid as the parameter space for s_Aa value. Format
#                        "(min,max,grid_num)".
#  --geom_s1_range GEOM_S1_RANGE
#                        Specify a geometrically-scaled grid as the parameter space for s_Aa value.
#                        Format "(min,max,grid_num)". A closest number to <grid_num> will be
#                        determined such that 10^(-n)-fold values of <min> or <max> are included and
#                        that the smallest absolute value is greater than 1/2Ne. When <min> and <max>
#                        are of opposite signs, s=0 will be included, and two grids with an equal
#                        number of grid points will be formed on either side of zero.
#  --fix_h H             Provide one or more fixed values for dominant coefficient h, separated by
#                        comma (without space).
#  --linear_h_range LIN_H_RANGE
#                        Specify a linear grid as the parameter space for dominant coefficient h.
#                        Format "(min,max,grid_num)".
#  --geom_h_range GEOM_H_RANGE
#                        Specify a geometrically-scaled grid as the parameter space for s_AA value.
#                        Format "(min,max,grid_num)". A closest number to <grid_num> will be
#                        determined such that 10^(-n)-fold values of <min> or <max> are included and
#                        that the smallest absolute value is greater than 1/2Ne. When <min> and <max>
#                        are of opposite signs, s=0 will be included, and two grids with an equal
#                        number of grid points will be formed on either side of zero.
#
#Output:
#  -o OUTPREFIX, --out_prefix OUTPREFIX
#                        Path and prefix of the output file.
#  --long_LL_table       Option to write out the computed loglikelihood surfaces to a plotting-
#                        friendly table with columns "locus, s1, s2, log-likelihood", instead of
#                        matrices.
#  --gzip_surface        Option to output the .gz file of computed likelihoods.
#  --get_on_grid_max     Option to output on-grid max log-likelihood and the matching (s1, s2) values
#                        among all values computed on each given locus, as `{PREFIX}_on-
#                        grid_maxLLs.txt`.
#  --get_off_grid_max    Option to interpolate the computed grid (must be at least 5x5) and infer
#                        local maximum likelihoods, written to file `{PREFIX}_off-grid_maxLLs.txt`.
#  --get_MLR             Option to write out the max-loglikelihood ratios relative to neutrality,
#                        i.e. 2(LL_opt - LL_0), for each locus, along side with the maximized log-
#                        likelihoods.
#  --get_chi2_pval       Option to write out the p-value of MLR in a standard chi-square distribution
#                        (df=1).
#
#HMM Parameters:
#  --numStates NUMSTATES
#                        Number of states across [0,1) for allele frequency. Default is 2000.
```

</details>


### Input
<a id="LL_input"> </a>

In order to compute likelihood for observing patterns of genetic variation in the given temporal samples, `DiploLocus likelihood` must be provided with: 
1) Allele counts (total number of observed alleles and its frequency in the sample), and 
2) Times when samples were taken.

To provide such information, the user has two options: 

- **Option 1**: **_VCF file_** (following `--vcf`) with accompanying _**sample information table**_ (`--info`, see [here](#input_samps) for detail),
  - The samples included in the VCF can be a mix of (pseudo-)haploid and diploid genotypes. Genotypes for variants in (pseudo-)haploid genomes must be coded as a single allele at each locus. 
  - The user can choose subsets of samples (`--inds`) or SNPs (`--snps`), either through comma-separated strings on the commandline, or through a plain text file.

- **Option 2**: **_Allele counts file_** (`--infile` or `-i`) with pooled batches of allele counts, and the **_sampling time_** of each pool _must_ be specified (with `--sample_times[_ago] <t1,t2,...>`). 
  - The values provided as sample times must be listed in the identical order as the sample batches listed in the input allele count file, where time is in units of 1 generation. User must specify the generation time with `--gen_time <# unit>` or `--gen_per_unit_time <# gen>` if it's not 1.

In cases when a likelihood surface has already been computed, it is possible to let the program reuse the pre-computed log-likelihoods to infer off-grid max likelihoods and parameter estimates on an 1D pre-computed parameter grid. The user can use `--read_LL_from <LLmatrix_file>` to achieve this.

The following subsections will describe in detail the requirements for either of the two options.

#### Sample allele counts from parsed text file or VCF file
 <a id="input_samps"> </a>

To present information on genetic variation of the samples, the user can provide either a VCF file or a tab-delimited plain text file. All the genetic variants are expected to be _bi-allelic_. A positive selection coefficient indicates that the alternative allele is beneficial (or the allele whose counts are given in the allele count file).

For VCF files (following `--vcf`), in addition to the [matching sample information table](#input_times), user should also pay attention to the format in which genotypes are presented. [Example 2](examples/README.md#example-2-selection-on-snps-in-lctmcm6-locus-within-uk-population-in-allen-acient-dna-resources-databse) shows an application of `DiploLocus likelihood` with VCF input.

When diploid genotype calls are specified in the data, but no heterozygote genotypes calls exist, the program will terminate, because this likely indicates a pseudo-haploid genotype that is incorrectly specified as diploid. In such cases, user must specify the sample IDs (as used in VCF column names) to indicate which samples represent haploids or diploids, using `--force_haps <id1,id2,...>` or `--force_dips <id1,id2,...>` tag, which will in turn override this termination.

- For the samples listed (comma-separated, without spaces) following `--force_haps`, the program will interpret haplotypes as-is while counting diploid genotypes as single alleles. That is, "0/0" will be intrpreted as "0", and "1/1" will be interpreted as "1". Note that the program will terminate with an error if any specified sample has a heterozygote genotype call.

- For the samples listed following `--force_dips`, the program will interpret diploid genotypes as-is and interpret haploid genotypes as double the allele count. That is, "0" will be interpreted as "0/0", and "1" as "1/1".

Alternatively, the user can supply a file where the temporal samples are grouped in discrete batches, with the observed count (dx) and the sample size (total number of alleles observed, nx) specified. This type of input files (`--infile` or `-i`) should be a _tab-delimited_ plain-text file with header. The first column has to be a unique identifier for the respective locus/dataset and has to have the header 'ID'. This unique identifier is used in the output. In addition, for _K_ batches of samples at different times, there have to be 2 x _K_ columns. Alternatingly, 'dx' for observed counts (could be either integers or fractions) and 'nx' for the sample size. Here x is an increasing index that starts at 1 to order the batches temporarily. Thus, the leftmost numbers belong to first sampling time and the rightmost to the last sampling time.

An example is the [input file for Example 1](examples/ex1_HorseLoci_count.txt):
```txt
# Data for loci related to horse coat coloration from Steinruecken & Song (2014)
ID	d1	n1	d2	n2	d3	n3	d4	n4	d5	n5	d6	n6
ASIP	0	10	1	22	15	20	12	20	15	36	18	38
MC1R	0	10	0	22	1	20	6	20	13	36	24	38
```

#### Sample information table for VCF input and sampling time points for allele count input
<a id="input_times"> </a>

If the user is providing samples using VCF (either `.vcf` or `.vcf.gz`) files, they must prepare an additional file with meta information for each of the samples to specify the respective sampling times. More specifically, this information table is expected to be a _**tab-delimited**_ plain text file with at least two columns, "ID" and "Time_Ago" (or "Gen_Ago", "Gen", or specified using `--time_col` or `--time_ago_col`), for each sample in the VCF file, _in the corresponding order_ that the samples are provided in the VCF file. One example of such file is [`examples/Ex2_v54_UK_1240K_noFam_strictPASS_from4500.info`](examples/Ex2_v54_UK_1240K_noFam_strictPASS_from4500.info), which specifies the mean dated time (years before 1950, listed under the `Time_Ago` column) for each sample in the [matching VCF file](examples/Ex2_UK_v54.1_all_chr2_135e6-137e6.vcf). 

For parsed allele count inputs, the user must provide additional information on timing with the argument `--sampling_times` or `--sampling_times_ago`, in the same order that sample batches are listed in the input file. That is, the times specified should be ascending when using the former (`--sampling_times`) or descending for the latter (`--sampling_times_ago`). 

##### Time units and generation times

When not specified, the program takes the numbers provided either in `--info` file or through `--sample_times[_ago]` as the numbers of generations. That is, the default time unit one generation.

In any case where the sampling times are not provided as the number of generations, the user should specify the generation time using either `--gen_time <# unit>` or `--gen_per_unit <# gen>`:

- `--gen_time <# unit>`: For organisms who live multiple time units (as used in the numbers provided), use this tag to indicate the average number of time units per generation. E.g., if the sampling time points are listed as "years ago" and the organism have a mean generation time of 5 years, then use `--gen_time 5`. If time is provided as the number of week and the organism has on average 12 weeks between generations, then use `--gen_time 12`
-  `--gen_per_unit <# gen>`: For organisms with short lifespans, sometimes it's more convenient to provide the number of generations within a given unit of time. For example, for a plant that flowers two times each year, one can use `--gen_per_unit 2` if the sampling times are presented in years.


### Initial condition
<a id="LL_init"> </a>
Initial condition refers to the initial probability distribution of the allele frequency assumed at time $t_0$ (specified by `--t0`), and selection starts at this time as well. Note that this time is interpreted the same way as the specified sampling times. Thus, if `--sample_times` or `--time_col` is used, this time is regular (forward), whereas with `--sample_times_ago` or `--time_ago_col`, this is time before present (backward).

There are three types of initial conditions available for users to specify (through `--init` flag):
* `--init uniform`: As all allele frequencies have equal probablity.
* `--init initFreq --initFreq <freq>`: This specifies the exact allele frequency at time $t_0$. This can be useful e.g. for experimental evolution data.
* `--init statBeta --init_u01 <mut rate>`: This option specifies a stationary Beta distribution with parameters $\alpha=2N_eu_{01}^{(init)}$ and $\beta=2N_eu_{10}^{(init)}$. If the mutation rate $u_{10}$ is not specified, it will be set to $u_{10}$.

Unless the exact selection scenario is known (see [Example 1](examples/README.md#example-1-log-likelihood-surfaces-for-ancient-horse-mc1r-and-asip-loci)) or the available sampling time points are too few, we recommend using `uniform` to minimize prior assumptions.
 
### Grid of values for selection coefficients
<a id="LL_grid"> </a>

`DiploLocus` can use the following two general diploid fitness schemes:

| Genotype (focal: A) | AA            | Aa            | aa  |
|---------------------|---------------|---------------|-----|
| Fitness scheme 1    | $s_2$         | $h\cdot s_2$  | 1   |
| Fitness scheme 2    | s<sub>2</sub> | s<sub>1</sub> | 1   | 

Here, $s_2$ and $s_1$ are are short for $s_{AA}$ and $s_{Aa}$, respectively. The user can use either fitness scheme to construct the parameter grid. In both schemes, `s2` has to be specified. The user can fix its value with `--fix_s2 [num]` or examine a range of values along a geometric (symmetric around 0; `--geom_s2_range="left,right,#step"`) or linear (`--linear_s1_range="left,right,#step"`) grid. Depending on the fitness scheme-of-choice, the user can also define a parameter space for `s1` or `h` in one of these two ways (`--geom_[x]_range` or `--linear_[x]_range`) or fix the parameter.

### Output

<a id="LL_output"></a>

In general, two types of output can be generated using the `likelihood` mode: 
* log-likelihood matrix computed under the given set of parameters, and 
* on-/off-grid maximum likelihood estimates on the specified 1-dimensional parameter grid (`--get_on_grid_max` or `--get_off_grid_max`) and the corresponding likelihood ratio statistics (`--get_MLR`). 

#### Compute Log-likelihood surfaces
<a id="LL_likelihood"></a>

By default, the program will output all log-likelihood values computed for each locus (say `L` loci) under each set of parameter values (say `P` pairs of selection coefficients) as an `L x P` matrix, saved to file `<outprefix>_LLsurfaces.table` or `<outprefix>_LLsurfaces.table.gz` (when using `--gzip_surface`).

Alternatively, with flag `--long_LL_table`, users can choose to have them written to a long table instead, for convenience of downstream processing. The file will be tab-delimited with header and four columns: `ID`, `s1`, `s2`, and `loglikelihood`. The `--gzip_surface` option can be used here as well.
 
#### On- and off-grid maximum log-likelihoods on 1D parameter space
<a id="interp"> </a>

With more than one set of selection parameters are considered, users can choose to output on-grid MLRs with the flag `--get_on_grid_max`, written to file `<outprefix>_on-grid_maxLLs.txt`, which is tab-delimited and has columns `ID`, `ongrid_s1hat`, `ongrid_s2hat`, and `ongrid_maxLogLikelihood`. Further, with `--get_MLR`, the output will include another `MLR` column with the likelihood ratio statistic, that is, two times the log ratio between the likelihood at the MLE and the likelihood for the same locus under neutrality (that is, $s_1=s_2=0$).

When only one of the selection parameter (out of $s_1$, $s_2$, and $h$) changes in value and more than five grid points are considered, users can choose to obtain interpolated off-grid MLR from the computed values (`--get_off_grid_max`). 

Either on- or off-grid, when the given parameters satisfy additive selection (dominance $h=0.5$, or homozygotes are twice as fit as heterozygote), the program can output the $p$ value assuming that the likelihood ratio statistics follow a standard $\chi^2$ distribution (degree of freedom =1) under the null hypothesis of no selection.

----------------------------------------------
 
## `simulate` Mode
<a id="sim_CLI"> </a>

### Quick Guide

<a id="quickSim"></a>

As with `likelihood`, indicating `simulate` to the main CLI tool is equivalent to calling `DiploLocus-simulate`:

```shell
DiploLocus simulate
DiploLocus-simulate
````
 Running one of the above commands will show

```shell
DiploLocus simulate
# usage: DiploLocus simulate [-h] --u01 U01 [--u10 U10] --Ne NE --s1 S1 --s2 S2
#                           [--selection_change_times SELECTIONCHANGETIMES] --sample_times SAMPLETIMES
#                           --sample_sizes SAMPLESIZES --num_rep NUMREP [--last_popfreq_nonzero]
#                           [--last_sample_nonzero] [--not_lost] [--not_fixed] [--minMAF MINMAF]
#                           --init {initFreq,statBeta,uniform,custom} [--initFreq INITFREQ]
#                           [--initDistn INITDISTNFILE] [--init_u01 INIT_U01] [--init_u10 INIT_U10] -o
#                           OUTPREFIX [--write_traj] [--gzip_output {s,t,st,ts}] [--plot_trajs]
#                           [--plot_samples] [--reps_to_plot REPS_TO_PLOT] [--deltaT DELTAT]
#                           [--seed SD]
#DiploLocus simulate: error: the following arguments are required: --u01, --Ne, --s1, --s2, --sample_times, --sample_sizes, --num_rep, --init, -o/--out_prefix
```
Likewise, the user can check out the full usage with `-h`

<details><summary>Click here to see the full usage</summary>

```shell
DiploLocus simulate -h
# usage: DiploLocus simulate [-h] --u01 U01 [--u10 U10] --Ne NE --s1 S1 --s2 S2
#                              [--selection_change_times SELECTIONCHANGETIMES] --sample_times
#                              SAMPLETIMES --sample_sizes SAMPLESIZES --num_rep NUMREP
#                              [--last_popfreq_nonzero] [--last_sample_nonzero] [--not_lost_ever]
#                              --init {initFreq,statBeta,uniform,custom} [--initFreq INITFREQ]
#                              [--initDistn INITDISTNFILE] [--init_u01 INIT_U01]
#                              [--init_u10 INIT_U10] [--deltaT DELTAT] [--seed SD] -o OUTPREFIX
#                              [--write_traj] [--gzip_output {s,t,st,ts}] [--plot_trajs]
#                              [--plot_samples] [--reps_to_plot REPS_TO_PLOT]
#
# optional arguments:
#  -h, --help            show this help message and exit
#
# Population Parameters:
#  --u01 U01             Mutation rate allele 0 to allele 1 (/site/generation).
#  --u10 U10             Mutation rate allele 1 to allele 0 (/site/generation). Default is u01.
#  --Ne NE               Effective diploid population size (i.e., total #(allele)=2Ne ).
#
# Selection Parameters:
#  --s1 S1               Value(s) of selection coefficient for heterozygotes. If more than one value
#                        is provided (separated by comma, no space), then user must also provide
#                        selection change time with "--selection_change_times".
#  --s2 S2               Value(s) of selection coefficient for homozygotes. If more than one value is
#                        provided (separated by comma, no space), then user must also provide
#                        selection change time with "--selection_change_times".
#  --selection_change_times SELECTIONCHANGETIMES
#                        Times (in forward generations) when either selection coefficient changes.
#                        Must be 1 fewer than the number of s values provided.
#
# Sample Parameters:
#  --sample_times SAMPLETIMES
#                        Numbers of generation times to take samples, separated by comma (no sapce).
#                        Must match the length of sample sizes.
#  --sample_sizes SAMPLESIZES
#                        Numbers of alleles to sample at each time points specified by "--
#                        sample_times". Must match the length of sample times.
#  --num_rep NUMREP      Number of replicates to generate.
#  --last_popfreq_nonzero
#                        Option to only count replicates with non-zero population frequency in the
#                        end.
#  --last_sample_nonzero
#                        Option to only count replicates where the last (i.e. most recent) sample is
#                        not zero.
#  --not_lost_ever       Option to only count replicates where the allele is not lost throughout the
#                        simulation.
#
# Initial Condition:
#  --init {initFreq,statBeta,uniform,custom}
#                        Specify initial probability distribution to generate the replicates
#  --initFreq INITFREQ   Required when --init="initFreq". Specify the allele frequency when selection
#                        started. If multiple values are provided (separated by comma), then length
#                        must match the number of replicates.
#  --initDistn INITDISTNFILE
#                        Path and name to the file specifying a particular frequency distribution as
#                        the initial condition.
#  --init_u01 INIT_U01   Specify mutation rate allele 0 to allele 1 for the initial distribution.
#                        Default is u01.
#  --init_u10 INIT_U10   Specify mutation rate allele 1 to allele 0 for the initial distribution.
#                        Default is u10.
#
# Output Options:
#  -o OUTPREFIX, --out_prefix OUTPREFIX
#                        Path and prefix of the output file.
#  --write_traj          Option to also write out the trajectories of replicates.
#  --gzip_output {none,s,t,st,ts}
#                        Option to write output in .gz format. Use "s" to indicate sample file, "t"
#                        to indicate trajectory file, and "st" (or "ts") for both. Use "none" to not
#                        zip either.
#  --plot_trajs          Option to plot the population allele frequency trajectories of up to 10
#                        replicates.
#  --plot_samples        Option to plot the sample allele frequency trajectories of up to 10
#                        replicates.
#  --reps_to_plot REPS_TO_PLOT
#                        Specify the replicates to plot. Default is the first 10 reps.
#
# Other Diffusion Parameters:
#  --deltaT DELTAT       Unit increment of time (in generations).
#  --seed SD             Specify random seed.
#
```

</details>

### Input and Output

<a id="sim_IO"> </a>

The input required for simulations are the population parameters, _i.e._ mutation rate (`u01`, `u10`) and effective population size (`Ne`), sample times, sample sizes, and the number of replicates. In addition, similar to `likelihood`, the user must also specify the initial condition with `--init`. It is highly recommended that the user also specify the random seed (through `--seed`). [Example 3](examples/README.md#example-3-simulating-temporal-allele-count-data-and-plot-their-trajectories) demonstrates in detail how to parameterize a set of simulations.

The table below lists the main input arguments:


|                              | <b>Required Args</b>                                                                   | <b>Optional Args</b>                                                                                                        |
|------------------------------|----------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------|
| <b>Population Parameters</b> | `--Ne [pop size]` <br> `--u01 [mut rate 0->1 (gen/nt)]`                             | `--u10 [mut rate 1->0 (gen/nt)]`                                                                                        |
| <b>Sample Information</b>    | ` --sample_sizes <n1,n2,...>` <br/> `--sample_times <t1,t2,...>` <br/> `--num_rep <N>` | `--seed [seed]`<br/> `--deltaT [dT]`                                                                                        |
| <b>Initial Condition</b>     | `--init {"uniform"/"initFreq"/"statBeta"}`                                             | `--initFreq [freq]`<br/> `--init_[u01/u10] [mut rate 0->1/1->0 (gen/nt)]`<br/> `--initDistn <distnFile>`               |
| <b>Selection Parameters</b>  | `--s1`, `--s2` <br/> followed by either a single value<br/>or comma-delimited values   | `--selection_change_times [t12,t23,...]`                                                                                    |
| <b>Output Options</b>        | `-o <prefix>` or `--out_prefix <prefix>`                                               | `--gzip_output {s,t,st,ts}` <br/>`--write_traj`, <br/> `--plot_trajs`, `--plot_samples` <br/>`--reps_to_plot [id1,id2,...]` |


The software writes out the simulated samples as an allele count table (with `_samples.count` suffix) in the same format expected by the `likelihood` mode. See e.g. the ouput of [Example 3](examples/README.md#example-3-simulating-temporal-allele-count-data-and-plot-their-trajectories).

When the flag `--write_trajs` is specified, an additional text file (with suffix `.traj`) will be generated to list out the allele frequency trajectories of all the replicates. More specifically, for `N` replicates simulated over a duration of `K` generations using `dT` time increments ("delta T"), the file would be a comma-delimited `N`x`K/dT` table of allele frequency at each iteration, with each row representing a replicate. 

When a large number of replicates or a large number of iterations are simulated, the user can choose to output simulated data in `.gz` format. To do so, use `--gzip_output {s|t}` with `s` and `t` for "samples" and "trajectories", respectively.

For the user's convenience, the software provides the option to plot the sample and/or population allele frequency trajectories at the end of the simulation. The flags `--plot_trajs` and `--plot_samples` allow the user to turn on these options. Without the `--reps_to_plot [id1,id2,...]` instruction, the software will only consider the first ten replicates.

To simulate discrete time-varying selection, the user can make use of `--selection_change_time` flag, pairing it with values provided to `--s1` and `--s2`. [Example 3](examples/README.md#example-3-simulating-temporal-allele-count-data-and-plot-their-trajectories) shows in detail how to simulate data under such selection models.







