#!/usr/bin/env python3
"""
Command-line interfacing Python script to compute likelihoods of allele frequency time series under general diploid selection in the Wright-Fisher diffusion.
"""
import sys, os, re, time
import gzip
import numpy as np
import pandas as pd

# local import
# from diplo_locus import diffusion_core
try:
    from diplo_locus import likelihood
    from diplo_locus.utility import _get_geom_grid, _reformat_LL_DF_to_matrix
    from diplo_locus.utility import _reformat_longLLtable_to_matrix
    from diplo_locus.utility import parse_allele_counts, parse_vcf_input, parse_ind_arg
    from diplo_locus.utility import get_onGrid_max_only, interpolate_offGrid_max
except ModuleNotFoundError or ImportError:
    # assume it's the repo directory
    from src import likelihood
    from src.utility import _get_geom_grid, _reformat_LL_DF_to_matrix
    from src.utility import _reformat_longLLtable_to_matrix
    from src.utility import parse_allele_counts, parse_vcf_input, parse_ind_arg
    from src.utility import get_onGrid_max_only, interpolate_offGrid_max


# new simplified way to parse sampling time
def parse_sampling_times(times: str, times_ago: str, gen_time=1, t0=None, force_t0=False):
    """Parse sampling times from command-line args. Return parsed sample times <: list>
     as in numbers of generations. Zero will be added at the beginning."""

    def _check_times(sampTimes: list, t_0=None):
        # sanity check: make sure it's ascending
        sampTimes_diff = np.diff(sampTimes)
        assert np.all(sampTimes_diff >= 0), \
            f'Please ensure \"--sample_times\" be given ascending numbers.\nnp.diff={sampTimes_diff}'
        # check t0
        if t_0 is not None:
            if t_0 > sampTimes[0]:
                if force_t0:
                    sample_cols_toKeep = [i for i, t in enumerate(sampTimes) if (t - t_0) >= 0]
                    print(f'Ignoring data sampled before t0={t_0} unit.')
                    sampTimes = [(t - t_0) / gen_time for t in sampTimes if (t - t_0) >= 0]
                else:
                    raise ValueError(f'Selection starting time \'t0\' is set to be more recent '
                                     f'than the oldest samples. Samples older than \'t0\' will not '
                                     f'be considered in the subsequent computations. To continue '
                                     f'despite discarding these data, please use \'--force_t0\' to force it.')
            else:
                sampTimes = [(t - t_0) / gen_time for t in sampTimes]
                sample_cols_toKeep = list(range(len(sampTimes)))
        else:
            t_0 = sampTimes[0]
            sampTimes = [(t - t_0) / gen_time for t in sampTimes]
            sample_cols_toKeep = list(range(len(sampTimes)))
        # sanity check in the end: all positive
        assert np.all(np.array(sampTimes) >= 0), f'Please make sure samples times are positive numbers.'
        return sampTimes, sample_cols_toKeep

    def _check_times_ago(Times_ago: list, t_0=None):
        # sanity check: make sure it's descending
        Times_ago_diff = np.diff(Times_ago)
        assert np.all(Times_ago_diff <= 0), \
            f'Please ensure \'--sample_times_ago\' be given descending numbers.\nnp.diff={Times_ago_diff}'
        # check t0
        if t_0 is not None:
            if t_0 < Times_ago[0]:
                if force_t0:
                    sample_cols_toKeep = [i for i, t in enumerate(Times_ago) if (t_0 - t) >= 0]
                    print(f'Ignoring data sampled before t0={t_0} unit ago.')
                    sampTimes = [(t_0 - t) / gen_time for t in Times_ago if (t_0 - t) >= 0]
                else:
                    raise ValueError(f'Selection starting time \'t0\' is set to be more recent '
                                     f'than the oldest samples. Samples older than \'t0\' will not '
                                     f'be considered in the subsequent computations. To continue '
                                     f'despite discarding these data, please use \'--force_t0\' to force it.')
            else:
                sampTimes = [(t - t_0) / gen_time for t in Times_ago]
                sample_cols_toKeep = list(range(len(sampTimes)))
        else:
            t_0 = Times_ago[0]
            sampTimes = [(t - t_0) / gen_time for t in Times_ago]
            sample_cols_toKeep = list(range(len(sampTimes)))
        sampTimes = np.array(sampTimes)
        # sanity check: ascending
        sampTimes_diff = np.diff(sampTimes)
        assert np.all(sampTimes_diff >= 0), f'sampTimes={sampTimes}'
        # sanity check in the end: all positive
        assert np.all(np.array(sampTimes) >= 0), f'Please make sure samples times are positive numbers.'
        return sampTimes, sample_cols_toKeep

    if times is not None:
        # parse the command-line argument first
        sampleTimes = np.array([float(x) for x in times.strip().split(',')])
        sampleTimes, sample_cols_toKeep = _check_times(sampleTimes, t0)
    elif times_ago is not None:
        # parse the command-line argument first
        sampleTimes_ago = np.array([float(x) for x in times_ago.strip().split(',')])
        sampleTimes, sample_cols_toKeep = _check_times_ago(sampleTimes_ago, t0)
    else:
        raise TypeError('Sampling times must be provided when input file is parsed allele counts.')
    # lastly: add zero
    sampleTimes = np.concatenate(([0], sampleTimes))
    # print(sampleTimes)
    return sampleTimes, sample_cols_toKeep


def parse_ind_arg(inds: str, type: str = "inds"):
    # the same function work for both inds and SNP IDs
    # check if it's a file
    if os.path.isfile(inds):
        print(f'Reading select {type} from {inds}')
        ## read in the file
        with open(inds, "r") as ind_file:
            ind_list = ind_file.read()
        ind_file.close()
        ## recognize separators: comma, tab, or newline
        ind_list = re.split(r'[,\r\n\t]{1}', ind_list)
    else:
        try:
            ind_list = inds.split(",")
        except Exception as e:
            print(e)
            return False
    print(f'{len(ind_list)} {type} in total.')
    return ind_list


ago_regex = re.compile(r'[A|a]go')


# file needs to be tab-delimited, must have header, with at least ID & Gen_Ago
# if Gen_Ago absent, must provide gen_time
def read_info_file(infofile: str, ID_colname: str = 'ID', time_colname=None, gen_time=1,
                   t0=None, force_t0: bool = False, inds=None):
    """Take in sample times from annotation file.
    Read `--info` file, return sample times (of length K+1) and a K-element list of
        sample ID lists that match the sampling times.
    -------------------------------------------
    :param infofile: str, filename for the tab-delimited info table. Must at least contain
        a column for sample IDs and a column for sampling times.
    :param ID_colname: Name of the column for sample IDs (as shown in the VCF).
    :param time_colname: Name of the column for sampling times. If the name contain \"[A|a]go\",
        then time values will be counted backward (increasing from present to past); otherwise,
        time will be presumed to go forward (increasing from past to present).
    :param gen_time: Number of time units per generation. Default is 1.
    :param t0: time when selection starts. Must be in the same unit as the info table
    :param force_t0: bool, decide whether to discard data (if any) collected before the `t0` provided. Default is False.
    :param inds: None or list, list of sample IDs (subset of all IDs in the ID column) to consider;
            if `None`, then all samples will be considered.
    :return:
    """
    infoTable = pd.read_csv(infofile, sep="\t")  # , comment="#", header=0
    assert ID_colname in infoTable.columns, ValueError(f"VCF Info file must have a column named \"{ID_colname}\".")
    if time_colname is not None:
        assert (time_colname in infoTable.columns) or (time_colname[:-4] in infoTable.columns), \
            ValueError(f"VCF Info file must have a column named \"{time_colname}\".")
    # only keep the individual in inds
    if inds is not None:
        assert isinstance(inds, list), f'type(inds) = {type(inds)}'
        assert all(isinstance(ind, str) for ind in inds), f'inds list: {inds}'
        infoTable = infoTable[infoTable[ID_colname].isin(inds)]

    # may be easier if we pre-define stuff
    def _count_time_forward(infotable, timeCol, t0):
        # forward time:
        print(f'Counting time from past to present, extracting info from \"{timeCol}\" column.')
        ## check for missing value
        if np.any(infotable[timeCol].isnull()):
            infotable = infotable[~infotable[timeCol].isnull()]
        if t0 is None:
            t0 = min(infoTable['Gen'])
            # t0 = 0
            print('Assuming selection starts at generation zero.')
        elif not force_t0:
            assert t0 <= min(infotable[
                                 "Gen"]), 'Please make sure selection start time `t0` is of the same unit as provided under the \"Gen\" column.\nIf selection is set to start after the oldest samples are collected, please use `--force_t0` to allow omitting data preceding the selection onset.'
        # infoTable['Gen'] = np.array(np.array((infoTable['Gen']) - t0), dtype=int)
        # convert to generation
        infotable['Gen'] = (infotable[timeCol] - t0) / gen_time
        # just in case
        if np.any(infotable["Gen"].isnull()):
            infotable = infotable[~infotable["Gen"].isnull()]
        infotable['Gen'] = infotable['Gen'].round().astype(int)
        # done processing
        return infotable

    def _count_time_backward(infotable, timeCol, t0):
        # backward time
        print(f'Counting time from present to past, extracting info from \"{timeCol}\" column.')
        ## check for missing value
        if np.any(infotable[timeCol].isnull()):
            infotable = infotable[~infotable[timeCol].isnull()]
        # t0 needs to be the same unit too
        if t0 is None:
            t0 = max(infotable[timeCol])
        elif not force_t0:
            assert t0 >= max(infotable[
                                 timeCol]), f'Please make sure selection start time `t0`={t0} is of the same unit as provided in the info table.\nIf selection is set to start after the oldest samples are collected, please use `--force_t0` to allow omitting data preceding the selection onset.'
        # convert to the number of generations
        infotable['Gen'] = (t0 - infotable[timeCol]) / gen_time
        # just in case
        if np.any(infotable["Gen"].isnull()):
            infotable = infotable[~infotable["Gen"].isnull()]
        # print(f't0={t0}, {infoTable["Time_Ago"].describe()}, gen_time={gen_time}') #\n{list(infoTable["Time_Ago"])}
        # algorithm only accept integer generations for now
        infotable['Gen'] = infotable['Gen'].round().astype(int)
        # done processing
        return infotable

    if time_colname is not None:
        find_ago = ago_regex.findall(time_colname)

        if len(find_ago) == 0:
            infoTable = _count_time_forward(infoTable, time_colname, t0)
        else:
            if len(find_ago) > 1:
                if time_colname.endswith('_Ago'):  # <-- remove the tag previously added
                    time_colname = time_colname[:-4]
            else:  # the OG name could be without "ago" too
                assert len(find_ago) == 1
                if time_colname not in infoTable.columns:
                    time_colname = time_colname[:-4]
                assert time_colname in infoTable.columns, ValueError(
                    f'Column \"{time_colname}\" missing in the info table.')
            infoTable = _count_time_backward(infoTable, time_colname, t0)
    # if any column has "ago" in it
    elif np.any(np.array([ago_regex.search(name) for name in infoTable.columns], dtype=bool)):
        assert 'Gen_Ago' in infoTable.columns or 'Time_Ago' in infoTable.columns, \
            f'Column for sampling time unspecified. Current column names: {infoTable.columns}, gen_time={gen_time}.'
        if 'Gen_Ago' in infoTable.columns and 'Time_Ago' not in infoTable.columns:
            timeColName = 'Gen_Ago'
        elif 'Time_Ago' in infoTable.columns and 'Gen_Ago' not in infoTable.columns:
            if gen_time == 1:
                print('Generation time not specified. Assume values under \"Time\" column are numbers of generations.')
            timeColName = 'Time_Ago'
        else:  # if both are in the table
            print('Both\"Time_Ago\" and \"Gen_Ago\" exist in table. Only reading numbers in \"Gen\" column.',
                  infoTable.columns)
            # shall we check whether time/gen_time == gen? <--- maybe not, in case gen_time changes over time
            timeColName = 'Gen_Ago'
        # count backward
        infoTable = _count_time_backward(infoTable, timeColName, t0)
    else:
        assert 'Gen' in infoTable.columns or 'Time' in infoTable.columns, \
            f'Column for sampling time unspecified. Current column names: {infoTable.columns}, gen_time={gen_time}.'
        if 'Gen' in infoTable.columns and 'Time' not in infoTable.columns:
            timeColName = 'Gen'
        elif 'Time' in infoTable.columns and 'Gen' not in infoTable.columns:
            if gen_time == 1:
                print('Generation time not specified. Assume values under \"Time\" column are numbers of generations.')
            timeColName = 'Time'
        else:  # if both are in the table
            print('Both\"Time\" and \"Gen\" exist in table. Only reading numbers in \"Gen\" column.')
            # shall we check whether time/gen_time == gen? <--- maybe not, in case gen_time changes over time
            timeColName = 'Gen'
        # count forward
        infoTable = _count_time_forward(infoTable, timeColName, t0)

    # now check whether force_f0 is needed.
    if np.any(np.array(infoTable["Gen"]) < 0):
        if force_t0:
            print(f'Ignoring data sampled before t0={t0}.')
            infoTable = infoTable[infoTable['Gen'] >= 0]
        else:
            print(infoTable[infoTable["Gen"] < 0])
            raise ValueError(f'Selection starting time \'t0\'={t0} is set to be more recent '
                             f'than the oldest samples. Samples older than \'t0\' will not '
                             f'be considered in the subsequent computations. To continue '
                             f'despite discarding these data, please use \'--force_t0\' .')
    # group by gen
    Pools = list(infoTable.groupby(by="Gen")[ID_colname])
    # pandas algorithm ensures that time_points are sorted
    time_points, sample_pools = map(list, zip(*Pools))
    sampleTimes = np.array([0] + time_points)
    samplePools = [list(pool) for pool in sample_pools]
    return sampleTimes, samplePools


def parse_grid_args(s2, lin_s2_range, geom_s2_range, s1, lin_s1_range, geom_s1_range, h, lin_h_range, geom_h_range,
                    Ne):
    """return an ordered list of (s1, s2) tuples"""
    if s2 is not None:
        s2_list = list(map(float, s2.split(',')))
        s2_list = sorted(list(set(s2_list)))
        s2_list = np.array(s2_list)
    elif lin_s2_range is not None:
        start, end, gridnum = map(float, lin_s2_range.split(','))
        s2_list = np.linspace(start, end, int(gridnum + 1))
    else:
        start, end, gridnum = map(float, geom_s2_range.split(','))
        s2_list = _get_geom_grid(start, end, gridnum, Ne)

    # same with s1 & h
    # reverting for-loops in list construction to match the order of np.meshgrid
    if s1 is not None:
        s1_list = list(map(float, s1.split(',')))
        s1_list = sorted(list(set(s1_list)))
        s_pairs = [(s1, s2) for s2 in s2_list for s1 in s1_list]
    elif lin_s1_range is not None:
        start, end, gridnum = map(float, lin_s1_range.split(','))
        s1_list = np.linspace(start, end, int(gridnum + 1))
        s_pairs = [(s1, s2) for s2 in s2_list for s1 in s1_list]
    elif geom_s1_range is not None:
        start, end, gridnum = map(float, geom_s1_range.split(','))
        s1_list = _get_geom_grid(start, end, gridnum, Ne)
        s_pairs = [(s1, s2) for s2 in s2_list for s1 in s1_list]
    # restraints on s2*h should be added here
    elif h is not None:
        h_list = list(map(float, h.split(',')))
        if len(h_list) > 1:
            h_list = np.array(sorted(list(set(h_list))))
        s_pairs = [(s2 * h, s2) for s2 in s2_list for h in h_list]
        s1_list = np.array(sorted(list(set(s1 for (s1, s2) in s_pairs))))
    elif lin_h_range is not None:
        start, end, gridnum = map(float, lin_h_range.split(','))
        h_list = np.linspace(start, end, int(gridnum + 1))
        s_pairs = [(s2 * h, s2) for s2 in s2_list for h in h_list]
        s1_list = np.array(sorted(list(set(s1 for (s1, s2) in s_pairs))))
    elif geom_h_range is not None:
        start, end, gridnum = map(float, geom_h_range.split(','))
        h_list = _get_geom_grid(start, end, gridnum, Ne)
        s_pairs = [(s2 * h, s2) for s2 in s2_list for h in h_list]
        s1_list = np.array(sorted(list(set(s1 for (s1, s2) in s_pairs))))
    else:
        raise ValueError('Must provide value grid for either s1 or h.')

    # s_pairs = itertools.product(s1_list.flatten(), s2_list)
    # make sure 1) no dups 2) sorted (if s1_list and s2_list doesn't have dup, pairs shouldn't have dup either
    num_pairs = len(s_pairs)
    print(f'Computing log-likelihoods for {num_pairs} pairs of (s1, s2) values.')
    return s1_list, s2_list, s_pairs, num_pairs


import argparse


def lik_args_parser(parser):
    # some pop gen parameters
    param = parser.add_argument_group('Population Parameters')
    param.add_argument('--u01', dest='u01', default=None, required=('--read_LL_from' not in sys.argv), type=float,
                       help="Forward mutation rate (/site/generation).")
    param.add_argument('--u10', dest='u10', default=None, type=float,
                       help='Backward mutation rate (/site/generation). Default to be equal to u01.')
    param.add_argument('--Ne', dest='Ne', default=None, required=('--read_LL_from' not in sys.argv), type=float,
                       help='Effective diploid population size (i.e., total #(allele)=2Ne ).')
    param.add_argument('--gen_time', dest='gen_time', default=1, type=float,
                       help='Average generation time in the same unit as provided to \'--sample_times\' or \'--sample_times_ago\' argument. Default is 1. ')

    inputs = parser.add_argument_group('Genetic Data Input (choose one)')
    inputFiles = inputs.add_mutually_exclusive_group(required=True)

    inputFiles.add_argument('-i', '--infile', dest='infile', default=None,
                            help='Path and name of the parsed input file.')
    inputFiles.add_argument('--vcf', dest='vcfFile', default=None, help='Path and name of the vcf file.')
    inputFiles.add_argument('--read_LL_from', dest='onGrid_LL_file', default=None,
                            # required=('--interpolate_max' in sys.argv),
                            help='Path and name to pre-computed log-likelihoods. File should be formatted the same as the on-grid log-likelihood output files generated by DiploLocus.')

    vcfs = parser.add_argument_group('For VCF input')
    vcfs.add_argument('--info', dest='infoFile', required=('--vcf' in sys.argv),
                      help='Path and name of the annotation file. Must be tab-delimited and '
                           'contain at least columns (\"ID\", \"Gen_Ago\") or (\"ID\", \"Time_Ago\") or (\"ID\", \"Gen\")')
    # vcfs.add_argument('--pseudo_hap', dest='pseudo', action='store_true', default=False,
    #                   help='Indicate whether the vcf files are for samples called as pseudo haplotypes.')
    # IDs_regex = re.compile(r',?(\w+),?')
    vcfs.add_argument('--ID_col', dest='IDcolName', default='ID',
                      help='Name of the column in info table for ID names of the sample (as shown in VCF). Default is \"ID\".')
    vcf_time = vcfs.add_mutually_exclusive_group()
    vcf_time.add_argument('--time_col', dest='TimeColName', default=None,
                          help='Name of the column in info table for each sample\'s sampling times (forward; ascending from past to present). Default name is \"Time\", unit is the number of generation unless specified with `--gen_time`.')
    vcf_time.add_argument('--time_ago_col', dest='TimeAgoColName', default=None,
                          help='Name of the column in info table for each sample\'s sampling times (backward; ascending from present to past). Default name is \"Time_Ago\", unit is the number of generation unless specified with `--gen_time`.')
    vcfs.add_argument('--inds', dest='ind_subset', default=None,
                      help='Comma-separated string or a plain-text file that lists a subset of sample IDs (among those included in the vcf) to consider.')
    vcfs.add_argument('--force_hap', dest='forced_haps', default='none',
                      # type=lambda x: x in {'all', 'none'} or len(IDs_regex.findall(x)) > 0,
                      help='Comma-separated string or a plain-text file that lists IDs (matching VCF column names) to be considered as haploids even though some may present diploid genotype calls. If \"all\", all samples will be considered as haploids, whose genotypes will be counted as matching haplotype (i.e. half the alleles). Force quit if any specified individual has heterozygote variant or any unmentioned diploid-formatted samples lack heterozygote variant calls. Default is \"none\".')
    vcfs.add_argument('--force_dip', dest='forced_dips', default='none',
                      # type=lambda x: x in {'all', 'none'} or len(IDs_regex.findall(x)) > 0,
                      help='Comma-separated string or a plain-text file that lists samples IDs (separated by comma) to be considered as diploids even though some may have loci with only haplotype calls. If \"all\", all samples will be considered as diploids, whose haploid GTs will be counted as matching homozygote diploid GTs (i.e. double-count the allele). Default is \"none\".')

    parsedInputs = parser.add_argument_group('For parsed allele counts (choose one)')
    sampTimeArgs = parsedInputs.add_mutually_exclusive_group(required=(('-i' in sys.argv) or ('--input' in sys.argv)))
    sampTimeArgs.add_argument('--sample_times', dest='sampleTimes', default=None,
                              help='Specify the time (ascending, from ancient to recent) of the same unit (e.g., generations, years, weeks, etc.) when each pool of samples were taken. Must be positive numbers. Format: \"<t1>,<t2>,...,<tk>\". Number of time points must match the number of sample pools.')
    sampTimeArgs.add_argument('--sample_times_ago', dest='sampleTimes_ago', default=None,
                              help='Specify the time before present (descending, from ancient to recent) of the same unit (e.g., generations, years, weeks, etc.) when each pool of samples were taken. Format: \"<t1>,<t2>,...,<tk>\". Number of time points must match the number of sample pools.')

    inputOpt = parser.add_argument_group('Input Options')
    inputOpt.add_argument('--snps', dest='snp_subset', default=None,
                          help='Comma-separated string or a plain-text file that lists a subset of variant IDs (among those included in the input file) to consider.')
    inputOpt.add_argument('--minMAF', dest='MAF', type=float, default=0.,
                          help='Minimum threshold (non-inclusive) for minor allele frequencies in the pooled samples. '
                               'SNPs with sub-threshold frequencies will not be considered for analyses.  Default is 0.')
    # initial conditions
    initConds = parser.add_argument_group('Initial Condition')
    initConds.add_argument('--init', dest='initCond', required=('--read_LL_from' not in sys.argv),
                           choices=['uniform', 'initFreq', 'statBeta'],  # , 'statBetaSel', 'custom'
                           help='Specify initial condition for computing likelihoods.')

    initConds.add_argument('--initFreq', dest='initFreq', type=float,
                           required=('initFreq' in sys.argv),
                           help='Required when --init=\"initFreq\". Specify the allele frequency when selection started.')

    initConds.add_argument('--t0', dest='t0', default=None, type=float,
                           help='The time point (in generations) when selection starts. Default is the earliest sampling time.')
    # 'Only applies when selection starts from de-novo mutations and initial condition is \"initFreq\". For selection on standing variation, consider `--piecewise` option.'
    initConds.add_argument('--force_t0', dest='force_t0', default=False, action='store_true',
                           help='Option to force implement a t0 later than the earliest samples. Data prior to t0 will be ignored.')
    initConds.add_argument('--init_u01', dest='init_u01', type=float,  # required=('statBetaSel' in sys.argv),
                           default=None,
                           help='Specify forward mutation rate for the stationary distribution \"statBeta\". Default to be identical to u01.')
    initConds.add_argument('--init_u10', dest='init_u10', type=float,  # required=('statBetaSel' in sys.argv),
                           default=None,
                           help='Specify backward mutation rate for the stationary distribution. Default to be identical to u10.')

    # piecewise
    # read the file that specify which one to optimize
    optional = parser.add_argument_group("Optional Mode(s)")
    optional.add_argument('--piecewise', dest='piecewise', action='store_true', default=False,
                          help='Option to assume changing selection coefficients and choose the time period when they vary (by the grid specified).\nIf selected, user must provide an additional text file, using \"--specify_piece\", to specify which period to vary [s].')
    optional.add_argument('--specify_piece', dest='specify_piece', default=None, required=('--piecewise' in sys.argv),
                          help='Path and name to the helper file when \"--piecewise\" is selected. It should be tab-delimited with 4 columns, using \"x\" to indicate the coefficients to vary in the matching time period:\n <from (forward gen#>\t<to (forward gen#>\t>\t<s1>\t<s2> ')

    # customize grids
    grid = parser.add_argument_group('Parameter Grids (specify s2 and [s1 | h])')
    s2_grid = grid.add_mutually_exclusive_group(required=('--read_LL_from' not in sys.argv))
    s2_grid.add_argument('--fix_s2', dest='s2',
                         help='Provide one or more fixed values for s_AA (fitness gain in homozygotes), separated by comma (without space).')
    s2_grid.add_argument('--linear_s2_range', dest='lin_s2_range',
                         help='Specify a linear grid as the parameter space for s_AA value. Format \"(min,max,grid_num)\".')
    s2_grid.add_argument('--geom_s2_range', dest='geom_s2_range',
                         help='Specify a geometrically-scaled grid as the parameter space for s_AA value. Format \"(min,max,grid_num)\". A closest number to <grid_num> will be determined such that 10^(-n)-fold values of <min> or <max> are included and that the smallest absolute value is greater than 1/2Ne. When <min> and <max> are of opposite signs, s=0 will be included, and two grids with an equal number of grid points will be formed on either side of zero.')
    s1h_grids = grid.add_mutually_exclusive_group(required=('--read_LL_from' not in sys.argv))
    s1_grid = s1h_grids.add_mutually_exclusive_group()
    s1_grid.add_argument('--fix_s1', dest='s1',
                         help='Provide one or more fixed values for s_Aa (fitness gain in heterozygotes), separated by comma (without space).')
    s1_grid.add_argument('--linear_s1_range', dest='lin_s1_range',
                         help='Specify a linear grid as the parameter space for s_Aa value. Format \"(min,max,grid_num)\".')
    s1_grid.add_argument('--geom_s1_range', dest='geom_s1_range',
                         help='Specify a geometrically-scaled grid as the parameter space for s_Aa value. Format \"(min,max,grid_num)\". A closest number to <grid_num> will be determined such that 10^(-n)-fold values of <min> or <max> are included and that the smallest absolute value is greater than 1/2Ne. When <min> and <max> are of opposite signs, s=0 will be included, and two grids with an equal number of grid points will be formed on either side of zero.')
    h_grid = s1h_grids.add_mutually_exclusive_group()
    h_grid.add_argument('--fix_h', dest='h',
                        help='Provide one or more fixed values for dominant coefficient h, separated by comma (without space).')
    h_grid.add_argument('--linear_h_range', dest='lin_h_range',
                        help='Specify a linear grid as the parameter space for dominant coefficient h. Format \"(min,max,grid_num)\".')
    h_grid.add_argument('--geom_h_range', dest='geom_h_range',
                        help='Specify a geometrically-scaled grid as the parameter space for s_AA value. Format \"(min,max,grid_num)\". A closest number to <grid_num> will be determined such that 10^(-n)-fold values of <min> or <max> are included and that the smallest absolute value is greater than 1/2Ne. When <min> and <max> are of opposite signs, s=0 will be included, and two grids with an equal number of grid points will be formed on either side of zero.')

    # output
    outputs = parser.add_argument_group('Output')
    outputs.add_argument('-o', '--out_prefix', required=True, dest='outPrefix',
                         help='Path and prefix of the output file.')
    outputs.add_argument('--long_LL_table', dest="long_table", action="store_true",
                         help="Option to write out the computed loglikelihood surfaces to a plotting-friendly table with columns \"locus, s1, s2, log-likelihood\", instead of matrices.")
    outputs.add_argument('--gzip_surface', dest='zip_output', action='store_true',
                         help='Option to output the .gz file of computed likelihoods.')
    outputs.add_argument('--get_on_grid_max', dest='onGrid_max', action="store_true", default=False,
                         help='Option to output on-grid max log-likelihood and the matching (s1, s2) values among all values computed on each given locus, as `{PREFIX}_on-grid_maxLLs.txt`.')
    outputs.add_argument('--get_off_grid_max', dest='interpolate_max', action='store_true', default=False,
                         help='Option to interpolate the computed grid (must be at least 5x5) and infer local maximum likelihoods, written to file `{PREFIX}_off-grid_maxLLs.txt`.')
    outputs.add_argument('--get_MLR', dest='get_MLR', default=False, action='store_true',
                         help='Option to write out the max-loglikelihood ratios relative to neutrality, i.e. 2(LL_opt - LL_0), for each locus, along side with the maximized log-likelihoods. ')
    outputs.add_argument('--get_chi2_pval', dest='get_pval', default=False, action='store_true',
                         help='Option to write out the p-value of MLR in a standard chi-square distribution (df=1).')

    # some HMM parameters
    HMMparam = parser.add_argument_group('HMM Parameters')
    HMMparam.add_argument('--numStates', dest='numStates', default=2000,
                          type=int, help='Number of states across [0,1) for allele frequency. Default is 2000.')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', default=False,
                        help='Option to write out selection coefficients when computing likelihoods.')

    return parser


def main(DL_args=None):
    # define parser
    parser = argparse.ArgumentParser(
        prog='DiploLocus_likelihood.py',
        # usage='%(prog)s [-h] (args for Population Parameter) (args for Input) (args for Grid of selection coefficients) (args for Output) (args for Initial Condition) [args for HMM parameters]',
        description='Light-weight toolkit for the inference of diffusion-based diploid selection on independent loci from time-stratified allele frequencies.',
        epilog='https://github.com/bioXiaoheng/diplo_locus'
    )

    lik_parser = lik_args_parser(parser)

    if DL_args is not None:
        # receive args passed from DL.py
        args = DL_args
    else:
        # receive commandline args directly w/ DL_lik.py
        # args = parser.parse_args(sys.argv[1:])
        args = lik_parser.parse_args(sys.argv[1:])

    # return usage if args blank
    if len(sys.argv[2:]) == 0:
        lik_parser.print_usage()
        sys.exit()

    # set default values
    if args.u10 is None:
        args.u10 = args.u01

    # initialize head notes
    notes = f'## DiploLocus-likelihood. u01 = {args.u01}, u10 = {args.u10}, Ne = {args.Ne}\n'

    # read pre-computed LLs if provided
    # header: locus s1 s2 LL
    if args.onGrid_LL_file is not None:
        print(args.onGrid_LL_file)
        notes += f'## Pre-computed log-likelihoods from {args.onGrid_LL_file}\n'
        assert (args.interpolate_max or args.onGrid_max) is True
        # Hey pd.read_csv is actually faster than np.loadtxt!
        onGrid_LLs = pd.read_csv(args.onGrid_LL_file, delimiter="\t", comment="#", header=[0])  # , index_col=0
        if "s1" in onGrid_LLs.columns:
            # sanity check the header for long (ggplot-fitting) tables
            assert np.all(onGrid_LLs.columns == ["ID", "s1", "s2", "loglikelihood"]), onGrid_LLs.head()
            LLmatrix, locus_names, s1_list_read, s2_list_read, s_pairs_read, num_pairs_read = _reformat_longLLtable_to_matrix(
                onGrid_LLs)
        else:
            # numpy-saved np.matrix
            # here, s_pairs_read is colname, locus_names is rowname
            LLmatrix, locus_names, s1_list_read, s2_list_read, s_pairs_read, num_pairs_read = _reformat_LL_DF_to_matrix(
                onGrid_LLs)
        numsites = len(locus_names)
        notes += f'## Recorded {numsites} loci, under {num_pairs_read} pairs of selection coefficients.\n'
        # if grid args are provided
        grid_args = (args.s2, args.lin_s2_range, args.geom_s2_range,
                     args.s1, args.lin_s1_range, args.geom_s1_range,
                     args.h, args.lin_h_range, args.geom_h_range)
        grid_args_set = set(grid_args)
        if len(grid_args_set) > 1 or grid_args_set != {None}:
            # generate grid from args
            s1_list, s2_list, s_pairs, num_pairs = parse_grid_args(args.s2, args.lin_s2_range, args.geom_s2_range,
                                                                   args.s1, args.lin_s1_range, args.geom_s1_range,
                                                                   args.h, args.lin_h_range, args.geom_h_range, args.Ne)
            try:
                # specified grids must match what's pre-computed
                assert np.all(np.isclose(s1_list, s1_list_read)) and np.all(np.isclose(s2_list,
                                                                                       s2_list_read)), "Grid points in the file do not match the grid defined by arguments."
                assert num_pairs == num_pairs_read, "Grid points in the file do not match the grid defined by arguments."
                # they must be of the same length past the previous checkpoint
                assert np.all([(np.isclose(s1i, s1j) and np.isclose(s2i, s2j)) for ((s1i, s2i), (s1j, s2j)) in
                               list(zip(s_pairs_read,
                                        s_pairs))]), "Grid points in the file do not match the grid defined by arguments."
            except AssertionError:
                # go through each pair & find match (need a more efficient alt)
                LL_subset_colids = []
                for s1, s2 in s_pairs:
                    # let's just look at .6f
                    find_pair = [(s1, s2, pair_idx) for pair_idx, (s1_read, s2_read) in enumerate(s_pairs_read)
                                 if np.isclose(s1_read, s1, atol=1e-6) and np.isclose(s2_read, s2, atol=1e-6)]
                    if len(find_pair) == 0:
                        raise ValueError(
                            f"Specified grid point ({s1}, {s2}) not in file {args.onGrid_LL_file}.")  # \n{s_pairs}\n{s_pairs_read}
                    else:
                        assert len(find_pair) == 1, find_pair
                        # (s1, s2, pair_idx) = find_pair
                        LL_subset_colids.append(find_pair[0][2])
                # after looping through, lengths should match
                assert np.all(
                    [(np.isclose(s_pairs_read[i][0], s1j) and np.isclose(s_pairs_read[i][1], s2j)) for (i, (s1j, s2j))
                     in list(zip(LL_subset_colids,
                                 s_pairs))]), "Specified grid is not a subset of the grid provided in the file."
                LLmatrix = LLmatrix[:, LL_subset_colids]
                assert LLmatrix.shape == (numsites, num_pairs)
        else:
            s1_list, s2_list, s_pairs, num_pairs = s1_list_read, s2_list_read, s_pairs_read, num_pairs_read
    else:
        # generate grid
        s1_list, s2_list, s_pairs, num_pairs = parse_grid_args(args.s2, args.lin_s2_range, args.geom_s2_range,
                                                               args.s1, args.lin_s1_range, args.geom_s1_range,
                                                               args.h, args.lin_h_range, args.geom_h_range, args.Ne)
        # print(f's1_list={s1_list}\ns2_list={s2_list}\ns_pairs={s_pairs}.')
        notes += f'## Analyze {num_pairs} pairs of sel. coef.\n'

        # get input
        if args.infile is not None:
            print(time.ctime(), f'\nLoading parsed allele counts from {args.infile}.')
            notes += f'## allele counts from {args.infile}\n'
            # parse snps
            if args.snp_subset is not None:
                snp_list = parse_ind_arg(args.snp_subset)
                notes += f'## {len(snp_list)} loci from `--snps`'
            else:
                snp_list = None
            locus_names, samples, sampleSizes = parse_allele_counts(args.infile, snp_list, minMAF=args.MAF)
            # sample times must be provided
            sampleTimes, sample_cols_toKeep = parse_sampling_times(args.sampleTimes, args.sampleTimes_ago,
                                                                   args.gen_time, args.t0, args.force_t0)
            num_time_points = len(sampleTimes) - 1
            # crop samples & sampleSizes too
            samples, sampleSizes = samples[:, sample_cols_toKeep], sampleSizes[:, sample_cols_toKeep]
        else:
            # sanity check
            assert args.vcfFile is not None, TypeError('Please provide an input file.')
            # check for info file
            assert args.infoFile is not None, TypeError('Please provide a matching info file for the VCF.')
            print(time.ctime(), f'Loading sample times from file {args.infoFile}.')

            # parse inds
            if args.ind_subset is not None:
                # or args.snp_subset is not None:
                ind_list = parse_ind_arg(args.ind_subset, "inds")
                notes += f'## {len(ind_list)} samples from `--inds`\n'
            else:
                ind_list = None

            # parse snps
            if args.snp_subset is not None:
                snp_list = parse_ind_arg(args.snp_subset, "snps")
                notes += f'## {len(snp_list)} loci from `--snps`\n'
            else:
                snp_list = None
            # notes += '## from\n'
            # to be followed by vcf file
            notes += f'## VCF: {args.vcfFile}\n## meta info: {args.infoFile}\n## generation time {args.gen_time} units,'

            # get info
            if args.TimeColName is None and args.TimeAgoColName is None:
                TimeColName = None
            elif args.TimeColName is None:
                TimeColName = args.TimeAgoColName
                find_ago = ago_regex.findall(TimeColName)
                if len(find_ago) == 0:
                    TimeColName = TimeColName + '_Ago'
                else:
                    try:
                        assert len(find_ago) == 1, f'Column name for sampling time (backward): {args.TimeAgoColName}'
                    except Exception as e:
                        print(e)
                        print(f'Column name for sampling time (backward): {args.TimeAgoColName}')
            elif args.TimeAgoColName is None:
                TimeColName = args.TimeColName
            else:
                raise ValueError('Please only specify one column to be used as sampling time.')

            sampleTimes, sampleIDPools = read_info_file(args.infoFile, ID_colname=args.IDcolName,
                                                        time_colname=TimeColName, t0=args.t0, gen_time=args.gen_time,
                                                        force_t0=args.force_t0, inds=ind_list)
            num_time_points = len(sampleIDPools)
            # sanity check
            assert (len(sampleTimes) - 1) == len(sampleIDPools), print(sampleTimes, sampleIDPools)
            # print(f'\tTaking samples at {len(sampleTimes[1:])} time points: {", ".join(map(str, sampleTimes[1:]))}.')

            # now check ploidy
            ## read all GTs as-is unless flags specify
            all_IDs = {ID for pool in sampleIDPools for ID in pool}
            if args.forced_haps != 'none' and args.forced_haps != "all":
                hap_IDs = parse_ind_arg(args.forced_haps)
                if len(hap_IDs) == 0:
                    print('Invalid input for `--force_hap`:', args.forced_haps)
                else:
                    hap_IDs = set(hap_IDs)
                    print(
                        f'{len(hap_IDs)} samples will be deliberately counted as haploids: {", ".join(list(hap_IDs))}')
            elif args.forced_haps == "all":
                hap_IDs = all_IDs
                print(f'All samples will be deliberately counted as haploids.')
            else:
                hap_IDs = {}

            if args.forced_dips != "none" and args.forced_dips != "all":
                dip_IDs = parse_ind_arg(args.forced_dips)
                if len(dip_IDs) == 0:
                    print('Invalid input for `--force_dip`:', args.forced_dips)
                else:
                    dip_IDs = set(dip_IDs)
                    print(
                        f'{len(dip_IDs)} samples will be deliberately counted as diploids: {", ".join(list(dip_IDs))}')
            elif args.forced_dips == "all":
                dip_IDs = all_IDs
                print(f'All samples will be deliberately counted as diploids.')
            else:
                dip_IDs = {}

            # sanity check
            if len(hap_IDs) > 0 and len(dip_IDs) > 0:
                assert not hap_IDs.isdisjoint(
                    dip_IDs), f'Sample(s) {hap_IDs.intersection(dip_IDs)} cannot be forced to be diploid and haploid at the same time.'
            # all IDs should be subset of all_IDs
            # print(hap_IDs, dip_IDs)
            assert set(hap_IDs).issubset(all_IDs) and set(dip_IDs).issubset(
                all_IDs), f'Sample IDs specified in `--force_hap` and `--force_dip` must already be included in the VCF. {hap_IDs - all_IDs} {dip_IDs - all_IDs} are not included.'

            # parse vcf
            print(f'Loading variant counts from {args.vcfFile}...')
            # the `trimmedSampleTimes` here already have 0 added too
            locus_names, samples, sampleSizes, trimmedSampleTimes, het_flags, trimmedSampleIndexPools = parse_vcf_input(
                args.vcfFile,
                sampleIDPools,
                sampleTimes=sampleTimes,
                snps=snp_list,
                forced_haps=hap_IDs,
                forced_dips=dip_IDs,
                minMAF=args.MAF)
            if len(trimmedSampleTimes) != len(sampleTimes):
                assert len(trimmedSampleTimes) < len(
                    sampleTimes), f'Before SNP-based trimming, sample times: {sampleTimes}\nAfter: {trimmedSampleTimes}'
                # sanity check
                trimmed = []
                kept = 0
                for k, time_point in enumerate(sampleTimes[1:]):
                    if time_point not in trimmedSampleTimes:
                        trimmed.append(k)
                    else:
                        kept += 1
                assert (
                                   len(trimmed) + kept) == num_time_points, f'len(trimmed)={len(trimmed)}, kept={kept}, num_time_points={num_time_points}'
                # matrix shape should match too
                assert samples.shape[1] == len(trimmedSampleTimes) - 1
                assert samples.shape == sampleSizes.shape
                # update other related variables
                num_time_points = len(trimmedSampleTimes) - 1
                sampleTimes = trimmedSampleTimes
                # sampleIDPools = [header[idx] for pool in trimmedSampleIndexPools for idx in pool] # <-- no use

        # decide emission type & prep relevant args
        samples_int = np.allclose(np.array(samples, dtype=float) - np.array(samples, dtype=int), 0)
        sampleSizes_int = np.allclose(np.array(sampleSizes, dtype=float) - np.array(sampleSizes, dtype=int), 0)
        # current implementation only supports int sample sizes:#
        if not sampleSizes_int:  # or sampleSizes.dtype not in (int, np.int)
            print(f'type(sampleSizes)={sampleSizes.dtype}. Enforcing integer sample sizes...')
            sampleSizes = np.array(sampleSizes, dtype=int)
        else:
            # sampleSizes is a list. Need to make sure to make it np.array
            sampleSizes = np.array(sampleSizes, dtype=int)
        # likewise, make sure to convert samples to np.array too
        if samples_int:
            emType = 'integer'
            # sampleSizeSet = set([tuple(locus) for locus in sampleSizes]) # <-- this is not what the algorithm expect
            sampleSizeSet = set(sampleSizes.flatten())
            samples = np.array(samples, dtype=int)
        else:
            emType = 'fractional'
            sampleSizeSet = None
            samples = np.array(samples)

        # for later
        numsites = len(locus_names)

        # double-checking
        assert samples.shape == sampleSizes.shape
        assert samples.shape[1] == len(sampleTimes) - 1, f'samples.shape={samples.shape}:{sampleSizes};\n ' \
                                                         f'len(sampleTimes)={len(sampleTimes)}: {sampleTimes}'
        maxSampleSizes = np.apply_along_axis(max, axis=0, arr=sampleSizes)
        print(
            f'Computing likelihoods for {numsites} loci at {num_time_points} time points: {", ".join(map(str, sampleTimes[1:]))}\n\tMax sample sizes at each time: {", ".join(map(str, maxSampleSizes))}')
        #: {sampleTimes[1:]} generations.
        notes += f'## Computing likelihoods for {numsites} loci at {num_time_points} time points: {", ".join(map(str, sampleTimes[1:]))}\n## \tMax sample sizes at each time: {", ".join(map(str, maxSampleSizes))}\n'

        # need to do this again now that sampleSizes updated
        sampleSizeSet = set(sampleSizes.flatten())

        # read piecewise info if needed
        if args.piecewise is False:
            notes += f'## Assume constant selection pressure throughout the duration considered\n'

            # obtain the wrapper function based on initCond args
            def _computeLL(sampleTimes, sampleSizes, samples, s1, s2):
                HMMcore = likelihood.SelHmm(args.Ne, s1, s2, args.u01, args.u10, initCond=args.initCond,
                                            initFreq=args.initFreq, sampleSizesSet=sampleSizeSet,
                                            emissionType=emType, transitionType='constant')
                LLs = HMMcore.computeLogLikelihood(sampleTimes, sampleSizes, samples)
                return LLs
        else:
            assert args.specify_piece is not None
            # read the file
            print(f'Loading stepwise time-varying selection intervals from {args.specify_piece}.')
            which_piece = np.loadtxt(args.specify_piece, delimiter='\t', comments="#",
                                     dtype={'names': ('start', 'end', 's1', 's2'),
                                            'formats': ('i', 'i', 'O', 'O')})
            which_piece = list(which_piece)
            num_pieces = len(which_piece)
            piece_starts, piece_ends, s1_pieces, s2_pieces = map(list, zip(*which_piece))
            selChangeTimes = [int(t) for t in piece_starts[1:]]
            # require that end_{t-1} == start_{t}
            assert np.allclose(piece_starts[1:], piece_ends[:-1]), ValueError(
                'Please make sure the starting time be equal to the previous ending time.')
            # we'll require both are "x" for now. In case they only specify one, the fact that the other one has a fixed value would be reflected in their grid parameters anyway.
            # check1: only one interval allowed
            assert s1_pieces.count("x") == 1 and s2_pieces.count("x") == 1, ValueError(
                f"Please choose only one interval to examine.")
            # check2: two x's at the same interval
            vary_i = s1_pieces.index("x")
            assert vary_i == s2_pieces.index("x"), ValueError(
                "Please use \"x\" for *both* s1 and s2 in the time interval to be examined.")

            # sanity check: last gen# match with sampleTimes
            assert which_piece[-1][1] >= sampleTimes[-1], ValueError(
                'Please make sure the specify_piece file covers every generation across all the sampling timepoints.')

            # assert num_pieces == len(selChangeTimes)+1 # this is already checked in likelihood.py
            print(
                f'\nDiploLocus will compute likelihoods for varying (s1,s2) values from generation {which_piece[vary_i][0]} to {which_piece[vary_i][1]}.\n')
            notes += f'## Assume distinct (constant) sel. coef. during time intervals {list(piece_starts) + [piece_ends[-1]]}: {zip(s1_pieces, s2_pieces)}\n'

            # obtain the wrapper function based on initCond and piecewise args
            def _computeLL(sampleTimes, sampleSizes, samples, s1, s2):
                s1s = [float(s1_pieces[j]) if j != vary_i else s1 for j in range(num_pieces)]
                s2s = [float(s2_pieces[j]) if j != vary_i else s2 for j in range(num_pieces)]
                # print(f'''s1s={s1s}, s2s={s2s},u01={args.u01}, u10={args.u10}, initCond={args.initCond}, initFreq={args.initFreq}, initMAlpha={args.init_u01}, initMBeta={args.init_u10}, initS1={args.init_s1,} initS2={args.init_s2}, emissionType={emType}, transitionType='piecewise', selectionChangeTimes={selChangeTimes}''')
                if emType == "integer":
                    sizeSet = set(np.array(sampleSizes).flatten())
                    HMMcore = likelihood.SelHmm(args.Ne, s1s, s2s, args.u01, args.u10, initCond=args.initCond,
                                                initFreq=args.initFreq, sampleSizesSet=sizeSet, emissionType="integer",
                                                transitionType='piecewise', selectionChangeTimes=selChangeTimes)
                    # initMAlpha=args.init_u01, initMBeta=args.init_u10, initS1=args.init_s1, initS2=args.init_s2,
                else:
                    HMMcore = likelihood.SelHmm(args.Ne, s1s, s2s, args.u01, args.u10, initCond=args.initCond,
                                                initFreq=args.initFreq, emissionType=emType, transitionType='piecewise',
                                                selectionChangeTimes=selChangeTimes)
                LL = HMMcore.computeLogLikelihood(sampleTimes, sampleSizes, samples)
                return LL

        # use the pre-defined wrapper function to compute LLs on the grid
        # axis-0 (row): position, axis-1 (col): s_pair
        LLmatrix = np.empty((numsites, num_pairs))
        for i, (s1, s2) in enumerate(s_pairs):
            if args.verbose and (i % 5 == 0): print(time.ctime(), f'Computing likelihood at ({s1}, {s2}).')
            LL = _computeLL(sampleTimes, sampleSizes, samples, s1, s2)
            LLmatrix[:, i] = LL

        # write output
        outputDF = pd.DataFrame(LLmatrix, columns=s_pairs)  # , index=locus_names
        outputDF['ID'] = locus_names  # outputDF.index
        if args.long_table:
            outputDF = outputDF.melt(id_vars='ID', var_name='s_pair', value_name='loglikelihood')
            outputDF['s1'], outputDF['s2'] = zip(*outputDF.s_pair)
            outputDF = outputDF[['ID', 's1', 's2', 'loglikelihood']]
            outfilename = args.outPrefix + '_longLLtable.tsv'
        else:
            outfilename = args.outPrefix + '_LLmatrices.table'

        if args.zip_output:
            outfilename += '.gz'
            import gzip
            with gzip.open(outfilename, "wb") as outhandle:
                outhandle.write(notes.encode())
            outhandle.close()
            outputDF.to_csv(outfilename, mode='ab', sep="\t", compression='gzip', index=False)  # , index_label='ID'
        else:
            with open(outfilename, "w") as outhandle:
                outhandle.write(notes)
            outhandle.close()
            outputDF.to_csv(outfilename, mode='a', sep="\t", index=False)  # , index_label='ID'

    # decide whether to interpolate & output stand-alone on-grid max
    write_onGrid_max = args.onGrid_max
    if args.get_pval:
        get_MLR = True
        from scipy.stats import chi2
    else:
        get_MLR = args.get_MLR

    # This only works for 1D for now
    if args.interpolate_max:
        # try:
        maxLLs, header = interpolate_offGrid_max(s1_list, s2_list, s_pairs, num_pairs, LLmatrix)
        # see if MLRs are needed
        if get_MLR:
            neutLL = LLmatrix[:, s_pairs.index((0, 0))]
            MLRs = (maxLLs[:, -1] - neutLL) * 2
            maxLLs = np.hstack((maxLLs, MLRs.reshape(numsites, 1)))
            header += '\tMLR'
            if args.get_pval:
                # this only works for h=0.5 (symmetry) for now
                ## hence the sanity check
                h_check = list(set([(s1 / s2) for (s1, s2) in s_pairs if s2 != 0]))
                if len(h_check) == 1 and np.isclose(h_check[0], 0.5):
                    pvals = chi2.sf(MLRs, df=1)
                    maxLLs = np.hstack((maxLLs, pvals.reshape(numsites, 1)))
                    header += '\tchi2_p'
                else:
                    print(f'h_check={h_check}, P-values from standard Chi-square only applies to additive (h=0.5)')
        # attach locus_names
        maxLLs = np.hstack((np.array(locus_names).reshape(numsites, 1), maxLLs))
        # write output
        maxLL_filename = args.outPrefix + '_off-grid_maxLLs.txt'
        np.savetxt(maxLL_filename, maxLLs, header=(notes + header), delimiter='\t', fmt='%s')
        # no need to output a separate file for it
        if write_onGrid_max:
            print('On-grid max are already included in the `_off-grid_maxLLs.txt` output file.')
            write_onGrid_max = False

    if write_onGrid_max:
        header = '\t'.join(['ID', 'ongrid_s1hat', 'ongrid_s2hat', 'ongrid_maxLogLikelihood'])
        maxLLs = get_onGrid_max_only(s_pairs, LLmatrix)
        # attach locus_names
        numsites = len(locus_names)
        # see if MLRs are needed
        if args.get_MLR:
            neutLL = LLmatrix[:, s_pairs.index((0, 0))]
            MLRs = (maxLLs[:, -1] - neutLL) * 2
            maxLLs = np.hstack((maxLLs, MLRs.reshape(numsites, 1)))
            header += '\tMLR'
        maxLLs = np.hstack((locus_names.reshape(numsites, 1), maxLLs))
        # write output
        maxLL_filename = args.outPrefix + '_on-grid_maxLLs.txt'
        np.savetxt(maxLL_filename, maxLLs, header=(notes + header), delimiter='\t', fmt='%s')

    print(time.ctime(), 'Done.')


if __name__ == '__main__':
    main()
