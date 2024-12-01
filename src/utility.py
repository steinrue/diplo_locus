import numpy
import pandas as pd
pd.set_option('display.max_colwidth', 255)
import re, os, sys
import scipy.interpolate
import scipy.optimize


def _get_geom_grid(left: float, right: float, grid_num, Ne):
    """Return a 1-D numpy.array of geometrically scaled grid. If [`left`, `right`] covers zero,
    either side will be (log10-)scaled with equal number of values, such that the total number
    of grid points is the closest to `grid_num`.
    Smallest absolute value be of the same (log10-)scale as 1/2Ne
    :param left: Left bound (inclusive; smallest value) of the grid
    :param right: Right bound (inclusive; largest value) of the grid
    :param grid_num: Number of grid points (or its closest approx. if both signs are considered)
    :param Ne: Effective diploid pop size
    :return: numpy.array(s_grid)
    """
    # assume symmetric. If not, make it so.
    if left < 0 and right > 0:
        bound = max(right, -1 * left)
        digits = round(numpy.log10(bound / (0.5 / Ne)))
        n = (grid_num / 2 - 1) // digits  # the # of grid pts bewteen any 10^k and 10^(k-1)
        pos_half = numpy.geomspace(right / (10 ** digits), right, int(n * digits + 1))
        neg_half = numpy.geomspace(left, left / (10 ** digits), int(n * digits + 1))
        s_range = numpy.hstack((neg_half, 0, pos_half))

    elif left * right > 0:  # 0 will be automatically included
        digits = round(numpy.log10(right / left))
        n = (grid_num - 1) // digits  # the # of grid pts bewteen any 10^k and 10^(k-1)
        s_range = numpy.sort(numpy.hstack((numpy.geomspace(left, right, int(n * digits + 1)), 0)))

    elif left * right == 0:
        bound = max(abs(left), abs(right))
        sign = ((left + right) > 0) * 2 - 1
        s_range = _get_geom_grid(sign * 1e-3 * bound, sign * bound, grid_num, Ne)

    else:
        print(f'Please double-check your input for geometric grid [{left}, {right}] grid_num = {grid_num}.')
        sys.exit()

    return s_range


# infile header needs to be fixed  ##, minMAF: float = 0.
def parse_allele_counts(infile: str, snp_list=None,  # minMAF=args.MAF,
                        chr_col: str=None, pos_col: str=None, snpID_col: str=None):
    """Parse the allele count file. Must have header.
    12/12/2023 Update: remove minMAF filter. Add column name option for chr, pos, and ID
    Return `numpy.array(locus_names)`, `numpy.array(samples)`, `numpy.array(sampleSizes)`
    """
    parsed = pd.read_csv(infile, sep="\t", header=0, comment="#")
    header = parsed.columns
    # retrieve snp names
    # pos_regex = re.compile(r'(locus|ch|pos|position|id|rep)')
    # id_cols = [colname for colname in header if pos_regex.search(colname.lower())]
    id_cols = [colname for colname in header if (colname == 'ID')]
    assert (len(id_cols) == 1), "Exactly one column 'ID' with unique identifier needed in input allele counts file."
    # see if any of the column names are provided:
    if chr_col is not None:
        assert chr_col in id_cols, f'Please provide `--chrom_col` column name that exactly match the input file. '
        if chr_col != 'Chr':
            parsed = parsed.rename({chr_col: 'Chr'}, axis='columns')
    if pos_col is not None:
        assert pos_col in id_cols, f'Please provide `--pos_col` column name that exactly match the input file. '
        if pos_col != 'physPos':
            parsed = parsed.rename({pos_col: 'physPos'}, axis='columns')
    if snpID_col is not None:
        assert snpID_col in id_cols, f'Please provide `--snp_col` column name that exactly match the input file. '
        if snpID_col != 'snpID':
            parsed = parsed.rename({chr_col: 'snpID'}, axis='columns')
    # re-collect id columns if anything changes:
    # if chr_col or pos_col or snpID_col:
    #     id_cols = [colname for colname in header if pos_regex.search(colname.lower())]
    # then decide what becomes `locus_names`
    ## prioritize SNP ids
    # if snpID_col is not None:
    #     if snp_list is not None:
    #         parsed = parsed[parsed.snpID.isin(snp_list)]
    #     locus_names = numpy.array(parsed.loc[:, 'snpID'])
    # elif 'locus' in id_cols:
    #     if snp_list is not None:
    #         parsed = parsed[parsed.locus.isin(snp_list)]
    #     locus_names = numpy.array(parsed.locus)
    # ## if both chr and pos are provided, then proceed with that
    # elif (chr_col is not None) and (pos_col is not None):
    #     parsed['locus_name'] = parsed.loc[:, ['Chr', 'physPos']].apply(lambda l: '_'.join([str(c) for c in l]), axis=1)
    #     # now down sample if needed
    #     if snp_list is not None:
    #         parsed = parsed[parsed.locus_name.isin(snp_list)]
    #     locus_names = numpy.array(parsed.locus_name)
    # ## in all other cases, concatenate all id columns:
    # else:
    parsed['locus_name'] = parsed.loc[:, id_cols].apply(lambda l: '_'.join([str(c) for c in l]), axis=1)
    # remove old ID column
    parsed = parsed.drop(columns=['ID'])
    # now down sample if needed
    if snp_list is not None:
        parsed = parsed[parsed.locus_name.isin(snp_list)]
    locus_names = numpy.array(parsed.locus_name)

    # just some sanity check
    if snp_list is not None:
        assert parsed.shape[1] > 0, f'No intersection between variants in allele count file and ' \
                                    f'those listed by `--snps`. Please double check your input, then' \
                                    f'specify the column name for variant IDs with `--snp_col`.'
        sys.exit()

    # more sanity checks
    # print (parsed.shape)
    # retrieve all xk column names
    # d_regex = re.compile(r'^[xd](\d+)$')
    d_regex = re.compile(r'^d(\d+)$')
    n_regex = re.compile(r'^n(\d+)$')
    d_cols = [colname for colname in header if d_regex.search(colname)]
    n_cols = [colname for colname in header if n_regex.search(colname)]
    # sanity check
    assert len(d_cols) == len(n_cols), f'Equal number of samples and sample sizes needed: d_cols={d_cols}; n_cols={n_cols}'
    assert (len(d_cols) >= 2), "Need at least two sampling batches in input file. Maybe file is not tab-seperated."
    # print (parsed.columns)
    assert (1 + len(d_cols) + len(d_cols) == parsed.shape[1]), f"Columns in input file that are not 'ID', samples sizes (nx), or num dervived alleles (dx). Check that file is tab-seperated, or remove before proceeding. ({1 + len(d_cols) + len(d_cols)}, {parsed.shape[1]})"
    # sanity check: the indices are monotonously ordered
    try:
        d_col_num = [int(d_regex.findall(colname)[0]) for colname in header if d_regex.search(colname)]
        n_col_num = [int(n_regex.findall(colname)[0]) for colname in header if n_regex.search(colname)]
    except Exception as e:
        print(e)
        print(f'Cannot parse column names: d_cols={d_cols}; n_cols={n_cols}')
        sys.exit(1)
    assert numpy.all(numpy.array(d_col_num) == numpy.array(n_col_num)), \
        f'Columns for samples and sample sizes must be ordered the same.'
    assert numpy.all(numpy.diff(d_col_num) >= 0) or numpy.all(numpy.diff(d_col_num) <= 0), \
        f'Please make sure the sample pools are ordered.'
    # assemble
    samples = numpy.array(parsed.loc[:, d_cols])
    sampleSizes = numpy.array(parsed.loc[:, n_cols])

    # down-sample by <del>pooled MAF</del>
    ## remove ones without observations first: (mostly to avoid zero division
    observed = (sampleSizes.sum(axis=1) > 0)
    samples = samples[observed, :]
    sampleSizes = sampleSizes[observed, :]
    locus_names = locus_names[observed]
    ## divide
    # freqs = samples.sum(axis=1) / sampleSizes.sum(axis=1)
    # ## fold
    # freqs = numpy.where(freqs < 0.5, freqs, 1 - freqs)
    # # filter
    # passed_rows = (freqs > minMAF)
    # down sample
    # samples = samples[passed_rows, :]
    # sampleSizes = sampleSizes[passed_rows, :]
    # locus_names = locus_names[passed_rows]

    return locus_names, samples, sampleSizes


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
    ## remove empty lines
    ind_list = [ind for ind in ind_list if len(ind.strip()) > 0]
    print(f'{len(ind_list)} {type} in total.')
    return ind_list


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
    :return: sampleTimes
    :return: samplePools

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
        if numpy.any(infotable[timeCol].isnull()):
            infotable = infotable[~infotable[timeCol].isnull()]
        if t0 is None:
            t0 = min(infoTable['Gen'])
            # t0 = 0
            print('Assuming selection starts at generation zero.')
        elif not force_t0:
            assert t0 <= min(infotable[
                                 "Gen"]), 'Please make sure selection start time `t0` is of the same unit as provided under the \"Gen\" column.\nIf selection is set to start after the oldest samples are collected, please use `--force_t0` to allow omitting data preceding the selection onset.'
        # infoTable['Gen'] = numpy.array(numpy.array((infoTable['Gen']) - t0), dtype=int)
        # convert to generation
        infotable['Gen'] = (infotable[timeCol] - t0) / gen_time
        # just in case
        if numpy.any(infotable["Gen"].isnull()):
            infotable = infotable[~infotable["Gen"].isnull()]
        infotable['Gen'] = infotable['Gen'].round().astype(int)
        # done processing
        return infotable

    def _count_time_backward(infotable, timeCol, t0):
        # backward time
        print(f'Counting time from present to past, extracting info from \"{timeCol}\" column.')
        ## check for missing value
        if numpy.any(infotable[timeCol].isnull()):
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
        if numpy.any(infotable["Gen"].isnull()):
            infotable = infotable[~infotable["Gen"].isnull()]
        # print(f't0={t0}, {infoTable["Time_Ago"].describe()}, gen_time={gen_time}') #\n{list(infoTable["Time_Ago"])}
        # algorithm only accept integer generations for now
        infotable['Gen'] = infotable['Gen'].round().astype(int)
        # done processing
        return infotable

    ago_regex = re.compile(r'([A|a][G|g][O|o])')
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
                assert time_colname in infoTable.columns, \
                    f'Column \"{time_colname}\" missing in the info table.'
            infoTable = _count_time_backward(infoTable, time_colname, t0)
    # if any column has "ago" in it
    elif numpy.any(numpy.array([ago_regex.search(name) for name in infoTable.columns], dtype=bool)):
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
    if numpy.any(numpy.array(infoTable["Gen"]) < 0):
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
    sampleTimes = numpy.array([0] + time_points)
    samplePools = [list(pool) for pool in sample_pools]
    return sampleTimes, samplePools


## looks like I need to write a function for flatten
def _flatten_list_of_lists(manyPools: list):
    flatPool = []
    for pool in manyPools:
        flatPool.extend(pool)
    return flatPool


def parse_forced_ploidy_arg(forced_ones, all_samples, Type):
    # start with sanity check
    # assert isinstance(forced_ones, str), \
    #     f'Unrecognized {Type} argument value type {type(forced_ones)}: {forced_ones}'
    if isinstance(forced_ones, str):
        # now they're all strings
        if forced_ones == "all":
            return set(all_samples)
        elif forced_ones == "none":
            return set([])
        else:
            return set(parse_ind_arg(forced_ones, Type))
    elif isinstance(forced_ones, set):  # or isinstance(forced_ones, dict):
        return forced_ones
    elif isinstance(forced_ones, dict):
        return set(forced_ones.keys())
    else:
        raise ValueError(f'Unrecognized `{Type}` argument value type {type(forced_ones)}: {forced_ones}')


import allel


def parse_vcf_input(vcffile: str, samplePools: list, sampleTimes=None, inds=None,
                    forced_haps: list or str = '', forced_dips: list or str = '',
                    snps=None, pos_range: str = '', minMAF: float = 0.):
    """Read in vcf-formatted input and generate `samples` & `sampleSizes` matrices
       for downstream LL computation. Return numpy.arrays of `locus_name`, `samples`,
       and `sampleSizes'.
    ---------------------------
    :param vcffile: str, VCF file name. Both `.vcf` and `.vcf.gz` work.
    :param samplePools: list, K-element list (K=#[time points]). The i-th element is
            the list of individual/sample IDs (as in the VCF) that constitute the sample
            pool at the i-th ancient sampling time point.
    :param snps: list of variant IDs (as a subset of all variants in the vcf) to consider.
            If not provided, all loci in the vcf will be included. Default is None.
    :param inds: list of sample IDs (as a subset of all sample IDs [as column names] in the vcf)
            to consider. If not provided, all loci in the vcf will be included. Default is None.
    :param sampleTimes: list of time points for each batch of samples
    :param forced_haps: list of sample IDs, as shown in VCF column names, to be counted as haploids
    :param forced_dips: list of sample IDs, as shown in VCF column names, to be double-counted as diploids if a GT is presented as haploid.
    :param pos_range: str, optional, argument to indicate the range of physical positions to consider.
    :param minMAF: float, *PHASING OUT* minimum threshold (non-inclusive) for minor allele frequencies in the pooled samples.
    -----------------------------
    :return: numpy.array(locus_names), numpy.array(samples), numpy.array(sampleSizes), numpy.array(sampleTimes), sampleIDpools
    """
    # setting up
    SNP = {'A', 'T', 'C', 'G'}

    # sanity check
    if sampleTimes is not None:
        assert len(sampleTimes) == (len(samplePools) + 1), \
            f'len(samplePools)={len(samplePools)}; len(sampleTimes)={len(sampleTimes)}'
    timePool_indice = range(1, len(samplePools) + 1)
    # if position range is defined
    if pos_range == '':
        pos_arg = None
    else:
        # can I just leave the parsing to `read_vcf()`?
        parts = re.split(r'\W+', pos_range)
        if len(parts) == 1:  # then this is chromosome name
            pos_arg = parts[0]
        elif len(parts) == 3:
            chrom = parts[0]
            from_pos, to_pos = map(float, parts[1:])
            pos_arg = f'{chrom}:{int(from_pos):d}-{int(to_pos):d}'
            print(pos_arg)
        else:  # let's just leave it like that. I doubt read_vcf takes open ranges
            raise ValueError(f'Cannot parse position description {pos_range}')

    # reading file
    if inds is not None and pos_arg != '':
        vcf = allel.read_vcf(vcffile, samples=inds, region=pos_arg)
    elif inds is not None:
        vcf = allel.read_vcf(vcffile, samples=inds)
    elif pos_arg != '':
        vcf = allel.read_vcf(vcffile, region=pos_arg)
        if vcf is None:
            raise ValueError(f'Please make sure chromosome name in pos_arg matches the vcf.')
    else:
        vcf = allel.read_vcf(vcffile)

    # locus_names = vcf['variants/ID']
    chrom, positions, snp_IDs = vcf['variants/CHROM'], vcf['variants/POS'], vcf['variants/ID']
    # locus_names = numpy.apply_along_axis("_".join, 1, list(zip(chrom, positions.astype(str), snp_IDs)))
    locus_names = ["_".join([str(chrom[i]), str(positions[i]), snp_IDs[i]]) for i in range(len(snp_IDs))]
    locus_names = numpy.array(locus_names)
    print(f'{len(locus_names)} snps read from VCF')
    GTs = vcf['calldata/GT']
    # downsample if specified
    if snps is not None:
        # snps_idx = numpy.where(locus_names in snps)[0]
        pick_snps = numpy.vectorize(lambda x: x in set(snps))
        # snps_included_bools = snp_IDs in snps
        snps_included_bools = pick_snps(snp_IDs)
        print(f'{sum(snps_included_bools)} snps out of {len(snp_IDs)} included after `numpy.vectorize(lambda x: x in set(snps))`.\n')
    else:  # <-- if None, then take all
        snps_included_bools = numpy.ones(len(snp_IDs)).astype(bool)
    # force bi-allelic
    # biallelic_idx = numpy.where((vcf['variants/ALT'] in SNP) & (vcf['variants/ALT'] in SNP))[0]
    is_SNP = lambda x: x in SNP
    # is_biallelic_bools = ((vcf['variants/ALT'] in SNP) & (vcf['variants/REF'] in SNP))
    is_biallelic_bools = numpy.array(list(map(is_SNP, vcf['variants/ALT'][:, 0])), dtype=bool) & \
                         numpy.all(vcf['variants/ALT'][:, 1:] == '', axis=1) & \
                         numpy.array(list(map(is_SNP, vcf['variants/REF'])), dtype=bool)
    print(f'{sum(is_biallelic_bools)} out of {len(snp_IDs)} snps are bi-allelic.\n')

    # update
    locus_names = locus_names[(snps_included_bools & is_biallelic_bools)]
    print(f'{locus_names.shape} snps after filtering for biallelic + overlaps between .snps & .vcf ')
    GTs = GTs[(snps_included_bools & is_biallelic_bools), :, :]

    ## sanity check: all GT calls are 1, 0, or -1
    all_GT_types = set(GTs.flatten())
    assert all_GT_types.issubset({0, 1, -1}), f'all_GT_types = {all_GT_types}'

    # now filter the samples if needed
    vcf_sample_pool = vcf['samples']
    ## get indices for each time pool
    ## samples in this pool are the intersection between vcf pool and info pool
    # TODO: maybe use np.where instead of list.index to make it faster? <-- after making sure all elements exist
    sampleIndexPools = [[list(vcf_sample_pool).index(sampID) for sampID in thisSamplePool if sampID in vcf_sample_pool]
                        for thisSamplePool in samplePools]
    ### these are the indice in OG vcf
    all_sample_idx = numpy.sort(_flatten_list_of_lists(sampleIndexPools)).astype(int)
    all_samples = _flatten_list_of_lists(samplePools)
    ### another sanity check
    assert len(all_sample_idx) <= len(all_samples), \
        f'len(all_sample_idx)={len(all_sample_idx)}; len(samplePools)={len(all_samples)}'

    # now parse forced_hap/dips
    ## these are sets of IDs
    forced_hap_IDs = parse_forced_ploidy_arg(forced_haps, all_samples, "--forced_hap")
    forced_dip_IDs = parse_forced_ploidy_arg(forced_dips, all_samples, "--forced_dip")
    all_forced_IDs = forced_hap_IDs.union(forced_dip_IDs)

    # reduce gt matrix & sample IDs to this subset
    GTs = GTs[:, all_sample_idx, :]
    all_sampleIDs = list(numpy.array(vcf_sample_pool)[all_sample_idx])
    # check for ploidy
    # hapSampleIdx = numpy.where( GTs[:,:,1] == -1, axis=0 )[0]
    hapSample_bools = numpy.all(GTs[:, :, 1] == -1, axis=0)  # <-- if True, then haploid
    hapSample_IDs = numpy.array(all_sampleIDs)[hapSample_bools]
    dipSample_IDs = numpy.array(all_sampleIDs)[~hapSample_bools]

    ## for each individual at each site, get the effective sample size,
    ## that is, the number of alleles that are not missing (either reference or alternative)
    n_alleleCalls = numpy.sum(GTs != -1, axis=2)  # <-- of shape (n_var, n_samp)
    assert (n_alleleCalls.min() >= 0) and (n_alleleCalls.max() <= 2), \
        f'unique values in n_alleleCalls={set(n_alleleCalls.flatten())}'
    # and the number of alternative alleles, that is, how many 1's
    n_alt = numpy.sum(GTs == 1, axis=2)
    assert (n_alt.min() >= 0) and (n_alt.max() <= 2), \
        f'unique values in n_alt={set(n_alt.flatten())}'
    assert numpy.all(n_alleleCalls >= n_alt)

    # find out which individuals have hets
    hasHet_bools = numpy.any(n_alt == 1, axis=0)
    # but the ones that are haploid to begin with don't count
    noHetDips_bools = ~hapSample_bools & ~hasHet_bools
    # print(len(noHetDips_bools), numpy.sum(noHetDips_bools))

    if noHetDips_bools.sum() > 0:
        noHetDips_IDs = numpy.array(all_sampleIDs)[noHetDips_bools]
        # if no instruction is given about the het-less samples, raise error
        if not set(noHetDips_IDs).issubset(all_forced_IDs):
            raise ValueError(f'Found diploid samples without heterozygotic sites:\n'
                             f'{set(noHetDips_IDs) - all_forced_IDs}\n'
                             f'\tTo continue after/without converting them to haploids, use '
                             f'\t`--force_hap/dip <comma-separated-IDs>` in the command.')
        else:
            noHetDips_asHaps = set(noHetDips_IDs).intersection(forced_hap_IDs)
            noHetDips_asDips = set(noHetDips_IDs).intersection(forced_dip_IDs)
            # get their indice
            noHetDips_asHaps_idx = [all_sampleIDs.index(ID) for ID in list(noHetDips_asHaps)]
            # noHetDips_asDips_idx = [all_sampleIDs.index(ID) for ID in list(noHetDips_asDips)] <-- no need to do anything to their counts
            print(f'Samples {noHetDips_asDips} will still be counted as diploids '
                  f'despite having no heterozygote GT calls.')
            # convert gts
            n_alleleCalls[:, noHetDips_asHaps_idx] = n_alleleCalls[:, noHetDips_asHaps_idx] // 2
            n_alt[:, noHetDips_asHaps_idx] = n_alt[:, noHetDips_asHaps_idx] // 2
            # maybe there are haps that user wants to count as dips?
            if len(set(forced_dips)) > len(noHetDips_asDips):
                asDips = set(forced_dips) - noHetDips_asDips
                hapsToDouble = asDips.intersection(set(hapSample_IDs))
                if len(hapsToDouble) > 0:
                    hapsToDouble_idx = [all_sampleIDs.index(ID) for ID in hapsToDouble]
                    # double
                    n_alleleCalls[:, hapsToDouble_idx] *= 2
                    n_alt[:, hapsToDouble_idx] *= 2
                elif len(asDips.intersection(set(all_sampleIDs))) == 0:
                    print(f'Sample(s) {asDips} are not in `all_sampleIDs`. '
                          f'In VCF? {asDips.issubset(set(vcf_sample_pool))}.')
            if len(set(forced_haps)) > len(noHetDips_asHaps):
                # maybe some will still be haps?
                asHaps = set(forced_haps) - noHetDips_asHaps
                # do not allow halving heterozygote GTs
                dipsToHalve = asHaps.intersection(set(dipSample_IDs))
                if len(dipsToHalve) > 0:
                    # do they all have hets?
                    dipsToHalve_idx = [all_sampleIDs.index(ID) for ID in dipsToHalve]
                    if numpy.any((GTs[:, dipsToHalve_idx, :])[:, :, 1] == 1, axis=2):
                        raise ValueError(f'Cannot force-halve diploid samples {dipsToHalve} '
                                         f'because they have heterozygotes')
                elif len(asHaps.intersection(set(all_sampleIDs))) == 0:
                    print(f'Sample(s) {asHaps} are not in `all_sampleIDs`. '
                          f'All in VCF? {asHaps.issubset(set(vcf_sample_pool))}.')

    # get sample IDs grouped by time/batch
    intsctSampleIDpools = [[vcf_sample_pool[idx] for idx in thisIndexPool]
                           for thisIndexPool in sampleIndexPools]
    # get their idx in the GT subset
    # newSampleIdxPools = [numpy.where(thisSamplePool in all_sampleIDs)[0] <-- doesn't work
    newSampleIdxPools = [[all_sampleIDs.index(sampID) for sampID in thisSamplePool]
                         for thisSamplePool in intsctSampleIDpools]
    # initiate
    samples = []
    sampleSizes = []
    kept_timePool_indice = []
    # in case some time points doesn't have any observations for certain samples
    for k, idx_pool in enumerate(newSampleIdxPools):
        if len(idx_pool) == 0:
            # let it pass
            print(f'No samples at sampling time {k} in vcf')
            continue
        else:
            # get counts
            this_pool_d = n_alt[:, idx_pool].sum(axis=1)
            this_pool_n = n_alleleCalls[:, idx_pool].sum(axis=1)
            samples.append(this_pool_d)
            sampleSizes.append(this_pool_n)
            # record time
            kept_timePool_indice.append(timePool_indice[k])
    # print(k+1, 'batches')
    # now transpose samples & sampleSizes
    samples = numpy.array(samples).transpose()
    sampleSizes = numpy.array(sampleSizes).transpose()

    if sampleTimes is not None:
        newSampleTimes = numpy.array([0] + [sampleTimes[k] for k in kept_timePool_indice])
    else:
        newSampleTimes = numpy.array([0] + kept_timePool_indice)

    # down-sample by pooled MAF  ## <-- *PHASING OUT* all minMAF filtering will be moved to `filter_snps`
    ## remove ones without observations first: (mostly to avoid zero division
    observed = (sampleSizes.sum(axis=1) > 0)
    samples = samples[observed, :]
    sampleSizes = sampleSizes[observed, :]
    locus_names = numpy.array(locus_names)[observed]
    # ## divide
    if minMAF > 0:
        print(f'All minMAF filters will be applied in `filter_snps()` function/step only.')
    # freqs = samples.sum(axis=1) / sampleSizes.sum(axis=1)
    # ## fold
    # freqs = numpy.where(freqs < 0.5, freqs, 1 - freqs)
    # # filter
    # passed_rows = (freqs > minMAF)
    # # down sample
    # samples = samples[passed_rows, :]
    # sampleSizes = sampleSizes[passed_rows, :]
    # locus_names = locus_names[passed_rows]
    # sanity check to make sure the dimensions match
    assert samples.shape == sampleSizes.shape, f'ValueError: samples.shape = {samples.shape} != sampleSizes.shape = {sampleSizes.shape}.'
    assert samples.shape[0] > 0, f'VCF file is empty. samples.shape = {samples.shape}.'
    # output
    return numpy.array(locus_names), numpy.array(samples), numpy.array(sampleSizes), newSampleTimes


def _format_to_dotless_e(number):
    # Format the number in "e" notation and convert it to a float
    formatted_number = "{:e}".format(float(number))
    # exponent = float(formatted_number.split("e")[1])
    e_part1, e_part2 = formatted_number.split("e")
    # remove zeros at the end:
    e_part2 = int(e_part2)
    if number >= 1:
        while not np.isclose(float(e_part1) - int(e_part1.split(".")[0]), 0):
            part1_update = float(e_part1) * 10
            e_part2 -= 1
            e_part1 = str(part1_update)
        # print(e_part1, e_part2)
        # only look at integer part
        while e_part1.endswith("0"):
            e_part1 = e_part1[:-1]
            e_part2 -= 1
            # print(e_part1, e_part2)
    else:
        # only look at decimal part
        while e_part1.startswith("0"):
            e_part1 = e_part1[1:]
            e_part2 += 1
            # print(e_part1, e_part2)
    # Convert the number to an integer to remove the decimal point and trailing zeros
    # formatted_integer = int(float(e_part1))
    return "{:d}e{:d}".format(int(float(e_part1)), int(e_part2) + 1)


def _parse_pos_arg(pos_arg: str):
    """Helper function to parse the string to indicate range of physical positions"""
    numbers = float_regex.findall(pos_arg)
    if len(numbers) == 2:
        left_pos, right_pos = sorted([float(p) for p in numbers])
        # pos_range = (left_pos, right_pos)
        pos_tag = f'_{_format_to_dotless_e(left_pos)}-{_format_to_dotless_e(right_pos)}'
    elif len(numbers) == 1:
        # need to tell which one
        assert "-" in pos_arg, pos_arg
        pos_range = pos_arg.split("-")
        if int_regex.search(pos_range[0]):
            left_pos = int(pos_range[0])
            right_pos = None
        elif int_regex.search(pos_range[1]):
            left_pos = None
            right_pos = int(pos_range[1])
        else:
            raise ValueError(f'cannot parse the position flag: {pos_arg}.')
            # sys.exit(1)
        # pos_range = (left_pos, right_pos)
    else:
        raise ValueError(f'cannot parse the position flag: {pos_arg}.')
        # sys.exit(1)
    return left_pos, right_pos, pos_tag


def _split_locus_name(locus_name):
    split_names = locus_name.split("_")
    if len(split_names) > 3:
        Chr, physPos = split_names[0], split_names[1]
        rsID = '_'.join(split_names[2:])
    elif len(split_names) == 3:
        Chr, physPos, rsID = split_names
    elif len(split_names) == 2:
        Chr, physPos = split_names
        rsID = locus_name
    else:
        return False
    return Chr, int(physPos), rsID


def filter_snps(Samples, K, n_max, minK, missingness: float, minMAF: float=0., chroms=None, pos_arg=None):
    """Remove loci whose data don't pass the filters. Input: allelic count matrix
    """
    sample_cols = [f'd{i}' for i in range(1, (K + 1))]
    size_cols = [f'n{j}' for j in range(1, (K + 1))]
    tag = ''

    ## parse chrom
    if chroms is not None:
        assert 'Chr' in Samples.columns, f'To specify chromosome, \"Chr\" column must exist in the input.' \
                                         f'Samples.columns = {Samples.clumns}'
        Chroms = chroms.split(',')
        Chroms = [str(c).strip() for c in Chroms if len(c) > 0]
        Samples = Samples[Samples.Chr.str.isin(Chroms)]
        print(f'Keeping {Samples.shape[0]} SNPs on chromosome(s) {Chroms}')
        tag += f'_chr{chroms}'
    else:
        Chroms = None
    ## then position
    if pos_arg is not None:
        if 'Chr' in Samples.columns:
            assert (len(Samples.Chr.values) == 1) or (Chroms is not None), \
                f'Chroms: {Chroms}, position range {pos_arg}'
        left_pos, right_pos, pos_tag = _parse_pos_arg(pos_arg)
        tag += pos_tag
        # now we filter
        before = Samples.shape[0]
        Samples = Samples[(Samples.physPos >= left_pos) & (Samples.physPos <= right_pos)]
        after = Samples.shape[0]
        print(f'Keeping {after} SNPs (from {before}) positioned between {int(left_pos)} to {int(right_pos)}')

    # another sanity check
    assert np.all(np.array(Samples[sample_cols]) <= np.array(Samples[size_cols])), \
        'allele counts must be no greater than respective samples sizes.'
    Samples.loc[:, 'pooled_d'] = Samples[sample_cols].sum(axis=1)
    Samples.loc[:, 'pooled_n'] = Samples[size_cols].sum(axis=1)

    ## now check out number of times observed
    if minK > 0:
        Samples['times_obs'] = (Samples[size_cols] > 0).sum(axis=1)
        before = Samples.shape[0]
        Samples = Samples[Samples.times_obs >= minK]
        after = Samples.shape[0]
        print(f'Remove {before - after} SNPs with <{minK} times/pools of observations. {after} SNPs left.')
        tag += f'_minK{minK:d}'

    # filter for missingness <--- float in (0, 1)
    before = Samples.shape[0]
    # print(Samples.pooled_n.dtype)
    Samples = Samples[(Samples.pooled_n > missingness * sum(n_max))]
    tag += f'_minObs{missingness:g}'
    after = Samples.shape[0]
    print(f'{after} SNPs (out of {before}) passed missingness filter (:= pooled_n > sum(n_max)*{missingness:g} ).')

    ## lastly: MAF
    assert 0 <= minMAF < 0.5, f'Invalid minor allele frequency cutoff: {minMAF}'
    tag += f'_MAF{str(minMAF)[1:]}'
    # get maf
    ## remove empty rows first
    before = Samples.shape[0]
    Samples = Samples[Samples.pooled_n > 0]
    after = Samples.shape[0]
    print(f'Remove {before - after} sites with zero samples')
    ## fold
    Samples.loc[:, 'pooled_mad'] = Samples.pooled_d.where(
        Samples.pooled_d.copy() <= 0.5 * Samples.pooled_n.copy(),
        (Samples.pooled_n.copy() - Samples.pooled_d.copy()))
    ## actual MAF filter
    before = Samples.shape[0]
    Samples['pooledMAF'] = Samples.pooled_mad / Samples.pooled_n  # .replace(0, np.nan) <-- there shouldn't be n=0 anymore
    Samples = Samples[Samples.pooledMAF > minMAF]
    after = Samples.shape[0]
    print(f' Remove {before - after} SNPs whose pooled MAF <= {minMAF}. {after} left.')
    # remove aux column?
    Samples = Samples.drop(columns=['pooled_d', 'pooled_mad', 'pooled_n'])

    # d columns with 0 counts throughout will not be removed to stay consistent with sample_times

    # done going through all filter args, return df
    ## exclude the "_" at the beginning of the tag
    return Samples, tag[1:]


def _reformat_LL_DF_to_matrix(onGrid_LLs):
    """Convert the DataFrame"""
    # assume the DF has both row names and column names
    loci_names = list(onGrid_LLs.ID)
    s_pairs = [colname for colname in onGrid_LLs.columns if colname != "ID"]
    # make sure they're numbers
    s_pairs = [list(map(float, pair.strip("()").split(','))) for pair in s_pairs]
    s1_list, s2_list = map(lambda x: numpy.array(sorted(list(set(x)))),
                           zip(*s_pairs))
    LLmatrix = numpy.array(onGrid_LLs.loc[:, onGrid_LLs.columns != "ID"], dtype=float)  # first column is ID
    # make sure s_pairs have tupples
    s_pairs = [tuple(pair) for pair in s_pairs]
    return numpy.array(LLmatrix), numpy.array(loci_names), s1_list, s2_list, s_pairs, len(s_pairs)


def _reformat_longLLtable_to_matrix(onGrid_LLs):
    onGrid_LLs['s_pair'] = pd.Series(zip(onGrid_LLs.s1, onGrid_LLs.s2))
    LLmatrix = onGrid_LLs.pivot(index="ID", columns="s_pair", values='loglikelihood')
    s_pairs = list(LLmatrix.columns)
    s1_list = sorted(list(set(onGrid_LLs.s1)))
    s2_list = sorted(list(set(onGrid_LLs.s2)))
    loci_names = numpy.array(LLmatrix.index)
    # sanity check
    assert LLmatrix.shape == (len(loci_names),
                              len(s_pairs)), f'LLmatrix.shape={LLmatrix.shape}, len(loci_names)={len(loci_names)}, len(s_pairs)={len(s_pairs)}'
    return numpy.array(LLmatrix), loci_names, numpy.array(s1_list), numpy.array(s2_list), s_pairs, len(s_pairs)


# only use when fix h or s
# start LLs_persite with row index
def _find_1D_Max_perSite(s_grid, ith_LLs_persite, bounds):
    idx, LLs_persite = ith_LLs_persite[0], ith_LLs_persite[1:]
    # remove inf
    inf_idx = numpy.where(numpy.isinf(LLs_persite))[0]
    # in-place update
    if len(inf_idx) > 0:
        # update LLs_persite last!
        s_grid = numpy.array(s_grid)[~numpy.isinf(LLs_persite)]
        LLs_persite = numpy.array(LLs_persite)[~numpy.isinf(LLs_persite)]
        bounds = (s_grid.min(), s_grid.max())
    # get on-grid max:
    ongrid_max_idx = numpy.argmax(LLs_persite)
    s_0 = s_grid[ongrid_max_idx]
    ongrid_max = LLs_persite[ongrid_max_idx]
    # interpolate
    LL = scipy.interpolate.interp1d(x=s_grid, y=-1 * LLs_persite, kind='cubic')
    # anyways, let's get some explicit 1d bounded brent optimizer going
    theOpt = scipy.optimize.minimize_scalar(LL, method='bounded', bounds=bounds)
    if theOpt.fun <= 0:
        # raise ValueError('Positive log-likelihood values during optimization')
        print(f'## POSITIVE LL!\toffGrid_max={-1 * theOpt.fun},')
        print(f'## ongrid_max_idx={ongrid_max_idx}; s_0={s_0}; ongrid max={LLs_persite[ongrid_max_idx]}')
        print(
            f'## LL(x={theOpt.x})={LL(x=theOpt.x)}, maxLL={-1 * theOpt.fun}; reported as x=-9, <row idx>={idx} in output\n')
        print('Positive log-likelihood values during optimization; use on-grid max instead')
        # return [s_0, ongrid_max, numpy.nan, numpy.nan]
        return [s_0, ongrid_max, -9, 0]
    # check for boundary
    elif numpy.any(numpy.isclose(bounds, theOpt.x)):
        print(f'## BOUNDARY LL!\t, <row idx>={idx}\toffGrid_max={-1 * theOpt.fun}, bounds={bounds}')
        print(f'## ongrid_max_idx={ongrid_max_idx}; s_0={s_0}; ongrid max={LLs_persite[ongrid_max_idx]}')
        print(f'## LL(x={theOpt.x})={LL(x=theOpt.x)}, maxLL={-1 * theOpt.fun}; reported in MLR output as x=-8 \n')
        return [s_0, ongrid_max, -8, -1 * theOpt.fun]
    else:
        return [s_0, ongrid_max, theOpt.x, -1 * theOpt.fun]


# pretty much the first part of `_find_1D_Max_perSite`
def _get_onGrid_Max_perSite(s_grid, ith_LLs_persite):
    idx, LLs_persite = ith_LLs_persite[0], ith_LLs_persite[1:]
    # get on-grid max:
    ongrid_max_idx = numpy.argmax(LLs_persite)
    s_0 = s_grid[ongrid_max_idx]
    ongrid_max = LLs_persite[ongrid_max_idx]
    return [s_0, ongrid_max, numpy.nan, numpy.nan]


"""CURRENTLY IN BETA. DO NOT USE
def _find_2D_Max_perSite(s1s, s2s, LLs_persite, bounds):
    '''
    Interpolate 2D LL surface and locate the local maxLL from the interpolated function.
    '''
    # interpolate first
    ## default is linear; it'd report positive LLs with cubic or quintic
    LL = scipy.interpolate.interp2d(x=s1s, y=s2s, z=LLs_persite, bounds_error=True)  # 'cubic', kind='quintic'
    # LL = scipy.interpolate.griddata((s1s, s2s), LLs_persite, (s1s.ravel(), s2s.ravel()))
    # get on-grid max:
    ongrid_max_idx = numpy.argmax(LLs_persite)
    s1_0 = s1s[ongrid_max_idx]
    s2_0 = s2s[ongrid_max_idx]
    ongrid_max = LLs_persite[ongrid_max_idx]
    # then maximize
    fLL = lambda x: -1 * LL(x=x[0], y=x[1])
    # bounds = 0.8*numpy.array(bounds)
    offGrid_max = scipy.optimize.minimize(fLL, x0=(s1_0, s2_0), bounds=bounds, method='Nelder-Mead')  # 'Powell
    s1hat, s2hat = offGrid_max.x
    try:
        maxLL = -1 * offGrid_max.fun[0]
    except:
        maxLL = -1 * offGrid_max.fun
    if maxLL >= 0:
        print(f'ongrid_max_idx={ongrid_max_idx}; s1_0={s1_0}; s2_0={s2_0}')
        print(f'ongrid max={LLs_persite[ongrid_max_idx]}')
        print(f'offGrid_max={-1 * offGrid_max.fun}')
        print(f'LL(x={s1hat}, y={s2hat})={LL(x=s1hat, y=s2hat)}, maxLL={-1 * offGrid_max.fun}\n')
        # sys.exit()
        # raise ValueError('Positive log-likelihood values during optimization')
        print('Positive log-likelihood values during optimization; use on-grid max instead')
        # s1hat, s2hat, maxLL = s1_0, s2_0, ongrid_max
        s1hat, s2hat, maxLL = numpy.nan, numpy.nan, numpy.nan
    return [s1_0, s2_0, ongrid_max, s1hat, s2hat, maxLL]
"""


def get_onGrid_max_only(s_pairs, LLmatrix):
    """Return the on-grid maxLL and its MLEs for the given numpy matrix of log-likelihoods.
    :param s_pairs: list of (s1, s2) tupples, matching the column names of `LLmatrix`
    :param LLmatrix: L x P numpy matrix of computed log-likelihoods. L and P are the
            numbers of loci (rows, axis 0) and (s1,s2)-pairs (cols, axis 1)
    :return: numpy.array([s1_hat, s2_hat, maxLL])
    """
    max_idx = numpy.argmax(LLmatrix, axis=1)
    s1_hat, s2_hat = zip(*[s_pairs[i] for i in max_idx])
    maxLLs = numpy.empty((LLmatrix.shape[0], 3))
    maxLLs[:, 0] = s1_hat
    maxLLs[:, 1] = s2_hat
    maxLLs[:, 2] = numpy.max(LLmatrix, axis=1)
    return maxLLs


def interpolate_offGrid_max(s1_list, s2_list, s_pairs, num_pairs, LLmatrix):
    """Interpolate off-grid local maxLL from *1D* surface
    Number of computed s_pairs must be no fewer than 4.
    :param s1_list: Ascending list/array of s1 grid points
    :param s2_list: Ascending list/array of s2 grid points
    :param s_pairs: list of (s1, s2) tupples, matching the column names of `LLmatrix`
    :param num_pairs: Total number of s pairs; len(s_pairs)
    :param LLmatrix: L x P numpy matrix of computed log-likelihoods. L and P are the
            numbers of loci (rows) and (s1,s2)-pairs (cols)
    :return: numpy.array(maxLL), a list of column names as the output header
    """
    # decide how many variables to optimize & define the function
    s1_gridnum = len(s1_list)
    s2_gridnum = len(s2_list)
    # check to see
    if min(s1_gridnum, s2_gridnum) > 3 and (s1_gridnum * s2_gridnum) == num_pairs:
        # 2d-grid
        print('WARNING: Do not interpolate off-grid max for 2D parameter grid. Only on-grid values will be reported.')
        header = '\t'.join(
            ['ID', 'ongrid_s1hat', 'ongrid_s2hat', 'ongrid_maxLogLikelihood'])

        # maxLLs is of dimension 3/
        maxLLs = get_onGrid_max_only(s_pairs, LLmatrix)

        # s1s, s2s = map(lambda x: numpy.array(list(x)), zip(*s_pairs))
        # bnds = ((min(s1s), max(s1s)), (min(s2s), max(s2s)))
        # header = '\t'.join(
        #     ['locus', 'ongrid_s1hat', 'ongrid_s2hat', 'ongrid_maxLogLikelihood', 's1hat', 's2hat', 'maxLogLikelihood'])
        #
        # def _get_max_persite(x):
        #     return _find_2D_Max_perSite(s1s, s2s, LLs_persite=x, bounds=bnds)
    # decide that it's 1d
    elif num_pairs in {s1_gridnum, s2_gridnum}:
        # fix h
        if s1_gridnum == num_pairs and s2_gridnum == num_pairs:
            var_list = list(zip(*s_pairs))[1]
            header = '\t'.join(['ID', 'ongrid_s2hat', 'ongrid_maxLogLikelihood', 's2hat', 'maxLogLikelihood'])
        # fix s1
        elif s1_gridnum == 1:
            var_list = list(zip(*s_pairs))[1]
            header = '\t'.join(['ID', 'ongrid_s2hat', 'ongrid_maxLogLikelihood', 's2hat', 'maxLogLikelihood'])
        # fix s2
        elif s2_gridnum == 1:
            var_list = list(zip(*s_pairs))[0]
            header = '\t'.join(['ID', 'ongrid_s1hat', 'ongrid_maxLogLikelihood', 's1hat', 'maxLogLikelihood'])
        else:
            print(f's1_gridnum={s1_gridnum}; s2_gridnum={s2_gridnum}; num_pairs={num_pairs}')
            print(f's1_list={s1_list}\ns2_list={s2_list}\ns_pairs={s_pairs}.')
            sys.exit()

        bnds = (min(var_list), max(var_list))

        def _get_max_persite(x):
            return _find_1D_Max_perSite(numpy.array(var_list), ith_LLs_persite=x, bounds=bnds)

        # add row index as a col
        iLLmatrix = numpy.insert(LLmatrix, 0, numpy.arange(1, LLmatrix.shape[0] + 1), axis=1)
        maxLLs = numpy.apply_along_axis(_get_max_persite, 1, iLLmatrix)
    else:
        print(f'''Insufficient grid points for effective interpolation. 
        At least 4 different values are needed for either `s1` or `s2`.
        s1_list={s1_list};\ns2_list={s2_list}\ns_pairs={s_pairs}.''')
        sys.exit()

    return maxLLs, header
