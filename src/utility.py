import numpy
import pandas as pd
import re, os, gzip, sys
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


# infile header needs to be fixed
def parse_allele_counts(infile: str, snp_list=None, minMAF: float=0.):
    """Parse the allele count file. Must have header.
    Return `numpy.array(locus_names)`, `numpy.array(samples)`, `numpy.array(sampleSizes)`
    """
    parsed = pd.read_csv(infile, sep="\t", header=0, comment="#")
    header = parsed.columns
    pos_regex = re.compile(r'^(locus|ch|pos|position|id)$')
    # sanity check
    # retrieve snp names + down-sample if needed
    id_cols = [colname for colname in header if pos_regex.search(colname.lower())]
    # locus_names = numpy.array(parsed.locus)
    if 'ID' in id_cols:
        if snp_list is not None:
            parsed = parsed[parsed.ID.isin(snp_list)]
        locus_names = numpy.array(parsed.ID)
    elif 'locus' in id_cols:
        if snp_list is not None:
            parsed = parsed[parsed.locus.isin(snp_list)]
        locus_names = numpy.array(parsed.locus)
    else:
        if snp_list is None:
            # no need to pick one, actually. just concatenate them
            locus_names = numpy.array(parsed.loc[:, id_cols])
            locus_names = numpy.apply_along_axis(lambda l: '_'.join(l), axis=1, arr=locus_names)
        else:
            id_colname = input('Please specify the ID column name used for `--snps` (Case-sensitive):')
            # down sample
            parsed = parsed[parsed[id_colname].isin(snp_list)]
            locus_names = numpy.array(parsed[id_colname])
    # just some sanity check
    if snp_list is not None:
        assert parsed.shape[1] > 0, f'No intersection between variants in allele count file and ' \
                                    f'those listed by `--snps`. Please double check your input.'
        sys.exit()
    # retrieve all xk column names
    d_regex = re.compile(r'^[xd](\d+)$')
    n_regex = re.compile(r'^n(\d+)$')
    d_cols = [colname for colname in header if d_regex.search(colname)]
    n_cols = [colname for colname in header if n_regex.search(colname)]
    # sanity check
    assert len(d_cols) == len(n_cols), f'd_cols={d_cols}; n_cols={n_cols}'
    # sanity check: the indices are monotonously ordered
    try:
        d_col_num = [int(d_regex.findall(colname)[0]) for colname in header if d_regex.search(colname)]
        n_col_num = [int(n_regex.findall(colname)[0]) for colname in header if n_regex.search(colname)]
    except Exception as e:
        print(e)
        print(f'Cannot parse column names: d_cols={d_cols}; n_cols={n_cols}')
        sys.exit(1)
    assert numpy.all(numpy.array(d_col_num) == numpy.array(n_col_num)), f'Columns for samples and sample sizes must be ordered the same.'
    assert numpy.all(numpy.diff(d_col_num) >= 0) or numpy.all(numpy.diff(d_col_num) <= 0), f'Please make sure the sample pools are ordered.'
    # assemble
    samples = numpy.array(parsed.loc[:, d_cols])
    sampleSizes = numpy.array(parsed.loc[:, n_cols])

    # down-sample by pooled MAF
    ## remove ones without observations first: (mostly to avoid zero division
    observed = (sampleSizes.sum(axis=1) > 0)
    samples = samples[observed,:]
    sampleSizes = sampleSizes[observed,:]
    locus_names = locus_names[observed]
    ## divide
    freqs = samples.sum(axis=1) / sampleSizes.sum(axis=1)
    ## fold
    freqs = numpy.where(freqs < 0.5, freqs, 1 - freqs)
    # filter
    passed_rows = (freqs > minMAF)
    # down sample
    samples = samples[passed_rows,:]
    sampleSizes = sampleSizes[passed_rows,:]
    locus_names = locus_names[passed_rows]

    return locus_names, samples, sampleSizes


def parse_ind_file(inds: str, type: str = "inds"):
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
        print(f'{len(ind_list)} {type} in total.')
        return ind_list
    else:
        print('Please provide a valid file name for the list of sample IDs.')
        return False


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
                if time_colname.endswith('_Ago'): # <-- remove the tag previously added
                    time_colname = time_colname[:-4]
            else: # the OG name could be without "ago" too
                assert len(find_ago) == 1
                if time_colname not in infoTable.columns:
                    time_colname = time_colname[:-4]
                assert time_colname in infoTable.columns, ValueError(f'Column \"{time_colname}\" missing in the info table.')
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
        else: # if both are in the table
            print('Both\"Time_Ago\" and \"Gen_Ago\" exist in table. Only reading numbers in \"Gen\" column.', infoTable.columns)
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
        else: # if both are in the table
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


def _get_GT_index(line):
    # get GT index
    GT_index = line[8].split(":").index("GT")
    return GT_index


def _get_ids_from_vcf_header(header, inds_list):
    # header = header.strip().split('\t') <-- already prepped before passing here
    # 0CHROM	1POS	2ID	3REF	4ALT	5QUAL	6FILTER	7INFO	8FORMAT	9:<individuals>
    # get ID list
    ## samples already stratified by time; all ids in the list are in the same pool. use set() to remove dups
    try:
        ID_list = set(map(header.index, inds_list))
    except ValueError: # some ids are not in the header
        # there's a certain individual not in list
        ID_list = [header.index(ind) for ind in inds_list if ind in header]
    # just in case
    if len(ID_list) != len(inds_list):
        assert len(ID_list) < len(inds_list), f'len(ID_list)={len(ID_list)}, len(inds_list)={len(inds_list)}.'
        print(f'{len(ID_list)} out of {len(inds_list)} in the inds_list exist/are included in the VCF file.')
    if len(ID_list) > 0:
        # sanity check (bc it's a vcf file
        assert min(ID_list) >= 9, f'ID_list={ID_list}, inds_list={inds_list}'

    return ID_list


# only consider bi-allelic locus here:
allele_regex = re.compile(r'[01]')
missing_regex = re.compile(r'(\.)')
GTseats_regex = re.compile(r'[01\.]')


def _GT_to_counts(vcfline, indexPools, het_flags: dict, force_hap_idx: list or set, force_dip_idx: list or set):
    """Turn a row in VCF file into counts of alleles among samples of the given col. indices"""
    GT_id = _get_GT_index(vcfline)

    # for each sample/col
    def _parseGT(idx):
        GT = vcfline[idx].split(":")[GT_id]
        alleles = allele_regex.findall(GT)
        assert len(set(alleles)) in {0, 1, 2}, f"Unconventional GT format: {GT}, alleles: {alleles}."
        # count all the "."
        missing = missing_regex.findall(GT)
        allSeats = GTseats_regex.findall(GT)
        # sanity check (make sure we counted all
        assert len(allSeats) == len(missing) + len(alleles), f"Unconventional GT format: {GT}, alleles: {alleles}."

        # convert GTs listed in `force_hap_idx` to haps
        if idx in force_hap_idx:
            if het_flags[idx]:
                print(f"Cannot force convert genotype {GT} to haplotype."
                    f" Variant{vcfline[2]}, column {idx + 1} GT: {GT}, len(alleles)={len(alleles)}.")
                het_flags[idx] = True
                # sys.exit(1)
                # ignore this GT
                xi, ni = 0, 0
            else:
                if len(alleles) == 0:
                    xi, ni = 0, 0
                else:
                    assert len(set(alleles)) == 1, f" Variant{vcfline[2]}, column {idx + 1} GT: {GT}, len(alleles)={len(alleles)}."
                    allele = int(alleles[0])
                    xi, ni = allele, 1
        # double the HTs listed in `force_hap_idx` to dips (but report dips as-is)
        elif idx in force_dip_idx:
            # if both alleles have calls
            if len(alleles) == 2:
                # just checking for hets?
                het_flags[idx] |= (alleles[0] != alleles[1])
                alleles = list(map(int, alleles))
                xi, ni = sum(alleles), 2
            # if it shows up as ./., [01]/., or [01]
            ## this is true haplotype [01]; double the counts
            elif len(alleles) == 1 and len(GT) == 1:
                allele = int(alleles[0])
                assert allele in {0, 1}, f'allele = {allele}; GT[idx={idx}]={GT}'
                xi, ni = 2*allele, 2
            ## this is [01]/.
            elif len(alleles) == 1 and len(GT) > 1:
                xi, ni = int(alleles[0]), 1
                try:
                    assert "." in GT, f'Unconventional genotype {GT}; ' \
                                      f'will not be converted as double-counted haplotype.'
                except Exception as e:
                    print(e)
                    print(f'will count sample[{idx}]\'s GT \"{GT}\" as xi={xi}, ni={ni}')
            # here len(alleles) must be 0, either ./. or .
            else:
                assert len(alleles) == 0
                xi, ni = 0, 0
        else:
            # check for hets
            check_het = (len(alleles) > 1 and alleles[0] != alleles[1])
            het_flags[idx] |= check_het
            # adding stuff up
            alleles = list(map(int, alleles))
            xi, ni = sum(alleles), len(alleles)

        return xi, ni

    # for each pool
    def _poolGTs(idx_pool):
        # n, x are tuples of alleles for each sample
        x, n = zip(*map(_parseGT, idx_pool))
        return sum(x), sum(n)

    # henceforth:
    samples, sampleSizes = zip(*map(_poolGTs, indexPools))
    samples = numpy.array(list(samples))
    sampleSizes = numpy.array(list(sampleSizes))

    return samples, sampleSizes, het_flags


def parse_vcf_input(vcffile: str, samplePools: list, sampleTimes=None, snps=None, inds=None,
                    forced_haps: list or set=[], forced_dips: list or set=[], minMAF: float=0.):
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
    :param sampleTimes: list of time points for each batch of samples
    :param forced_haps: list of sample IDs, as shown in VCF columns, to be counted as haploids
    :param forced_dips: list of sample IDs, as shown in VCF columns, to be double-counted as diploids if a GT is presented as haploid.
    :param minMAF: float, minimum threshold (non-inclusive) for minor allele frequencies in the pooled samples.
    -----------------------------
    :return: numpy.array(locus_names), numpy.array(samples), numpy.array(sampleSizes), numpy.array(sampleTimes), sampleIDpools
    """
    # setting up
    twohash_regex = re.compile(r'^##')
    onehash_regex = re.compile(r'^#[^#](.*)')
    SNP = {'A', 'T', 'C', 'G'}

    # reading file
    if vcffile.endswith('.vcf.gz'):
        vcf = gzip.open(vcffile, 'rt')
    else:
        assert vcffile.endswith('vcf')
        vcf = open(vcffile, 'r')
    # l = vcf.readline()
    locus_names = []
    samples = []
    sampleSizes = []
    sampleIndexPools = []
    # while l != "":
    with vcf:
        for l in vcf:
            if twohash_regex.search(l):
                continue
            # retrieve headerline
            elif onehash_regex.search(l):
                header = l.strip().split("\t")
                # samplePools is a list of lists, each includes IDs of the same time point
                # if data only a subset of all data, some time points may have empty index lists.
                if inds is None:
                    sampleIndexPools = [_get_ids_from_vcf_header(header, inds_pool) for inds_pool in samplePools]
                else:
                    smaller_samplePools = [[ID for ID in pool if ID in inds] for pool in samplePools]
                    sampleIndexPools = [_get_ids_from_vcf_header(header, inds_pool) for inds_pool in smaller_samplePools]
                # list of indice (in the vcf file) for samples to be forced haps
                force_hap_idx = _get_ids_from_vcf_header(header, forced_haps)
                force_dip_idx = _get_ids_from_vcf_header(header, forced_dips)
                ## indicator for existence of hets
                het_flags = {idx: False for pools in sampleIndexPools for idx in pools}
                # record all IDs, just in case
                sample_IDs = {idx: header[idx] for idx in het_flags.keys()}
                # in case some time points doesn't have any observations for certain samples
                if [] in sampleIndexPools:
                    assert sampleTimes is not None, 'VCF doesn\'t have all listed individuals, need `sampleTimes` arg too'
                    # remove empty lists from OG sampleTimes
                    trimmedTimePoints = [sampleTimes[k+1] for k, id_pool in enumerate(sampleIndexPools) if len(id_pool) > 0]
                    trimmedIndexPools = [id_pool for id_pool in sampleIndexPools if len(id_pool) > 0]
                    # update
                    sampleTimes = [0] + trimmedTimePoints
                    sampleIndexPools = trimmedIndexPools
                    samplePools = [[header[idx] for idx in pool] for pool in sampleIndexPools]
            else:
                l = l.strip().split("\t")
                assert len(l) > 9, ValueError('Please check the vcf file format.')
                assert sampleIndexPools != [], f'samplePools = {samplePools},\nheader = {header}'
                # pass if not in the snp list
                if snps is not None:
                    assert isinstance(snps, list), f'type(snps)={type(snps)}'
                    if l[2] not in snps:
                        continue
                # do we condition on SNPs? <-- yes.
                # but at least it should be bi-allelic; can check if alt has ","?
                if l[3] in SNP and l[4] in SNP:
                    locus_names.append('_'.join(l[:3]))
                    # het_flags is a dict storing whether a sample has het GT calls
                    Xs, Ns, het_flags = _GT_to_counts(l, sampleIndexPools, het_flags, force_hap_idx, force_dip_idx)
                    samples.append(Xs)
                    sampleSizes.append(Ns)
        # sanity check after reading the whole VCF
        if sum(het_flags.values()) > 0:
            # assert these individuals are not forced to be haps
            het_idx = {idx for idx in het_flags.keys() if het_flags[idx]}
            assert het_idx.isdisjoint(set(force_hap_idx)), f'Sample(s) {", ".join([sample_IDs[idx] for idx in het_idx])} have heterozygote GT call(s) despite being commanded to be haploid(s).'
        # alert user if no hets observed
        else:
            option = input("No heterozygote genotypes found in this VCF. Is this to be expected? [y|n] ")
            if "y" in option.lower():
                print('Program will continue...')
            elif "n" in option.lower():
                sys.exit()
            else:
                print("Invalid response.")
                sys.exit()
        # l = vcf.readline()
    vcf.close()
    # make sure to output the right format
    samples, sampleSizes = numpy.array(samples), numpy.array(sampleSizes)
    # sanity check to make sure the dimensions match (ie `sampleTimes` be updated if needed)
    assert samples.shape == sampleSizes.shape, f'ValueError: samples.shape = {samples.shape} != sampleSizes.shape = {sampleSizes.shape}.'
    assert samples.shape[0] > 0, f'VCF file is empty. samples.shape = {samples.shape}.'

    # down-sample by pooled MAF
    ## remove ones without observations first: (mostly to avoid zero division
    observed = (sampleSizes.sum(axis=1) > 0)
    samples = samples[observed,:]
    sampleSizes = sampleSizes[observed,:]
    locus_names = numpy.array(locus_names)[observed]
    ## divide
    freqs = samples.sum(axis=1) / sampleSizes.sum(axis=1)
    ## fold
    freqs = numpy.where(freqs < 0.5, freqs, 1 - freqs)
    # filter
    passed_rows = (freqs > minMAF)
    # down sample
    samples = samples[passed_rows,:]
    sampleSizes = sampleSizes[passed_rows,:]
    locus_names = locus_names[passed_rows]

    # output
    return numpy.array(locus_names), numpy.array(samples), numpy.array(sampleSizes), sampleTimes, het_flags, sampleIndexPools


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
            ['locus', 'ongrid_s1hat', 'ongrid_s2hat', 'ongrid_maxLogLikelihood'])

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
            header = '\t'.join(['locus', 'ongrid_s2hat', 'ongrid_maxLogLikelihood', 's2hat', 'maxLogLikelihood'])
        # fix s1
        elif s1_gridnum == 1:
            var_list = list(zip(*s_pairs))[1]
            header = '\t'.join(['locus', 'ongrid_s2hat', 'ongrid_maxLogLikelihood', 's2hat', 'maxLogLikelihood'])
        # fix s2
        elif s2_gridnum == 1:
            var_list = list(zip(*s_pairs))[0]
            header = '\t'.join(['locus', 'ongrid_s1hat', 'ongrid_maxLogLikelihood', 's1hat', 'maxLogLikelihood'])
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


