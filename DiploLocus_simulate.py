#!/usr/bin/env python3
"""Command-line interface for simulating independent replicates of temporal samples with given parameters."""
import sys, time, gzip
import numpy as np

# local import
try:
    from diplo_locus import simulate
except ModuleNotFoundError:
    from src import simulate


def load_custom_distn(filename):
    distn = np.loadtxt(filename, delimiter="\t")
    initValues = distn[:, 0]
    initProbs = distn[:, 1]
    return (initValues, initProbs)


def write_output(samples, trajs, sampleSizes, notes, outPrefix, write_traj, gzip_output=""):
    """Write out samples and, if choose to, trajectories, of simulated replicates"""
    samples_outfile = f'{outPrefix}_samples.count'
    if "s" in gzip_output:
        samples_outfile += ".gz"
        samples_handle = gzip.open(samples_outfile, "wt")
    else:
        samples_handle = open(samples_outfile, "w")
    print(f'Saving sample counts to {samples_outfile}...')
    # annotate:
    samples_handle.write(notes)
    # samples_handle.write(f'##Sample.Times.gen: {", ".join(list(map(str, sampleTimes[1:])))}\n')

    num_reps, num_samps = samples.shape
    # write header:
    header = "ID\t" + "\t".join([f'd{k}\tn{k}' for k in range(1, num_samps + 1)]) + "\n"
    samples_handle.write(header)
    # sanity check
    assert num_samps == len(sampleSizes)
    for row_i in range(num_reps):
        counts = "\t".join(["%g\t%g" % (samples[row_i, k], sampleSizes[k]) for k in range(num_samps)])
        line = f'rep{row_i}\t{counts}\n'
        samples_handle.write(line)
    # end writing. Close file.
    samples_handle.close()
    # see if traj needs to be written
    if write_traj:
        traj_outfile = f'{outPrefix}.traj'
        if "t" in gzip_output:
            traj_outfile += ".gz"
        # no need to add rep id
        print(f'Saving trajectories to {traj_outfile}...')
        np.savetxt(traj_outfile, trajs, delimiter=", ")


import matplotlib.pyplot as plt


# import pandas as pd
# import seaborn as sns
def plot_reps(samples, trajectories, sampleTimes, notes, figname, reps=range(10), sampleSizes=None, deltaT=None,
              include_traj=False, include_samples=False):
    """Plot the sample frequency (solid) and/or population frequency (dashed) of (up to 10) simulated replicates."""
    # built-in args & dimension check
    assert include_traj or include_samples
    if include_traj:
        assert deltaT is not None
        assert (trajectories.shape[1] - 2) * deltaT <= max(
            sampleTimes), f'traj.shape={trajectories.shape}, deltaT={deltaT},' \
                          f'sample_times={sampleTimes}.'
    if include_samples:
        assert sampleSizes is not None
        assert sampleTimes is not None
        # sampleSizes must be array-like object that has a dimension of (K,)
        # only support uniform sample sizes at the moment
        assert len(sampleSizes) == samples.shape[1]
        # sampleTimes will not have the zero at the beginning
        assert len(sampleSizes) == len(sampleTimes)

    # prune reps
    if len(reps) > 10:
        print('Only the first 10 replicates will be considered.')
        reps = reps[:10]
    elif len(reps) > samples.shape[0]:
        reps = np.arange(samples.shape[0])

    # extract data of reps
    if include_samples:
        samp_freqs = samples[reps, :] / np.array(sampleSizes)
    if include_traj:
        traj_to_plot = trajectories[reps, :]
        x_values = deltaT * np.arange(trajectories.shape[1])

    # assign color to reps
    rep_colors = plt.cm.tab10(range(len(reps)))
    assert len(rep_colors) == len(reps), f'rep_colors.shape={len(rep_colors)}, len(reps)={len(reps)}\n{rep_colors}'

    # now we plot
    fig, ax = plt.subplot_mosaic([[0], [1]], gridspec_kw={'height_ratios': [3, 1]})
    if include_traj and include_samples:
        # assert traj_to_plot is not None and x_values is not None
        plot_obj = "samples and trajectories"
        for i in range(len(reps)):
            ax[0].plot(x_values, traj_to_plot[i, :], color=rep_colors[i], label=reps[i], ls='-')
            # ax[0].plot(sampleTimes, samp_freqs[i,:], color=rep_colors[i], ls=':', alpha=0.3)
            ax[0].plot(sampleTimes, samp_freqs[i, :], color=rep_colors[i], ls='-.', marker='x', alpha=0.5)
    elif include_samples:
        plot_obj = "samples"
        for i in range(len(reps)):
            # ax[0].plot(sampleTimes, samp_freqs[i,:], color=rep_colors[i], ls='-')#, label=reps[i]
            ax[0].plot(sampleTimes, samp_freqs[i, :], color=rep_colors[i], label=reps[i], ls='-', marker='x')
    elif include_traj:
        plot_obj = "trajectories"
        for i in range(len(reps)):
            ax[0].plot(x_values, traj_to_plot[i, :], color=rep_colors[i], label=reps[i], ls='-')
    else:
        return False

    if include_samples:
        ax[0].legend(ncol=int(round(len(reps) / 2)), fontsize='x-small', loc='best',
                     framealpha=0.2, title='Replicate', title_fontsize='small')
    elif include_traj:
        ax[0].legend(ncol=int(round(len(reps) / 2)), fontsize='x-small', loc='best',
                     framealpha=0.2, title='Replicate', title_fontsize='small')

    # aesthetic stuff
    ax[0].set_ylim(0, 1)
    ax[0].set_ylabel('Allele Frequency')
    ax[0].set_xlabel('Generation')
    # write out info
    notes = notes.replace("## ", "")
    ax[1].text(0, 0, notes, ha='left', va='center')
    ax[1].set_axis_off()

    plt.rcParams['axes.unicode_minus'] = False
    print(f'Plotting {plot_obj} to {figname}...')
    plt.tight_layout()
    plt.savefig(figname, dpi=300)


import argparse


def sim_args_parser(parser):
    # some pop gen parameters
    param = parser.add_argument_group('Population Parameters')
    param.add_argument('--u01', dest='u01', default=None, required=True, type=float,
                       help="Mutation rate allele 0 to allele 1 (/site/generation).")
    param.add_argument('--u10', dest='u10', default=None, type=float,  # required=True,
                       help='Mutation rate allele 1 to allele 0 (/site/generation). Default is u01.')
    param.add_argument('--Ne', dest='Ne', default=None, required=True, type=float,
                       help='Effective diploid population size (i.e., total #(allele)=2Ne ).')

    sel = parser.add_argument_group('Selection Parameters')
    sel.add_argument('--s1', dest='S1', required=True,
                     help='Value(s) of selection coefficient for heterozygotes. If more than one value is provided \
                        (separated by comma, no space), then user must also provide selection change time with \"--selection_change_times\".')
    sel.add_argument('--s2', dest='S2', required=True,
                     help='Value(s) of selection coefficient for homozygotes. If more than one value is provided \
                        (separated by comma, no space), then user must also provide selection change time with \"--selection_change_times\".')
    sel.add_argument('--selection_change_times', dest='selectionChangeTimes', default=None,
                     # required=("," in sys.argv[sys.argv.index('--s')+1]),
                     help='Times (in forward generations) when either selection coefficient changes. Must be 1 fewer than the number of s values provided.')

    samp = parser.add_argument_group('Sample Parameters')
    samp.add_argument('--sample_times', dest="sampleTimes", required=True,
                      help='Numbers of generation times to take samples, separated by comma (no sapce). Must match the length of sample sizes.')
    samp.add_argument('--sample_sizes', dest='sampleSizes', required=True,
                      help='Numbers of alleles to sample at each time points specified by \"--sample_times\". '
                           'Must match the length of sample times.')
    samp.add_argument('--num_rep', dest='numRep', required=True, default=1, type=int,
                      help='Number of replicates to generate.')
    samp.add_argument('--last_popfreq_nonzero', dest='last_popfreq_nonzero', default=False, action='store_true',
                      help='Option to only count replicates with non-zero population frequency in the end.')
    samp.add_argument('--last_sample_nonzero', dest='last_sample_nonzero', default=False, action='store_true',
                      help='Option to only count replicates where the last (i.e. most recent) sample is not zero.')
    samp.add_argument('--not_lost', dest='not_lost', default=False, action='store_true',
                      help='Option to only count replicates whose selected alleles are not lost throughout the simulation.')
    samp.add_argument('--not_fixed', dest='not_fixed', default=False, action='store_true',
                      help='Option to only count replicates whose selected alleles are not fixed throughout the simulation.')
    samp.add_argument('--minMAF', dest='minMAF', default=0., type=float,
                      help='Minimal minor allele frequency (not included) in the meta-sample pooled from all sampling time. Default value is 0, where replicates with all empty samples or all full samples are excluded.')

    init = parser.add_argument_group('Initial Condition')
    init.add_argument('--init', dest='initCond', required=True, choices=['initFreq', 'statBeta', 'uniform', 'custom'],
                      help='Specify initial probability distribution to generate the replicates')
    init.add_argument('--initFreq', dest='initFreq', type=float,
                      required=('initFreq' in sys.argv), default=None,
                      help='Required when --init=\"initFreq\". Specify the allele frequency when selection started.\nIf \
                       multiple values are provided (separated by comma), then length must match the number of replicates.')
    init.add_argument('--initDistn', dest='initDistnFile', required=('custom' in sys.argv),
                      help='Path and name to the file specifying a particular frequency distribution as the initial condition.')
    init.add_argument('--init_u01', dest='init_u01', type=float, default=None,
                      help='Specify mutation rate allele 0 to allele 1 for the initial distribution. Default to be identical to u01.')
    init.add_argument('--init_u10', dest='init_u10', type=float, default=None,
                      help='Specify mutation rate allele 1 to allele 0 for the initial distribution. Default to be identical to u10.')

    out = parser.add_argument_group('Output Options')
    out.add_argument('-o', '--out_prefix', required=True, dest='outPrefix', help='Path and prefix of the output file.')
    out.add_argument('--write_traj', dest='write_traj', default=False, action='store_true',
                     help='Option to also write out the trajectories of replicates.')
    out.add_argument('--gzip_output', dest='gzip_output', choices=['none', 's', 't', 'st', 'ts'], default='t',
                     help='Option to write output in .gz format. Use \"s\" to indicate sample file, \"t\" to indicate trajectory file, and \"st\" (or \"ts\") for both. Use \"none\" to not zip either.')
    out.add_argument('--plot_trajs', dest='plot_trajs', default=None, action='store_true',
                     help='Option to plot the population allele frequency trajectories of up to 10 replicates.')
    out.add_argument('--plot_samples', dest='plot_samples', default=None, action='store_true',
                     help='Option to plot the sample allele frequency trajectories of up to 10 replicates.')
    out.add_argument('--reps_to_plot', dest='reps_to_plot', default=None,
                     help='Specify the replicates to plot. Default is the first 10 reps.')

    diffParam = parser.add_argument_group('Other Diffusion Parameters')
    diffParam.add_argument('--deltaT', dest='deltaT', type=float, default=0.05,
                           help='Time increment used for simulating under the diffusion (in generations). Default is 0.05.')
    diffParam.add_argument('--seed', dest='sd', type=int, default=None, help='Specify random seed.')

    return parser


def main(DL_args=None):
    # initiate parser
    parser = argparse.ArgumentParser(
        prog='DiploLocus_simulate.py',
        description='Light-weight tool to simulate independent replicates of temporal samples under given parameters.'
        )

    sim_parser = sim_args_parser(parser)

    # parse args
    if DL_args is not None:
        # receive args passed from DL.py
        args = DL_args
    else:
        # args = parser.parse_args(sys.argv[1:])
        args = sim_parser.parse_args(sys.argv[1:])

    # return usage if args blank
    if len(sys.argv[2:]) == 0:
        sim_parser.print_usage()
        sys.exit()

    # now let's go through CL arg groups one by one
    # pop parameters:
    if args.u10 is None:
        args.u10 = args.u01
    # record it for later use
    notes = f'## Parameters: Ne={args.Ne}, u01={args.u01}, u10={args.u10};\n'

    # sample parameters
    sampleTimes = list(map(float, args.sampleTimes.split(",")))
    # make sure it's ascending
    # this condition works when len(sampleTimes)==1 too
    assert np.all(np.array(sampleTimes) >= 0) and all(
        np.diff(np.array(sampleTimes)) > 0), 'Sampling times must be in ascending order.'
    if ',' in args.sampleSizes:
        sampleSizes = list(map(float, args.sampleSizes.split(",")))
        # length must match
        assert len(sampleSizes) == len(
            sampleTimes), 'When more than one sample sizes are provided, their number should match the number of sampling time points.'
    else:  # when only provide a single value
        size = float(args.sampleSizes)
        assert size > 0
        sampleSizes = [size] * len(sampleTimes)
    # take notes
    notes += f'## Sampling times: {", ".join(list(map(lambda x: "%g" % x, sampleTimes)))};\n'
    notes += f'## Sample sizes:  {", ".join(list(map(lambda x: "%g" % x, sampleSizes)))};\n'
    # add 0 in front of sample times
    sampleTimes = [0.] + list(sampleTimes)

    # initial conditions
    condInitSeg = True
    if args.initCond == 'initFreq':
        initCond = args.initCond
        assert args.initFreq is not None
        initValues, initProbs = None, None
        try:
            initFreq = float(args.initFreq)
        except Exception as e:
            print(e)
            assert "," in args.initFreq
            assert args.initFreq.count(",") == (
                    args.numRep - 1), 'Length of initial frequency provided should either be one or match the number of replicates.'
            initFreq = np.array(list(map(float, args.initFreq.split(","))))

        if type(initFreq) is float:
            assert 0 <= initFreq <= 1
            # take notes
            notes += f'## Initial condition: fixed frequency at {initFreq};\n'
        else:
            assert isinstance(initFreq, (list, np.ndarray))
            assert all(initFreq <= 0) and all(initFreq >= 1)
            # take notes
            notes += f'## Initial condition: fixed freq. specified for each replicates;\n'
    elif args.initCond == "statBeta":
        initCond = 'contBeta'
        initValues, initProbs = None, None
        initFreq = None
        condInitSeg = False
        if args.init_u01 is None:
            args.init_u01 = args.u01
        if args.init_u10 is None:
            args.init_u10 = args.init_u01
        notes += f'## Initial condition: stationary Beta distribution with alpha={args.init_u10}, beta={args.init_u10};\n'
    elif args.initCond == 'uniform':
        initCond = 'contBeta'
        initValues, initProbs = None, None
        initFreq = None
        condInitSeg = False
        args.init_u01 = 1 / (4 * args.Ne)
        args.init_u10 = 1 / (4 * args.Ne)
    elif args.initCond == "custom":
        initCond = 'choice'
        assert args.initDistnFile is not None, 'Please provide the file name of the initial condition probability distribution.'
        initValues, initProbs = load_custom_distn(args.initDistnFile)
        initFreq = None
        condInitSeg = False  # this won't be used anyways (for now)
        # take notes
        notes += f'## Initial condition: loaded from {args.initDistnFile};\n'
    else:
        print(f'\"--init\" does not take \"{args.initCond}\"')
        parser.print_usage()
        sys.exit()

    # get seed, if any
    if args.sd is not None:
        seed = args.sd
    else:
        seed = np.random.get_state()[1][0]
    notes = f'## Seed {seed};\n' + notes
    print(time.ctime(), f'Simulating {args.numRep} independent replicates with seed {seed}...\n')
    rng = np.random.RandomState(seed)
    # np.random.seed(seed)

    # selection parameters:
    if "," in args.S1 or "," in args.S2:
        # must provide sel change time
        assert args.selectionChangeTimes is not None, 'Must provide selection change times when more than one value is provided for either selection coefficient.'
        # length must match
        if args.S1.count(",") != args.S2.count(","):
            assert int("," in args.S1) * int(
                "," in args.S2) == 0, 'Values for selection coefficients must have the same length or have either s1 or s2 be constant.'
            # make them the same length
            if "," not in args.S1:
                s2s = list(map(float, args.S2.split(",")))
                num_pieces = args.S2.count(",")
                assert args.selectionChangeTimes.count(",") == (
                        num_pieces - 1), 'Length of selection change times must match the number of s values provided.'
                s1s = [float(args.S1)] * num_pieces
            # elif "," not in args.S2:
            else:
                s1s = list(map(float, args.S1.split(",")))
                num_pieces = args.S1.count(",")
                s2s = [float(args.S2)] * num_pieces
        else:
            s2s = list(map(float, args.S2.split(",")))
            s1s = list(map(float, args.S1.split(",")))
            assert len(s1s) == len(
                s2s), 'Values for selection coefficients must have the same length or have either s1 or s2 be constant.'
            num_pieces = len(s1s)
        # parse change times
        selChangeTimes = list(map(float, args.selectionChangeTimes.split(",")))
        # sanity check
        assert len(selChangeTimes) == (
                num_pieces - 1), 'Length of selection change times must match the number of s values provided.'
        # make sure it's ascending
        assert all(np.array(selChangeTimes) > 0) and all(
            np.diff(np.array(selChangeTimes)) > 0), 'Selection change times must be in ascending order.'

        # make wrapper
        def _batch_simulate(numRep):
            seed_here = rng.randint(1e12, size=1)[0]
            tempHandle = simulate.simulateSamples(args.Ne, s1s, s2s, mAlpha=args.u01, mBeta=args.u10,
                                                  times=sampleTimes, sampleSizes=sampleSizes, seed=seed_here,
                                                  initCond=initCond, initFreq=initFreq,
                                                  numReplicates=numRep, deltaT=args.deltaT,
                                                  condInitSeg=condInitSeg, initProbs=initProbs, initValues=initValues,
                                                  initMAlpha=args.init_u01, initMBeta=args.init_u10,
                                                  selectionChangeTimes=selChangeTimes)
            return tempHandle

        # take notes
        s_pair_strings = [f'({s1}, {s2})' for s1, s2 in zip(s1s, s2s)]
        time_string = ", ".join(list(map(lambda x: "%g" % x, (selChangeTimes + [sampleTimes[-1]]))))
        notes += f'## Selection coefficients on [0, {time_string}] time intervals:\n##     (s1, s2)[t] = {", ".join(s_pair_strings)}.\n'
    else:
        # both have single value
        s1, s2 = float(args.S1), float(args.S2)

        # now run stuff
        def _batch_simulate(numRep):
            seed_here = rng.randint(1e12, size=1)[0]
            tempHandle = simulate.simulateSamples(args.Ne, s1, s2, mAlpha=args.u01, mBeta=args.u10,
                                                  times=sampleTimes, sampleSizes=sampleSizes, seed=seed_here,
                                                  initCond=initCond, initFreq=initFreq,
                                                  numReplicates=numRep, deltaT=args.deltaT,
                                                  condInitSeg=condInitSeg,
                                                  initProbs=initProbs, initValues=initValues,
                                                  initMAlpha=args.init_u01, initMBeta=args.init_u10)
            return tempHandle

        # take notes
        notes += f'## Selection coefficients: ({s1}, {s2}), constant.\n'

    def _iter_mulate(_batch_sim, numReps, not_lost=False, not_fixed=False,
                     last_sample_nonzero=False, last_popfreq_nonzero=False):
        r = 0
        samples = np.empty((numReps, len(sampleSizes)))
        trajectories = np.empty((numReps, 2 + int((sampleTimes[-1] - sampleTimes[0]) / args.deltaT)))
        numThrownOut = 0
        while r < numReps:
            batch_samples, batch_trajs = _batch_sim(numReps - r)
            # toKeep = [True] * (numReps - r)
            toKeep = np.ones(numReps - r).astype(bool)
            if last_popfreq_nonzero:
                toKeep &= (batch_trajs[:, -1] != 0)
            if last_sample_nonzero:
                toKeep &= (batch_samples[:, -1] != 0)
            if not_lost:
                toKeep &= np.all((batch_trajs[:, 2:] > 0.5 / (2 * args.Ne)), axis=1)
            if not_fixed:
                toKeep &= np.all(batch_trajs[:, 1:] < 1. - 0.5 / (2 * args.Ne), axis=1)
            # check for MAF
            ## sample freqs
            pooledMAF = np.sum(batch_samples, axis=1) / np.sum(sampleSizes)
            ## fold
            pooledMAF = np.where(pooledMAF > 0.5, 1. - pooledMAF, pooledMAF)
            ## filter
            toKeep &= (pooledMAF > args.minMAF)

            numToKeep = sum(toKeep)
            assert numToKeep <= numReps - r
            # apply filter
            passed_samples = batch_samples[toKeep, :]
            passed_trajs = batch_trajs[toKeep, :]
            # stack on
            samples[r:(r + numToKeep), :] = passed_samples
            trajectories[r:(r + numToKeep), :] = passed_trajs
            numThrownOut += (numReps - r - numToKeep)
            # update
            r += numToKeep
        # end of while-loop
        return samples, trajectories, numThrownOut

    # now let's actually simulate
    samples, trajectories, numThrownOut = _iter_mulate(_batch_simulate, args.numRep, args.not_lost, args.not_fixed,
                                                       args.last_sample_nonzero, args.last_popfreq_nonzero)
    print(time.ctime(), f'Simulation finished. Threw away {numThrownOut} replicates.\n')
    # write output
    outPrefix = f'{args.outPrefix}_seed{seed}'
    ## note the filters
    if args.last_sample_nonzero or args.last_popfreq_nonzero or args.not_lost or args.not_fixed:
        notes += '## Only count replicates '
        filters = []
        if args.last_popfreq_nonzero:
            filters.append('with non-zero population allele freq. in the end')
            outPrefix += '_existNow'
        if args.last_sample_nonzero:
            filters.append('with non-zero count in the last (most recent) sample')
            outPrefix += '_nonzeroLastSamples'
        if args.not_lost:
            filters.append('whose population allele freq. never drop below 1/(4Ne),\n##        i.e. not ever lost')
            outPrefix += '_notLost'
        if args.not_fixed:
            filters.append('whose population allele freq. never exceeds 1 - 1/(4Ne),\n##        i.e. not ever fixed')
            outPrefix += '_notFixed'
        notes += ', and\n##    '.join(filters)
    notes += '.\n'
    print(notes)
    ## write tables
    write_output(samples, trajectories, sampleSizes, notes, outPrefix, args.write_traj, args.gzip_output)

    # see if we need to plot
    if args.plot_trajs is not None or args.plot_samples is not None:
        figname = outPrefix
        if args.plot_trajs:
            figname += '_popFreq'
        if args.plot_samples:
            figname += '_sampleFreq'
        figname += '.png'
        # parse the reps to plot
        if args.reps_to_plot is not None:
            reps_to_plot = list(map(int, args.reps_to_plot.split(',')))
        else:
            reps_to_plot = range(min([10, args.numRep]))
        plot_reps(samples, trajectories, sampleTimes[1:], notes, figname, reps=reps_to_plot,
                  sampleSizes=sampleSizes, deltaT=args.deltaT, include_traj=args.plot_trajs,
                  include_samples=args.plot_samples)

    print(time.ctime(), "Done.")


if __name__ == '__main__':
    main()
