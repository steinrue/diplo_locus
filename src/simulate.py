# Python library for functions to perform simulations under a Wright-Fisher diffusion and taking samples at given times.
import numpy

# local imports
# from . import diffusion_core
# the python package management system is a pure delight
import pathlib

module_path = pathlib.Path(__file__).parent.resolve().__str__()
import sys

sys.path.insert(0, module_path)
import diffusion_core


# now do some WF model simulations
def simulateSamples(Ne, allS1, allS2, mAlpha, mBeta, times, sampleSizes, seed, initCond=None, initFreq=None,
                    numReplicates=1, deltaT=0.05, condInitSeg=True, initGridResolution=None, initProbs=None,
                    initValues=None, initMAlpha=None, initMBeta=None, selectionChangeTimes=None):
    r"""Simulate temporal samples from a population evolving according to the Wright-Fisher diffusion with general diploid selection.

    Parameters
    ----------
    Ne : int or float
        Effective population size (assuming diploid individuals).

    allS1 : int, float, numpy numbers, or array_like
        Selection coefficient of the heterozygote. One value if selection is constant. When simulating piecewise
        time-varying selection, the length of `list` or `numpy.ndarray` should match that of `selectionChangeTimes`
        so that ``len(allS1) = len(selectionChangeTimes) + 1``.

    allS2 : int, float, numpy numbers, list, or array_like
        Selection coefficient of the homozygote. Requirements are the same as `allS1`.

    mAlpha : float
        Per-site per-generation forward mutation rate.

    mBeta : float
        Per-site per-generation backward mutation rate.

    times : array_like
        Numbers in (forward) generations when samples were taken.
        Must start with zero to specify the time at which the initial condition is assumed. The remaining
        numbers are the times when samples were taken. Thus``len(times) = len(sampleSizes) + 1``.

    sampleSizes: array_like
        The total number of sampled alleles at each sampling time point.
        Must satisfy ``len(times) = len(sampleSizes) + 1``.

    seed: int
        Integer to set as a seed for the random number generator used for the simulation.

    initCond : {"initFreq", "contBeta", "choice"}
        Indicate the initial condition (at generation 0).
        - "initFreq" :
            Each replicate starts with a fixed given frequency. Either all replicates same frequency
            or each with a specific one. The parameter `initFreq` specifies the initial frequency(ies).
        - "contBeta" :
            Initial frequency for each replicate will be drawn from a Beta distribution (stationary distribution
            under recurrent mutation). Use `initAlpha` and `initBeta` the mutation rates used to determine the parameters
            of the Beta distribution. If not specified explicitly, the rates `mAlpha` and `mBeta` are used.
        - "choice" :
            Draw initial frequency from a distribution given by the parameters `initValues` and `initProbs`.

    numReplicates : int, optional, default=1
        Number of independent replicate loci to simulate.


    Other Parameters
    ----------------
    initFreq : float, optional
        Required when ``initCond="initFreq"``. Must be between 0 and 1.

    initAlpha : float, optional, default=`mAlpha`
        Per-site per-generation forward mutation rate underlying the initial Beta distribution.

    initBeta : float, optional, default=`mBeta`
        Per-site per-generation backward mutation rate underlying the initial Beta distribution.

    deltaT : float, optional, default=0.05
        Unit increment of time (in generations).

    condInitSeq : bool, default=False
        Experimental, do not change this value.

    initGridResolution : int, optional
        Experimental, do not change this value.

    initProbs : array_like, optional
        Required when initial condition is set to be "choice".
        This is an array of allele frequencies for which
        `initValue` specifies their corresponding probability densities.

    initValues : numpy.ndarray, optional
        Required when initial condition is set to be "choice".

    selectionChangeTimes : int or array_like, optional
        Set the generation time when selection coefficients change. Only needs to be specified if multiple
        selection coefficients are given. Must be one more than the length of `allS1` and `allS2`.

    Returns
    -------
    samples : numpy.ndarray
        A `numReplicates` x len(`sampleSizes`) matrix of simulated samples.
        samples[i,j] records the number of derived alleles observed
        on replicate i at sampling time point j.

    wfDiffReplicates : numpy.ndarray
        A `numReplicates` x steps matrix of simulated allele frequency trajectories,
        with steps = #(generations spanned) / `deltaT`.
    """
    # make sure all parameters are good
    assert (len(times) - 1 == len(sampleSizes))

    # see about times
    assert (all([x <= y for (x, y) in zip(times[:-1], times[1:])]))
    # let's fix the first one to be at zero
    assert (numpy.isclose(times[0], 0))
    # just the last one
    numGenerations = times[-1]
    # so how many steps do we have
    numSteps = int(numGenerations / deltaT) + 1

    # see about rng
    assert (numpy.issubdtype(type(seed), numpy.integer)), f'seed={seed}, type(seed)={type(seed)}'
    thisRNG = numpy.random.default_rng(seed)

    numSampleTimes = len(times) - 1
    # we also need the sampling indeces (omit first one)
    sampleTicks = [int(x / deltaT) for x in times[1:]]
    assert (all([(0 <= x) and (x < numSteps) for x in sampleTicks]))
    # and a map for convenience
    sampleIdx = {}
    for i in range(len(sampleTicks)):
        sampleIdx[sampleTicks[i]] = i
    # just some checks to make sure nothing went wrong
    idxList = list(sampleIdx.values())
    assert (numpy.min(idxList) == 0)
    assert (numpy.max(idxList) == numSampleTimes - 1)
    assert (len(set(idxList)) == len(idxList))

    # see about maybe multiple selection stuff
    if selectionChangeTimes is None:
        # only one coefficient
        assert (isinstance(allS1, (int, numpy.integer, float, numpy.floating))), type(allS1)
        assert (isinstance(allS2, (int, numpy.integer, float, numpy.floating))), type(allS2)
        initS1 = allS1
        initS2 = allS2

        # make a vector of only this coefficients
        s1s = numpy.ones(numSteps) * initS1
        s2s = numpy.ones(numSteps) * initS2
    else:
        # list of coefficients
        assert (type(allS1) in [list, numpy.ndarray])
        assert (type(allS2) in [list, numpy.ndarray])
        assert (len(selectionChangeTimes) == len(allS1) - 1)
        assert (len(selectionChangeTimes) == len(allS2) - 1)
        assert (all([x <= y for (x, y) in zip(selectionChangeTimes[:-1], selectionChangeTimes[1:])]))
        # should be more general now
        # assert (selectionChangeTimes[0] >= times[0])
        # assert (selectionChangeTimes[-1] <= times[-1])

        # build the list of coefficients
        s1s = numpy.ones(numSteps) * allS1[0]
        s2s = numpy.zeros(numSteps) * allS2[0]
        for tIdx in numpy.arange(len(selectionChangeTimes)):
            realTime = int(selectionChangeTimes[tIdx] / deltaT)
            realTime = numpy.clip(realTime, 0, len(s1s))
            s1s[realTime:] = allS1[tIdx + 1]
            s2s[realTime:] = allS2[tIdx + 1]
        initS1 = s1s[0]
        initS2 = s2s[0]
        # should be all good

    # see about mutation rates for initial condition
    if initCond in ["discBeta", "discBetaSel", "contBeta"]:
        if initMAlpha == None:
            initMAlpha = mAlpha
        if initMBeta == None:
            initMBeta = mBeta
        assert (initMAlpha > 0), "'initMAlpha' has to be greater than zero. Either specify or adjust mutation rate."
        assert (initMBeta > 0), "'initMBeta' has to be greater than zero. Either specify or adjust mutation rate."

    # initialize some storage
    wfDiffReplicates = numpy.zeros((numReplicates, numSteps + 1))
    samples = numpy.zeros((numReplicates, numSampleTimes)).astype(int)

    # initialize a vector,  cause we will do all replicates in parallel =).
    if initCond == "initFreq":

        assert (initFreq is not None)

        # do we have one frequency or many frequencies
        if (type(initFreq) == float) or (len(initFreq) == 1):
            # one 
            if type(initFreq) != float:
                initFreq = initFreq[0]
            assert (initFreq >= 0)
            assert (initFreq <= 1)

            # have given init frequency
            initFrequencies = numpy.ones(numReplicates) * initFreq
        else:
            # many
            assert (isinstance(initFreq, (list, numpy.ndarray)))
            assert (len(initFreq) > 1)
            initFreq = numpy.array(initFreq)
            assert (len(initFreq) == numReplicates), "Need as many initial frequencies as replicates."
            assert ((initFreq >= 0).all())
            assert ((initFreq <= 1).all())

            # already good
            initFrequencies = numpy.array(initFreq)

    elif initCond == "contBeta":

        assert (
            not condInitSeg), "Cannot condition on segregating with continuous initial distribution. Set 'condInitSeg = False'"
        # sample from beta distribution
        alpha = 4 * Ne * initMAlpha
        beta = 4 * Ne * initMBeta
        initFrequencies = thisRNG.beta(alpha, beta, size=numReplicates)

    # elif initCond in ["discBeta", "discBetaSel", "choice"]:
    elif initCond == "choice":

        # if initCond in ["discBeta", "discBetaSel"]:
        #     # here we need grids
        #     assert (initGridResolution is not None), "'initGridResolution' needed, but not given."
        #     yGrid = diffusion_core.getLinearGrid(initGridResolution)
        #     weights = diffusion_core.getWeights(diffusion_core.getLinearBounds(yGrid))
        #     # if we use this function, we now that lowest is 0 and highest is 1

        # TODO put these back in
        # if (initCond == "statNeutralPRF"):
        #     samplingProbs = Stationary.stationaryNeutralPRF (Ne, m, yGrid, weights)
        # elif (initCond == "statSelPRF"):
        #     samplingProbs = Stationary.stationarySelectionPRF (Ne, m, initS1, initS2, yGrid, weights)
        # if initCond == "discBeta":
        #     samplingProbs = diffusion_core.Stationary.stationaryNeutralBeta(Ne, initMAlpha, initMBeta, yGrid, weights)
        # elif initCond == "discBetaSel":
        #     samplingProbs = diffusion_core.Stationary.stationarySelectionBeta(Ne, initMAlpha, initMBeta, initS1, initS2,
        #                                                                       yGrid, weights)
        # elif initCond == "choice":
        #     assert (not initProbs is None)
        #     samplingProbs = initProbs
        # else:
        #     assert False
        assert (not initProbs is None)
        samplingProbs = numpy.array(initProbs)

        # # see if we need to modify
        # if initCond in ["discBeta", "discBetaSel"]:
        #     if condInitSeg:
        #         samplingProbs = samplingProbs[1:-1]
        #         samplingProbs /= numpy.sum(samplingProbs)
        #         daGrid = yGrid[1:-1]
        #     else:
        #         daGrid = yGrid
        # elif initCond == "choice":
        #     assert (not initValues is None)
        #     daGrid = initValues
        # else:
        #     assert False
        assert (not initValues is None)
        daGrid = numpy.array(initValues)

        assert (len(daGrid) == len(samplingProbs))

        idxs = thisRNG.choice(numpy.arange(len(samplingProbs)), size=numReplicates, replace=True, p=samplingProbs)
        initFrequencies = daGrid[idxs]
    else:
        assert False, f"Unknown initial condition: {initCond}"

    # so I think that should give us the right initial conditions

    # go through time
    for t in numpy.arange(0, wfDiffReplicates.shape[1]):

        # are we at the beginning
        if t == 0:
            wfDiffReplicates[:, t] = initFrequencies
        else:
            # get the parameters for the gaussian
            # mu
            daMu = diffusion_core.mu(wfDiffReplicates[:, t - 1], s1s[t - 1], s2s[t - 1], mAlpha, mBeta)
            # some safety check for now that hopefully catches stuff we can't do yet
            # should be fine, but check anyways
            assert (numpy.logical_not(numpy.isin(daMu, [float("-inf"), float("inf")])).all())
            # sigma square
            daSigmaSq = diffusion_core.sigmaSq(wfDiffReplicates[:, t - 1], Ne)

            # get loc and scale
            (daLoc, daScale) = diffusion_core.diffusionLocScale(wfDiffReplicates[:, t - 1], daMu, daSigmaSq, deltaT)

            # and do the gaussian steps
            wfDiffReplicates[:, t] = thisRNG.normal(daLoc, daScale)
            wfDiffReplicates[:, t] = numpy.clip(wfDiffReplicates[:, t], 0, 1)

        # sampling if we have to
        if t in sampleTicks:
            # get the idx
            thisIdx = sampleIdx[t]

            # and then do binomial sampling
            samples[:, thisIdx] = thisRNG.binomial(sampleSizes[thisIdx], wfDiffReplicates[:, t])

    # if (t > 20):
    #     return (None, None)

    return samples, wfDiffReplicates
