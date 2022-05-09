# `diplo_locus`:

Python modules for simulating and computing log-likelihoods for time series genetic data based on the diploid Wright-Fisher diffusion

This is the beta version. For testing only.

# Table of Contents

* [`likelihood` module](#likelihood)
  * [class `SelHmm`](#likelihood.SelHmm)
    * [Constructors](#SelHmm.constructors) 
    * [Functions](#SelHmm.functions)
* [`simulate` module](#simulate)


<a id="likelihood"></a>

## `likelihood` module

Python module for HMM class to compute likelihood of temporal data under the Wright-Fisher diffusion and general diploid selection.

<a id="likelihood.SelHmm"></a>

### *class* `SelHmm`

Object to compute likelihood of observing temporal sampled allele frequency data under the given parameters. Likelihood is computed assuming that underlying population allele frequency evolves according to Wright-Fisher diffusion.

<a id="SelHmm.constructors"></a>

### Constructors


#### `SelHmm (Ne, s1, s2, mAlpha, mBeta, initCond, initFreq=None, initMAlpha=None, initMBeta=None, initS1=None, initS2=None, sampleSizesSet=None, numStates=1001, deltaT=1, emissionType="integer", transitionType="constant", selectionChangeTimes=None, selectionGridObject=None)`


<details>
<summary>Click to see parameter details</summary>

- **`Ne`** `float`

   Effective population size (assuming diploid individuals).

- **`s1`** `float or array_like`

   Selection coefficient(s) of the heterozygote. One value if selection is constant. When simulating piecewise time-varying selection, the number of values should be `selectionChangeTimes + 1`

- **`s2`** `float or array_like`

   Selection coefficient(s) of the homozygote. One value if selection is constant. When simulating piecewise time-varying selection, the number of values should be `selectionChangeTimes + 1`

- **`mAlpha`** `float`

   Per-generation mutation probability towards the focal allele.

- **`mBeta`** `float`

   Per-generation mutation probability towards the non-focal allele.

- **`initCond`** `{"uniform"` `"initFreq"` `"statBeta"}`

   Specify the initial condition for the HMM at generation zero, which is before the first sampling time.
   - "uniform" :
         Distribution of inital frequency in HMM is distribution on [0,1].
   - "initFreq" :
         Specfiy a fixed frequency as the initial condition for the HMM. Paramter `initFreq` specifies the initial frequency(ies). Either one frequency for the analysis of all datasets or one frequency for each dataset (has to be compatible with the number of datasets given).
   - "statBeta" :
         Uses the stationary Beta distribution as the initial distribution of the HMM. The mutation rates to specify this distributoon are `mAlpha` and `mBeta`, ot they can be speficied separately using `initMAlpha` and `initMBeta`.

- **`initFreq`** `float or array_like` `optional`

   Required when ``initCond="initFreq"``. Specifies intial frequencies for HMM. Either one frequency for all datasets, or a list of frequencies, one for each dataset. Must be between 0 and 1.

- **`initAlpha`** `float` `default=`mAlpha``

   Per-generation mutation probability towards the focal allele used for the initial Beta distribution, if specified by ``initCond="statBeta"``.

- **`initBeta`** `float` `default=`mBeta``

   Per-generation mutation probability towards the non-focal allele used for the initial Beta distribution, if specified by ``initCond="statBeta"``.

- **`initS1`** `float` `optional`

   Experimental, do not change this value.

- **`initS2`** `float` `optional`

   Experimental, do not change this value.

- **`sampleSizesSet`** `set of int`

   Required if ``emissionType="integer"``. Set of sample sizes that will be in the given data. Used to initialize the emission probabilities.

- **`numStates`** `int` `default=1001`

   Number of states to discretize the allele frequency space [0,1]. Currently implemented as equidistant.

- **`deltaT`** `float` `default=1`

   Experimental. Unit increment of time (in generations). Values unequal to 1 have not been tested.

- **`emissionType`** `{"integer"` `"fractional"}` `default="integer"`

   Type of emission data, that is, sample sizes and number of focal alleles in data.
   - "integer" :
         Data is assumed to have integer samples sizes and allele counts.
   - "fractional" :
         Data can have (positve) fractional samples sizes and allele counts. This is useful in situations where allele count is estimated and can assume non-integer values.

- **`transitionType`** `{"constant"` `"piecewise"}` `default="constant"`

   Type of transitions to consider in the HMM:
   - "constant" :
         Selection coefficients stay constant throughout the entire duration considered. `s1` and `s2` must be floats (not `array_like`) for this option.
   - "piecewise" :
         This has to be chosen if selection coefficients change over time, in a piece-wise manner. Must also specify `selectionChangeTimes` for this option, and `s1` and `s2` have to specify the correct number of coefficients.

- **`selectionChangeTimes`** `array_like` `optional`

   The generation times when selection coefficients change. Must match the length of `s1` and `s2`.

- **`selectionGridObject`** `default=None`

   Experimental, do not change this value.

</details>

<a id="SelHmm.functions"></a>

### Functions

#### `computeLogLikelihood (times, samplesSizes, samples)`


After specifying the underlying parameters to construct a `SelHmm` object, this function compute the log-likelihoods of observing the given samples at the sampling times `times`. The locus is assumed to be bi-allelic. If multiple datasets are given, then they are assumed independent.

<details>
<summary>Click to see parameter details</summary>

- **`times`** `array_like`

   The generation times when the samples were taken. Must start with zero and ascend from past to present. The first value specfies the time of the iniital condition, and the second value the time of the first sample. This enables specifying a given condition before the first sample is taken. If the first sample is taken at the tome of the initial condition, the first two values can be specified both as 0.

- **`samplesSizes`** `array_like`

   An N by K matrix of the total numbers of the sample sizes for each dataset and each time point. N --> number of datsets ;  K --> number of sampling times. ``sampleSizes[i,j]`` records the sample size of dataset ``i`` at time point ``j``.

- **`samples`** `array_like`

   An N by K matrix of the numbers of focal alleles observed in the respective dataset at the respective time. N --> number of datasets ;  K --> number of sampling times. ``samples[i,j]`` records the number of observed focal alleles in dataset ``i`` at time point ``j``.

Returns
-------
   numpy.array of log-likelihoods for each dataset.

</details>

<a id="simulate"></a>

## `simulate` module

Python module with functions to perform simulations under a Wright-Fisher diffusion and taking samples at given times.

#### `simulateSamples(Ne, allS1, allS2, mAlpha, mBeta, times, sampleSizes, initCond=None, initFreq=None, numReplicates=1, deltaT=0.05, condInitSeg=True, initGridResolution=None, initMAlpha=None, initMBeta=None, selectionChangeTimes=None)`

Simulate temporal samples from a population evolving according to the Wright-Fisher diffusion with general diploid selection.

<details>
<summary>Click to see parameter details</summary>

- **`Ne`** `float`

   Effective population size (assuming diploid individuals).

- **`allS1`** `float or array_like`

   Selection coefficient(s) of the heterozygote. One value if selection is constant. When simulating piecewise time-varying selection, the number of values should be `selectionChangeTimes + 1`

- **`allS2`** `float or array_like`

   Selection coefficient(s) of the homozygote. Requirements are the same as `allS1`.

- **`mAlpha`** `float`

   Per-generation mutation probability towards the focal allele.

- **`mBeta`** `float`

   per-generation mutation probability towards the non-focal allele.

- **`times`** `array_like`

   Generation times. Must start with zero to specify the time at which the initial condition is assumed. The remaining generations are the times when samples were taken. Thus ``len(times) = len(sampleSizes) + 1`` must hold.

- **`sampleSizes`** `array_like`

   The total number of sampled alleles at each sampling time point. Must satisfy ``len(times) = len(sampleSizes) + 1``.

- **`initCond`** `{"initFreq"` `"contBeta"` `"choice"}`

   Indicate the initial condition (at genration 0).
   - "initFreq" :
         Each replicate starts with a fixed given frequency. Either all replicates same frequency or each with a specific on. Paramter `initFreq` specifies the initial frequency(ies).
   - "contBeta" :
         Initial frequency for each relpicate will be drawn from a Beta distribution (stationary distribution under reccurrent mutation). Use `initAlpha` and `initBeta` the mutation rates used to determine the parameters of the Beta distribution. If not specified explicitly, the rates `mAlpha` and `mBeta` are used.
   - "choice" :
         Draw initial frequency from a distrubtion given by the parameter `initProbs`.

- **`initFreq`** `float or array_like` `optional`

   Required when ``initCond="initFreq"`` is set. Either one frequency for all replicates or one frequency per replicate. Values must be between 0 and 1.

- **`numReplicates`** `int` `default=1`

   Number of independent replicates to simulate.

- **`deltaT`** `float` `default=0.05`

   Specify the time increment for each step of the Normal approximation to the Wright-Fisher diffusion (in generations).

- **`condInitSeq`** `bool` `default=False`

   Experimental, do not change this value.

- **`initGridResolution`** `int` `optional`

   Experimental, do not change this value.

- **`initProbs`** `array_like` `optional`

   Required when initial condition is set to be "choice". This is an array of probabilities that is used to choose the initial frequencies for each replicate. `initValues` specifies the initial frequency that each probability corresponds to.

- **`initValues`** `numpy.ndarray` `optional`

   Required when initial condition is set to be "choice". This is an array of frequencies that is used to choose the initial frequencies for each replicate. `initProbs` specifies the probabilities that a certain frequency is chosen as initial frequency.

- **`initAlpha`** `float` `default=`mAlpha``

   Per-generation mutation probability towards the focal allele used to specify the initial Beta distribution.

- **`initBeta`** `float` `default=`mBeta``

   Per-generation mutation probability towards the non-focal allele used to specify the initial Beta distribution.

- **`selectionChangeTimes`** `array_like` `optional`

   Set the generation times when selection coefficients change. Only needs to be spedified if multiple selection coefficients are given. Must be one more than the length of `allS1` and `allS2`.

Returns
-------
   - **`samples`** `numpy.ndarray`

A `numReplicates` x len(`sampleSizes`) matrix of simulated samples. samples[i,j] records the number of derived alleles observed for replicate i at sampling time point j.

- **`wfDiffReplicates`** `numpy.ndarray`

   A `numReplicates` x steps matrix of simulated allele frequency trajectories, with number of steps = #(generations spanned) / `deltaT`.
 </details>


