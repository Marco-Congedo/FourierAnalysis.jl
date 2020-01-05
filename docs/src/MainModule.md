# MainModule

This is the main unit containing the **PosDefManifold** *module* (FourierAnalysis.jl).

## dependencies

| standard Julia packages |     external packages    |
|:-----------------------:|:-----------------------:|
| [LinearAlgebra](https://bit.ly/2W5Wq8W) | [FFTW](https://github.com/JuliaMath/FFTW.jl) |
| [Statistics](https://bit.ly/2Oem3li) | [AbstractFFTs](https://github.com/JuliaMath/AbstractFFTs.jl) |
|  |[RecipesBase](https://github.com/JuliaPlots/RecipesBase.jl) |
|  | [DSP](https://github.com/JuliaDSP/DSP.jl) |


The main module does not contains functions, but it declares all
**types** and objects holding data (named **data objects**) used in all units.

| Contents  |
|:----------:|
| [types](@ref) |
| [data objects](@ref) |
| [tips & tricks](@ref) |


## types

### IntOrReal
```
IntOrReal = Union{Int, Real}
```

This type is used internally and not exported.

### fInterval
```
fInterval = Union{IntOrReal,Tuple{IntOrReal, IntOrReal}, Colon}
```

In several functions you are allowed to select a **frequency range** (in Hz).
This can be
- a single frequency (in Hz), given as an integer or a real number, e.g., 10,
- a 2-tuple of integers or reals holding the lower and upper frequency limit (in Hz), e.g. (8, 10.5),
- a colon, to indicate as usual in Julia all available frequencies, e.g., :

!!! note "'frequencies' in different domains"
    Frequency ranges are used as argument in the functions [`mean`](@ref), [`extract`](@ref), [`meanAmplitude`](@ref), [`concentration`](@ref), [`meanDirection`](@ref), [`comodulation`](@ref) and [`coherence`](@ref). They apply both to frequency domain and to time-frequency domain data. In the frequency domain, the frequencies are the Fourier discrete frequencies with resolution ``sr/wl``, with or without the DC level in the first position, depending on how the object has been constructed. In the time-frequency domain, the frequencies actually are the center frequencies of the [filter bank](https://en.wikipedia.org/wiki/Filter_bank) used for constructing the object. Thus, The actual frequencies actually contained in a given position depends on the `bandwidth` argument used to construct the oject. See [`filterbank`](@ref).

### tInterval
```
tInterval = Union{Int, Tuple{Int, Int}, Colon}
```

In several functions you are allowed to select a **time range** (in samples).
This can be
- a single sample, given as an integer, e.g., 16,
- a 2-tuple of integers holding the lower and upper time limit (in sampes), e.g. (1, 120),
- a colon, to indicate as usual in Julia all available samples, e.g., :

These ranges are used as argument in the same functions where [fInterval](@ref)
ranges are used, however they apply only to time-frequency domain data.

### Smoother
```
@enum Smoother begin
    noSmoother       = 1
    hannSmoother     = 2
    hammingSmoother  = 3
    blackmanSmoother = 4
end
```

An instance of this type is requested by the [`smooth`](@ref) function,
and as an optional keyword argument by several constructors. It apply both to objects created in the frequency domain and in the time-frequency domain.

- `noSmoother` corresponds to no *smoothing*.
- `hannSmoother` is the *Hann* smoothing window (3-points)
- `hammingSmoother` is the *Hamming* smoothing window (3-points)
- `blackmanSmoother` is the *Blackman* smoothing window (5-points)

See [`smooth`](@ref) for details on these windows.

### SpectraVector

A vector of [Spectra](@ref) objects.

### CrossSpectraVector

A vector of [CrossSpectra](@ref) objects.

### CoherenceVector

A vector of [Coherence](@ref) objects.

**CoherenceVector₂**

A vector of [CoherenceVector](@ref) objects.

### TFAnalyticSignalVector

A vector of [TFAnalyticSignal](@ref) objects.

### TFAmplitudeVector

A vector of [TFAmplitude](@ref) objects.

### TFPhaseVector

A vector of [TFPhase](@ref) objects.

### FDobjects

An object in the frequency domain (FD), that is, the union of [Spectra](@ref),
[CrossSpectra](@ref) and [Coherence](@ref) objects.

### FDobjectsVector

A vector of objects in the frequency domain (FD), that is, the union of [SpectraVector](@ref),
[CrossSpectraVector](@ref) and [CoherenceVector](@ref) types.

### TFobjects

An object in the time-frequency (TF) domain, that is, the union of [TFAnalyticSignal](@ref), [TFAmplitude](@ref), [TFPhase](@ref) objects.

### TFobjectsVector

A vector of objects in the time-frequency (TF) domain, that is, the union of [TFAnalyticSignalVector](@ref), [TFAmplitudeVector](@ref) and [TFPhaseVector](@ref) types.

## data objects

*FourierAnalysis* creates an *operator object*, the [`Planner`](@ref), which is used to create FFTW plans for FFT computations, and several *data objects*. All of them are Julia [structures](https://docs.julialang.org/en/v1/base/base/#struct)). Data objects all have a corresponding vector type:

|      data objects        |     domain     |         vector type        |
|:------------------------:|:--------------:|:--------------------------:|
| [Taper](@ref)            |      time      | *none*                     |
| [Spectra](@ref)          |   frequency    | [SpectraVector](@ref)      |
| [CrossSpectra](@ref)     |   frequency    | [CrossSpectraVector](@ref) |
| [Coherence](@ref)        |   frequency    | [CoherenceVector](@ref)    |
| [TFAnalyticSignal](@ref) | time-frequency | [TFAnalyticSignalVector](@ref) |
| [TFAmplitude](@ref)      | time-frequency | [TFAmplitudeVector](@ref)  |
| [TFPhase](@ref)          | time-frequency | [TFPhaseVector](@ref)      |

For all data objects, *FourierAnalysis* overwrites the Julia [Base.show](https://docs.julialang.org/en/v1/base/io-network/#Base.show-Tuple{Any,Any,Any})
function to display in the REPL relevant information about the
name and value of the struct's fields in an easily-readable tabular form.

The field holding the data of all data objects is consistently named `.y`. As a summary, this is what `.y` holds in all data objects:

- [Taper](@ref): for all simple tapers, a real ``t``-vector holding the tapering window. For Slepians tapers (multi-tapering), a real ``t``x``h`` matrix where each column holds the tapering window for one of the ``h`` Slepian tapering windows.

- [Spectra](@ref): a real ``f``x``n`` matrix, where ``f`` is the number of Fourier discrete frequencies and ``n`` the number of time-series. Each column of `y` holds the spectrum for the corresponding time-series. In the case of one time-series only, `y` is a vector.

- [CrossSpectra](@ref): a ``f``-vector of complex ``n``x``n`` matrices  where ``f`` is the number of Fourier discrete frequencies and ``n>1`` the number of time-series that have generated the cross-spectra. Each matrix is the cross-spectral matrix for the corresponding frequency.

- [Coherence](@ref): a ``f``-vector of real ``n``x``n`` matrices  where ``f`` is the number of Fourier discrete frequencies and ``n>1`` the number of series that have generated the coherence. Each matrix is the coherence matrix for the corresponding frequency.

- [TFAnalyticSignal](@ref): a complex ``f``x``t`` matrix, where ``f`` is the number of band-pass regions of the filter bank used to generate the analytic signal and ``t`` is the number of time samples (time-frequency representation).

- [TFAmplitude](@ref): a real ``f``x``t`` matrix, where ``f`` is the number of band-pass regions of the filter bank used to generate the amplitude and ``t`` is the number of time samples (time-frequency representation).

- [TFPhase](@ref): a real ``f``x``t`` matrix, where ``f`` is the number of band-pass regins  of the filter bank used to generate the phase and ``t`` is the number of time samples (time-frequency representation).

## tips & tricks

See the [recipes.jl](@ref) for tips on how to plot tapering windows,
spectra and time-frequency data.

By convention, the *frequency* dimension of data arrays in
[data objects](@ref)
is always the first dimension. For instance,
the data of multivariate [Spectra](@ref) and of all time-frequency objects
([TFobjects](@ref)), that is, [TFAnalyticSignal](@ref),
[TFAmplitude](@ref) and [TFPhase](@ref), is a matrix in which the first
dimension unrolls frequencies.
This is the case also for [CrossSpectra](@ref) and [Coherence](@ref) objects,
which data is a vector of cross-spectral or coherence matrices across
frequencies.

#### window length in FFTW

For effective use of the FFTW package, keep in mind that
in FFTW the window length ``wl`` does not need to be a power of two,
however FFTW is best for window lengths of the form
``2^a, 3^b, 5^c, 7^d, 11^e,  13^f``, where ``e+f`` is either ``0`` or ``1``,  and the other exponents are arbitrary.

#### derive your own time-frequency measures

All functions implemented in unit [timefrequencyuni.jl](@ref) and
[timefrequencybi.jl](@ref) allow to compute time-frequency measures
averaging across several analytic signal objects. For each analytic signal
you can extract a time-frequency region passing with argument `mode` the
[`extract`](@ref) function (default) or the average in such a region passing
with argument `mode` the [`mean`](@ref) function.

Argument `func` allows you to apply any function to the extracted regions
or means. Using this you can create a hell of new averaging
procedures. As an example, suppose that you want to use the
[`meanAmplitude`](@ref) of a [TFAnalyticSignalVector](@ref) object,
but instead of the arithmetic mean of the amplitude you want to compute
the geometric mean of the power. You will then pass as argument to
the [`meanAmplitude`](@ref) function `func=x->log.(x.^2)` and take
the `exp.` function of the output. It is as simple as that and works
both using the `mode=extract` argument (default) or the `mode=mean`
argument.


#### Threads
Some functions in *FourierAnalysis* calls BLAS routine implicitly
via Julia. You can set the number of threads
the BLAS library should use by:

```
using LinearAlgebra
BLAS.set_num_threads(n)
```

where `n` is the number of threads.
By default, *FourierAnalysis* reserves to BLAS
all CPU threads available on your computer (given by the output of `Sys.CPU_THREADS`). The number threads used by Julia
for multi-threaded computations is given by the output of `Threads.nthreads()`.
In Windows this latter number of threads is set to half the available threads.
In Linux and OSX defaults to one and is controlled by an environment variable, i.e.,

```
export JULIA_NUM_THREADS=4
```

In Linux, working with the Atom IDE, you also have to
set to `global` the field found in Atom under
`Settings(or Preferences)/julia-client/Settings/Julia Options/Number of Threads`.

In Windows, set the desired number of threads in the settings
of the julia-client Juno package.

See this [post](https://discourse.julialang.org/t/issue-number-of-threads/14593), this [post](https://discourse.julialang.org/t/customize-number-of-threads-interactively/11574/2) and julia
[doc on threads](https://docs.julialang.org/en/v1/manual/parallel-computing/#Multi-Threading-(Experimental)-1).

Notice that *FourierAnalysis* features many multi-threaded functions and these
may allow a gain in computation time only if Julia is instructed to use
at least two threads.
