# timefrequency.jl

This unit implements time-frequency representations based
on analytic signal estimations. It uses the [filters.jl](@ref)
unit for filtering the signals in a
[filter bank](https://en.wikipedia.org/wiki/Filter_bank)
and the [hilbert.jl](@ref) unit for computing the
[Hilbert transform](https://en.wikipedia.org/wiki/Hilbert_transform).

The main object in the time-frequency domain is the [TFAnalyticSignal](@ref), from which the [TFAmplitude](@ref) and [TFPhase](@ref) objects are derived.
Taken together these three objects are referred to as [TFobjects](@ref) and vectors of them are a [TFobjectsVector](@ref) type.

|         Object          |         Vector type          |  Constructors from data  |
|:-----------------------:|:----------------------------:|:------------------------:|
| [TFAnalyticSignal](@ref)|[TFAnalyticSignalVector](@ref)|[`TFanalyticsignal`](@ref)|
|   [TFAmplitude](@ref)   |[TFAmplitudeVector](@ref)     |[`TFamplitude`](@ref)     |
|    [TFPhase](@ref)      |[TFPhaseVector](@ref)         |[`TFphase`](@ref)         |


## TFAnalyticSignal

**Categories**: [data objects](@ref), [TFobjects](@ref).

An analytic signal object has the following structure:

```
struct TFAnalyticSignal
    y          :: Matrix{T} where T<:Complex
    bandwidht  :: IntOrReal
    flabels    :: Vector{S} where S<:Real
    nonlinear  :: Bool
    fsmoothing :: Smoother
    tsmoothing :: Smoother
end
```

**Fields**:

`y`: a *complex* matrix holding the *analytic signal* in the time-frequency domain, with frequency band-pass regions (in Hz) in rows and time (samples) in columns.

`bandwidht`: the *bandwidht* (in Hz) of the filter bank band-pass regions.
See constructor [`TFanalyticsignal`](@ref) for details.
It can be an integer or a real number.

`flabels`: a vector holding all center frequencies (in Hz) of the filter bank
band-pass regions. Those are the *frequency labels* for the rows of `y`.

`non-linear`: a flag indicating whether the analytic signal has been normalized
so as to have amplitude``1.0`` at all points.
Such normalization allows non-linear univariate and bivariate estimations
(see [timefrequencyuni.jl](@ref) and [timefrequencybi.jl](@ref)).

`fsmoothing`: a flag of [Smoother::Enumerated Type](@ref) indicating
whether the analytic signal has been smoothed across adjacent
frequency band-pass regions. If no frequency smoothing has been applied,
it is equal to `noSmoother`, which is the default in all constructors.

`tsmoothing`: a flag of [Smoother::Enumerated Type](@ref) indicating
whether the analytic signal has been smoothed across adjacent
samples (time). If no time smoothing has been applied,
it is equal to `noSmoother`, which is the default in all constructors.

**Note**: In Julia the fields are accessed by the usual dot notation, e.g., you may verify that for *TFAnalyticSignal* object Y,
length(Y.flabels) == dim(Y.z, 1).

A vector of *TFAnalyticSignal* objects is of type [TFAnalyticSignalVector](@ref).

**Methods for TFAnalyticSignal and TFAnalyticSignalVector objects**

|        method          | TFAnalyticSignal | TFAnalyticSignalVector |
|:-----------------------|:---------------:|:---------------------:|
| [`amplitude`](@ref)    |       ✔        |           ✔           |
| [`phase`](@ref)        |       ✔        |           ✔           |
| [`polar`](@ref)        |       ✔        |                        |
| [`extract`](@ref)      |       ✔        |           ✔           |
| [`mean`](@ref)         |       ✔        |           ✔           |
| [`smooth`](@ref)       |       ✔        |           ✔           |
| [`sameParams`](@ref)   |                |           ✔           |
| [`isLinear`](@ref)     |                |           ✔           |
| [`isNonLinear`](@ref)  |                |           ✔           |
| [`meanAmplitude`](@ref)|                |           ✔           |
| [`concentration`](@ref)|                |           ✔           |
| [`comodulation`](@ref) |                |           ✔           |
| [`meanDirection`](@ref)|                |           ✔           |
| [`coherence`](@ref)    |                |           ✔           |


**Generic Constructors**:

In order to construct a *TFAnalyticSignal* object from univariate
data, *FourierAnalysis* provides two [`TFanalyticsignal`](@ref) constuctors,
which is what you will use in practice most of the time.

Manual constructors are also possible, for which you have to provide
appropriate arguments. The default manual constructor of *TFAnalyticSignal*
objects is

```
TFAnalyticSignal(y, bandwidht, flabels, nonlinear, fsmoothing, tsmoothing).
```

No other generic constructor is provided for this object.


## TFAmplitude

**Categories**: [data objects](@ref), [TFobjects](@ref).

An amplitude object has the following structure:

```
struct TFAmplitude
    y          :: Matrix{T} where T<:Real
    bandwidht  :: IntOrReal
    flabels    :: Vector{S} where S<:Real
    fsmoothing :: Smoother
    tsmoothing :: Smoother
    func       :: Function
end
```

**Fields**:

`y`: a *real* matrix holding the *amplitude* (modulus) of an analytic signal
in the time-frequency domain, often referred to as the *envelope*,
with frequency band-pass regions (in Hz) in rows and time (samples) in columns.

`bandwidht`: the *bandwidht* (in Hz) of the filter bank band-pass regions.
See constructor [`TFanalyticsignal`](@ref) for details.
It can be an integer or a real number.

`flabels`: a vector holding all center frequencies (in Hz) of the filter bank
band-pass regions. Those are the *frequency labels* for the rows of `y`.

`fsmoothing`: a flag of [Smoother::Enumerated Type](@ref) indicating
whether the amplitude has been smoothed across adjacent
frequency band-pass regions. If no frequency smoothing has been applied,
it is equal to `noSmoother`, which is the default in all constructors.

`tsmoothing`: a flag of [Smoother::Enumerated Type](@ref) indicating
whether the amplitude has been smoothed across adjacent
samples (time). If no time smoothing has been applied,
it is equal to `noSmoother`, which is the default in all constructors.

**Note**: Smoothing flags in this object indicate that the amplitude has
been smoothed, whereas in the [TFAnalyticSignal](@ref) object indicate
that the analytic signal has been smoothed. Note that it is not equivalent
to obtain amplitude from smoothed analytic signal
(e.g., using the [`amplitude`](@ref) function) or to smooth the amplitude
of analytic signal, e.g., using the [`TFamplitude`](@ref) constructor.

`func`: A name of a function that has been applied element-wise to the
matrix `.y` holding the amplitude. All constructors from data by default
use the `identity` (do nothing) function.

A vector of *TFAmplitude* objects is of type [TFAmplitudeVector](@ref).

**Methods for TFAmplitude and TFAmplitudeVector objects**

|        method          |   TFAmplitude   |   TFAmplitudeVector   |
|:-----------------------|:---------------:|:---------------------:|
| [`extract`](@ref)      |       ✔        |           ✔           |
| [`mean`](@ref)         |       ✔        |           ✔           |
| [`smooth`](@ref)       |       ✔        |           ✔           |
| [`sameParams`](@ref)   |                |           ✔           |
| [`isLinear`](@ref)     |                |           ✔           |
| [`isNonLinear`](@ref)  |                |           ✔           |
| [`meanAmplitude`](@ref)|                |           ✔           |
| [`comodulation`](@ref) |                |           ✔           |

**Generic Constructors**:

In order to construct a *TFAmplitude* object from univariate
data, *FourierAnalysis* provides four [`TFamplitude`](@ref) constuctors,
which is what you will use in practice most of the time.

Manual constructors are also possible, for which you have to provide
appropriate arguments. The default manual constructor of *TFAmplitude*
objects is

```
TFAmplitude(y, bandwidht, flabels, fsmoothing, tsmoothing, func).
```

Other generic constructors are also provided:

```
TFAmplitude(y, bandwidht, flabels, fsmoothing, tsmoothing)
```
enables construction giving only `y`, `bandwidth`, `flabels`, `fsmoothing`
and `tsmoothing`. `func` is set automatically to `identity`;

```
TFAmplitude(y, bandwidht, flabels)
```
acts like the constructor above, but sets by default also both `fsmoothing`
and `tsmoothing` to `noSmoother`.


## TFPhase

**Categories**: [data objects](@ref), [TFobjects](@ref).

A phase object has the following structure:

```
struct TFPhase
    y          :: Matrix{T} where T<:Real
    bandwidht  :: IntOrReal
    flabels    :: Vector{S} where S<:Real
    nonlinear  :: Bool
    fsmoothing :: Smoother
    tsmoothing :: Smoother
    unwrapped  :: Bool
    func       :: Function
end
```

**Fields**:

`y`: a *real* matrix holding the *phase* (argument) of an analytic signal
in the time-frequency domain, with frequency band-pass regions (in Hz)
in rows and time (samples) in columns. By default the phase is
represented in ``[−π, π]``.

`bandwidht`: the *bandwidht* (in Hz) of the filter bank band-pass regions.
See constructor [`TFanalyticsignal`](@ref) for details.
It can be an integer or a real number.

`flabels`: a vector holding all center frequencies (in Hz) of the filter bank
band-pass regions. Those are the *frequency labels* for the rows of `y`.

`non-linear`: a flag indicating whether the phase has been estimated from
analytic signal normalized so as to have amplitude``1.0`` at all points.
Such normalization allows non-linear univariate and bivariate estimations
(see [timefrequencyuni.jl](@ref) and [timefrequencybi.jl](@ref)).

`fsmoothing`: a flag of [Smoother::Enumerated Type](@ref) indicating
whether the phase has been smoothed across adjacent
frequency band-pass regions. If no frequency smoothing has been applied,
it is equal to `noSmoother`, which is the default in all constructors.

`tsmoothing`: a flag of [Smoother::Enumerated Type](@ref) indicating
whether the phase has been smoothed across adjacent
samples (time). If no time smoothing has been applied,
it is equal to `noSmoother`, which is the default in all constructors.

**Note**: Smoothing raw phase estimates is unappropriate since the phase
is a discontinous function, however it makes sense to smooth phase
if the phase is unwrapped (see below).

**Note**: Smoothing flags in this object indicate that the (unwrapped) phase has
been smoothed, whereas in the [TFAnalyticSignal](@ref) object indicate
that the analytic signal has been smoothed. Note that it is not equivalent to
obtain unwrapped phase from smoothed analytic signal
(e.g., using the [`phase`](@ref) function) or to smooth the unwrapped phase
of analytic signal, e.g., using the [`TFphase`](@ref) constructor.

`unwrapped`: a flag indicating if the phase has been unwrapped.
The unwrapped phase is defined as the cumulative sum of the phase
along the time dimension once this is represented in ``[0, 2π]``.

`func`: A name of a function that has been applied element-wise to the
matrix `.y` holding the phase. All constructors from data by default
use the `identity` (do nothing) function.
- passing `func=x->x+π` will return the phase in [0, 2π],
- passing `func=x->x/π` will return the phase in [-1, 1],
- passing `func=sin` will return the sine of the phase.

A vector of *TFPhase* objects is of type [TFPhaseVector](@ref).

**Methods for TFPhase and TFPhaseVector objects**

|        method          |   TFAmplitude   |   TFAmplitudeVector   |
|:-----------------------|:---------------:|:---------------------:|
| [`unwrapPhase`](@ref)  |       ✔        |           ✔           |
| [`isUnwrapped`](@ref)  |       ✔        |           ✔           |
| [`extract`](@ref)      |       ✔        |           ✔           |
| [`mean`](@ref)         |       ✔        |           ✔           |
| [`smooth`](@ref)       |       ✔        |           ✔           |
| [`sameParams`](@ref)   |                |           ✔           |
| [`isLinear`](@ref)     |                |           ✔           |
| [`isNonLinear`](@ref)  |                |           ✔           |

**Generic Constructors**:

In order to construct a *TFPhase* object from univariate
data, *FourierAnalysis* provides four [`TFphase`](@ref) constuctors,
which is what you will use in practice most of the time.

Manual constructors are also possible, for which you have to provide
appropriate arguments. The default manual constructor of *TFPhase*
objects is

```
TFPhase(y, bandwidht, flabels, nonlinear,
        fsmoothing, tsmoothing, unwrapped, func).
```

Other generic constructors are also provided:

```
TFPhase(y, bandwidht, flabels, nonlinear, fsmoothing, tsmoothing)
```
enables construction giving only `y`, `bandwidht`, `flabels`, `fsmoothing`
and `tsmoothing`. `unwrapped` is set to `false` and `func` is set
to `identity`;

```
TFPhase(y, bandwidht, flabels)
```
acts like the constructor above, but sets by default also `nonlinear` to `true`
and both `fsmoothing` and `tsmoothing` to `noSmoother`.

**Constructors from data**:

```@docs
TFanalyticsignal
TFamplitude
TFphase
```
