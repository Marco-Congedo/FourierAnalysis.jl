#   Unit "timefrequencyuni" of the FourierAnalysis Package for julia language
#   v 0.0.1 - last update 5th of September 2019
#
#   MIT License
#   Copyright (c) 2019, Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENTS :
#   This unit implements time-frequency univariate measure based
#   on timefrequency.jl and tools.jl


################### UNIVARIATE MEASURE ##########################
# All measures are weighted versions of measures described in:
# Congedo(2018): https://hal.archives-ouvertes.fr/hal-01868538/document

# Mean Amplitude (MAmp), Eq. 0.1 therein: <|z|> = <|re^i𝜑|> = <r> if w=[].
# If w is a positive weight vector a weighted version <wr>/<w> is computed.
# NB: w must be normalized so that the denominator is 1.0
"""
```
(1)
function meanAmplitude( 𝐀      :: TFAmplitudeVector,
                        frange :: fInterval,
                        trange :: tInterval;
                    mode   :: Function = extract,
                    func   :: Function = identity,
                    w      :: Vector   = [],
                    check  :: Bool     = true)

(2)
function meanAmplitude( 𝐙 :: TFAnalyticSignalVector,
                    < same arguments as method (1) >

(3)
function meanAmplitude( 𝐱      :: Vector{Vector{T}},
                        sr     :: Int,
                        wl     :: Int,
                        frange :: fInterval,
                        trange :: tInterval,
                     bandwidht :: IntOrReal    = 2;
                 mode          :: Function     = extract,
                 func          :: Function     = identity,
                 w             :: Vector       = [],
                 check         :: Bool         = true,
                 filtkind      :: FilterDesign = Butterworth(2),
                 fmin          :: IntOrReal    = bandwidht,
                 fmax          :: IntOrReal    = sr÷2,
                 fsmoothing    :: Smoother     = noSmoother,
                 tsmoothing    :: Smoother     = noSmoother,
                 planner       :: Planner      = getplanner,
                 ⏩           :: Bool         = true) where T<:Real

```
**alias**: mamp

(1)
Given a [TFAmplitudeVector](@ref) object, estimate the
[(weighted) mean amplitude](@ref) measure across those objects.
The time-frequency planes of all the objects in `𝐀` should be congruent.

**arguments**:

`frange` and `trange` define a time-frequency region on which the estimation
is sought. `frange` is a [fInterval](@ref) type
and delimits center frequencies of the filter bank. `trange` is a
[tInterval](@ref) type and delimits time samples. To obtain the estimation
in the whole time-frequency plane use the column (:) sign for both arguments.

**optional keyword arguments**

If `mode=extract` is passed (default), the measure will be computed
for all points in the chosen time-frequency region. If `mode=mean` is passed,
it will be computed on the mean of these points (grand-average).
The [`extract`](@ref) and [`mean`](@ref) functions are
[generic methods](@ref) of *FourierAnalysis*.

**Note**: with `mode=mean` the output of the function is always a real number,
whereas with `mode=extract` the output may be a real number, a real
row or column vector or a real matrix, depending on the shape of the
chosen time-frequency region.

Passing a function with the `func` argument you can
[derive your own time-frequency measures](@ref).

`w` may be a vector of non-negative real weights associated to
each input object. By default the unweighted version of the measure
is computed.

If `check` is true (default), check that if the column sign is passed
- as `frange` argument, all input objects have the same number of rows (center frequencies);
- as `trange` argument, all input objects have the same number of columns (time samples).
If either check fails, print an error message and return `Nothing`.
No other range checks are performed.

(2)
Given a [TFAnalyticSignalVector](@ref) object, compute the amplitude
of all objects in `𝐙` and estimate the [(weighted) mean amplitude](@ref)
measure across those objects as per method (1). In addition,
since using this method all [TFAnalyticSignal](@ref) in `𝐙` must be `linear`,
if `check` is true (default) and this is not the case,
print an error and return `Nothing`. The checks on `frange` and `trange`
performed by method (1) are also performed by this method.

(3)
Estimate the amplitude of all data vectors in `𝐱` calling the
[`TFamplitude`](@ref) constructor and then estimate the
[(weighted) mean amplitude](@ref) measure across the constructed
amplitude objects as per method (1).

`frange`, `trange`, `mode`, `func`, `w` and `check` have the same meaning as
in method (1). The other arguments are passed to the
[`TFamplitude`](@ref) constructor, to which the reader is referred
for their meaning.

**See also**: [`concentration`](@ref), [`meanDirection`](@ref),
[timefrequencybi.jl](@ref).

**Examples**:
```
using Plots, FourierAnalysis

# generate 100 vectors of data
sr, t, bandwidht=128, 512, 2
h=taper(harris4, t)
𝐱=[sinusoidal(2, 10, sr, t, 0).*h.y+randn(t) for i=1:100]

𝐘=TFanalyticsignal(𝐱, sr, t, bandwidht; fmax=32)
𝐀=TFamplitude(𝐘)

# mean amplitude in a TF region from a TFAnalyticSignalVector object
MAmp=meanAmplitude(𝐘, (4, 16), :)
heatmap(MAmp; c=:fire) # y axis labels are not correct

# mean amplitude in a TF region from a TFAmplitudeVector object
MAmp=meanAmplitude(𝐀, (4, 16), :)

# mean amplitude in a TF region directly from data
MAmp=meanAmplitude(𝐱, sr, t, (4, 16), :, bandwidht)

# NB: in the above, the analytic signal objects must all
# be linear, since meanAmplitude is computed from amplitude
# and the amplitude of non-linear analytic signal is uniformy equal to 1.

# All these computations can be obtained averaging in a TF region, e.g.,
MAmp=meanAmplitude(𝐘, (4, 16), :; mode=mean) # output a real number

# and can be obtained on smoothed Amplitude, e.g.,
MAmp=meanAmplitude(𝐱, sr, t, (4, 16), :;
                   fsmoothing=blackmanSmoother,
                   tsmoothing=blackmanSmoother)
# or, equivalently, and using the alias `mamp`,
MAmp=mamp(smooth(blackmanSmoother, blackmanSmoother, 𝐀), (4, 16), :)

# A similar syntax is used for the other univariate measures, e.g.,
# concentration averaging in a TF region from a TFAnalyticSignalVector object
ConM=concentration(𝐘, (4, 16), (128, 384); mode=mean)

# concentration in a TF region directly from data (using the alias `con`)
ConE=con(𝐱, sr, t, (4, 16), (128, 384), bandwidht; mode=extract)
heatmap(Con; c=:fire) # y axis labels are not correct

NB: ConM is not at all equivalent to mean(ConE) !

# mean direction averaging in a TF region directly from data
MDir=meanDirection(𝐱, sr, t, (4, 16), :, bandwidht; mode=mean)

# mean direction in a TF region from a TFAnalyticSignalVector object
MDir=meanDirection(𝐘, (4, 16), :)

# and for the non-linear counterpart:
# phase concentration in a TF region directly from data
Con=concentration(𝐱, sr, t, (8, 12), :; nonlinear=true)

# phase concentration at a single TF point
Con=concentration(𝐱, sr, t, 10, 256; nonlinear=true)

# phase mean direction averaging in a TF region directly from data
# and using the alias `mdir`
MDir=mdir(𝐱, sr, t, (8, 12), :; mode=mean, nonlinear=true)

# If you try to compute a non-linear measure from a linear
# AnalyticSignal object you will get en error (see the REPL), e.g.,
Con=con(𝐘, (8, 12), (1, 512); mode=mean, nonlinear=true)

# In order to compute non-linear measures from analytic signal objects
# first we need to compute non-linear analytic signal objects:
𝐘=TFanalyticsignal(𝐱, sr, t, bandwidht; fmax=32, nonlinear=true)

# then, we can obtain for example the phase concentration
Con=con(𝐘, (8, 12), :; mode=mean, nonlinear=true)

# and the phase mean direction
MDir=meanDirection(𝐘, (8, 12), :; nonlinear=true)
```
"""
meanAmplitude( 𝐀      :: TFAmplitudeVector,
               frange :: fInterval,
               trange :: tInterval;
           mode   :: Function = extract,
           func   :: Function = identity,
           w      :: Vector   = [],
           check  :: Bool     = true) =
    if !check || _checkUniData(𝐀, frange, trange, false, true, "meanAmplitude")
        return mean(func(mode(𝐀, frange, trange; w=w, check=check)))
    end


meanAmplitude( 𝐙      :: TFAnalyticSignalVector,
               frange :: fInterval,
               trange :: tInterval;
           mode   :: Function = extract,
           func   :: Function = identity,
           w      :: Vector   = [],
           check  :: Bool     = true) =
    if !check || _checkUniData(𝐙, frange, trange, false, true, "meanAmplitude")
        return mean(func(mode(TFamplitude(𝐙), frange, trange; w=w, check=check)))
    end


# NB  if `t` > 2^14 then `t` is set to 2^10.
meanAmplitude( 𝐱      :: Vector{Vector{T}},
               sr     :: Int,
               wl     :: Int,
               frange :: fInterval,
               trange :: tInterval,
            bandwidht :: IntOrReal    = 2;  # <>
        mode          :: Function     = extract,
        func          :: Function     = identity,
        w             :: Vector       = [],
        check         :: Bool         = true,
        filtkind      :: FilterDesign = Butterworth(2),
        fmin          :: IntOrReal    = bandwidht,
        fmax          :: IntOrReal    = sr÷2,
        fsmoothing    :: Smoother     = noSmoother,
        tsmoothing    :: Smoother     = noSmoother,
        planner       :: Planner      = getplanner,
        ⏩           :: Bool         = true) where T<:Real =

    meanAmplitude(TFamplitude(𝐱,
                              sr,
                              wl,
                              bandwidht;
                        fmin       = fmin,
                        fmax       = fmax,
                        fsmoothing = fsmoothing,
                        tsmoothing = tsmoothing,
                        filtkind   = filtkind,
                        planner    = planner,
                        ⏩        = ⏩),
                  frange,
                  trange;
            mode  = mode,
            func  = func,
            w     = w,
            check = check)


mamp = meanAmplitude



"""
```
(1)
function concentration( 𝐙       :: TFAnalyticSignalVector,
                        frange  :: fInterval,
                        trange  :: tInterval;
                    nonlinear :: Bool     = false,
                    mode      :: Function = extract,
                    func      :: Function = identity,
                    w         :: Vector   = [],
                    check     :: Bool     = true)

(2)
function concentration( 𝐱      :: Vector{Vector{T}},
                        sr     :: Int,
                        wl     :: Int,
                        frange :: fInterval,
                        trange :: tInterval,
                     bandwidht :: IntOrReal    = 2;
                 nonlinear  :: Bool         = false,
                 mode       :: Function     = extract,
                 func       :: Function     = identity,
                 w          :: Vector       = [],
                 check         :: Bool      = true,
                 filtkind   :: FilterDesign = Butterworth(2),
                 fmin       :: IntOrReal    = bandwidht,
                 fmax       :: IntOrReal    = sr÷2,
                 fsmoothing :: Smoother     = noSmoother,
                 tsmoothing :: Smoother     = noSmoother,
                 planner    :: Planner      = getplanner,
                 ⏩        :: Bool         = true) where T<:Real

```

**alias**: con

If optional keyword parameter `nonlinear` is false (default),
estimate the [(weighted) concentration](@ref) measure,
otherwise estimate the [(weighted) phase concentration](@ref) measure.

(1)
The desired measure is obtained averaging across the
[TFAnalyticSignal](@ref) objects in `𝐙`. Since this method uses
pre-computed analytic signal objects, their `.nonlinear` field
must agree with the `nonlinear` argument passed to
this function.

`frange`, `trange`, `w`, `mode` and `func` have the same meaning as
in the [`meanAmplitude`](@ref) function, however keep in mind that
the two possible `mode` functions, i.e.,
[`extract`](@ref) and [`mean`](@ref), in this function operate on
complex numbers.

The checks performed in the [`meanAmplitude`](@ref) function
are performed here too. In addition, if `check` is true, also check that
- if `nonlinear` is true, all objects in `𝐙` are nonlinear;
- if `nonlinear` is false, all objects in `𝐙` are linear.
If either check fails, print an error message and return `Nothing`.

(2)
Estimate the analytic signal of all data vectors in `𝐱` calling the
[`TFanalyticsignal`](@ref) constructor and then use method (1)
to obtained the desired measure.

`frange`, `trange`, `mode`, `func`, `w` and `check` have the same meaning as
in the [`meanAmplitude`](@ref) function.
The other arguments are passed to the [`TFanalyticsignal`](@ref) constructor,
to which the reader is referred for understanding their action.

**See also**: [`meanAmplitude`](@ref), [`meanDirection`](@ref),
[timefrequencybi.jl](@ref).

**Examples**: see the examples of [`meanAmplitude`](@ref) function.
"""
concentration( 𝐙         :: TFAnalyticSignalVector,
               frange    :: fInterval,
               trange    :: tInterval;
               nonlinear :: Bool     = false,
               mode      :: Function = extract,
               func      :: Function = identity,
               w         :: Vector   = [],
               check     :: Bool     = true) =
    if !check || _checkUniData(𝐙, frange, trange, nonlinear, !nonlinear, "concentration")
        return amplitude(mean(func(mode(𝐙, frange, trange; w=w, check=check))))
    end


# NB  if `t` > 2^14 then `t` is set to 2^10.
concentration( 𝐱      :: Vector{Vector{T}},
               sr     :: Int,
               wl     :: Int,
               frange :: fInterval,
               trange :: tInterval,
            bandwidht :: IntOrReal    = 2;
        nonlinear     :: Bool         = false,
        mode          :: Function     = extract,
        func          :: Function     = identity,
        w             :: Vector       = [],
        check         :: Bool         = true,
        filtkind      :: FilterDesign = Butterworth(2),
        fmin          :: IntOrReal    = bandwidht,
        fmax          :: IntOrReal    = sr÷2,
        fsmoothing    :: Smoother     = noSmoother,
        tsmoothing    :: Smoother     = noSmoother,
        planner       :: Planner      = getplanner,
        ⏩           :: Bool         = true) where T<:Real =

    concentration(TFanalyticsignal(𝐱,
                                   sr,
                                   wl,
                                   bandwidht;
                              fmin       = fmin,
                              fmax       = fmax,
                              fsmoothing = fsmoothing,
                              tsmoothing = tsmoothing,
                              filtkind   = filtkind,
                              nonlinear  = nonlinear,
                              planner    = planner,
                              ⏩        = ⏩),
                  frange,
                  trange;
              nonlinear = nonlinear,
              mode      = mode,
              func      = func,
              w         = w,
              check     = check)

con = concentration


"""
```
(1)
function meanDirection( 𝐙         :: TFAnalyticSignalVector,
                        frange    :: fInterval,
                        trange    :: tInterval;
                    nonlinear :: Bool       = false,
                    mode      :: Function   = extract,
                    func      :: Function   = identity,
                    w         :: Vector     = [],
                    check     :: Bool       = true)

(2)
function meanDirection( 𝐱      :: Vector{Vector{T}},
                        sr     :: Int,
                        wl     :: Int,
                        frange :: fInterval,
                        trange :: tInterval,
                     bandwidht :: IntOrReal    = 2;
                 nonlinear     :: Bool         = false,
                 mode          :: Function     = extract,
                 func          :: Function     = identity,
                 w             :: Vector       = [],
                 check         :: Bool         = true,
                 filtkind      :: FilterDesign = Butterworth(2),
                 fmin          :: IntOrReal    = bandwidht,
                 fmax          :: IntOrReal    = sr÷2,
                 fsmoothing    :: Smoother     = noSmoother,
                 tsmoothing    :: Smoother     = noSmoother,
                 planner       :: Planner      = getplanner,
                 ⏩           :: Bool         = true) where T<:Real
```
This function features two methods that use exactly the same syntax
as the two corresponding methods of the [`concentration`](@ref) function.
All arguements have exactly the same meaning as therein.
Only the output differs:

If optional keyword parameter `nonlinear` is false (default),
estimate the [(weighted) mean direction](@ref) measure,
otherwise estimate the [(weighted) phase mean direction](@ref) measure.

**alias**: mdir

**See also**: [`meanAmplitude`](@ref), [`concentration`](@ref),
[timefrequencybi.jl](@ref).

**Examples**: see the examples of [`meanAmplitude`](@ref).
"""
meanDirection( 𝐙         :: TFAnalyticSignalVector,
               frange    :: fInterval,
               trange    :: tInterval;
           nonlinear :: Bool       = false,
           mode      :: Function   = extract,
           func      :: Function   = identity,
           w         :: Vector     = [],
           check     :: Bool       = true) =
    if !check || _checkUniData(𝐙, frange, trange, nonlinear, !nonlinear, "meanDirection")
        return phase(mean(func(mode(𝐙, frange, trange; w=w, check=check))))
    end


# NB  if `t` > 2^14 then `t` is set to 2^10.
meanDirection( 𝐱      :: Vector{Vector{T}},
               sr     :: Int,
               wl     :: Int,
               frange :: fInterval,
               trange :: tInterval,
            bandwidht :: IntOrReal    = 2;
        nonlinear     :: Bool         = false,
        mode          :: Function     = extract,
        func          :: Function     = identity,
        w             :: Vector       = [],
        check         :: Bool         = true,
        filtkind      :: FilterDesign = Butterworth(2),
        fmin          :: IntOrReal    = bandwidht,
        fmax          :: IntOrReal    = sr÷2,
        fsmoothing    :: Smoother     = noSmoother,
        tsmoothing    :: Smoother     = noSmoother,
        planner       :: Planner      = getplanner,
        ⏩           :: Bool         = true) where T<:Real =

    meanDirection(TFanalyticsignal(𝐱,
                                   sr,
                                   wl,
                                   bandwidht;
                               fmin       = fmin,
                               fmax       = fmax,
                               fsmoothing = fsmoothing,
                               tsmoothing = tsmoothing,
                               filtkind   = filtkind,
                               nonlinear  = nonlinear,
                               planner    = planner,
                               ⏩        = ⏩),
                  frange,
                  trange;
              nonlinear = nonlinear,
              mode      = mode,
              func      = func,
              w         = w,
              check     = check)

mdir = meanDirection


# Check that:
# - if `frange` is a Colon (:), then the data in all objects in `𝒀` have the same number of rows,
# - if `trange` is a Colon (:), then the data in all objects in `𝒀` have the same number of columns,
# - if `mustBeAllNonLinear` is `true`, then the `nonlinear` field in all objects in `𝒀` is set to `true`,
# - if `mustBeAllLinear` is `true`, then the `nonlinear` field in all objects in `𝒀` is set to `false`.
# Return `true` if all checks are positive, `false` otherwise.
function _checkUniData(𝒀::TFobjectsVector, frange::fInterval, trange::tInterval,
                       mustBeAllNonLinear::Bool, mustBeAllLinear::Bool, funcname::String)
    OK = true
    frange isa Colon && !_allsame([size(Y.y, 1) for Y ∈ 𝒀]) ?
        (@error 📌*", "*funcname*" function: the number of filter banks is not the same in all input TF objects. Don't use a colon for the `frange` or make sure the input TF objects are of similar size"; OK = false) : nothing
    trange isa Colon && !_allsame([size(Y.y, 2) for Y ∈ 𝒀]) ?
        (@error 📌*", "*funcname*" function: the number of samples is not the same in all input TF objects. Don't use a colon for the `trange` or make sure the input TF objects are of similar size"; OK = false) : nothing
    mustBeAllNonLinear && !isNonLinear(𝒀) ?
        (@error 📌*", all input TF objects must be nonlinear in order to compute "*funcname; OK = false) : nothing
    mustBeAllLinear && !isLinear(𝒀) ?
        (@error 📌*", all input TF objects must be linear in order to compute "*funcname; OK = false) : nothing
    return OK
end
