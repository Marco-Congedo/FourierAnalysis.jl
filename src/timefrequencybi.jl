#   Unit "timefrequencybi" of the FourierAnalysis Package for julia language
#   v 0.0.1 - last update 5th of September 2019
#
#   MIT License
#   Copyright (c) 2019, Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENTS :
#   This unit implements time-frequency bivariate measure based
#   on timefrequency.jl and tools.jl


################### BIVARIATE MEASURE ##########################
# see Congedo, 2018: https://hal.archives-ouvertes.fr/hal-01868538v2/document

"""
```
(1)
function comodulation( 𝐀₁     :: TFAnalyticSignalVector,
                       𝐀₂     :: TFAnalyticSignalVector,
                       frange :: fInterval,
                       trange :: tInterval;
                  mode  :: Function = extract,
                  func  :: Function = identity,
                  w     :: Vector = [],
                  check :: Bool   = true) =

(2)
function comodulation( 𝐙₁     :: TFAnalyticSignalVector,
                       𝐙₂     :: TFAnalyticSignalVector,
    < arguments frange, trange, mode, func, w and check as in method (1) >

(3)
function comodulation(𝐱₁        :: Vector{Vector{T}},
                      𝐱₂        :: Vector{Vector{T}},
                      sr        :: Int,
                      wl        :: Int,
                      frange    :: fInterval,
                      trange    :: tInterval,
                      bandwidht :: IntOrReal    = 2;
                mode            :: Function     = extract,
                func            :: Function     = identity,
                w               :: Vector       = [],
                filtkind        :: FilterDesign = Butterworth(2),
                fmin            :: IntOrReal    = bandwidht,
                fmax            :: IntOrReal    = sr÷2,
                fsmoothing      :: Smoother     = noSmoother,
                tsmoothing      :: Smoother     = noSmoother,
                planner         :: Planner      = getplanner,
                ⏩             :: Bool         = true) where T<:Real

```
**alias**: com

(1)

Given a pair of [TFAmplitudeVector](@ref) objects, estimate their
[(weighted) comodulation](@ref) measure. `𝐀₁` and `𝐀₂` must hold the
same number of objects and the time-frequency planes
of all the objects should be congruent.

**arguments**:

`frange` and `trange` have the same meaning
as in the [`meanAmplitude`](@ref) method.

**optional keyword arguments**

`mode`, `func` and `check` have the same meaning
as in the [`meanAmplitude`](@ref) method.

`w` may be a vector of non-negative real weights associated to
each pair of input objects. By default the unweighted version of the measure
is computed.

(2)

Given a pair of [TFAnalyticSignalVector](@ref) object,
compute the amplitude of all objects and estimate the
[(weighted) comodulation](@ref) as per method (1).
Like in method (1), `𝐙₁` and `𝐙₂` must hold the same number of objects and
the time-frequency planes of all the objects should be congruent.
In addition, since using this method all [TFAnalyticSignal](@ref) objects
in `𝐙₁` and `𝐙₂` must be `linear`,
if `check` is true (default) and this is not the case,
print an error and return `Nothing`. The checks on `frange` and `trange`
performed by method (1) are also performed by this method.

(3)

Estimate the amplitude of all data vectors in `𝐱₁` and `𝐱₂` calling the
[`TFamplitude`](@ref) constructor and then estimate the
[(weighted) comodulation](@ref) measure across the constructed
amplitude objects as per method (1).

`frange`, `trange`, `mode`, `func`, `w` and `check` have the same meaning as
in method (1). The other arguments are passed to the
[`TFamplitude`](@ref) constructors, to which the reader is referred
for their meaning.

If a `planner` for FFT computations is not provided, it is computed once
and applied for all amplitude estimations.

**See also**: [`coherence`](@ref), [timefrequencybi.jl](@ref).

**Examples**:
```
using Plots, FourierAnalysis

# generate 100 pairs of data vectors
sr, t, bandwidht=128, 512, 2
h=taper(harris4, t)
𝐱₁=[sinusoidal(2, 10, sr, t, 0).*h.y+randn(t) for i=1:100]
𝐱₂=[sinusoidal(2, 10, sr, t, 0).*h.y+randn(t) for i=1:100]

# compute their (linear) analytic signal
𝐘₁=TFanalyticsignal(𝐱₁, sr, wl, bandwidht; fmax=32, nonlinear=false)
𝐘₂=TFanalyticsignal(𝐱₂, sr, wl, bandwidht; fmax=32, nonlinear=false)

# compute their amplitude
𝐀₁=TFamplitude(𝐘₁)
𝐀₂=TFamplitude(𝐘₂)

# compute the Com averaging in a TF region from TFAnalyticSignal objects
# 𝐘₁ and 𝐘₂ must be linear
Com=comodulation(𝐘₁, 𝐘₂, (8, 12), :; mode=mean)

# compute the Com averaging in a TF region from TFAmplitudeVector objects
# 𝐀₁ and 𝐀₂ must be linear
Com=comodulation(𝐀₁, 𝐀₂, (8, 12), :; mode=mean)

# compute the Com averaging in a TF region directly from data
# In this case you don't have to worry about linearity
Com=comodulation(𝐱₁, 𝐱₂, sr, wl, (8, 12), :, bandwidht; mode=mean)

# compute comodulation from smoothed amplitude:
Com=comodulation(𝐱₁, 𝐱₂, sr, wl, (8, 12), :, bandwidht;
                 mode=mean,
                 fsmoothing=blackmanSmoother,
                 tsmoothing=blackmanSmoother)

# you can go faster pre-computing a FFTW plan.
# This is useful when you have to call the comodulation function several times
plan=Planner(plan_patient, 5, wl, Float64, true)
Com=comodulation(𝐱₁, 𝐱₂, sr, wl, (8, 12), :, bandwidht; mode=mean, planner=plan)

# compute the Com in a TF region from TFAnalyticSignalVector objects
Com=comodulation(𝐘₁, 𝐘₂, (8, 12), :; mode=extract)

# compute the Com in a TF region from TFAmplitudeVector objects
Com=comodulation(𝐀₁, 𝐀₂, (8, 12), :; mode=extract)

# compute the Com in a TF region directly from data
Com=comodulation(𝐱₁, 𝐱₂, sr, wl, (8, 12), :, bandwidht; mode=extract)

# All these operations can be done also for coherence measures, for example
Coh=coherence(𝐘₁, 𝐘₂, (8, 12), :; mode=mean)

Coh=coherence(𝐘₁, 𝐘₂, (8, 12), :; mode=extract)

# Compute all 5 coherence types
Coh=coherence(𝐘₁, 𝐘₂, (8, 12), :; mode=extract, allkinds=true)

# phase coherence (phase-locking value)
# we obtain this measure from non-linear TFAnalyticSignalVector objects
𝐘₁=TFanalyticsignal(𝐱₁, sr, wl, bandwidht; fmax=32, nonlinear=true)
𝐘₂=TFanalyticsignal(𝐱₂, sr, wl, bandwidht; fmax=32, nonlinear=true)

Coh=coherence(𝐘₁, 𝐘₂, (8, 12), :; mode=mean, nonlinear=true)

# or directly from data (no need to worry about non-linearity in this case)
Coh=coherence(𝐱₁, 𝐱₂, sr, wl, (8, 12), :, bandwidht; mode=mean, nonlinear=true)

```
"""
comodulation( 𝐀₁     :: TFAmplitudeVector,
              𝐀₂     :: TFAmplitudeVector,
              frange :: fInterval,
              trange :: tInterval;
          mode       :: Function = extract,
          func       :: Function = identity,
          w          :: Vector = [],
          check      :: Bool   = true) =

    if !check || _checkBiData(𝐀₁, 𝐀₂, frange, trange, false, true, "comodulation")
        b(f::Function, g::Function)=
            mean(g(f(𝐀₁*𝐀₂, frange, trange; w=w, check=check)))
        u(f::Function, g::Function, 𝐀::TFAmplitudeVector)=
            (mean(g(f(sqr(𝐀), frange, trange; w=w, check=check))))
        return b(mode, func) ./ sqrt.(u(mode, func, 𝐀₁) .* u(mode, func, 𝐀₂))
    end


comodulation( 𝐙₁     :: TFAnalyticSignalVector,
              𝐙₂     :: TFAnalyticSignalVector,
              frange :: fInterval,
              trange :: tInterval;
          mode       :: Function = extract,
          func       :: Function = identity,
          w          :: Vector = [],
          check      :: Bool   = true) =

    if !check || _checkBiData(𝐙₁, 𝐙₂, frange, trange, false, true, "comodulation")
        b(f::Function, g::Function)=
            mean(g(f(TFamplitude(𝐙₁)*TFamplitude(𝐙₂), frange, trange; w=w, check=check)))
        u(f::Function, g::Function, 𝐙::TFAnalyticSignalVector)=
            (mean(g(f(TFpower(𝐙), frange, trange; w=w, check=check))))
        return b(mode, func) ./ sqrt.(u(mode, func, 𝐙₁) .* u(mode, func, 𝐙₂))
    end


# NB  if `t` > 2^14 then `t` is set to 2^10.
function comodulation(𝐱₁        :: Vector{Vector{T}},
                      𝐱₂        :: Vector{Vector{T}},
                      sr        :: Int,
                      wl        :: Int,
                      frange    :: fInterval,
                      trange    :: tInterval,
                      bandwidht :: IntOrReal    = 2;
                mode            :: Function     = extract,
                func            :: Function     = identity,
                w               :: Vector       = [],
                filtkind        :: FilterDesign = Butterworth(2),
                fmin            :: IntOrReal    = bandwidht,
                fmax            :: IntOrReal    = sr÷2,
                fsmoothing      :: Smoother     = noSmoother,
                tsmoothing      :: Smoother     = noSmoother,
                planner         :: Planner      = getplanner,
                ⏩             :: Bool         = true) where T<:Real

    planner ≠ getplanner ? plan=planner : plan=Planner(plan_estimate, -1.0, wl, T, true)

    _TFamplitude(𝐱::Vector{Vector{T}}) =
        TFamplitude(𝐱,
                    sr,
                    wl,
                    bandwidht;
                filtkind   = filtkind,
                fmin       = fmin,
                fmax       = fmax,
                fsmoothing = fsmoothing,
                tsmoothing = tsmoothing,
                planner    = plan,
                ⏩        = ⏩)

    comodulation(_TFamplitude(𝐱₁),
                 _TFamplitude(𝐱₂),
                 frange,
                 trange;
             mode  = mode,
             func  = func,
             w     = w,
             check = false)
end

com = comodulation



"""
```
(5)
function coherence(𝐙₁     :: TFAnalyticSignalVector,
                   𝐙₂     :: TFAnalyticSignalVector,
                   frange :: fInterval,
                   trange :: tInterval;
              nonlinear :: Bool     = false,
              allkinds  :: Bool     = false,
              mode      :: Function = extract,
              func      :: Function = identity,
              w         :: Vector   = [],
              check     :: Bool     = true)

(6)
function coherence(𝐱₁        :: Vector{Vector{T}},
                   𝐱₂        :: Vector{Vector{T}},
                   sr        :: Int,
                   wl        :: Int,
                   frange    :: fInterval,
                   trange    :: tInterval,
                   bandwidht :: IntOrReal = 2;
              nonlinear  :: Bool         = false,
              allkinds   :: Bool         = false,
              mode       :: Function     = extract,
              func       :: Function     = identity,
              w          :: Vector       = [],
              filtkind   :: FilterDesign = Butterworth(2),
              fmin       :: IntOrReal    = bandwidht,
              fmax       :: IntOrReal    = sr÷2,
              fsmoothing :: Smoother     = noSmoother,
              tsmoothing :: Smoother     = noSmoother,
              planner    :: Planner      = getplanner,
              ⏩        :: Bool         = true) where T<:Real

```
**alias**: coh

If optional keyword parameter `nonlinear` is false (default),
estimate linear [(weighted) coherence](@ref) measure,
otherwise estimate the the [(weighted) phase concentration](@ref) measure.

If optional keyword argument `allkinds` is true all five
[kinds of coherence](@ref) are returned. In this case the output
is a 5-tuple of `Coherence` matrices, in the order:
- *total* coherence,
- *real* coherence,
- *instantaneous* coherence
- *imaginary* coherence,
- *lagged* coherence.

If `allkinds` is false (default) only the *total* (classical) coherence
is returned as a single `Coherence` matrix.

(5)
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

(6)
Estimate the analytic signal of all data vectors in `𝐱` calling the
[`TFanalyticsignal`](@ref) constructor and then use method (5)
to obtained the desired measure.

`frange`, `trange`, `mode`, `func`, `w` and `check` have the same meaning as
in the [`meanAmplitude`](@ref) function.
The other arguments are passed to the [`TFanalyticsignal`](@ref) constructor,
to which the reader is referred for understanding their action.

**See also**: [`meanAmplitude`](@ref), [`meanDirection`](@ref),
[timefrequencybi.jl](@ref).

**Examples**: see the examples of [`comodulation`](@ref) function.
"""
coherence(𝐙₁     :: TFAnalyticSignalVector,
          𝐙₂     :: TFAnalyticSignalVector,
          frange :: fInterval,
          trange :: tInterval;
      nonlinear :: Bool     = false,
      allkinds  :: Bool     = false,
      mode      :: Function = extract,
      func      :: Function = identity,
      w         :: Vector   = [],
      check     :: Bool     = true) =

    if !check || _checkBiData(𝐙₁, 𝐙₂, frange, trange, nonlinear, !nonlinear, "coherence with argument nonlinear=$nonlinear")
        𝐙₁₂=𝐙₁*conj(𝐙₂)
        bRe(f::Function, g::Function)=
            power(mean(g(f(real(𝐙₁₂), frange, trange; w=w, check=check))))
        bIm(f::Function, g::Function)=
            power(mean(g(f(imag(𝐙₁₂), frange, trange; w=w, check=check))))
        u(f::Function, g::Function, 𝐙::TFAnalyticSignalVector)=
            mean(g(f(TFpower(𝐙), frange, trange; w=w, check=check)))
        reGij=bRe(mode, func)
        imGij=bIm(mode, func)
        if allkinds
            if nonlinear
                reGij isa Number ? uno=1 : uno=ones(eltype(reGij), size(reGij))
                NormC=reGij+imGij
                RealC=reGij
                ImagC=imGij
                InstC=reGij ./ (uno-imGij)
                LaggC=imGij ./ (uno-reGij)
                return [NormC, RealC, ImagC, InstC, LaggC]
            else # linear
                GiGj=u(mode, func, 𝐙₁) .* u(mode, func, 𝐙₂)
                NormC=(reGij+imGij) ./ GiGj
                RealC=reGij ./ GiGj
                ImagC=imGij ./ GiGj
                InstC=reGij ./ (GiGj-imGij)
                LaggC=imGij ./ (GiGj-reGij)
                return [NormC, RealC, ImagC, InstC, LaggC]
            end
        else # only classical
            return nonlinear ? reGij+imGij :
                              (reGij+imGij) ./ ( u(mode, func, 𝐙₁) .* u(mode, func, 𝐙₂) )
        end # if allkinds
    end


function coherence(𝐱₁        :: Vector{Vector{T}},
                   𝐱₂        :: Vector{Vector{T}},
                   sr        :: Int,
                   wl        :: Int,
                   frange    :: fInterval,
                   trange    :: tInterval,
                   bandwidht :: IntOrReal = 2;
             nonlinear  :: Bool         = false,
             allkinds   :: Bool         = false,
             mode       :: Function     = extract,
             func       :: Function     = identity,
             w          :: Vector       = [],
             filtkind   :: FilterDesign = Butterworth(2),
             fmin       :: IntOrReal    = bandwidht,
             fmax       :: IntOrReal    = sr÷2,
             fsmoothing :: Smoother     = noSmoother,
             tsmoothing :: Smoother     = noSmoother,
             planner    :: Planner      = getplanner,
             ⏩        :: Bool         = true) where T<:Real

    planner ≠ getplanner ? plan=planner : plan=Planner(plan_estimate, -1.0, wl, T, true)

    _TFanalyticsignal(𝐱::Vector{Vector{T}}) =
        TFanalyticsignal(𝐱,
                         sr,
                         wl,
                         bandwidht;
                     filtkind   = filtkind,
                     fmin       = fmin,
                     fmax       = fmax,
                     fsmoothing = fsmoothing,
                     tsmoothing = tsmoothing,
                     nonlinear  = nonlinear,
                     planner    = plan,
                    ⏩        = ⏩)

    coherence(_TFanalyticsignal(𝐱₁),
              _TFanalyticsignal(𝐱₂),
              frange,
              trange;
          nonlinear = nonlinear,
          allkinds  = allkinds,
          mode      = mode,
          func      = func,
          w         = w,
          check     = false)
end

coh=coherence



# Check that:
# - if `frange` is a Colon (:), then the data in all objects in `𝒀` and `𝐙` have the same number of rows,
# - if `trange` is a Colon (:), then the data in all objects in `𝒀` and `𝐙` have the same number of columns,
# - if `mustBeAllNonLinear` is `true`, then the `nonlinear` field in all objects in `𝒀` and `𝐙` is set to `true`,
# - if `mustBeAllLinear` is `true`, then the `nonlinear` field in all objects in `𝒀` and `𝐙` is set to `false`.
# Return `true` if all checks are positive, `false` otherwise.
function _checkBiData(𝒀::TFobjectsVector, 𝐙::TFobjectsVector, frange::fInterval, trange::tInterval,
                       mustBeAllNonLinear::Bool, mustBeAllLinear::Bool, funcname::String)
    OK = true
    frange isa Colon && !_allsame(vcat([size(Y.y, 1) for Y ∈ 𝒀], [size(Z.y, 1) for Z ∈ 𝐙])) ?
        (@error 📌*", "*funcname*" function: the number of filter banks is not the same in all input TF objects. Don't use a colon for the `frange` or make sure the input TF objects are of similar size"; OK = false) : nothing
    trange isa Colon && !_allsame(vcat([size(Y.y, 2) for Y ∈ 𝒀], [size(Z.y, 2) for Z ∈ 𝐙])) ?
        (@error 📌*", "*funcname*" function: the number of samples is not the same in all input TF objects. Don't use a colon for the `trange` or make sure the input TF objects are of similar size"; OK = false) : nothing
    mustBeAllNonLinear && !isNonLinear(𝒀) && !isNonLinear(𝐙) ?
        (@error 📌*", all input TF objects must be nonlinear in order to compute "*funcname; OK = false) : nothing
    mustBeAllLinear && !isLinear(𝒀) && !isLinear(𝐙) ?
        (@error 📌*", all input TF objects must be linear in order to compute "*funcname; OK = false) : nothing
    return OK
end


#######################  internal utilities  ##############################
*(Z₁::TFAnalyticSignal, Z₂::TFAnalyticSignal) =
    TFAnalyticSignal(Z₁.y.*Z₂.y, Z₁.bandwidht, Z₁.flabels, Z₁.nonlinear, Z₁.fsmoothing, Z₁.tsmoothing)

*(𝐙₁::TFAnalyticSignalVector, 𝐙₂::TFAnalyticSignalVector) =
    TFAnalyticSignalVector([Z₁*Z₂ for (Z₁, Z₂) in zip(𝐙₁, 𝐙₂)])

*(A₁::TFAmplitude, A₂::TFAmplitude) =
    TFAmplitude(A₁.y.*A₂.y, A₁.bandwidht, A₁.flabels, A₁.fsmoothing, A₁.tsmoothing, A₁.func)

*(𝐀₁::TFAmplitudeVector, 𝐀₂::TFAmplitudeVector) =
    TFAmplitudeVector([A₁*A₂ for (A₁, A₂) in zip(𝐀₁, 𝐀₂)])

conj(Z::TFAnalyticSignal) =
    TFAnalyticSignal(conj(Z.y), Z.bandwidht, Z.flabels, Z.nonlinear, Z.fsmoothing, Z.tsmoothing)

conj(𝐙::TFAnalyticSignalVector) = TFAnalyticSignalVector([conj(Z) for Z ∈ 𝐙])

real(Z::TFAnalyticSignal) =
    TFAmplitude(real(Z.y), Z.bandwidht, Z.flabels, Z.fsmoothing, Z.tsmoothing, identity)

real(𝐙::TFAnalyticSignalVector) = TFAmplitudeVector([real(Z) for Z ∈ 𝐙])

imag(Z::TFAnalyticSignal) =
    TFAmplitude(imag(Z.y), Z.bandwidht, Z.flabels, Z.fsmoothing, Z.tsmoothing, identity)

imag(𝐙::TFAnalyticSignalVector) = TFAmplitudeVector([imag(Z) for Z ∈ 𝐙])

#dummy function for internal use only, temporary way to flag that it is not known what function has been applied
unknown()=Nothing

sqr(a)=a^2
sqr(A::TFAmplitude) = TFAmplitude(A.y.^2, A.bandwidht, A.flabels, A.fsmoothing, A.tsmoothing, A.func==identity ? sqr : unknown)
sqr(𝐀::TFAmplitudeVector) =  TFAmplitudeVector([sqr(A) for A ∈ 𝐀])

power(r::Real)=abs2(r)
power(c::Complex) = abs2(c)
power(Z::AbstractArray{T}) where T<:Real = abs2.(Z)
power(Z::AbstractArray{T}) where T<:Complex = abs2.(Z)

TFpower(Z::TFAnalyticSignal) =
    TFAmplitude(abs2.(Z.y), Z.bandwidht, Z.flabels, Z.fsmoothing, Z.tsmoothing, identity)

TFpower(𝐙::TFAnalyticSignalVector) =
    TFAmplitudeVector([TFpower(Z) for Z ∈ 𝐙])
######################  end internal utilities  ##############################
