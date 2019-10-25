#   Unit "timefrequency" of the FourierAnalysis Package for julia language
#   v 0.2.0 - last update 20th of October 2019
#
#   MIT License
#   Copyright (c) 2019, Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENTS :
#   This unit implements time-frequency representations based
#   on analytic signal estimations using FFTW and DSP filters.


#### Time-Frequency Analytic Signal ####

"""
```
(1)
function TFanalyticsignal(x         :: Vector{T},
                          sr        :: Int,
                          wl        :: Int          = 0,
                          bandwidth :: IntOrReal    = 2;
                      fmin          :: IntOrReal    = bandwidth,
                      fmax          :: IntOrReal    = sr÷2,
                      filtkind      :: FilterDesign = Butterworth(2),
                      nonlinear     :: Bool         = false,
                      fsmoothing    :: Smoother     = noSmoother,
                      tsmoothing    :: Smoother     = noSmoother,
                      planner       :: Planner      = getplanner,
                      ⏩           :: Bool         = true) where T<:Real

(2)
function TFanalyticsignal(𝐱         :: Vector{Vector{T}},
                      <same arguments as method (1)>

```
(1)

Given sampling rate `sr`, construct a [TFAnalyticSignal](@ref) object
from univariate data `x`, that is, a time-frequency representation of real
vector `x`.

Call three functions in series performing the three following operations:

i. pass the data throught a bank of pass-band filters calling function
[`filterbank`](@ref),

ii. compute the analytic signal calling function [`analyticsignal`](@ref)
(method (1) therein),

iii. If requested, smooth the analytic signal along the time and/or frequency
dimension calling function [`smooth`](@ref).

The arguments passed to each of these functions are:

|   filterbank   | analyticsignal |   smooth      |
|:--------------:|:--------------:|:-------------:|
|      `x`       |     `wl`       | `fsmoothing`  |
|     `sr`       |  `nonlinear`   | `tsmoothing`  |
|  `bandwidth`   |   `planner`    |               |
|  `filtkind`    |    `⏩`       |               |
|    `fmin`      |                |               |
|    `fmax`      |                |               |
|     `⏩`      |                |               |

Refer to the documentation of each function to learn the meaning
of each argument.

By default the `wl` argument is set to `length(x)`.

If `⏩` is true (default), filtering and computation fo the analytic signal
are multi-threaded across band-pass regions as long as the number of the regions
is at least twice the number of threads Julia is instructed to use.
See [Threads](@ref).

(2)

Given sampling rate `sr`, construct a
[TFAnalyticSignalVector](@ref) object from a vector of univariate data `𝐱`,
that is, from the time-frequency representations of the vectors in `𝐱`.

This method operates as method (1) with the following exceptions:

By default the `wl` argument is set to `length(x)` for each vectors in `𝐱`.
If another values is given, it will be used for all of them.

If `⏩` is true (default), the method is run in multi-threaded mode across the
vectors in `𝐱` as long as their number is at least twice the number of
threads Julia is instructed to use, otherwise this method attempts to run
each analytic signal estimation in multi-threaded mode like in method (1).
See [Threads](@ref).

If a `Planner` is not explicitly passed as an argument,
the FFTW plan is computed once and applied for all analytic signal
estimations.

**See**: [`filterbank`](@ref), [`analyticsignal`](@ref), [`smooth`](@ref).

**Examples**:
```
using Plots, FourierAnalysis

# generate some data
sr, t, bandwidth=128, 512, 2
h=taper(harris4, t)
x1=sinusoidal(10, 8, sr, t, 0)
x2=sinusoidal(10, 19, sr, t, 0)
x=Vector((x1+x2).*h.y+randn(t))
y1=sinusoidal(10, 6, sr, t, 0)
y2=sinusoidal(10, 16, sr, t, 0)
y=Vector((y1+y2).*h.y+randn(t))
plot([x, y])

# TFAnalyticSignal object for x (method (1))
Y = TFanalyticsignal(x, sr, sr*4; fmax=32)

# vector of TFAnalyticSignal objects for x and y (method (2))
𝒀 = TFanalyticsignal([x, y], sr, sr*4; fmax=32)

# (for shortening comments: TF=time-frequency, AS=analytic signal)

# mean AS in a TF region (8:12Hz and samples 1:128) for the first object
m=mean(𝒀[1], (8, 12), (1, 128)) # Output a complex number

# extract the AS in a TF region (8:12Hz and samples 1:128) for the first object
E=extract(𝒀[1], (8, 12), (1, 128)) # Output a complex matrix

# mean AS in a TF region (8:12Hz and samples 1:128) for the two objects
𝐦=mean(𝒀, (8, 12), (1, 128)) # Output a vector of complex numbers
# same computation without checking homogeneity of the two objects in 𝒀 (faster)
𝐦=mean(𝒀, (8, 12), (1, 128); check=false)

# extract the AS in a TF region (8:12Hz and samples 1:128) for the two objects
𝐄=extract(𝒀, (8, 12), (8, 12)) # Output a vector of complex matrices

# plot the real part of the AS of x (see unit recipes.jl)

# gather first useful attributes for the plot
tfArgs=(right_margin = 2mm,
        top_margin = 2mm,
        xtickfont = font(10, "Times"),
        ytickfont = font(10, "Times"))

heatmap(tfAxes(Y)..., real(Y.y);
        c=:pu_or, tfArgs...)

# ...the imaginary part
heatmap(tfAxes(Y)..., imag(Y.y);
        c=:bluesreds, tfArgs...)

# ...the amplitude
heatmap(tfAxes(Y)..., amplitude(Y.y);
        c=:amp, tfArgs...)

# ...the amplitude of the AS smoothed in frequency and time
heatmap(tfAxes(Y)..., amplitude(smooth(hannSmoother, hannSmoother, Y).y);
        c=:fire, tfArgs...)

# ...the phase
heatmap(tfAxes(Y)..., phase(Y.y);
        c=:bluesreds, tfArgs...)

# or generate a TFAmplitude object from the AS
A=TFamplitude(Y)
# and plot it (with different colors)
heatmap(tfAxes(A)..., A.y;
        c=:fire, tfArgs...)

# generate a TFPhase object
ϴ=TFphase(Y)
# and plot it (with custom colors)
heatmap(tfAxes(ϴ)..., ϴ.y;
        c=:pu_or, tfArgs...)

# compute and plot phase in [0, 2π]
heatmap(tfAxes(Y)..., TFphase(Y; func=x->x+π).y;
        c=:amp, tfArgs...)

# compute and plot unwrapped phase
heatmap(tfAxes(Y)..., TFphase(Y; unwrapped=true).y;
        c=:amp, tfArgs...)

# smooth time-frequency AS: smooth frequency
Z=smooth(blackmanSmoother, noSmoother, Y)

# plot amplitude of smoothed analytic signal
heatmap(tfAxes(Z)..., amplitude(Z.y);
        c=:amp, tfArgs...)

# not equivalently (!), create an amplitude object and smooth it:
# in this case the amplitude is smoothed, not the AS
A=smooth(blackmanSmoother, noSmoother, TFamplitude(Y))
heatmap(tfAxes(A)..., A.y;
        c=:fire, tfArgs...)

# Smoothing raw phase estimates is unappropriate
# since the phase is a discontinous function, however it makes sense
# to smooth phase if the phase is unwrapped.
heatmap(tfAxes(Y)...,
        smooth(blackmanSmoother, noSmoother, TFphase(Y; unwrapped=true)).y;
        c=:amp, tfArgs...)

# smooth AS: smooth both frequency and time
E=smooth(blackmanSmoother, blackmanSmoother, Y)

# plot amplitude of smoothed analytic signal
heatmap(tfAxes(E)..., amplitude(E.y);
        c=:fire, tfArgs...)

# plot phase of smoothed analytic signal
heatmap(tfAxes(E)..., phase(E.y);
        c=:bluesreds, tfArgs...)

# not equivalently (!), create amplitude and phase objects and smooth them
A=smooth(blackmanSmoother, blackmanSmoother, TFamplitude(Y))
heatmap(tfAxes(A)..., A.y;
        c=:fire, tfArgs...)

ϴ=smooth(blackmanSmoother, blackmanSmoother, TFphase(Y, unwrapped=true))
heatmap(tfAxes(ϴ)..., ϴ.y;
        c=:pu_or, tfArgs...)

# smooth again
ϴ=smooth(blackmanSmoother, blackmanSmoother, ϴ)
heatmap(tfAxes(ϴ)..., ϴ.y;
        c=:pu_or, tfArgs...)
# and again ...

# create directly smoothed AS
Y=TFanalyticsignal(x, sr, t, bandwidth;
                   fmax=32,
                   fsmoothing=hannSmoother,
                   tsmoothing=hannSmoother)

# plot amplitude of smoothed analytic signal
heatmap(tfAxes(Y)..., amplitude(Y.y);
        c=:amp, tfArgs...)

# create directly smoothed Amplitude
A=TFamplitude(x, sr, t, bandwidth;
              fmax=32,
              fsmoothing=hannSmoother,
              tsmoothing=hannSmoother)

# plot smoothed amplitude
heatmap(tfAxes(A)..., A.y;
        c=:amp, tfArgs...)

# compute a TFAnalyticSignal object with non-linear AS
Y=TFanalyticsignal(x, sr, t, bandwidth; fmax=32, nonlinear=true)

# check that it is non-linear
Y.nonlinear

# check that the amplitude is now 1.0 everywhere
norm(amplitude(Y.y)-ones(eltype(Y.y), size(Y.y))) # must be zero

# plot non-linear phase
heatmap(tfAxes(Y)..., phase(Y.y);
        c=:bkr, tfArgs...)

# get the center frequencies of TFAmplitude object A
A.flabels

# extract the amplitude in a time-frequency region
extract(A, (2, 10), (1, 256)) # output a real matrix

# extract the amplitude in a time-frequency region at only one frequency
extract(A, 10, (1, 256)) # output a row vector

# extract the amplitude at one temporal sample at one frequency
extract(A, 10, 12) # or extract(A, 10.0, 12)

# extract amplitude at one temporal sample in a frequency region
extract(A, (10, 12), 12) # or extract(A, (10.0, 12.0), 12)

# extract amplitude at one temporal sample and all frequencies
extract(A, :, 12) # output a (column) vector

# compute the mean in a time-frequency region:
mean(A, (2, 10), (1, 256)) # output a real number
# is equivalent to (but may be less efficient than)
mean(extract(A, (2, 10), (1, 256)))

# using column sign for extracting all time samples
extract(A, (2, 10), :)

# This:
extract(A, :, :)
# is equivalent to this:
A.y
# but if you don't need to extract all frequencies,
# use the extract function to control what frequencies will be extracted:
# This
extract(A, (4, 5), 10)
# is not equivalent to this
A.y[4:5, 10]
# since the `extract` function finds the correct rows corresponding
# to the sought frequencies (in Hz), while A.y[4:5, 10]
# just returns the elements [4:5, 10] in the TF amplitude object

# Although the first center frequency in A is 2Hz, its
# band-pass region is 1-3Hz, therefore the frequency range 1:10 is accepted
mean(A, (1, 10), (1, 256))
# but this result in an error (see the REPL) since 0.5 is out of range:
mean(A, (0.5, 10), (1, 256))

# using a colon sign for time range
a=mean(A, (1, 10), :)
# using an integer for time range (indicates one specific sample)
a=mean(A, (1, 10), 16)

# using a colon sign for frequency range
a=mean(A, :, (1, 16))
# using a real number for frequency range
a=mean(A, 8.5, (1, 16))

# NB: the `extract` and `mean` functions work with the same syntax
#     for objects TFAnayticSignal, TFAmplitude and TFPhase.
```
"""
function TFanalyticsignal(x         :: Vector{T},
                          sr        :: Int,
                          wl        :: Int          = 0,
                          bandwidth :: IntOrReal    = 2;
                      fmin          :: IntOrReal    = bandwidth,
                      fmax          :: IntOrReal    = sr÷2,
                      filtkind      :: FilterDesign = Butterworth(2),
                      nonlinear     :: Bool         = false,
                      fsmoothing    :: Smoother     = noSmoother,
                      tsmoothing    :: Smoother     = noSmoother,
                      planner       :: Planner      = getplanner,
                      ⏩           :: Bool         = true) where T<:Real
    # NB  if `t` > 2^14 then `t` is set to 2^10 by `anaytic signal function`.
    f, F=filterbank(x, sr, bandwidth;
                filtkind=filtkind, fmin=fmin, fmax=fmax, ⏩=⏩)

    smooth(fsmoothing,
           tsmoothing,
           TFAnalyticSignal(Matrix(transpose(analyticsignal(F,
                                                            wl==0 ? length(x) : wl;
                                                        nonlinear=nonlinear,
                                                        planner=planner,
                                                        ⏩=⏩))),
                            bandwidth,
                            f,
                            nonlinear,
                            noSmoother,
                            noSmoother),
            "TFanalyticsignal")
end


function TFanalyticsignal(𝐱         :: Vector{Vector{T}},
                          sr        :: Int,
                          wl        :: Int,
                          bandwidth :: IntOrReal    = 2;
                      filtkind      :: FilterDesign = Butterworth(2),
                      fmin          :: IntOrReal    = bandwidth,
                      fmax          :: IntOrReal    = sr÷2,
                      nonlinear     :: Bool         = false,
                      fsmoothing    :: Smoother     = noSmoother,
                      tsmoothing    :: Smoother     = noSmoother,
                      planner       :: Planner      = getplanner,
                      ⏩           :: Bool         = true) where T<:Real
    # NB  if `t` > 2^14 then `t` is set to 2^10 by `anaytic signal function`.
    # differently from previous method, wl here must be explicitly provided
    k = length(𝐱)
    𝐙 = TFAnalyticSignalVector(undef, k)

    planner ≠ getplanner ? plan=planner : plan=Planner(plan_estimate, -1.0, wl, T, true)

    function TF!(i::Int, threaded::Bool)
        f, F = filterbank(𝐱[i],
                          sr,
                          bandwidth;
                      filtkind=filtkind,
                      fmin=fmin,
                      fmax=fmax,
                      ⏩=threaded)

        𝐙[i] = TFAnalyticSignal(Matrix(transpose(analyticsignal(F,
                                                                wl;
                                                            nonlinear=nonlinear,
                                                            planner=plan,
                                                            ⏩=threaded))),
                                bandwidth,
                                f,
                                nonlinear,
                                noSmoother,
                                noSmoother)
    end

    _thread(⏩, k) ? (@threads for i=1:k TF!(i, false) end) : (for i=1:k TF!(i, true) end)
    return smooth(fsmoothing, tsmoothing, 𝐙, "TFanalyticsignal")
end



# ++++++++++++++++++++  Show override  +++++++++++++++++++ # (REPL output)
function Base.show(io::IO, ::MIME{Symbol("text/plain")}, Y::TFAnalyticSignal)
    r=size(Y.y, 1)
    c=size(Y.y, 2)
    l=length(Y.flabels)
    println(io, titleFont, "▤ TFAnayticSignal type; $r freq. x $c samples")
    #println(io, "□  □    □      □        □           □", defaultFont)
    println(io, separatorFont, "⭒  ⭒    ⭒      ⭒        ⭒           ⭒", defaultFont)
    println(io, "non-linear(.nonlinear): $(Y.nonlinear)")
    println(io, "freq. sm.(.fsmoothing): ", string(Y.fsmoothing))
    println(io, "time  sm.(.tsmoothing): ", string(Y.tsmoothing))
    println(io, "bandwidth             : $(Y.bandwidth) Hz")
    println(io, "freq. lab.  (.flabels): $(l)-", typeof(Y.flabels))
    println(io, "data              (.y): $(r)x$(c)-", typeof(Y.y))
    r≠l && @warn "number of frequency labels does not match the data matrix size" l r
end

function Base.show(io::IO, ::MIME{Symbol("text/plain")}, 𝐘::TFAnalyticSignalVector)
    println(io, titleFont, "▤ ⋯ ▤ TFAnalyticSignalVector Type")
    println(io, separatorFont, "⭒  ⭒    ⭒      ⭒        ⭒           ⭒", defaultFont)
    println(io, "$(length(𝐘))-element Vector{TFAnalyticSignal}")
end
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #


#### Time-Frequency Amplitude ####

#############################################################################
# Generic constructors of TFAmplitude objects:
##############################################
# Enable construction giving only `y`, `bandwidth`, `flabels`, `fsmoothing`
# and `tsmoothing`. `func` is set to `identity`.
TFAmplitude(y, bandwidth, flabels, fsmoothing, tsmoothing) =
    TFAmplitude(y, bandwidth, flabels, fsmoothing, tsmoothing, identity)
#
# As above, but setting by default also both
# `fsmoothing` and `tsmoothing` to `noSmoother`.
TFAmplitude(y, bandwidth, flabels) =
    TFAmplitude(y, bandwidth, flabels, noSmoother, noSmoother, identity)

#############################################################################
# Constructors from Data
##############################################
"""
```
(1)
function TFamplitude(Z::TFAnalyticSignal;
                func::Function=identity)

(2)
function TFamplitude(𝐙::TFAnalyticSignalVector;
                func::Function=identity)

(3)
function TFamplitude(x         :: Vector{T},
                     sr        :: Int,
                     wl        :: Int,
                     bandwidth :: IntOrReal = 2;
                func       :: Function     = identity,
                filtkind   :: FilterDesign = Butterworth(2),
                fmin       :: IntOrReal    = bandwidth,
                fmax       :: IntOrReal    = sr÷2,
                fsmoothing :: Smoother     = noSmoother,
                tsmoothing :: Smoother     = noSmoother,
                planner    :: Planner      = getplanner,
                ⏩        :: Bool         = true) where T<:Real =

(4)
function TFamplitude(𝐱 :: Vector{Vector{T}},
                < same arguments as method (3) >
```

(1)

Construct a [TFAmplitude](@ref) object computing the amplitude of
[TFAnalyticSignal](@ref) object `Z`.
Optional keyword argument `func` is a function to be applied
element-wise to the data matrix of the output. By default,
the `identity` (do nothing) function is applied.

(2)

Construct a [TFAmplitudeVector](@ref) object from a
[TFAnalyticSignalVector](@ref) object executing method (1)
for all [TFAnalyticSignal](@ref) objects in `𝐙`

(3)

Call [`TFanalyticsignal`](@ref) to obtain the time-frequency
analytic signal of real signal vector `x` and construct a [TFAmplitude](@ref)
object holding the time-frequency amplitude (the mudulus, often
referred to as the *envelope*) of `x`.

All arguments are used for regulating the estimation of the analytic signal,
with the exception of `func`, `fsmoothing` and `fsmoothing`:

`func` is an optional function to be applied to the amplitude data matrix output.

In order to estimate the analytic signal in the time-frequency domain
this function calls the [`TFanalyticsignal`](@ref) constructor
(method (1) therein), with both `fsmoothing` and `tsmoothing` arguments
set to `noSmoother`.
Arguments `fsmoothing` and `fsmoothing` are then used to smooth the amplitude.

In order to obtain amplitude estimations on smoothed analytic signal instead,
create a [TFAnalyticSignal](@ref) object passing a
[Smoother](@ref) to the [`TFanalyticsignal`](@ref)
constructor and then use method (1) to obtain the amplitude.
Such amplitude estimation can be further smoothed
using the [`smooth`](@ref) function, as shown in the examples.

For the meaning of all other arguments, which are passed to function
[`TFanalyticsignal`](@ref), see the documentation therein.

(4)

Construct a [TFAmplitudeVector](@ref) object from a
vector of real signal vectors `𝐱`, executing method (3) for all of them.
In order to estimate the time-frequency analytic signal for a
vector of signals, method (2) of [`TFanalyticsignal`](@ref) is called.

**See**: [`TFanalyticsignal`](@ref), [TFAmplitude](@ref).

**See also**: [`amplitude`](@ref).

**Examples**: see the examples of [`TFanalyticsignal`](@ref).
"""
TFamplitude(Z::TFAnalyticSignal; func::Function=identity) =
    TFAmplitude(amplitude(Z.y; func=func), Z.bandwidth, Z.flabels, Z.fsmoothing, Z.tsmoothing, func)

TFamplitude(𝐙::TFAnalyticSignalVector; func::Function=identity) =
    TFAmplitudeVector([TFAmplitude(amplitude(Z.y; func=func), Z.bandwidth, Z.flabels, Z.fsmoothing, Z.tsmoothing, func) for Z ∈ 𝐙])

TFamplitude(𝐱         :: Union{Vector{T}, Vector{Vector{T}}},
            sr        :: Int,
            wl        :: Int,
            bandwidth :: IntOrReal    = 2;
        func          :: Function     = identity,
        filtkind      :: FilterDesign = Butterworth(2),
        fmin          :: IntOrReal    = bandwidth,
        fmax          :: IntOrReal    = sr÷2,
        fsmoothing    :: Smoother     = noSmoother,
        tsmoothing    :: Smoother     = noSmoother,
        planner       :: Planner      = getplanner,
        ⏩           :: Bool         = true) where T<:Real =

    smooth(fsmoothing,
           tsmoothing,
           TFamplitude(TFanalyticsignal(𝐱,
                                        sr,
                                        wl,
                                        bandwidth;
                                    filtkind  = filtkind,
                                    fmin      = fmin,
                                    fmax      = fmax,
                                    planner   = planner,
                                    ⏩       = ⏩);
                        func=func),
           "TFamplitude")

# ++++++++++++++++++++  Show override  +++++++++++++++++++ # (REPL output)
function Base.show(io::IO, ::MIME{Symbol("text/plain")}, Y::TFAmplitude)
   r=size(Y.y, 1)
   c=size(Y.y, 2)
   l=length(Y.flabels)
   println(io, titleFont, "▤ TFAmplitude type; $r freq. x $c samples")
   #println(io, "□  □    □      □        □           □", defaultFont)
   println(io, separatorFont, "⭒  ⭒    ⭒      ⭒        ⭒           ⭒", defaultFont)
   println(io, "freq. sm.(.fsmoothing): ", string(Y.fsmoothing))
   println(io, "time  sm.(.tsmoothing): ", string(Y.tsmoothing))
   println(io, "bandwidth             : $(Y.bandwidth) Hz")
   println(io, "function     (.func)  : ", string(Y.func))
   println(io, "freq. lab.(.flabels)  : $(l)-", typeof(Y.flabels))
   println(io, "data            (.y)  : $(r)x$(c)-", typeof(Y.y))
   r≠l && @warn "number of frequency labels does not match the data matrix size" l r
end


function Base.show(io::IO, ::MIME{Symbol("text/plain")}, 𝐘::TFAmplitudeVector)
   println(io, titleFont, "▤ ⋯ ▤ TFAmplitudeVector Type")
   println(io, separatorFont, "⭒  ⭒    ⭒      ⭒        ⭒           ⭒", defaultFont)
   println(io, "$(length(𝐘))-element Vector{TFAmplitude}")
end
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #


#### Time-Frequency Phase ####

#############################################################################
# Generic constructors of TFAmplitude objects:
##############################################
#
# Enable construction giving only `y`, `bandwidth`, `flabels`, `fsmoothing`
# and `tsmoothing`. `unwrapped` is set to false and `func` is set to `identity`.
TFPhase(y, bandwidth, flabels, nonlinear, fsmoothing, tsmoothing) =
    TFPhase(y, bandwidth, flabels, nonlinear, fsmoothing, tsmoothing, false, identity)
#
# As above, but setting by default also `nonlinear` to true and both
# `fsmoothing` and `tsmoothing` to `noSmoother`.
TFPhase(y, bandwidth, flabels) =
    TFPhase(y, bandwidth, flabels, true, noSmoother, noSmoother, false, identity)

#############################################################################
# Constructors from Data
##############################################
"""
```
(1)
function TFphase(Z :: TFAnalyticSignal;
            func      :: Function = identity,
            unwrapped :: Bool     = false)

(2)
function TFphase(𝐙 :: TFAnalyticSignalVector;
            func      ::Function = identity,
            unwrapped ::Bool     = false)

(3)
function TFphase(x         :: Vector{T},
                 sr        :: Int,
                 wl        :: Int,
                 bandwidth :: IntOrReal = 2;
           unwrapped  :: Bool         = false,
           func       :: Function     = identity,
           filtkind   :: FilterDesign = Butterworth(2),
           fmin       :: IntOrReal    = bandwidth,
           fmax       :: IntOrReal    = sr÷2,
           nonlinear  :: Bool         = false,
           fsmoothing :: Smoother     = noSmoother,
           tsmoothing :: Smoother     = noSmoother,
           planner    :: Planner      = getplanner,
           ⏩        :: Bool         = true) where T<:Real

(4)
function TFphase(𝐱 :: Vector{Vector{T}},
                < same arguments as method (3) >
```

(1)

Construct a [TFPhase](@ref) object computing the phase of
[TFAnalyticSignal](@ref) object `Z`. By default the phase is
represented in ``[−π, π]``.

If optional keyword argument `unwrapped` is true (false by defalut),
the phase is unwrapped, that is, it holds the cumulative sum of the phase
along the time dimension once this is represented in ``[0, 2π]``.

Optional keyword argument `func` is a function to be applied
element-wise to the data matrix of the output. By default,
the `identity` (do nothing) function is applied. If `unwrapped`
is true, the function is applied on the unwrapped phase.

(2)

Construct a [TFPhaseVector](@ref) object from a
[TFAnalyticSignalVector](@ref) object executing method (1)
for all [TFAnalyticSignal](@ref) objects in `𝐙`

(3)

Call [`TFanalyticsignal`](@ref) to obtain the time-frequency
analytic signal of real signal vector `x` and construct a [TFPhase](@ref)
object holding the time-frequency phase (argument) of `x`.

All arguments are used for regulating the estimation of the analytic signal,
with the exception of `unwrapped`, `func`, `fsmoothing` and `fsmoothing`.

`unwrapped` has the same meaning as in method (1) and (2).

`func` is an optional function to be applied to the phase data matrix output.
If `unwrapped` is true, it is applied to the unwrapped phase.

In order to estimate the analytic signal in the time-frequency domain
this function calls the [`TFanalyticsignal`](@ref) constructor
(method (1) therein),
with both `fsmoothing` and `tsmoothing` arguments set to `noSmoother`.
`fsmoothing` and `fsmoothing` are then used to smooth the phase if `unwrapped`
is true.

In order to obtain phase estimations on smoothed analytic signal instead,
create a [TFAnalyticSignal](@ref) object passing a
[Smoother](@ref) to the [`TFanalyticsignal`](@ref)
constructor and then use method (1) to obtain the phase.

For the meaning of all other arguments, which are passed to function
[`TFanalyticsignal`](@ref), see the documentation therein.

(4)

Construct a [TFPhaseVector](@ref) object from a
vector of real signal vectors `𝐱`, executing method (3) for all of them.
In order to estimate the time-frequency analytic signal for a
vector of signals, method (2) of [`TFanalyticsignal`](@ref) is called.

**See**: [`TFanalyticsignal`](@ref), [TFPhase](@ref).

**See also**: [`phase`](@ref), [`unwrapPhase`](@ref), [`polar`](@ref).

**Examples**: see the examples of [`TFanalyticsignal`](@ref).
"""
TFphase(Z::TFAnalyticSignal; func::Function=identity, unwrapped::Bool=false) =
    TFPhase(phase(Z.y; unwrapdims=2*unwrapped, func=func), Z.bandwidth, Z.flabels, Z.nonlinear, Z.fsmoothing, Z.tsmoothing, unwrapped, func)

TFphase(𝐙::TFAnalyticSignalVector; func::Function=identity, unwrapped::Bool=false) =
    TFPhaseVector([TFphase(Z; func=func, unwrapped=unwrapped) for Z in 𝐙])

function TFphase(𝐱         :: Union{Vector{T}, Vector{Vector{T}}},
                 sr        :: Int,
                 wl        :: Int,
                 bandwidth :: IntOrReal    = 2;
           unwrapped  :: Bool         = false,
           func       :: Function     = identity,
           filtkind   :: FilterDesign = Butterworth(2),
           fmin       :: IntOrReal    = bandwidth,
           fmax       :: IntOrReal    = sr÷2,
           nonlinear  :: Bool         = false,
           fsmoothing :: Smoother     = noSmoother,
           tsmoothing :: Smoother     = noSmoother,
           planner    :: Planner      = getplanner,
           ⏩        :: Bool         = true) where T<:Real

   TFP = TFphase(TFanalyticsignal(𝐱,
                                 sr,
                                 wl,
                                 bandwidth;
                                 filtkind  = filtkind,
                                 fmin      = fmin,
                                 fmax      = fmax,
                                 nonlinear = nonlinear,
                                 planner   = planner,
                                 ⏩       = ⏩);
                         unwrapped = unwrapped,
                         func      = func)

   !unwrapped && (fsmoothing≠noSmoother || tsmoothing≠noSmoother) &&
   @warn 📌*", call to TFPhase constructor; smoothing is meaningful only on unwrapped phase. Smoothing has been disabled.)" unwrapped

   return unwrapped ? smooth(fsmoothing,
                             tsmoothing,
                             TFP,
                             "TFphase") :
                      TFP
end

# ++++++++++++++++++++  Show override  +++++++++++++++++++ # (REPL output)
function Base.show(io::IO, ::MIME{Symbol("text/plain")}, Y::TFPhase)
    r=size(Y.y, 1)
    c=size(Y.y, 2)
    l=length(Y.flabels)
    println(io, titleFont, "▤ TFPhase type; $r freq. x $c samples")
    #println(io, "□  □    □      □        □           □", defaultFont)
    println(io, separatorFont, "⭒  ⭒    ⭒      ⭒        ⭒           ⭒", defaultFont)
    println(io, "non-linear(.nonlinear): $(Y.nonlinear)")
    println(io, "freq. sm.(.fsmoothing): ", string(Y.fsmoothing))
    println(io, "time  sm.(.tsmoothing): ", string(Y.tsmoothing))
    println(io, "bandwidth             : $(Y.bandwidth) Hz")
    println(io, "unwrapped             : $(Y.unwrapped)")
    println(io, "function     (.func)  : ", string(Y.func))
    println(io, "freq. lab.(.flabels)  : $(l)-", typeof(Y.flabels))
    println(io, "data            (.y)  : $(r)x$(c)-", typeof(Y.y))
    r≠l && @warn "number of frequency labels does not match the data matrix size" l r
end

function Base.show(io::IO, ::MIME{Symbol("text/plain")}, 𝐘::TFPhaseVector)
    println(io, titleFont, "▤ ⋯ ▤ TFPhaseVector Type")
    println(io, separatorFont, "⭒  ⭒    ⭒      ⭒        ⭒           ⭒", defaultFont)
    println(io, "$(length(𝐘))-element Vector{TFPhase}")
end
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
