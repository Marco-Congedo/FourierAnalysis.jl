#   Unit "plots" of the FourierAnalysis Package for julia language
#   v 0.0.1 - last update 7th of August 2019
#
#   MIT License
#   Copyright (c) 2019, Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENTS :
#   This unit implements utilities fr plotting data produced by
#   Fourier Analysis

########################################################

"""
```
plot( h      :: Taper,
    showkind :: Bool = false;
   args...)
```

Plot the data of a [Taper](@ref) object, that is, the tapering
window(s) hold in its `.y` field.
This is a single tapering window for all kind of tapering windows,
but for Slepian's discrete prolate spheroid sequences,
which are multi-tapering.

If optional keyword argument `showkind` is true, the output of the
[`taperinfo`](@ref) function is printed as title of the plot.
By defaut `showkind` is false.

The function allows to pass any other optional keyword argument
of the standard Julia plot function.

**Examples**:
```
using FourierAnalysis, Plots
sr, t = 128, 128
H=taper(hamming, t)

# plot tapering window
plot(H; title="my own title") # title is a Plot.jl argument
plot(H; showkind=true)

# plot Slepian's dpss
plot(slepians(sr, t*4, 2); showkind=true)
plot(slepians(sr, t*16, 1); title="my own title")
```
"""
plot( h      :: Taper;
    showkind :: Bool = false,
   args...)  =
  plot(h.y;
    title=showkind ? taperinfo(h) : "", gridcolor=:white, labels="", args...)


"""
```
plot(S      :: Spectra;
    maxf    :: Union{Real, Int} = 0,
    space   :: Int              = 4,
    ylabel  :: String = "Power (\\muV²)",
   args...)
```
Produce the line plot of the data in a [Spectra](@ref) object,
that is, plot all its spectra with separate lines.

Optional keyword argument `maxf` limits the plot on the right
to the given frequency (in Hz).

Optional keyword argument `space` (in Hz) determines the spacing between
frequency ticks and labels (x-axis).

The function allows to pass any other optional keyword argument
of the standard Julia plot function.

**See also**: [`bar`](@ref)

**Examples**:
```
using FourierAnalysis, Plots
sr, t, f, a = 128, 128, 16, 0.5
# create a sinusoidal at 16Hz superimposed to white noise
v=sinusoidal(a, f, sr, t*16, 0) + randn(t*16)
# create a data matrix
X=broadcast(+, v, randn(t*16, 3)) * randn(3, 3)
# compute spectra using hamming tapering window
H=taper(hamming, t)
S=spectra(X, sr, t; tapering=H)

# plot spectra
plot(S)
plot(S; maxf=32)
plot(S; maxf=32, space=2)

# compute and plot amplitude spectra
S=spectra(X, sr, t; tapering=H, func=sqrt)
plot(S; maxf=32, space=2, ylabel="Amplitude (\\muV)")

```
"""
plot(S      :: Spectra;
    maxf    :: Union{Real, Int} = 0,
    space   :: Int              = 4,
    ylabel  :: String = "Power (\\muV²)",
   args...) =
  plot(S.y[1:(maxf==0 ? f2b(S.sr/2, S.sr, S.wl; DC=S.DC) : f2b(maxf, S.sr, S.wl; DC=S.DC)), :];
     _spectraplotargs(S;
                  maxf   = maxf,
                  space  = space,
                  ylabel = ylabel)...,
     args...)


 """
 ```
 bar(S      :: Spectra;
     maxf    :: Union{Real, Int} = 0,
     space   :: Int              = 4,
     ylabel  :: String = "Power (\\muV²)",
    args...)
```

 Produce the bar plot of the data in a [Spectra](@ref) object,
 that is, plot all its spectra with separate bars.
 Use the same syntax as the [`plot`](@ref) fuction,
 which produces a line plot instead.
 Bar plots for spectra should be used insted of the line plots
 only for univariate spectra.

 **See**: [`plot`](@ref)
 """
bar( S      :: Spectra;
    maxf    :: Union{Real, Int} = 0,
    space   :: Int              = 4,
    ylabel  :: String = "Power (\\muV²)",
   args...) =
  bar(S.y[1:(maxf==0 ? f2b(S.sr/2, S.sr, S.wl; DC=S.DC) : maxf), :];
      _spectraplotargs(S;
                   maxf   = maxf,
                   space  = space,
                   ylabel = ylabel)...,
      args...)



# internal function:
# generate appropriate arguments for plotting spectra using the PLOTS package,
# given a Spectra object and optional keyword argumants max frequency `maxf`,
# `space` and `ylabel`. This function is called by the plot and bar
# functions declared in this unit.
# It needs unpacking (...), for example:
# plot(S, _spectraplotargs(128, 512, 32; space=4)...)
function _spectraplotargs( S  :: Spectra;
                  maxf   :: Union{Real, Int} = 0,
                  space  :: Int              = 4,
                  ylabel :: String = "Power (\\muV²)")
  if maxf==0 maxf=S.sr/2 end
  b=f2b(maxf, S.sr, S.wl)
  xticstring=["$(space*i)" for i=1:maxf÷space]
  pushfirst!(xticstring, "$(fres(S.sr, S.wl))")
  xticrange=[1:S.wl/S.sr*space:b+1;]
  for i=2:length(xticrange) xticrange[i]=xticrange[i]-1 end
  return (labels        = "",
          dpi           = 300, #linecolor=:grey,
          xaxis         = ("Frequency (Hz)", font(11, "Courier New")),
          xticks        = (xticrange, xticstring),
          yaxis         = (ylabel, font(11, "Courier New")),
          left_margin   = 2mm,
          bottom_margin = 2mm)
end



"""
```
function heatmap(Y    :: TFAnalyticSignal,
                 func :: Function;
            c :: Symbol = :default,
            args...)
```

Generate heatmaps from a [TFAnalyticSignal](@ref) object.
Since the analytic signal (AS) is complex, a function `func`
must be specified and must transform the AS data
matrix into a real matrix. Typical functions are:

- `real`: heatmap of the real part of the AS
- `imag`: heatmap of the imaginary part of the AS
- `amplitude`: heatmap of the AS [`amplitude`](@ref)
- `phase`: heatmap of the AS [`phase`](@ref)

Note:

Appropriate colors for plotting analytic signal and phase are `:pu_or`, `:bluesreds`.
Appropriate colors for plotting amplitude are `:amp`, `:fire`, `:dimgray`, `:gwv`.

**Examples**:
```
using FourierAnalysis, Plots
sr, t, f, a = 128, 128, 8, 2
# create a sinusoidal at 8Hz
wave=sinusoidal(a, f, sr, t*4, 0)
# create a modulating sinusoidal at 0.25Hz
modu=sinusoidal(a, 0.25, sr, t*4, 0)
# create modulated wave and superimpose white noise
v=wave .* modu + randn(t*4)
plot(v)

# compute analytic signal (AS) of vector v
Y=TFanalyticsignal(v, sr, t; fmax=32)

# plot the real part of the AS
heatmap(Y, real)

# ...the imaginary part of the AS
heatmap(Y, imag)

# ...the amplitude of the AS
heatmap(Y, amplitude)

# ...the amplitude of the AS smoothed in the freq. dim.
heatmap(smooth(hannSmoother, noSmoother, Y), amplitude)

# ...the amplitude of the AS smoothed in freq. and time
heatmap(smooth(hannSmoother, hannSmoother, Y), amplitude)

# ...the phase of the AS and use custom colors
heatmap(Y, phase; c=:pu_or)
```
"""
function heatmap(Y    :: TFAnalyticSignal,
                 func :: Function;
            c :: Symbol = :amp,
            args...)
   if c==:default
      (func==real || func==imag || func==phase) ? c=:bluesreds :
      c=:amp
   end
   heatmap(_tfaxis(Y)..., func(Y.y), yflip = false, c=c, args...)
end

"""
```
function heatmap(Y :: TFAmplitude;
        c :: Symbol = :default,
        args...) =
```

Generate heatmaps of the data in a [TFAmplitude](@ref) object.

Note: appropriate colors for plotting amplitude are `:amp`, `:fire`, `:dimgray`, `:gwv`.

**Examples**:
```
using FourierAnalysis, Plots
sr, t, f, a = 128, 128, 8, 2
# create a sinusoidal at 8Hz
wave=sinusoidal(a, f, sr, t*4, 0)
# create a modulating sinusoidal at 0.25Hz
modu=sinusoidal(a, 0.25, sr, t*4, 0)
# create modulated wave and superimpose white noise
v=wave .* modu + randn(t*4)
plot(v)

# compute amplitude of vector v
A=TFamplitude(v, sr, t; fmax=32)

heatmap(A)
# or, quicly: heatmap(TFamplitude(v, sr, t; fmax=32))

heatmap(A; c=:fire)
```
"""
heatmap(Y :: TFAmplitude;
  c :: Symbol = :amp,
  args...) =
heatmap(_tfaxis(Y)..., Y.y, yflip = false, c=c, args...)



"""
```
function heatmap(Y :: TFPhase;
              c :: Symbol = :default,
              args...)
```
Generate heatmaps of the data in a [TFPhase](@ref) object.
Default colors for the heatmap are `c=:amp` if the phase object is unwrapped,
`c=:bluesreds` otherwise.

Note:

Appropriate colors for plotting phase: `:pu_or`, `:bluesreds`
Appropriate colors for plotting unwrapped phase are `:amp`, `:fire`, `:dimgray`, `:gwv`.

**See also**: [`isUnwrapped`](@ref)

**Examples**:
```
using FourierAnalysis, Plots
sr, t, f, a = 128, 128, 8, 2
# create a sinusoidal at 8Hz
wave=sinusoidal(a, f, sr, t*4, 0)
# create a modulating sinusoidal at 0.25Hz
modu=sinusoidal(a, 0.25, sr, t*4, 0)
# create modulated wave and superimpose white noise
v=wave .* modu + randn(t*4)
plot(v)

# heatmap of the phase of vector v
heatmap(TFphase(v, sr, t; fmax=32))

# heatmap of the unwrapped phase of vector v
heatmap(unwrapPhase(TFphase(v, sr, t; fmax=32)); c=:fire)
```
"""
function heatmap(Y :: TFPhase;
              c :: Symbol = :default,
              args...)
  c==:default ? (Y.unwrapped ? c=:amp : c=:bluesreds) : nothing
  heatmap(_tfaxis(Y)..., Y.y, yflip = false, c=c, args...)
end



#Generate labels for the axes to be used in heatmap plots of [TFobjects](@ref).
#On the y-axis the labels correspon to the center frequencies of the band-pass
#filter bank used to obtain the time-frequency representation, that is,
#to the .flabels field of the object.

## Examples:
#    using FourierAnalysis, Plots
#    sr, t, f, a = 128, 128, 8, 2
    # create a sinusoidal at 8Hz
#    wave=sinusoidal(a, f, sr, t*4, 0)
    # create a modulating sinusoidal at 0.25Hz
#    modu=sinusoidal(a, 0.25, sr, t*4, 0)
    # create modulated wave and superimpose white noise
#    v=wave .* modu + randn(t*4)

    # compute analytic signal of vector v
#    Y=TFanalyticsignal(v, sr, t; fmax=32)
    # plot amplitude
#    heatmap(_tfaxis(Y)..., amplitude(Y.y), yflip=false, c=:amp)
    # or create amplitude object
#    A=TFamplitude(Y)
    # and plot it
#    heatmap(_tfaxis(A)..., A.y, c=::pu_or)

    # plot phase
#    heatmap(_tfaxis(Y)..., phase(Y.y), c=:bluesreds)
    # plot phase in [0, 2π]
#    heatmap(_tfaxis(Y)..., phase(Y.y, func=x->x+π), c=:amp)
    # unwrapped phase
#    heatmap(_tfaxis(Y)..., phase(Y; unwrapped=true), c=:amp)
    # or create phase object
#    ϴ=TFphase(Y)
    # and plot it
#    heatmap(_tfaxis(ϴ)..., ϴ.y, c=:bluesreds)

# Be careful with phase! Smoothing phase is unappropriate
# since the phase is a discontinous function, however it makes sense to smooth
# phase if the phase is unwrapped.
# Now, this will 'illegally' plot the smoothed phase nonetheless:
## heatmap(tfaxis(Z)..., phase(Z.y), c=:bkr) # bluesreds
# The line here below would not create anything and return an error
# ϴ=smooth(blackmanSmoother, noSmoother, TFphase(Y))
# The line here below instead will smooth the unwrapped phase (in frequencies)
## ϴ=smooth(blackmanSmoother, noSmoother, TFphase(Y, unwrapped=true))
# plot it
## heatmap(tfaxis(ϴ)..., ϴ.y, c=:amp) # bluesreds

#Note:
#appropriate colors for plotting analytic signal and phase: :pu_or, :bluesreds
#appropriate colors for plotting amplitude: :amp, :fire, :dimgray, ;gwv
function _tfaxis(Y::TFobjects)
    function formatf(r::Real)
        s="$r"
        ss=split(s, ".")
        length(ss)>1 && parse(Float64, ss[2])==0.0 ? (return ss[1]) : (return s)
    end
    xs = [string(i) for i = 1:size(Y.y, 2)]
    ys = Vector{String}([formatf(f) for f in Y.flabels])
    return xs, ys
end
