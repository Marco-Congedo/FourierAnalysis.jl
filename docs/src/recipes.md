# recipes.jl

This unit implements plot recipes for plotting data produced by
**Fourier Analysis**.

In the following recipes, operator `-->` sets the default plot attribute,
while operator `:=` forces the plot attribute to a specific value.

!!! note "Nota Bene"
    Working with plot recipes *you must* put the `;` symbol in your
    argument list to separate optional keyword arguments.

## plot tapering windows

```
@recipe function f(::Type{Taper}, taper::Taper)
    title               --> taperinfo(taper)
    background_color    --> :white
    gridcolor           --> :white
    legend              --> false
    dpi                 --> 300
    taper.y
end
```
Plot the data of a [Taper](@ref) object, that is, the tapering
window(s) hold in its `.y` field.
This is a single tapering window for all kind of tapering windows,
but for Slepian's discrete prolate spheroidal sequences,
which are multi-tapering.

The recipes allows to pass any other optional keyword argument
of the standard Julia plot package.

**Examples**:

```
using FourierAnalysis, Plots

plot(taper(parzen, 256))
plot!(taper(harris4, 256))

# discrete prolate spheroidal sequences
H=taper(slepian, 256, α=4, n=7)
plot(H)
```

## plot spectra

```
@recipe function f(::Type{Spectra}, S::Spectra;
            fmax   = S.sr/2,
            xspace = 4,
            ytitle = "Power (\\muV²)")

    b=f2b(fmax, S.sr, S.wl)
    xticstring=["$(xspace*i)" for i=1:fmax÷xspace]
    pushfirst!(xticstring, "$(fres(S.sr, S.wl))")
    xticrange=[1:S.wl/S.sr*xspace:b+1;]
    for i=2:length(xticrange) xticrange[i]=xticrange[i]-1 end

    legend         --> false
    dpi            --> 300
    xlabel         --> "Frequency (Hz)"
    xticks         --> (xticrange, xticstring)
    ylabel         --> ytitle
    delete!(plotattributes, :fmax)
    delete!(plotattributes, :xspace)
    delete!(plotattributes, :ytitle)
    S.y[1:(f2b(fmax, S.sr, S.wl; DC=S.DC)), :]
end
```

Plot the data in a [Spectra](@ref) object,
*e.g.*, plot all its spectra with separate lines.

Optional keyword argument `fmax` limits the plot on the right
to the given frequency (in Hz).

Optional keyword argument `xspace` (in Hz) determines the spacing between frequency ticks (x-axis).

Optional keyword argument `ytitle` specifies a title for the y axis. The default is "Power (\\muV²)"

The function allows to pass any other optional keyword argument of the standard Julia plot function.

Useful attributes for plotting spectra are:

- ```xtickfont      = font(11, "Times"))```
- ```ytickfont      = font(11, "Times"))```
- ```left_margin    = 2mm```
- ```bottom_margin  = 2mm```
- ```xaxis          = ("My x axis title", font(11, "Courier New"))```
- ```yaxis          = ("My y axis title", font(11, "Courier New"))```
- ```linecolor      = :grey```

!!! note "Nota Bene"
    In order to use the `mm` measure for margins you need to add
    `using Plots.Measures` in your script.

**Examples**:

```
using FourierAnalysis, Plots, Plots.Measures

sr, t, f, a = 128, 128, 16, 0.5
# create a sinusoidal at 16Hz superimposed to white noise
v=sinusoidal(a, f, sr, t*16, 0) + randn(t*16)
# create a data matrix
X=broadcast(+, v, randn(t*16, 3)) * randn(3, 3)

# compute spectra using hamming tapering window
S=spectra(X, sr, t; tapering=hann)

# plot spectra

plot(S) # quick and dirt

# gather useful attributes to make nicer plots
spectraArgs=(fmax=32,
             xspace=4,
             left_margin = 2mm,
             bottom_margin = 2mm,
             xtickfont = font(11, "Times"),
             ytickfont = font(11, "Times"))

plot(S; spectraArgs...)

# compute and plot amplitude spectra
S=spectra(X, sr, t; tapering=riesz, func=√)
plot(S; ytitle="Amplitude (\\muV)", spectraArgs...)

```

## plot time-frequency objects

Plots of data held by [TFAnalyticSignal](@ref), [TFAmplitude](@ref)
and [TFPhase](@ref) objects can be obtained with a simple call to the `heatmap` function of the
[Plots.jl](https://github.com/JuliaPlots/Plots.jl) package.
The [`tfAxes`](@ref) function provides automatic labeling of the axes
and is passed as the first argument.

The second argument is the data matrix to be plotted. The matrix held in [TFAmplitude](@ref)
and [TFPhase](@ref) objects is real, therefore it can be passed directly as argument. The matrix held in [TFAnalyticSignal](@ref) objects holds the analytic signal (AS), therefore it is complex.
A function should be specified so as to transform it into a real matrix. Typical functions are:

- `real`: return the real part of the AS
- `imag`: return the imaginary part of the AS
- `amplitude`: return the AS [`amplitude`](@ref)
- `phase`: return the AS [`phase`](@ref)

**Note**:

Appropriate colors for plotting the real and imaginary part of the analytic signal, as well as the phase, are `:pu_or`, `:bluesreds`.

Appropriate colors for plotting amplitude and unwraped phase are `:amp`, `:fire`, `:dimgray`, `:gwv`.

Useful attributes for plotting time-frequency objects are:
- ```xtickfont      = font(10, "Times"))```
- ```ytickfont      = font(10, "Times"))```
- ```right_margin   = 2mm```
- ```top_margin     = 2mm```

**Examples**:
```
using FourierAnalysis, Plots, Plots.PlotMeasures
sr, t, f, a, s = 128, 128, 8, 2.5, sinusoidal
# create a sinusoidal at 8Hz modulated by a sinusoidal at 0.125Hz
# and superimpose white noise
v=s(a, f, sr, t*4, 0) .* s(a, 0.125, sr, t*4, 0) + randn(t*4)
plot(v)

# generate a times series with two frequencies
# and a clear envelope
sr, t, bandwidth=128, 512, 2
h=taper(harris4, t)
x1=sinusoidal(10, 8, sr, t, 0)
x2=sinusoidal(10, 19, sr, t, 0)
v=Vector((x1+x2).*h.y+randn(t))
plot(v)

# gather useful attributes to obtain nice heatmpap plots
tfArgs=(right_margin = 2mm,
        top_margin = 2mm,
        xtickfont = font(10, "Times"),
        ytickfont = font(10, "Times"))

# compute the analytic signal (AS) of vector v
Y=TFanalyticsignal(v, sr, t; fmax=32)

# plot the real part of the AS
heatmap(tfAxes(Y)..., real(Y.y); c=:pu_or, tfArgs...)

# ...the imaginary part of the AS
heatmap(tfAxes(Y)..., imag(Y.y); c=:bluesreds, tfArgs...)

# ...the amplitude of the AS
heatmap(tfAxes(Y)..., amplitude(Y.y); c=:amp, tfArgs...)

# ...the amplitude of the AS smoothed in the freq. dim.
heatmap(tfAxes(Y)...,
        amplitude(smooth(hannSmoother, noSmoother, Y).y);
        c=:amp, tfArgs...)

# ...the amplitude of the AS smoothed in freq. and time
heatmap(tfAxes(Y)...,
        amplitude(smooth(hannSmoother, hannSmoother, Y).y);
        c=:amp, tfArgs...)

# ...the phase of the AS
heatmap(tfAxes(Y)..., phase(Y.y); c=:bluesreds, tfArgs...)

# ...the phase of the AS weighted by the amplitude
heatmap(tfAxes(Y)..., phase(Y.y).*amplitude(Y.y);
        c=:bluesreds, tfArgs...)

# compute a TF Amplitude obejct
A=TFamplitude(v, sr, t; fmax=32)

# plot the amplitude
heatmap(tfAxes(A)..., A.y; c=:amp, tfArgs...)

# compute a TF Phase obejct
ϴ=TFphase(v, sr, t; fmax=32)

# plot the phase
heatmap(tfAxes(ϴ)..., ϴ.y; c=:pu_or, tfArgs...)

# compute a TF Phase obejct with the phase unwrapped
uwϴ=TFphase(v, sr, t; unwrapped=true)

# plot the unwrapped phase
heatmap(tfAxes(uwϴ)..., uwϴ.y; c=:fire, tfArgs...)
```

```@docs
tfAxes
```
