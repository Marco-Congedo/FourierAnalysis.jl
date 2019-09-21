# FourierAnalysis Documentation

## Requirements

**Julia** version ≥ 1.1.1

**Packages:**
[FFTW](https://github.com/JuliaMath/FFTW.jl),
[DSP](https://github.com/JuliaDSP/DSP.jl),
[Plots](https://github.com/JuliaPlots/Plots.jl).

## Installation

Execute the following command in Julia's REPL:

    ]add FourierAnalysis

To obtain the latest development version execute instead

    ]add FourierAnalysis#master

## Disclaimer

This package is still in a preliminary stage.
It needs throughout testing.
Any independent reviewer is welcome.

## About the Author

[Marco Congedo](https://sites.google.com/site/marcocongedo) is
a research scientist of CNRS (Centre National de la Recherche Scientifique), working in Grenoble, France.

## Overview

*FourierAnalysis* allows the analysis of
- sets of *real multivariate time-series* in the *frequency* domain
- sets of *real time-series* in the *time-frequency* domain.

In the frequency domain *FourierAnalysis* computes **spectra**, linear and
non-linear **cross-spectral matrices** and several linear and non-linear **coherence matrices** using the sliding windows Welch method.

For those estimations, keep in mind that given sampling rate `sr`
and window length ``wl``, the **discrete Fourier frequencies** are
``0, 1r, 2r,..., qr``, where ``r=sr/wl`` is the **frequency resolution**
and ``q=wl÷2`` (integer division).
The 0 (zero) frequency corresponds to the DC level.

Time-frequency (TF) representations are obtained applying a *pass-band filter-bank* and the *Hilber transform*. This way *FourierAnalysis* computes the **analytic signal**, from which the **instantaneous amplitude** (envelope) and **instantaneous phase** are obtained, along with several popular *linear* and *non-linear*, *weighted*, *univariate* and *bivariate* statistics, such as the **mean amplitude**, **mean direction**, **phase concentration** (the non-linear version is a directional statistic known as **circular mean resultant length**), **amplitude co-modulation**, **coherence** (the non-linear version is a synchronization statistic known as **phase-locking values** or **phase coherence**), etc.

A large panel of measures are provided in a simple and unified fashion,
following the approach illustrated in
[Congedo(2018)](https://hal.archives-ouvertes.fr/hal-01868538/document)
in the context of electroencephalography (EEG), for which all default settings
have been tailored. The package has been written with the "do-it-with-one-line" spirit
and with the aim of allowing full control over relevat options for the
Fourier analysis of multivariate time-series.

*FourierAnalysis* is based on packages [FFTW](https://github.com/JuliaMath/FFTW.jl),
[AbstractFFTs](https://github.com/JuliaMath/AbstractFFTs.jl) and
[DSP](https://github.com/JuliaDSP/DSP.jl) providing a simple interface to the
parts of them it uses.

For starting using this package, browse the code units listed here below and
execute the many **code examples** you will find therein or execute
the example .jl units collected in the "example" folder distributed
in the github repository.

## Code units

*FourierAnalysis* includes fourteen code units (.jl files):

| Main API Units   | Description |
|:----------|:----------|
| [MainModule](@ref) | (FourierAnalysis.jl) constants, types, structs, aliases, tips & tricks |
| [tapers.jl](@ref) | tapering windows for spectral, cross-spectral and coherence analysis |
| [spectra.jl](@ref) | spectra of a series set or of a multivariate time series set |
| [crossspectra.jl](@ref) | cross-spectral matrices of a multivariate time series set |
| [coherence.jl](@ref) | coherence matrices of a multivariate time series set |
| [timefrequency.jl](@ref) | analytic signal, instantaneous amplitude and phase of a series set |
| [timefrequencyuni.jl](@ref) | univariate measures of a series set |
| [timefrequencybi.jl](@ref) | bivariate measures of a series set |
| [plots.jl](@ref) | simple interface to Plots.jl to plot spectra and time-frequency representations |
| [tools.jl](@ref) | collecton of useful functions |

| Other API Units  | Description |
|:----------|:----------|
| [fftw.jl](@ref) | simple interface to FFTW.jl to obtain forward and backward FFT planners |
| [goertzel.jl](@ref) | Goertzel's algorithm for estimating a single DFT coefficient |
| [filters.jl](@ref) | simple interface to DSP.jl to apply band-pass filter-banks |
| [hilbert.jl](@ref) | computation of the analytic signal via Hilbert transform |


## Contents

```@contents
Pages = [       "index.md",
                "MainModule.md",
                "tapers.md",
                "spectra.md",
                "crossspectra.md",
                "coherence.md",
                "timefrequency.md",
                "timefrequencyuni.md",
                "timefrequencybi.md",
                "plots.md",
                "tools.md",
                "fftw.md",
                "goertzel.md",
                "filters.md",
                "hilbert.md"]
Depth = 1
```

## Index

```@index
```
