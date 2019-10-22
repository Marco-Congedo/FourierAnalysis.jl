# FourierAnalysis Documentation

## Requirements & Installation

**Julia** version ≥ 1.1.1

**Packages:** see the [dependencies](@ref) of the main module.

The package is still not registered. To install it
execute the following command in Julia's REPL:

    ]add https://github.com/Marco-Congedo/FourierAnalysis.jl

## Disclaimer

This package is still in a preliminary stage.
It needs throughout testing.
Any independent reviewer for both the code and the documentation is welcome.

## About the Author

[Marco Congedo](https://sites.google.com/site/marcocongedo) is
a research scientist of [CNRS](http://www.cnrs.fr/en) (Centre National de la Recherche Scientifique), working in [UGA](https://www.univ-grenoble-alpes.fr/english/) (University of Grenoble Alpes), in Grenoble, France, the city where [Jean-Baptiste Joseph Fourier](https://en.wikipedia.org/wiki/Joseph_Fourier) has served as a Governor. **Contact**: first name dot last name at gmail dot com

## Overview

*FourierAnalysis* allows the analysis of
- sets of *real multivariate time-series* in the *frequency* domain,
- sets of *real univariate time-series* in the *time-frequency* domain.

The implemented tools may suit also other kind of data.

Frequency-domain representations include **spectra**, linear and
non-linear **cross-spectral matrices** and several linear and non-linear **coherence matrices**, which are all estimated using the sliding-windows (Welch) method.

For those estimations, keep in mind that given sampling rate ``sr``
and window length ``wl``, the **discrete Fourier frequencies** are
``0, 1r, 2r,..., qr``, where ``r=sr/wl`` is the **frequency resolution**
and ``q=wl÷2`` (integer division).
The 0 (zero) frequency corresponds to the **DC level**, which estimation in *FourierAnalysis* is always given as an option.

Time-frequency (TF) representations are obtained applying a
[filter-bank](https://en.wikipedia.org/wiki/Filter_bank) and the
[Hilbert  transform](https://en.wikipedia.org/wiki/Hilbert_transform). This way *FourierAnalysis* computes the **analytic signal**, from which the **instantaneous amplitude** (envelope) and **instantaneous phase** are obtained, along with several popular *linear* and *non-linear*, *weighted*, *univariate* and *bivariate* statistics, such as the **mean amplitude**, **mean direction**, **phase concentration**, the non-linear version of which is a directional statistic known as **circular mean resultant length**, **amplitude co-modulation**, **coherence**, the non-linear version of which is a synchronization measure known as **phase-locking values** or **phase coherence**, etc.

Such a large panel of measures is provided in a simple and unified fashion,
following the approach illustrated in
[Congedo(2018)](https://hal.archives-ouvertes.fr/hal-01868538v2/document)
in the context of electroencephalography (EEG), for which all default settings have been tailored. The package has been written with the *do-it-with-one-line*
spirit and with the aim of allowing full control over relevant options for the Fourier analysis of multivariate time-series.

For starting using this package, browse the code units listed here below and
execute the many **code examples** you will find therein or execute
the 'example.jl' units collected in the "example" folder distributed
in the github repository.

## Code units

*FourierAnalysis* includes fourteen code units (.jl files):

| Main API Units   | Description |
|:----------|:----------|
| [MainModule](@ref) | (FourierAnalysis.jl) constants, types, some structs, tips & tricks |
| [tapers.jl](@ref) | tapering windows for spectral, cross-spectral and coherence analysis |
| [spectra.jl](@ref) | spectra of a time-series (set) or of a multivariate time series (set) |
| [crossspectra.jl](@ref) | cross-spectral matrices of a multivariate time series (set) |
| [coherence.jl](@ref) | coherence matrices of a multivariate time series (set) |
| [goertzel.jl](@ref) | Goertzel's algorithms for estimating a single DFT coefficient |
| [timefrequency.jl](@ref) | analytic signal, instantaneous amplitude and phase of a time-series set |
| [timefrequencyuni.jl](@ref) | univariate measures of a time-series set |
| [timefrequencybi.jl](@ref) | bivariate measures of a time-series set |
| [recipes.jl](@ref) | plot recipes and tips to plot data created by *Fourier Analysis* |
| [tools.jl](@ref) | collection of useful functions |

| Other API Units  | Description |
|:----------|:----------|
| [fftw.jl](@ref) | interface to FFTW.jl to obtain forward and backward FFT planners |
| [filters.jl](@ref) | interface to DSP.jl to apply filter-banks |
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
                "recipes.md",
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
