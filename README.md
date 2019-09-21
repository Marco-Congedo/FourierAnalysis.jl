# FourierAnalysis

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://Marco-Congedo.github.io/FourierAnalysis.jl/latest)

**FourierAnalysis** is a [**Julia**](https://julialang.org/) package for
performing the analysis of *real multivariate data* (e.g., time series)
in the *frequency* domain and in the *time-frequency* domain.

In the frequency domain *FourierAnalysis* computes **spectra**, *linear* and
*non-linear* **cross-spectral matrices** and several *linear* and *non-linear* **coherence matrices** using the sliding-windows [Welch method](https://en.wikipedia.org/wiki/Welch%27s_method).

Time-frequency representations are obtained applying a
[filter-bank](https://en.wikipedia.org/wiki/Filter_bank) and the
[Hilber transform](https://en.wikipedia.org/wiki/Hilbert_transform).
This way *FourierAnalysis* computes the **analytic signal**, from which the **instantaneous amplitude** (envelope) and **instantaneous phase** are obtained, along with several popular *linear* and *non-linear*, *weighted*, *univariate* and *bivariate* statistics, such as the **mean amplitude**, **mean direction**, **phase concentration**, the non-linear version of which is a directional statistic known as **circular mean resultant length**, **comodulation** and **coherence**, the non-linear version of which is a synchronization statistic known as **phase-locking values** or **phase coherence**.

Such large panel of measures is provided in a simple and unified fashion,
following the conceptual approach illustrated in
[Congedo(2018)](https://hal.archives-ouvertes.fr/hal-01868538/document)
in the context of electroencephalography (EEG), for which all default settings
have been tailored. The package has been written with the "do-it-with-one-line" spirit and with the aim of allowing full control over relevant options.

*FourierAnalysis* is based on packages [FFTW](https://github.com/JuliaMath/FFTW.jl),
[AbstractFFTs](https://github.com/JuliaMath/AbstractFFTs.jl) and
[DSP](https://github.com/JuliaDSP/DSP.jl), providing a simple interface to the parts of them it uses.

## Installation

The package is till under preliminary testing and is not registered.
Execute the following command in Julia's REPL:

    ]add https://github.com/Marco-Congedo/FourierAnalysis

To obtain the latest development version execute instead

    ]add https://github.com/Marco-Congedo/FourierAnalysis#master

## Disclaimer

This package is still in a preliminary stage.
It needs throughout testing.
Independent reviewers are more then welcome.

## About the Author

[Marco Congedo](https://sites.google.com/site/marcocongedo) is
a research scientist of CNRS (Centre National de la Recherche Scientifique), working in Grenoble, France. Contact: first name *dot* last name at gmail *dot* com
