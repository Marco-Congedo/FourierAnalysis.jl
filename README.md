| **Documentation**  | 
|:---------------------------------------:|
| [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://Marco-Congedo.github.io/FourierAnalysis.jl/dev) |

![](/docs/src/assets/Fig1.jpg)

**FourierAnalysis** is a signal-processing [**Julia**](https://julialang.org/) package for
performing the analysis of *real multivariate data* (e.g., multivariate time series)
in the *frequency* domain and in the *time-frequency* domain. It is based upon the
[DSP.jl](https://github.com/JuliaDSP/DSP.jl), [FFTW.jl](https://github.com/JuliaMath/FFTW.jl) and [AbstractFFTs.jl](https://github.com/JuliaMath/AbstractFFTs.jl) packages.

In the frequency domain *FourierAnalysis* computes **spectra**, *linear* and
*non-linear* **cross-spectral matrices** and several *linear* and *non-linear* **coherence matrices** using the sliding-windows [Welch method](https://en.wikipedia.org/wiki/Welch%27s_method).

Time-frequency representations are obtained applying a
[filter-bank](https://en.wikipedia.org/wiki/Filter_bank) and the
[Hilber transform](https://en.wikipedia.org/wiki/Hilbert_transform).
This way *FourierAnalysis* computes the **analytic signal**, from which the **instantaneous amplitude** (envelope) and **instantaneous phase** are obtained, along with several popular *linear* and *non-linear*, *weighted*, *univariate* and *bivariate* statistics, such as
- **mean amplitude** 
- **mean direction** 
- **phase concentration** (the non-linear version of which is the directional statistic *circular mean resultant length*)
- **comodulation**
- **coherence** (the non-linear version of which is known as *phase-locking values* or *phase coherence*).

All these measures are provided in a simple and unified fashion, following the conceptual approach previously illustrated in
in the context of electroencephalography ([Congedo, 2018](https://hal.archives-ouvertes.fr/hal-01868538v2/document)), for which all default settings have been tailored. The package has been written with the *do-it-with-one-line* spirit, but without sacrificing full control over relevant options.

## Installation

The package is not registered yet.
Execute the following command in Julia's REPL:

    ]add https://github.com/Marco-Congedo/FourierAnalysis

## Disclaimer

This package is still in a pre-release stage. It needs throughout testing.
Independent reviewers are more then welcome.

## About the Author

[Marco Congedo](https://sites.google.com/site/marcocongedo) is
a research scientist of [CNRS](http://www.cnrs.fr/en) (Centre National de la Recherche Scientifique), working at [UGA](https://www.univ-grenoble-alpes.fr/english/) (University of Grenoble Alpes), in Grenoble (France), the city where Jean-Baptiste Joseph Fourier has served as a Prefect:).

## Contact 
first name *dot* last name *at* gmail *dot* com

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://Marco-Congedo.github.io/FourierAnalysis.jl/dev)
