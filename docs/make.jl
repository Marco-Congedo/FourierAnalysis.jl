push!(LOAD_PATH,"../src/")
push!(LOAD_PATH, @__DIR__)
using Documenter
using FourierAnalysis

makedocs(
   sitename="FourierAnalysis",
   format = Documenter.HTML(),
   authors="Marco Congedo, CNRS, France",
   modules=[FourierAnalysis],
   pages =  [
      "index.md",    
      "Main Module" => "MainModule.md",    
      "Tapering Windows" => "tapers.md",
      "frequency domain" => Any[
                           "Spectral Estimations" => "spectra.md",
                           "Cross-Spectral Matrices" => "crossspectra.md",
                           "Coherence Matrices" => "coherence.md",
                           "Goertzel's Algorithms" => "goertzel.md"
      ],
      "time-frequency(TF) domain" => Any[
                           "TF Representations" => "timefrequency.md",
                           "TF Univariate Measures" => "timefrequencyuni.md",
                           "TF Bivariate Measures " => "timefrequencybi.md"
      ],
      "utilities" => Any[
                        "Plots" => "recipes.md",
                        "Tools" => "tools.md",
                        "FFTW planners"  => "fftw.md",
                        "Filter Banks" => "filters.md",
                        "Hilbert Transform" => "hilbert.md"
      ]
   ]
)

deploydocs(
    repo = "github.com/Marco-Congedo/FourierAnalysis.jl.git",
    target = "build",
    devurl = "dev"
)