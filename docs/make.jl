push!(LOAD_PATH,"../src/")
using Documenter, FourierAnalysis

makedocs(
   sitename="FourierAnalysis",
   modules=[FourierAnalysis],
   pages =  [
      "index.md",
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
      "hilbert.md",
   ]
)

#deploydocs(
#    repo = "github.com/Marco-Congedo/PosDefManifold.jl.git",
#    target = "build",
#    devurl = "dev",
#)
