push!(LOAD_PATH,"../src/")
push!(LOAD_PATH, @__DIR__)
using Documenter, DocumenterTools, DocumenterCitations, DocumenterInterLinks
using FourierAnalysis

makedocs(
   sitename = "FourierAnalysis",
   format = Documenter.HTML(),
   authors = "Marco Congedo, CNRS, France",
   modules = [FourierAnalysis],
   pages =  ["index.md"]
)

deploydocs(
    repo = "github.com/Marco-Congedo/FourierAnalysis.jl.git",
    target = "build",
    devurl = "dev",
)