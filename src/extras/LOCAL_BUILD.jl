#   This script is not part of the FourierAnalysis.jl package.
#   It allows to build the package locally from the source code,
#   without actually installing the package.
#   You won't need this script for using the package.
#
#   MIT License
#   Copyright (c) 2019-2025, Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home
#
#   DIRECTIONS:
#   1) If you have installed the FourierAnalysis.jl from github or Julia registry, uninstall it.
#   3) Run this block (With VS code, click anywhere here and hit ALT+Enter)
#
begin
  push!(LOAD_PATH, abspath(@__DIR__, ".."))
  using FourierAnalysis
end

