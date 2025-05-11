#   Unit examples of the FourierAnalysis Package for julia language
#
#   MIT License
#   Copyright (c) 2019-2022,
#   Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENTS :
#   This example shows how to compute coherence
#   and how to extract inforation from them.

#   ~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~  #
#                                                                             #
#   ~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~  #


using FourierAnalysis, FFTW, LinearAlgebra, Statistics, Plots

# add module for reading the two EEG text files to be used ater
push!(LOAD_PATH, @__DIR__)
using IOtxt

# Get EEG file names with complete path (they have extension .txt)
S=getFilesInDir(@__DIR__; ext=(".txt",))

# read the two EEG data files and put them in a Matrix object
X1=readEEG(S[1])
X2=readEEG(S[2])

#############################################################
# test coherence using sinusoids with same frequency and different phase.
using FourierAnalysis, FFTW, LinearAlgebra, Statistics, Plots
t=32
sr=32
f=2
noise=1000
tapering=rectangular
v=sinusoidal(10, f, sr, t*4, 0).+[randn()/noise for i=1:t*4]
w=sinusoidal(10, f, sr, t*4, π/4).+[randn()/noise for i=1:t*4]
y=sinusoidal(10, f, sr, t*4, π/2).+[randn()/noise for i=1:t*4]
z=sinusoidal(10, f, sr, t*4, π).+[randn()/noise for i=1:t*4]
X=[v w y z]
# as compared to v, w is out-of-phase, y is in between and z is at phase opposition
plot(v)
plot!(w)
plot!(y)
plot!(z)

# in order to avoid conflict with the DSP package also exporting 'coherence'
𝘾=FourierAnalysis.coherence(X, sr, t; tapering=tapering, tril=true)
𝘾=FourierAnalysis.coherence(X, sr, t; tapering=tapering, tril=false)
# get all five kinds of coherences
𝘾₁, 𝘾₂, 𝘾₃, 𝘾₄, 𝘾₅=FourierAnalysis.coherence(X, sr, t; tapering=tapering, tril=true, allkinds=true)
𝘾₁, 𝘾₂, 𝘾₃, 𝘾₄, 𝘾₅=FourierAnalysis.coherence(X, sr, t; tapering=tapering, tril=false, allkinds=true)
# check result against cross-spectra
𝐒=crossSpectra(X, sr, t; tapering=tapering, tril=true)
# compute coherence from a CrossSpectra object (lower triangles only)
𝘾=FourierAnalysis.coherence(𝐒)

𝐒=crossSpectra(X, sr, t; tapering=tapering, tril=false)
# compute coherence from a CrossSpectra object (full matrices)
𝘾=FourierAnalysis.coherence(𝐒)

# compute all five kinds of coherence from a CrossSpectra object
𝘾₁, 𝘾₂, 𝘾₃, 𝘾₄, 𝘾₅=FourierAnalysis.coherence(𝐒, allkinds=true)

# average full coherence in range 4:8Hz
Y=mean(𝘾, (4, 8))

# full coherence matrices in range 4:8Hz
E=extract(𝘾, (4, 8))

# average coherence in 1Hz-bands
Z=bands(𝘾, 1)

# smooth coherence a-posteriori
𝘾2=smooth(blackmanSmoother, 𝘾)
C=𝘾.y
D=𝘾2.y

# or get directly smoothed coherence
𝘾3=FourierAnalysis.coherence(X, sr, t; tapering=tapering, smoothing=hannSmoother, tril=true)

############################################
# coherence of several data matrices at once

# Compute the coherence altogether
𝓒=FourierAnalysis.coherence([X1, X2], sr, t; tapering=slepians(sr, t))

𝐗=[X1, X2]

# get all five kinds of coherences
𝓒₁, 𝓒₂, 𝓒₃, 𝓒₄, 𝓒₅=FourierAnalysis.coherence(𝐗, sr, t; tapering=slepians(sr, t), allkinds=true)
