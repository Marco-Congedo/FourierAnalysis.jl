#   Unit examples of the FourierAnalysis Package for julia language
#   v 0.0.1 - last update 24th of September 2019
#
#   MIT License
#   Copyright (c) 2019, Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENTS :
#   This example shows how to compute the analytic signal via
#   Hilbert transform and how to extract inforation from them.

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

######################################################################
# Test Hilbert Transform Analytic Signal

using FourierAnalysis, FFTW, LinearAlgebra, Statistics, Plots, DSP
t=128; lab=["x", "real(y)", "imag(y)"]

# Analytic signal of one vector
x=sinusoidal(10, 2, 128, t, π/2; DC=10) # cosine
y=analyticsignal(x)
# note that real(y) is x without the DC level, i.e., x=real(y)+DC
plot([x, real(y), imag(y)]; labels=lab)

# make a check
s=sinusoidal(10, 2, 128, t, 0) # sine
norm(s-imag(y)) # should be zero

# Analytic Signal by DSP.jl
y2=hilbert(x)
norm(s-imag(y2)) # should be zero
# DSP.jl does not remove the DC level
# thus x=real(y2) in this case
plot([x, real(y2), imag(y2)]; labels=lab)

# Analytic signal of multiple vectors
x=hcat(x, sinusoidal(10, 3, 128, t, π/2; DC=10))
y=analyticsignal(x)

# welch-like analytic signal of one vector
# (note edge effects)
x=sinusoidal(10, 2, 128, t*4, π/2; DC=0)
y=analyticsignal(x, t)
plot([x, real(y), imag(y)]; labels=lab)

# Welch-like analytic signal of multiple vectors
x=hcat(x, sinusoidal(10, 3, 128, t*4, π/2; DC=0))
y=analyticsignal(x, t)


t=256
lim=512
Y=analyticsignal(X1[1:lim, :], lim)
plot([X1[1:lim, 1], real(Y[:, 1]), imag(Y[:, 1])])

# compute Welch-like Analytis Signal altogether for X1 and X2
𝒀 = analyticsignal([X1[1:lim, :], X2[1:lim, :]], lim)

## test amplitude, phase and polar functions

x=sinusoidal(10, 2, 128, t*4, 0).*sinusoidal(10, 1, 128, t*4, 0)

# amplitude and phase of a vector using analytic signal standard method
y=analyticsignal(x)
a=amplitude(y)
ϕ=phase(y, func=x->(x+π)/2π*50)
plot(x; labels="signal")
plot!(a; labels="amplitude")
plot!(ϕ; labels="phase")

# see what happen if `x` contains energy in frequencies below sr/wl Hz
# (see documentation of `analyticSignal` function)
y=analyticsignal(x, 64)
a=amplitude(y)
ϕ=phase(y, func=x->(x+π)/2π*50)
plot(x; labels="signal")
plot!(a; labels="amplitude")
plot!(ϕ; labels="phase")

# get Amplitude from analytic Signal of a data matrix holding multiple time-series
A=amplitude(Y)
plot(A[:, 1:2])

# get Phase from analytic Signal
𝛷=phase(Y)
plot(𝛷[:, 1:1])

# get Phase from analytic Signal and transform it in [-1, 1]
𝛷=phase(Y, func=x->(x+π)/2π)
plot(𝛷[:, 1:1])

# get Phase from analytic Signal and transform it to the sine of the phase
𝛷=phase(Y, func=x->sin(x))
plot(𝛷[:, 1:1])

# get Amplitude and phase from analytic Signal
A, 𝛷=polar(Y)
