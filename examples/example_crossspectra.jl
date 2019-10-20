#   Unit examples of the FourierAnalysis Package for julia language
#   v 0.0.1 - last update 24th of September 2019
#
#   MIT License
#   Copyright (c) 2019, Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENTS :
#   This example shows how to compute cross-spectra
#   and how to extract inforation from them.

#   ~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~  #
#                                                                             #
#   ~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~  #


using FourierAnalysis, FFTW, LinearAlgebra, Statistics, Plots, Plots.Measures

# add module for reading the two EEG text files to be used ater
push!(LOAD_PATH, @__DIR__)
using IOtxt

# Get EEG file names with complete path (they have extension .txt)
S=getFilesInDir(@__DIR__; ext=(".txt",))

# read the two EEG data files and put them in a Matrix object
X1=readEEG(S[1])
X2=readEEG(S[2])


# Cross-Spectra of EEG data
##########################################

t, sr, slide, tapering = 512, 128, 64, harris4

# gather some attributes to obtain nice spectra plots
spectraArgs=(left_margin   = 2mm,
             bottom_margin = 2mm,
             xtickfont     = font(10, "Times"),
             ytickfont     = font(10, "Times"))

# spectra
S=spectra(X1, sr, t; tapering=tapering, func=√)
plot(S; ytitle="Amplitude (\\muV)", spectraArgs...)
# smoothed spectra
S2=spectra(X1, sr, t; tapering=tapering, smoothing=blackmanSmoother, func=√)
plot(S2; ytitle="Amplitude (\\muV)", spectraArgs...)

#cross-spectra
𝙎=crossSpectra(X1, sr, t; tapering=tapering, tril=true)
# smooth a-posteriori the cross-spectra
𝙎2=smooth(blackmanSmoother, 𝙎)

# mean cross-spectra in 8Hz-12Hz range
alpha=mean(𝙎, (8, 12))

# extract all cross-spectra in 8Hz-12Hz range
E=extract(𝙎, (8, 12))

# mean smoothed cross-spectra in 8Hz-12Hz range
alpha=mean(𝙎2, (8, 12))

# extract smoothed cross-spectra in 8Hz-12Hz range
E=extract(𝙎2, (8, 12))

# get amplitude spectra from cross-spectra and compare with S (the long way)
S__=sqrt.(Real.([𝙎.y[i][j, j] for i=1:length(𝙎.y), j=1:size(X1, 2)]))
# the short way using FourierAnalysis
S_=Spectra(𝙎; func=√)

# check they are the same
norm(S.y-S_.y)

# cross-spectra of several data matrix at once
t, sr, slide, tapering = 1024, 128, 512, harris4
𝐗=[X1, X2]

# Compute the cross-spectra altogether
𝓢=crossSpectra(𝐗, sr, t; tapering=slepians(sr, t))

# mean cross-spectrum in 8Hz-12Hz range for each cross-spectra in 𝓢
alpha=mean(𝓢, (8, 12))
# mean of the above cross-spectra in one pass
alpha=mean(mean(𝓢, (8, 12)))
# the above computation without checking for homogeneity of elements in 𝓢 (faster)
alpha=mean(mean(𝓢, (8, 12), check=false))

# extract cross-spectra in 8Hz-12Hz range for all cross-spectra in 𝓢
E=extract(𝓢, (8, 12))

# average cross-spectra in 8Hz-12Hz range for 𝓢[1] only
Y=mean(𝓢[1], (8, 12))

# This (8 is the frequency in Hz)
Y=mean(𝓢[1], 8)
# is equivalent to this (64 is the bin where discrete Fourier Frequency 8Hz is)
𝓢[1].y[64]
# and equivalent to this
Y=extract(𝓢[1], 8)
# and equivalent to this
Y=extract(𝓢[1], (8, 8))

# Compute the lower-trianguler part of the cross-spectra
𝓢=crossSpectra(𝐗, sr, t; tapering=slepians(sr, t), tril=true)

# average cross-spectra in 8Hz-12Hz range for all 𝕃Vector in 𝓢
Y=mean(𝓢, (8, 12))

# extract cross-spectra in 8Hz-12Hz range for all 𝕃Vector in 𝓢
E=extract(𝓢, (8, 12))

# cross-spectra 𝓢[1] averaged in 1Hz band-pass regions
Z=bands(𝓢[1], 1)

# all cross-spectra in 𝓢 averaged in 1Hz band-pass regions
Z=bands(𝓢, 1)

# get amplitude spectra from cross-spectra (the long way)
S_=sqrt.(Real.([𝓢[1].y[i][j, j] for i=1:length(𝓢[1].y), j=1:size(X1, 2)]))*1.5687746
plot(S_)

# the short way using FourierAnalysis
S_sl=Spectra(𝓢[1]; func=x->sqrt(x)*1.5687746)
norm(S_sl.y-S_) # must be zero

# Compute non-linear cross-spectra
𝓢=crossSpectra(𝐗, sr, t; tapering=slepians(sr, t), tril=true, nonlinear=true)
