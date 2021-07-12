#   Unit examples of the FourierAnalysis Package for julia language
#
#   MIT License
#   Copyright (c) 2019-2021,
#   Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENTS :
#   This example shows how to compute univariate measures in
#   the time-frequency domain.

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

#############################################################################
# Univariate measures
# see Congedo, 2018: https://hal.archives-ouvertes.fr/hal-01868538/document

sr, wl, bandwidth=128, 512, 2
Pz=15
𝐱=[X1[:, Pz], X2[:, Pz]] # get the two times-series at electrode Pz
𝐘=TFanalyticsignal(𝐱, sr, wl, bandwidth; fmax=32, nonlinear=false)
𝐀=TFamplitude(𝐘)

# Mean Amplitude Eq. 0.1 (MAmp)
# compute the MAmp averaging in a TF region from a TFAnalyticSignalVector object
MAmp=meanAmplitude(𝐘, (8, 12), (1, 512); mode=mean)
# compute the MAmp averaging in a TF region from a TFAmplitudeVector object
MAmp=meanAmplitude(𝐀, (8, 12), (1, 512); mode=mean)
# compute the MAmp averaging in a TF region directly from data
MAmp=meanAmplitude(𝐱, sr, wl, (8, 12), (1, 512), bandwidth; mode=mean)

# compute the MAmp in a TF region from a TFAnalyticSignalVector object
MAmp=meanAmplitude(𝐘, (8, 12), (1, 512); mode=extract)
# compute the MAmp in a TF region from a TFAmplitudeVector object
MAmp=meanAmplitude(𝐀, (8, 12), (1, 512); mode=extract)
# compute the MAmp in a TF region directly from data
MAmp=meanAmplitude(𝐱, sr, wl, (8, 12), (1, 512), bandwidth; mode=extract)
# NB the Analytic Signal or the amplitude objects must be linear (note the nonlinear=false above)

# All these operation can be obtained on smoothed Amplitude, e.g.,
MAmp=meanAmplitude(𝐱, sr, wl, (8, 12), (1, 512), bandwidth;
                   mode=extract,
                   fsmoothing=hannSmoother,
                   tsmoothing=hannSmoother)

# The same operations can be done with the other univariate measures. For example:
# Concentration Eq. 0.2 (Con)
# compute the Con averaging in a TF region from a TFAnalyticSignalVector object
Con=concentration(𝐘, (8, 12), (1, 512); mode=mean)
# compute the Con in a TF region directly from data
Con=concentration(𝐱, sr, wl, (8, 12), (1, 512), bandwidth; mode=extract)
# Mean Direction Eq. 0.3 (MDir)
# compute the MDir averaging in a TF region directly from data
MDir=meanDirection(𝐱, sr, wl, (8, 12), (1, 512), bandwidth; mode=mean)
# compute the MDir in a TF region from a TFAnalyticSignalVector object
MDir=meanDirection(𝐘, (8, 12), (1, 512); mode=extract)

# and for the non-linear counterpart:
# Phase Concentration Eq. 0.4 (PCon)
# compute the Con in a TF region directly from data
Con=concentration(𝐱, sr, wl, (8, 12), (1, 512), bandwidth; mode=extract, nonlinear=true)
# Phase Mean Direction Eq. 0.5 (PMDir)
# compute the MDir averaging in a TF region directly from data
MDir=meanDirection(𝐱, sr, wl, (8, 12), (1, 512), bandwidth; mode=mean, nonlinear=true)

# If you try to compute a non-linear measure from a linear AnalyticSignal object
# you will get en error (see the REPL)
Con=concentration(𝐘, (8, 12), (1, 512); mode=mean, nonlinear=true)

# In order to compute those quantities from Analytic Signal objects
# first we need to compute a non-linear Analytic Signal objects:
𝐘=TFanalyticsignal(𝐱, sr, wl, bandwidth; fmax=32, nonlinear=true)
# then
# Phase Concentration Eq. 0.4 (PCon)
Con=concentration(𝐘, (8, 12), (1, 512); mode=mean, nonlinear=true)
# compute the MDir in a TF region from a TFAnalyticSignalVector object
MDir=meanDirection(𝐘, (8, 12), (1, 512); mode=extract, nonlinear=true)
# NB the Analytic Signal or the amplitude objects must be non-linear (note the nonlinear=true above)

############################################################################
