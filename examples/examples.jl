#   Unit examples of the FourierAnalysis Package for julia language
#   v 0.0.1 - last update 10th of July 2019
#
#   MIT License
#   Copyright (c) 2019, Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENTS :
#   This unit provides examples to learn how to use the FourierAnalysis package.

#   ~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~  #
#                                                                             #
#   ~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~  #


########## Test goertzel algorithms #############
using FourierAnalysis, FFTW, LinearAlgebra, Statistics, Plots
sr, t, f, a = 128, 128, 5, 10
v=sinusoidal(a, f, sr, t, 0)
c=goertzel(v, f, sr, t) # should be 0+aim
d=goertzel2(v, f, sr, t) # should be 0+aim
abs2(c) # power, should be a²
abs2(d) # power, should be a²

ftrue=f+0.15
v=sinusoidal(a, ftrue, sr, t, 0 )
c=goertzel(v, f, sr, t) # should be 0+aim
d=goertzel2(v, ftrue, sr, t) # should be 0+aim
abs2(c) # power, should be a²
abs2(d) # power, should be a²
#################################################


############ plot tapering windows #############
using FourierAnalysis, FFTW, LinearAlgebra, Statistics, Plots

# using the standard plot function
tapers=[TaperKind(i) for i=1:8]
X=zeros(t, 8)
for i=1:8 X[:, i] = taper(tapers[i], t).y end
labels=[string(tapers[i]) for i=1:8]
plot(X; labels=labels)

begin
    # discrete prolate spheroid sequences
    H=taper(slepian, t, α=4, n=7)
    # This is a special plot function declared in the `plots.jl unit`
    plot(H, title="Slepian multi-tapering (dpss)")
end


################################################

####### Plot Sinusoidals at all positive Fourier Discrete Frequencies ########
using FourierAnalysis, FFTW, LinearAlgebra, Statistics, Plots
begin
    sr, t, = 128, 64
    t¼ = t÷4
    V=Matrix{Float64}(undef, t, t¼)
    for i=1:t¼ V[:, i]=sinusoidal(10, b2f(i, sr, t), sr, t, 1) end

    # plot the sinusoids in V
    Z=deepcopy(V)
    for i=2:size(Z, 2), j=1:size(Z, 1) Z[j, i]+=(25.0*i) end
    plot(Z, label="")
end
#################################################################


######## compute the FFT of one epoch of a time series ########
using FourierAnalysis, FFTW, LinearAlgebra, Statistics, Plots
sr, t, f, a = 128, 128, 8, 2
v=sinusoidal(a, f, sr, t, 0.6)
plot(v) # plot the time series
tapering=harris4

# The long way using FFTW.jl only:

# get a plan.  The normalization 2/t outputs correct peak amplitude
P=plan_rfft(v)*(2/t)

# get a tapering window
a=taper(tapering, t)

# plot the tapered signal
plot!(v.*a.y)

# do the FFT on tapered data
Ψ=P*(v.*a.y)

# get the power spectrum
Σ=@. sqrt(abs2(Ψ))
bar!(Σ[2:end])

# check the power at frequency f and neighbooring bins
pos=f2b(f, sr, t, DC=true)

Σ[pos-1]
Σ[pos]
Σ[pos+1]

# The fast way using FourierAnalysis.jl:
Σ2=spectra(v, sr, t; tapering=tapering, func=sqrt)
bar!(Σ2.y)

# equivalently:
Σ2=spectra(v, sr, t; tapering=taper(harris4, t), func=sqrt)
bar!(Σ2.y)

####################################################

# test spectra with odd samples and with DC
using FourierAnalysis, FFTW, LinearAlgebra, Statistics, Plots

sr=64; t=64
f, a = fres(sr, t)*4, 2
v=sinusoidal(a, f, sr, t, 0.6; DC=3)
plot(v) # plot the time series
Σ2=spectra(v, sr, t; tapering=rectangular, DC=true, func=√)
bar!(Σ2.y)

####################################################
### Check amplitude spectra at all Fourier Frequency (rectangular taper) ###
# the amplitudes are in increment of 10 along frequencies
# NB when t is even, correct amplitude for the last frequency is obtained
# only if the sinusoidal has a phase=π/6
sr, t, = 16, 32
t½ = t÷2
V=Matrix{Float64}(undef, t, t½)
for i=1:t½ V[:, i]=sinusoidal(10*i, b2f(i, sr, t), sr, t, π/6) end

# plot the sinusoids in V
Z=deepcopy(V)
for i=2:size(Z, 2), j=1:size(Z, 1) Z[j, i]+=(10.0*i) end
plot(Z, label="")

# long way using FFTW only.jl
P=plan_rfft(V, 1)*(2/t)
Ψ=P*V
Ψ[end, :]*=0.5*2^0.25 # correction for Nyquist frequency
Σ=abs2.(Ψ)
bar(sqrt.(Σ[brange(t, DC=true), :]), labels="")

# fast way using FourierAnalysis.jl
Σ2=spectra(V, sr, t; tapering=rectangular, func=√, DC=true)
bar(Σ2.y[brange(t, DC=true), :], labels="")

#############################################################################


### Check amplitude spectra on long signals obtained by welch methods
# one sinusoidal is at an exact discrete Fourier Frequency and the other not
# Rectangular window
sr, t, f, a = 128, 128, 10, 0.5
v=sinusoidal(a, f, sr, t*16)+sinusoidal(a, f*3.5+0.5, sr, t*16)+randn(t*16)
Σ=spectra(v, sr, t; tapering=rectangular, func=√)
bar(Σ.y, labels="raw")

# harris4 window (default)
sr, t, f, a = 128, 128, 10, 0.5
v=sinusoidal(a, f, sr, t*16)+sinusoidal(a, f*3.5+0.5, sr, t*16)+randn(t*16)
Σ=spectra(v, sr, t; func=√)
bar(Σ.y, labels="raw")

#smooth spectra
Σ2=smooth(blackmanSmoother, Σ)
bar!(Σ2.y, labels="smoothed")



######## Check amplitude spectrum obtained with dpss (slepian tapers) ########
sr, t = 128, 128
t½ = t÷2
V=Matrix{Float64}(undef, t, t½)
for i=1:t½ V[:, i]=sinusoidal(1, b2f(i, sr, t), sr, t, 0 ) end
# plot half of the the first 7 sinusoids in V and -V
plot(V[1:t½+1, 1:7], color=:gray, legend=false)
# plot!(-V[1:t½+1, 1:7], color=:gray, legend=false)

# the correction 1.5687746 on amplitude spectra for slepian return approx. correct amp.
Σ=spectra(V, sr, t; tapering=taper(slepian, t; α=2, n=3), func=x->sqrt(x)*1.5687746)
plot(Σ.y, labels="")
#############################################################

#### Compute Welch amplitude spectra of EEG data ####
using FourierAnalysis, FFTW, LinearAlgebra, Statistics, Plots, IOtxt

dir="C:\\Users\\congedom\\Documents\\My Data\\EEG data\\NTE 84 Norms"
# Gat all file names with complete path
S=getFilesInDir(dir)
# read some files of NTE database and put it in a Matrix object
X1=readEEG(S[1])
X2=readEEG(S[2])
X3=readEEG(S[3])

t, sr, slide, tapering = 1024, 128, 512, harris4
S=spectra(X1, sr, t; tapering=tapering, slide=sr, func=sqrt)
plot(S)

#smooth spectra
S2=smooth(blackmanSmoother, S)
plot(S2)

# extract spectra in alpha range (8Hz to 12Hz) at all electrodes
e=extract(S, (8, 12))

# You can use any combination of integer and real numbers
e=extract(S, (8, 12.5))
e=extract(S, (8.5, 12))

# the following equivalently extract the spectra at 10Hz only (1 bin)
e=extract(S, 10)
e=extract(S, (10))
e=extract(S, (10, 10))

# extract spectra at all frequencies. Equivalent to S.y
e=extract(S, :)

# All these ways to indicate frequency ranges work with the mean function as well

# compute average spectra in alpha range for each time-series
bar(mean(S, (8, 12)))
# or
bar(mean(S, (8.0,12.0)))

# this actually extract the power at 10Hz for each time series
# and is equivalent to extract(S, 10)
mean(S, 10)
# thus this computes the average power at 10Hz across time-series
mean(extract(S, 10)) # or mean(mean(S, 10))
# and this computes the mean across time-series and across frequencies in range [8Hz, 12Hz]
mean(extract(S, (8, 12))) # or equivalently mean(mean(S, (8, 12)))

# average spectrum across all frequencies for each time-series
bar(mean(S, :))

# plot average spectra in 2Hz-band-pass regions for all time-series
plot(bands(S, 2))

# plot average spectra in in 2Hz-band-pass regions for series 15-19
plot(bands(S, 2)[:, 15:19])

# plot average spectra in in 2Hz-band-pass regions for series 19
plot(bands(S, 2)[:, 19])

# plot average spectra in alpha range for series 15-19
bar(mean(S, (8,12))[15:19])

# extract spectra in alpha range for series 15-19
e=extract(S, (8, 12))[:, 15:19]


# Compute the spectra altogether
𝐗=[X1, X2, X3]
𝐒=spectra(𝐗, sr, t; tapering=tapering, slide=sr, func=sqrt)
plot(𝐒[1])
plot(𝐒[2])
plot(𝐒[3])

# do the same thing using a fast FFTW plan (wait up to 10s for computing the plan)
plan=Planner(plan_exhaustive, 10.0, t, eltype(𝐗[1]))
𝐒_fast=spectra(𝐗, sr, t; planner=plan, tapering=tapering, slide=sr, func=sqrt)

# compute the average spectra in the alpha range for all time-series and all subjects
mean(𝐒, (8,12))

# compute the average spectra in the alpha range across subjects for all time-series
mean(mean(𝐒, (8,12)))
# compute and plot it
plot(mean(mean(𝐒, (8,12))))

# extract spectra in alpha range for all time-series and all subjects
extract(𝐒, (8, 12))

# if you enter en illegal range, nothing will be done and you will get
# an error in the REPL explaining what is wrong in your argument
extract(𝐒, (0, 12))
mean(𝐒, (0, 128))

# extract 4Hz-band-pass average spectra for all electrodes and all subjects
bands(𝐒, 4)

# Apply smoothing in the above spectra computations
𝐒=spectra(𝐗, sr, t; tapering=tapering, smoothing=blackmanSmoother, slide=sr, func=sqrt)
plot(𝐒[1])
plot(𝐒[2])
plot(𝐒[3])

# plot the average spectrum across all electrodes for the three subjects
# using Julia standard mean function
plot(mean(𝐒[1].y[:, i] for i=1:size(𝐒[1].y, 2)))
plot!(mean(𝐒[2].y[:, i] for i=1:size(𝐒[2].y, 2)))
plot!(mean(𝐒[3].y[:, i] for i=1:size(𝐒[3].y, 2)))

# plot spectra in in 1Hz band-pass regions for all electrodes in 𝐒[1]
plot(bands(𝐒[1], 1))

# use slepian multi-tapering
S_sl=spectra(X1, sr, t; tapering=slepians(sr, t, 1.5), func=x->sqrt(x)*1.5687746)
plot(S_sl.y)


# Cross-Spectra
##########################################
using FourierAnalysis, FFTW, LinearAlgebra, Statistics, Plots, IOtxt
dir="C:\\Users\\congedom\\Documents\\My Data\\EEG data\\NTE 84 Norms"

# Gat all file names with complete path
S=getFilesInDir(dir)

# read some files of NTE database and put it in a Matrix object
X1=readEEG(S[1])

t, sr, slide, tapering = 512, 128, 64, harris4
# spectra
S=spectra(X1, sr, t; tapering=tapering, func=sqrt)
plot(S)
# smoothed spectra
S2=spectra(X1, sr, t; tapering=tapering, smoothing=blackmanSmoother, func=sqrt)
plot(S2)

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
S_=Spectra(𝙎; func=sqrt)

plot(S.y-S_.y)
norm(S.y-S_.y)

# get amplitude spectra from smoothed cross-spectra and compare with S2
S2_=Spectra(𝙎2, func=sqrt)
plot(S2.y-S2_.y)
norm(S2.y-S2_.y)



# cross-spectra
using FourierAnalysis, FFTW, LinearAlgebra, Statistics, Plots, IOtxt
t, sr, slide, tapering = 1024, 128, 512, harris4
dir="C:\\Users\\congedom\\Documents\\My Data\\EEG data\\NTE 84 Norms"
# Gat all file names with complete path
S=getFilesInDir(dir)
# read some files of NTE database and put it in a Matrix object
X1=readEEG(S[1])
X2=readEEG(S[2])
X3=readEEG(S[3])
𝐗=[X1, X2, X3]

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
# you can check that real
plot(v)
plot!(w)
plot!(y)
plot!(z)

𝘾=coherence(X, sr, t; tapering=tapering, tril=true)
𝘾=coherence(X, sr, t; tapering=tapering, tril=false)
# get all five kinds of coherences
𝘾₁, 𝘾₂, 𝘾₃, 𝘾₄, 𝘾₅=coherence(X, sr, t; tapering=tapering, tril=true, allkinds=true)
𝘾₁, 𝘾₂, 𝘾₃, 𝘾₄, 𝘾₅=coherence(X, sr, t; tapering=tapering, tril=false, allkinds=true)
# check result against cross-spectra
𝐒=crossSpectra(X, sr, t; tapering=tapering, tril=true)
# compute coherence from a CrossSpectra object (lower triangles only)
𝘾=coherence(𝐒)

𝐒=crossSpectra(X, sr, t; tapering=tapering, tril=false)
# compute coherence from a CrossSpectra object (full matrices)
𝘾=coherence(𝐒)

# compute all five kinds of coherence from a CrossSpectra object
𝘾₁, 𝘾₂, 𝘾₃, 𝘾₄, 𝘾₅=coherence(𝐒, allkinds=true)

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
𝘾3=coherence(X, sr, t; tapering=tapering, smoothing=hannSmoother, tril=true)

############################################
# coherence of several data matrices at once
using FourierAnalysis, FFTW, LinearAlgebra, Statistics, Plots, IOtxt
t, sr, slide, tapering = 1024, 128, 512, harris4
dir="C:\\Users\\congedom\\Documents\\My Data\\EEG data\\NTE 84 Norms"
# Gat all file names with complete path
S=getFilesInDir(dir)
# read some files of NTE database and put it in a Matrix object
X1=readEEG(S[1])
X2=readEEG(S[2])
X3=readEEG(S[3])
𝐗=[X1, X2, X3]

# Compute the coherence altogether
𝓒=coherence(𝐗, sr, t; tapering=slepians(sr, t))

# get all five kinds of coherences
𝓒₁, 𝓒₂, 𝓒₃, 𝓒₄, 𝓒₅=coherence(𝐗, sr, t; tapering=slepians(sr, t), allkinds=true)

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



using FourierAnalysis, FFTW, LinearAlgebra, Statistics, Plots, IOtxt, DSP
dir="C:\\Users\\congedom\\Documents\\My Data\\EEG data\\NTE 84 Norms"

# Gat all file names with complete path
S=getFilesInDir(dir)

# read some files of NTE database and put it in a Matrix object
X1=readEEG(S[1])
X2=readEEG(S[2])

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


####################
# test Filters

using Plots, DSP, FourierAnalysis
x=sinusoidal(2, 8, 128, 512, 0)
f, Y=filterbank(x, 128)
plot(x)
plot!((Y[:, 7]))
plot!((x-Y[:, 7]))

f, Y=filterbank(x, 128; filtkind=Butterworth(4))

x2=randn(512)+sinusoidal(2, 8, 128, 512, 0)
f, Y=filterbank(x2, 128)
plot(x)
plot!(x2)
plot!((Y[:, 7]))
plot!((x-Y[:, 7]))


####################
# test timefrequency

using Plots, FourierAnalysis
sr, t, bandwidht=128, 512, 2
h=taper(harris4, t)
x1=sinusoidal(10, 8, sr, t, 0)
x2=sinusoidal(10, 19, sr, t, 0)
x=Vector((x1+x2).*h.y+randn(t))
y1=sinusoidal(10, 6, sr, t, 0)
y2=sinusoidal(10, 16, sr, t, 0)
y=Vector((y1+y2).*h.y+randn(t))

# Vector of Time-Frequency Object for x and y
𝒀 = TFanalyticsignal([x, y], sr, sr*4)

# compute the mean in a TF region (8:12Hz and samples 1:128) for the first object
y=mean(𝒀[1], (8, 12), (1, 128))

# extract a TF region (8:12Hz and samples 1:128) for the first object
Z=extract(𝒀[1], (8, 12), (1, 128))

# compute the mean in a TF region (8:12Hz and samples 1:128) for the two objects
Y=mean(𝒀, (8, 12), (1, 128))
# do the same computation without checking homogeneity of the two objects in 𝒀
Y=mean(𝒀, (8, 12), (1, 128); check=false)

# extract the data in a TF region (8:12Hz and samples 1:128) for the two objects
E=extract(𝒀, (8, 12), (8, 12))
# do the same computation without checking homogeneity of the two objects in 𝒀
E=extract(𝒀, (8, 12), (8, 12); check=false)

# compute time-frequency object for vector x
Y=TFanalyticsignal(x, sr, t, bandwidht; fmax=32)

# plot the real part of the AS
heatmap(Y, real)

# ...the imaginary part of the AS
heatmap(Y, imag)

# ...the amplitude of the AS
heatmap(Y, amplitude)

# or generate a TFAmplitude object
A=TFamplitude(Y)
# and plot it (with different colors)
heatmap(A; c=:fire)

# ...the amplitude of the AS smoothed in the freq. dim.
heatmap(smooth(hannSmoother, noSmoother, Y), amplitude)

# ...the amplitude of the AS smoothed in freq. and time
heatmap(smooth(hannSmoother, hannSmoother, Y), amplitude)

# ...the phase of the AS and use custom colors
heatmap(Y, phase)

#or, generate a TFPhase object
ϴ=TFphase(Y)
# and plot it (with custom colors)
heatmap(ϴ; c=:pu_or)

# compute and plot phase in [0, 2π]
heatmap(TFphase(Y; func=x->x+π); c=:amp)

# compute and plot unwrapped phase
heatmap(TFphase(Y; unwrapped=true); c=:amp)

# smooth time-frequency analytic signal: smooth frequency
Z=smooth(blackmanSmoother, noSmoother, Y)
# plot amplitude of smoothed analytic signal
heatmap(Z, amplitude)
# not equivalently (!), you can create an amplitude object and smooth them:
# in this case the amplitude is smoothed, not the analutic signal
A=smooth(blackmanSmoother, noSmoother, TFamplitude(Y))
heatmap(A) # bluesreds

# Smoothing raw phase estimates is unappropriate
# since the phase is a discontinous function, however it makes sense to smooth
# phase if the phase is unwrapped.
heatmap(smooth(blackmanSmoother, noSmoother, TFphase(Y; unwrapped=true)); c=:amp)

# smooth analytic signal:smooth frequency and time
E=smooth(blackmanSmoother, blackmanSmoother, Y)
# plot amplitude of smoothed analytic signal
heatmap(E, amplitude)
# plot phase of smoothed analytic signal
heatmap(E, phase) # bluesreds
# not equivalently (!), you can create amplitude and phase object and smooth them
A=smooth(blackmanSmoother, blackmanSmoother, TFamplitude(Y))
heatmap(A)
ϴ=smooth(blackmanSmoother, blackmanSmoother, TFphase(Y, unwrapped=true))
heatmap(ϴ; c=:amp) # bluesreds
# smooth again
ϴ=smooth(blackmanSmoother, blackmanSmoother, ϴ)
heatmap(ϴ; c=:amp)
# and again ...
ϴ=smooth(blackmanSmoother, blackmanSmoother, ϴ)
heatmap(ϴ; c=:amp)

# you may also create all these objects already smoothed, for example
Y=TFanalyticsignal(x, sr, t, bandwidht; fmax=32, fsmoothing=hannSmoother, tsmoothing=hannSmoother)
# plot amplitude of smoothed analytic signal
heatmap(Y, amplitude)

A=TFamplitude(x, sr, t, bandwidht; fmax=32, fsmoothing=hannSmoother, tsmoothing=hannSmoother)
# plot smoothed amplitude
heatmap(A)

# compute time-frequency object with non-linear analytic signal
Y=TFanalyticsignal(x, sr, t, bandwidht; fmax=32, nonlinear=true)
# check that the amplitude is now 1.0 everywhere
amplitude(Y.y)
# plot phase
heatmap(Y, phase; c=:bkr)

# get the center frequencies of TF Amplitude A
A.flabels

# extract the amplitude in a time-frequency region
extract(A, (2, 10), (1, 256))
# extract the amplitude in a time-frequency region at only one frequency
extract(A, 10, (1, 256))
# extract the amplitude at one temporal sample at one frequency
extract(A, 10, 12)
# or
extract(A, 10.0, 12)
# extract amplitude at one temporal sample in a frequency band-pass
extract(A, (10, 12), 12)
#or
extract(A, (10.0, 12.0), 12)
# extract amplitude at one temporal sample and all frequencies
extract(A, :, 12)
# NB: the extract function work in the same way
#     for objects TFAnayticSignal and TFPhase


# compute the mean in a time-frequency region:
mean(A, (2, 10), (1, 256))
# is equivalent to (but may be less efficient than)
mean(extract(A, (2, 10), (1, 256)))

# using column sign for extracting all time samples
extract(A, (2, 10), :)
# This:
extract(A, :, :)
# is equivalent to this:
A.y
# but if you don't need to extract all frequencies,
# use the extract function to control what frequencies will be extracted:
# This
extract(A, (4, 5), 10)
# is not equivalent to this
A.y[4:5, 10]
# since the `extract` function finds the correct rows corresponding
# to the sought frequencies (in Hz), while A.y[4:5, 10]
# just returns the element [4:5, 10] in the TF amplitude object


mean(A, 10, (1, 256))
mean(A, 10, 12)
# this is not equivalent to the following
A.y[10, 12]
# since the `mean` function also finds the correct rows corresponding
# to the sought frequencies (in Hz), while A.y[10, 12]
# just returns the element [10, 12] in the TF amplitude object

# Although the first center frequency in A is 2Hz, its
# band-pass region is 1-3Hz, therefore the frequency range 1:10 is accepted
mean(A, (1, 10), (1, 256))
# but this result in an error (see the REPL) since 0.5 is out of range
mean(A, (0.5, 10), (1, 256))

# using a colon sign for time range
a=mean(A, (1, 10), :)
# using an integer for time range (indicates one specific sample)
a=mean(A, (1, 10), 16)

# using a colon sign for frequency range
a=mean(A, :, (1, 16))
# using a real number for frequency range
a=mean(A, 8.5, (1, 16))

# This
a=mean(A, :, :)
# is equivalent to this
a=mean(A.y)


###############################################################

# multiple time-frequency at the same time
using FourierAnalysis, FFTW, LinearAlgebra, Statistics, Plots, IOtxt, DSP
dir="C:\\Users\\congedom\\Documents\\My Data\\EEG data\\NTE 84 Norms"

# Gat all file names with complete path
S=getFilesInDir(dir)

# read some files of NTE database and put it in a Matrix object
X1=readEEG(S[1])
X2=readEEG(S[2])
X3=readEEG(S[3])
Pz=15
𝐱=[X1[:, Pz], X2[:, Pz], X3[:, Pz]] # get the three times-series at electrode Pz
𝐘=TFanalyticsignal(𝐱, sr, t, bandwidht; fmax=32, nonlinear=false)
𝐀=TFamplitude(𝐘)
heatmap(𝐀[1])
heatmap(𝐀[2])
heatmap(𝐀[3]; c=:pu_or) # c=:pu_or

# plot the power over time from instantaneous amplitude
plot([sum(𝐀[3].y[:, t]) for t=1:size(𝐀[3].y, 2)])
# plot the smoothed power over time from the instantaneous amplitude
A=smooth(noSmoother, blackmanSmoother, 𝐀[3])
plot([sum(A.y[:, t]) for t=1:size(A.y, 2)])

# 𝐀 holds 3 TFAmplitude objects;
# The following will create a 3-Vector of matrices
# holding the amplitude in the TF region defined by frequency 8:12 and samples 1:128:
extract(𝐀, (8, 12), (1, 128))

# use of a column is allowed for the frequency range as all the amplitudes have the same number of frequencies:
extract(𝐀, :, (1, 128))

# however, the use of a column is not allowed for the time range as the number
# of samples is different in the three TFAmplitude objects
extract(𝐀, (8, 12), :)

# the same goes for the mean function
mean(𝐀, (8, 12), :)

# to obtain the grand-average across all objects:
a=mean(mean(𝐀, (8, 12), (1, 128)))

#############################################################################
# Univariate measures
# see Congedo, 2018: https://hal.archives-ouvertes.fr/hal-01868538/document

using FourierAnalysis, FFTW, LinearAlgebra, Statistics, Plots, IOtxt, DSP
dir="C:\\Users\\congedom\\Documents\\My Data\\EEG data\\NTE 84 Norms"
sr, wl, bandwidht=128, 512, 2
# Gat all file names with complete path
S=getFilesInDir(dir)

# read some files of NTE database and put it in a Matrix object
X1=readEEG(S[1])
X2=readEEG(S[2])
X3=readEEG(S[3])
Pz=15
𝐱=[X1[:, Pz], X2[:, Pz], X3[:, Pz]] # get the three times-series at electrode Pz
𝐘=TFanalyticsignal(𝐱, sr, wl, bandwidht; fmax=32, nonlinear=false)
𝐀=TFamplitude(𝐘)

# Mean Amplitude Eq. 0.1 (MAmp)
# compute the MAmp averaging in a TF region from a TFAnalyticSignalVector object
MAmp=meanAmplitude(𝐘, (8, 12), (1, 512); mode=mean)
# compute the MAmp averaging in a TF region from a TFAmplitudeVector object
MAmp=meanAmplitude(𝐀, (8, 12), (1, 512); mode=mean)
# compute the MAmp averaging in a TF region directly from data
MAmp=meanAmplitude(𝐱, sr, wl, (8, 12), (1, 512), bandwidht; mode=mean)

# compute the MAmp in a TF region from a TFAnalyticSignalVector object
MAmp=meanAmplitude(𝐘, (8, 12), (1, 512); mode=extract)
# compute the MAmp in a TF region from a TFAmplitudeVector object
MAmp=meanAmplitude(𝐀, (8, 12), (1, 512); mode=extract)
# compute the MAmp in a TF region directly from data
MAmp=meanAmplitude(𝐱, sr, wl, (8, 12), (1, 512), bandwidht; mode=extract)
# NB the Analytic Signal or the amplitude objects must be linear (note the nonlinear=false above)

# All these operation can be obtained on smoothed Amplitude, e.g.,
MAmp=meanAmplitude(𝐱, sr, wl, (8, 12), (1, 512), bandwidht;
                   mode=extract,
                   fsmoothing=hannSmoother,
                   tsmoothing=hannSmoother)

# The same operations can be done with the other univariate measures. For example:
# Concentration Eq. 0.2 (Con)
# compute the Con averaging in a TF region from a TFAnalyticSignalVector object
Con=concentration(𝐘, (8, 12), (1, 512); mode=mean)
# compute the Con in a TF region directly from data
Con=concentration(𝐱, sr, wl, (8, 12), (1, 512), bandwidht; mode=extract)
# Mean Direction Eq. 0.3 (MDir)
# compute the MDir averaging in a TF region directly from data
MDir=meanDirection(𝐱, sr, wl, (8, 12), (1, 512), bandwidht; mode=mean)
# compute the MDir in a TF region from a TFAnalyticSignalVector object
MDir=meanDirection(𝐘, (8, 12), (1, 512); mode=extract)

# and for the non-linear counterpart:
# Phase Concentration Eq. 0.4 (PCon)
# compute the Con in a TF region directly from data
Con=concentration(𝐱, sr, wl, (8, 12), (1, 512), bandwidht; mode=extract, nonlinear=true)
# Phase Mean Direction Eq. 0.5 (PMDir)
# compute the MDir averaging in a TF region directly from data
MDir=meanDirection(𝐱, sr, wl, (8, 12), (1, 512), bandwidht; mode=mean, nonlinear=true)

# If you try to compute a non-linear measure from a linear AnalyticSignal object
# you will get en error (see the REPL)
Con=concentration(𝐘, (8, 12), (1, 512); mode=mean, nonlinear=true)

# In order to compute those quantities from Analytic Signal objects
# first we need to compute a non-linear Analytic Signal objects:
𝐘=TFanalyticsignal(𝐱, sr, wl, bandwidht; fmax=32, nonlinear=true)
# then
# Phase Concentration Eq. 0.4 (PCon)
Con=concentration(𝐘, (8, 12), (1, 512); mode=mean, nonlinear=true)
# compute the MDir in a TF region from a TFAnalyticSignalVector object
MDir=meanDirection(𝐘, (8, 12), (1, 512); mode=extract, nonlinear=true)
# NB the Analytic Signal or the amplitude objects must be non-linear (note the nonlinear=true above)


############################################################################




#############################################################################
# Bivariate measures
# see Congedo, 2018: https://hal.archives-ouvertes.fr/hal-01868538/document

using FourierAnalysis, FFTW, LinearAlgebra, Statistics, Plots, IOtxt, DSP
dir="C:\\Users\\congedom\\Documents\\My Data\\EEG data\\NTE 84 Norms"
sr, wl, bandwidht=128, 512, 2
# Gat all file names with complete path
S=getFilesInDir(dir)

# read some files of NTE database and put it in a Matrix object
X1=readEEG(S[1])
X2=readEEG(S[2])
X3=readEEG(S[3])
Pz=15
Fz=5
𝐱₁=[X1[:, Pz], X2[:, Pz], X3[:, Pz]] # get the three times-series at electrode Pz
𝐱₂=[X1[:, Fz], X2[:, Fz], X3[:, Fz]] # get the three times-series at electrode Fz
𝐘₁=TFanalyticsignal(𝐱₁, sr, wl, bandwidht; fmax=32, nonlinear=false)
𝐘₂=TFanalyticsignal(𝐱₂, sr, wl, bandwidht; fmax=32, nonlinear=false)
𝐀₁=TFamplitude(𝐘₁)
𝐀₂=TFamplitude(𝐘₂)


# Comodulation Eq. 0.11 (Com)
# compute the Com averaging in a TF region from a TFAnalyticSignalVector object
# NB If you compute it from Analytic Signal or Amplitude objects, those
# must be linear (note the nonlinear=false above)
Com=comodulation(𝐘₁, 𝐘₂, (8, 12), (1, 512); mode=mean)
# compute the Com averaging in a TF region from a TFAmplitudeVector object
Com=comodulation(𝐀₁, 𝐀₂, (8, 12), (1, 512); mode=mean)
# compute the Com averaging in a TF region directly from data
# In this care you don't have to worry about linearity of the analytic signal
Com=comodulation(𝐱₁, 𝐱₂, sr, wl, (8, 12), (1, 512), bandwidht; mode=mean)
# You can compute comodulation from smoothed amplitude:
Com=comodulation(𝐱₁, 𝐱₂, sr, wl, (8, 12), (1, 512), bandwidht;
                 mode=mean,
                 fsmoothing=blackmanSmoother,
                 tsmoothing=blackmanSmoother)


# you can go faster pre-computing a FFTW plan.
# This is aslo useful when you have to call the comodulation function several times
plan=Planner(plan_patient, 5, wl, Float64, true)
Com=comodulation(𝐱₁, 𝐱₂, sr, wl, (8, 12), (1, 512), bandwidht; mode=mean, planner=plan)


# compute the Com in a TF region from a TFAnalyticSignalVector object
Com=comodulation(𝐘₁, 𝐘₂, (8, 12), (1, 512); mode=extract)
# compute the Com in a TF region from a TFAmplitudeVector object
Com=comodulation(𝐀₁, 𝐀₂, (8, 12), (1, 512); mode=extract)
# compute the Com in a TF region directly from data
Com=comodulation(𝐱₁, 𝐱₂, sr, wl, (8, 12), (1, 512), bandwidht; mode=extract)

# All these operations can be done also for coherence measures, for example
Coh=coherence(𝐘₁, 𝐘₂, (8, 12), (1, 512); mode=mean)
Coh=coherence(𝐘₁, 𝐘₂, (8, 12), (1, 512); mode=extract)
# Compute all 5 coherence types
Coh=coherence(𝐘₁, 𝐘₂, (8, 12), (1, 512); mode=extract, allkinds=true)


# phase coherence (phase-locking value)
𝐘₁=TFanalyticsignal(𝐱₁, sr, wl, bandwidht; fmax=32, nonlinear=true)
𝐘₂=TFanalyticsignal(𝐱₂, sr, wl, bandwidht; fmax=32, nonlinear=true)
Coh=coherence(𝐘₁, 𝐘₂, (8, 12), (1, 512); mode=mean, nonlinear=true)

# or directly from data (no need to compute non-linear analytic signal in this case)
Coh=coherence(𝐱₁, 𝐱₂, sr, wl, (8, 12), (1, 512), bandwidht; mode=mean, nonlinear=true)

# and also for non-linear meausures
# compute non-linear analyticSignal
𝐘₁=TFanalyticsignal(𝐱₁, sr, wl, bandwidht; fmax=32, nonlinear=true)
𝐘₂=TFanalyticsignal(𝐱₂, sr, wl, bandwidht; fmax=32, nonlinear=true)
# although you are allowed to compute the amplitude of non-linear
# analytic sygnal this way, it does not make sense as the amplitude
# is 1.0 everywhere, as you can check
𝐀₁=TFamplitude(𝐘₁)
𝐀₂=TFamplitude(𝐘₂)

##################################################################




struct ShortFun{T<:Function}
    f::T
    expr
end
Base.show(io::IO, sf::ShortFun) = print(io, sf.expr)
(sf::ShortFun)() = sf.f()
macro shortfun(expr)
    esc(:(ShortFun($expr, $(Expr(:quote, expr)))))
end

# does not work, it does not accept the function to be passed as argument
f = string(@shortfun sqrt)
# or f = string(@shortfun x->x+2)

f1=findfirst("#=", f)
f2=findfirst("=#", f)
if f1≠nothing && f2≠nothing
    f_=f[1:first(f1)-1]*f[last(f2)+1: end]
    f_=replace(f_, "begin" => "")
    f_=replace(f_, "end" => "")
    f_=replace(f_, "\n" => "")
    f_=replace(f_, " " => "")
else f_=f
end

g=@shortfun sqrt



a=x->x^2


# temp
######## compute the FFT of one epoch of a time series ########
using FourierAnalysis, FFTW, LinearAlgebra, Statistics, Plots
sr, t, a = 12, 13, 2
f=sr/t*6
v=sinusoidal(a, f, sr, t, π/6; DC=3)
v
plot(v) # plot the time series
P=plan_rfft(v)
Ψ=P*v*(2/t)
Σ=abs.(Ψ)
bar(Σ)
S=spectra(v, sr, t; tapering=rectangular, DC=true, func=sqrt)
bar!(S.y)
