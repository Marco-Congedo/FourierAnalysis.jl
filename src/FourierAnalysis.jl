#    Main Module of the FourierAnalysis Package for julia language
#    v 0.2.0 - last update 20th of October 2019
#
#    MIT License
#    Copyright (c) 2019, Marco Congedo, CNRS, Grenobe, France:
#    https://sites.google.com/site/marcocongedo/home

# __precompile__()

module  FourierAnalysis

using   Base.Threads,
        LinearAlgebra,
        Statistics,
        FFTW,
        AbstractFFTs,
        DSP,
        RecipesBase


# # # # # Special instructions and variables # # # # #
# N.B. multi-threaded planners are not yet exported by FFTW
BLAS.set_num_threads(Sys.CPU_THREADS-Threads.nthreads())


# # # # # constants # # # # #
const 📌            = "FourierAnalysis"
#const qroot2        = 2^0.25
const titleFont     = "\x1b[36m"
const separatorFont = "\x1b[33;1m"
const defaultFont   = "\x1b[0m"
const greyFont      = "\x1b[90m"

# # # # # types # # # # #

IntOrReal    = Union{Int, Real}

# Frequency ranges (in Hz) for the mean function applying to (Cross)Spcetra objects
fInterval    = Union{IntOrReal, Tuple{IntOrReal, IntOrReal}, Colon}

# time ranges (in samples) for the mean function applying to TFx objects
tInterval    = Union{Int, Tuple{Int, Int}, Colon}


# smoother kind; see tools.jl
@enum Smoother begin
    noSmoother       = 1
    hannSmoother     = 2
    hammingSmoother  = 3
    blackmanSmoother = 4
end

# see spectra.jl
SpectraData      = AbstractArray{T} where T<:Real
struct Spectra
    y           :: SpectraData  # spectra
    sr          :: Int          # sampling rate
    wl          :: Int          # FFT window length
    DC          :: Bool         # DC component in first bin if true, not present otherwise
    taper       :: String       # time-domain tapering type as a string (with parameters for Slepians)
    flabels     :: Vector{T} where T<:Union{Real, Int} # Fourier Discrete Frequencies
    func        :: Function
    smoothing   :: Smoother # enum type, see above
end
SpectraVector   = Vector{Spectra}

# see crossspectra.jl
CrossSpectraData = Union{Vector{LowerTriangular}, Vector{Hermitian}}
struct CrossSpectra
    y           :: CrossSpectraData
    sr          :: Int          # sampling rate
    wl          :: Int          # FFT window length
    DC          :: Bool         # DC component in first bin if true, not present otherwise
    taper       :: String       # time-domain tapering type as a string (with parameters for Slepians)
    flabels     :: Vector{T} where T<:Union{Real, Int} # Fourier Discrete Frequencies
    nonlinear   :: Bool
    smoothing   :: Smoother
    tril        :: Bool
end
CrossSpectraVector = Vector{CrossSpectra}


# see coherence.jl
CoherenceData = Union{Vector{LowerTriangular}, Vector{Hermitian}}
struct Coherence
    y           :: CoherenceData
    sr          :: Int          # sampling rate
    wl          :: Int          # FFT window length
    DC          :: Bool         # DC component in first bin if true, not present otherwise
    taper       :: String       # time-domain tapering type as a string (with parameters for Slepians)
    flabels     :: Vector{T} where T<:Union{Real, Int} # Fourier Discrete Frequencies
    nonlinear   :: Bool
    smoothing   :: Smoother
    tril        :: Bool
end
CoherenceVector = Vector{Coherence}
CoherenceVector₂ = Vector{CoherenceVector}

# NB data field `y` in CrossSpectra and Coherence structures is one of these types:
# Array{LowerTriangular,1}, Array{Hermitian,1}
# CrossSpectraVector and CoherenceVector are of one of these types:
# Array{Array{LowerTriangular,1},1}, Array{Array{Hermitian,1},1}
# All these types are straightforwardly compatible with package PosDefMaifold.


# see timefrequency.jl
TFAnalyticSignalData = Matrix{T} where T<:Complex
struct TFAnalyticSignal
    y          :: TFAnalyticSignalData # Analytic Signal, dim1=freq, dim2=time
    bandwidht  :: IntOrReal
    flabels    :: Vector{S} where S<:Real # center frequency axis labels (in Hz)
    nonlinear  :: Bool # true if nonlinear, false otherwise
    fsmoothing :: Smoother
    tsmoothing :: Smoother
end
TFAnalyticSignalVector = Vector{TFAnalyticSignal}


TFAmplitudeData = Matrix{T} where T<:Real
struct TFAmplitude
    y          :: TFAmplitudeData # Analytic (Instantaneous) Amplitude, dim1=freq, dim2=time
    bandwidht  :: IntOrReal
    flabels    :: Vector{S} where S<:Real # center frequency axis labels (in Hz)
    fsmoothing :: Smoother
    tsmoothing :: Smoother
    func       :: Function
end
TFAmplitudeVector = Vector{TFAmplitude}

TFPhaseData = Matrix{T} where T<:Real
struct TFPhase
    y          :: TFPhaseData # Analytic (Instantaneous) Phase, dim1=freq, dim2=time
    bandwidht  :: IntOrReal
    flabels    :: Vector{S} where S<:Real # center frequency axis labels (in Hz)
    nonlinear  :: Bool # true if nonlinear, false otherwise
    fsmoothing :: Smoother
    tsmoothing :: Smoother
    unwrapped  :: Bool
    func       :: Function
end
TFPhaseVector = Vector{TFPhase}

# all frequency domain objects
FDobjects=Union{Spectra, CrossSpectra, Coherence}
FDobjectsVector=Union{SpectraVector, CrossSpectraVector, CoherenceVector}
# all Time-frequency domain objects
TFobjects=Union{TFAnalyticSignal, TFAmplitude, TFPhase}
TFobjectsVector=Union{TFAnalyticSignalVector, TFAmplitudeVector, TFPhaseVector}

#import DSP: dpss
import Statistics.mean
import Base: show, *, conj, real, imag #,
import DSP: dB

export

# From this module
fInterval,
tInterval,
Smoother,
    noSmoother,
    hannSmoother,
    hammingSmoother,
    blackmanSmoother,

SpectraType,            Spectra,            SpectraVector,
CrossSpectraType,       CrossSpectra,       CrossSpectraVector,
CoherenceType,          Coherence,          CoherenceVector,
TFAnalyticSignalType,   TFAnalyticSignal,   TFAnalyticSignalVector,
TFAmplitudeType,        TFAmplitude,        TFAmplitudeVector,
TFPhaseType,            TFPhase,            TFPhaseVector,

FDobjects, FDobjectsVector,
TFobjects, TFobjectsVector,

# from fftw.jl
plan_estimate,
plan_measure,
plan_patient,
plan_exhaustive,
plan_conserve_memory,
plan_wisdom_only,
plan_unaligned,
plan_preserve_input,
Planner,

# From tapers.jl
TaperKind,
    rectangular,
    triangular,
    hann,
    hamming,
    blackman,
    harris4,
    riesz,
    parzen,
    slepian,
taper,
slepians,
taperinfo,

# from tools.jl
sinusoidal,
fres,
f2b,
b2f,
fdf,
brange,
bbands,
fbands,
dB,
amplitude,
phase,
polar,
unwrapPhase,
smooth,
extract, extr,
mean,
bands,
sameParams,
isLinear,
isNonLinear,
isUnwrapped,

# from goertzel.jl
goertzel,
goertzel_fast,
goertzel2,

# from spectra.jl
spectra,

# from crossspectra.jl
crossSpectra,

# from coherence.jl
coherence,

# from hilbert.jl
analyticsignal,

# from filters.jl
FilterDesign,
filterbank,

# from timefrequency.jl
TFanalyticsignal,
TFamplitude,
TFphase,

# from timefrequencyuni.jl
meanAmplitude, mamp,
concentration, con,
meanDirection, mdir,

# from timefrequencybi.jl
comodulation, com,
coherence, coh,

# from recipes.jl
tfAxes

include("fftw.jl")
include("tapers.jl")
include("tools.jl")
include("goertzel.jl")
include("spectra.jl")
include("crossspectra.jl")
include("coherence.jl")
include("hilbert.jl")
include("filters.jl")
include("timefrequency.jl")
include("timefrequencyuni.jl")
include("timefrequencybi.jl")
include("recipes.jl")

# welcome message
println("\n⭐ "," Welcome to the","\x1b[36m"," FourierAnalysis ","\x1b[0m","package", " ⭐\n")
@info(" ")
println(" Your Machine `",gethostname(),"` (",Sys.MACHINE, ")")
println(" runs on kernel ",Sys.KERNEL," with word size ",Sys.WORD_SIZE,".")
println(" CPU  Threads: ",Sys.CPU_THREADS)
# Sys.BINDIR # julia bin directory
println(" Base.Threads: ", "$(Threads.nthreads())")
println(" BLAS Threads: ", "$(Sys.CPU_THREADS-Threads.nthreads())", "\n")

end # module
