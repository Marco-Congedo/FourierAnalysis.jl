module  IOtxt

# types
StringVector = Vector{String}

export
    getFilesInDir,
    readEEG


# read EEG data from a .txt file in LORETA format and put it in a matrix
# of dimension txn, where n=#electrodes and t=#samples.
# if optional keyword argument `msg` is not empty, print `msg` on exit
function readEEG(filename::String; msg::String="")
    S=readlines(filename) # read the lines of the file as a vector of strings
    t=length(S) # number of samples
    n=length(split(S[1])) # get the number of electrodes
    X=Matrix{Float64}(undef, t, n) # declare the X Matrix
    for j=1:t
        x=split(S[j]) # this get the n potentials from a string
        for i=1:n
            X[j, i]=parse(Float64, x[i])
        end
    end
    if !isempty(msg) println(msg) end
    return X
end

# Read several EEG data from .txt files in LORETA format given in `filenames`
# ( a Vector of strings) and put them in a vector of matrices object
# print "read file "*[filenumber]*": "*[filename] after each file has been read.
readEEG(filenames::StringVector, skip::Vector{Int}=[]) =
        [readEEG(filenames[f]; msg="read file $f: "*basename(filenames[f])) for f in eachindex(filenames) if f ∉ skip]


# get the complete path of all files in `dir` as a vector of strings
# `ext` is an optional tuple of file extensions given strings.
# If it is provided, only files with those extensions will be included
# in the returned vector.
# the extensions must be entered in lowercase.
## Examples
# S=getFilesInDir(@__DIR__)
# S=getFilesInDir(@__DIR__; ext=(".txt", ))
# S=getFilesInDir(@__DIR__; ext=(".txt", ".jl"))
function getFilesInDir(dir::String; ext::Tuple=())
    if !isdir(dir) @error "Function `getFilesInDir`: input directory is incorrect!"
    else
        S=[]
        for (root, dirs, files) in walkdir(dir)
            if root==dir
                for file in files
                    if ext==() || ( lowercase(string(splitext(file)[2])) ∈ ext )
                       push!(S, joinpath(root, file)) # complete path and file name
                    end
                end
            end
        end
        isempty(S) && @warn "Function `getFilesInDir`: input directory does not contain any files"
        return StringVector(S)
    end
end

end # module

# Example
# dir="C:\\Users\\congedom\\Documents\\My Data\\EEG data\\NTE 84 Norms"

# Gat all file names with complete path
# S=getFilesInDir(dir)

# Gat all file names with complete path with extension ".txt"
# S=getFilesInDir(@__DIR__; ext=(".txt", ))

# read one file of NTE database and put it in a Matrix object
# X=readEEG(S[1])

# read all files of NTE database and put them in a vector of matrix
# 𝐗=readEEG(S)
