#   This script allows to precompile the FourierAnalysis package
#   and `using` it locally without installin the package.
#   v 0.O.1 - last update 10th of July  2019
#
#   MIT License
#   Copyright (c) 2019, Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home
#
#   DIRECTIONS:
#   1) If you have installed the FourierAnalysis from github or Julia registry,
#      uninstall it.
#   2) Change the `juliaCodeDir` path here below to the path
#           where the FourierAnalysis folder is located on your computer.
#   3) Under Linux, replace all '\\' with `/`
#   4) Put the cursor in this unit and hit SHIFT+CTRL+ENTER
#
#   Nota Bene: all you need is actually the 'push' lines and
#   the 'using' line. You can safely delete the rest once
#   you have identified the 'srcDir' to be used in the push command.

begin
    # change the 'juliaCodeDir' path to the folder where your projects are
    juliaCodeDir= homedir()*"\\Documents\\Code\\julia\\"

    # push!(LOAD_PATH, juliaCodeDir*"PosDefManifold"*"\\src\\")
    using Documenter, DocumenterTools, BenchmarkTools, Revise

    # add simple modules to be used
    # push!(LOAD_PATH, juliaCodeDir*"Modules")
    # using IOtxt

    projectName = "FourierAnalysis"
    srcDir      = juliaCodeDir*projectName*"\\src\\"
    docsDir     = juliaCodeDir*projectName*"\\docs\\"
    push!(LOAD_PATH, srcDir)
    using FourierAnalysis

    # for compiling the documentation
    cd(docsDir)
    clipboard("""makedocs(sitename="$projectName", modules=[$projectName])""")
    @info("\nhit CTRL+V+ENTER on the REPL for building the documentation.");
end
