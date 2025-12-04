# Not used (runs via adding export in the run-pluto-slider.sh script)
import PlutoSliderServer
using Base.Filesystem: mv
using Base: cp

function create_static_notebook(
    nbname::String; 
    nbpath = joinpath(@__DIR__, "..", "lecturenotes"),
    staticpath = joinpath(@__DIR__, "..", "static")
    )
    nbfile = joinpath(nbpath, nbname)
    htmlname = replace(nbname, ".jl" => ".html")
    mktempdir() do folder
        cp(nbfile, joinpath(folder, nbname))
        PlutoSliderServer.export_notebook(joinpath(folder, nbname))
        if !isfile(joinpath(folder, htmlname))
            @warn "No html-file generated: $htmlname"
            println("Files in temporary directory: ", folder)
            foreach(println, readdir(folder))
        else
            cp(joinpath(folder, htmlname), joinpath(staticpath, htmlname); force = true)
        end
    end
end

function create_static_notebooks(; nbpath = joinpath(@__DIR__, "..", "lecturenotes"), kwargs...)
    for nbname in readdir(nbpath)
        if endswith(nbname, ".jl")
            create_static_notebook(nbname; nbpath, kwargs...)
        end
    end
end

create_static_notebooks()