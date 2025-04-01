#=
*
*              jlEGADS --- Julia version of EGADS API
*
*
*      Copyright 2011-2025, Massachusetts Institute of Technology
*      Licensed under The GNU Lesser General Public License, version 2.1
*      See http://www.opensource.org/licenses/lgpl-2.1.php
*      Written by: Julia Docampo Sanchez
*      Email: julia.docampo@bsc.es
*
*
=#

__precompile__()
module egadslite

function __init__()
    if "egads" in string.(Base.loaded_modules_array())
        throw(ErrorException("egads has already been loaded! egads and egadslite cannot be loaded at the same time." ))
    end
end

using Libdl
# Dependencies
const EGADS_COMMON = joinpath(@__DIR__, "..", "..", "egadscommon", "egadscommon.jl")
if !isfile(EGADS_COMMON)
    throw(ErrorException("egadscommon.jl file not found !! got $EGADS_COMMON"))
end

include(EGADS_COMMON)

#Set library path

const egadsLIB = Base.Sys.isapple()   ? "libegadslite.dylib" :
                 Base.Sys.islinux()   ? "libegadslite.so"    :
                 Base.Sys.iswindows() ? "egadslite.dll"      :
                 error("Unsupported operating system!")

C_egadslib = Libdl.find_library(egadsLIB, [ESP_ROOT_LIB])
if isempty(C_egadslib)
    throw(ErrorException(" EGADSLITE library $egadsLIB not found !!"))
end


"""
Imports a model from a bytestream

Returns
-------

stream: the model as a Ego
"""
function importModel(context::Context,stream::Vector{UInt8})
    ptr = Ref{ego}()
    nbytes = length(stream)
    ccall((:EG_importModel, C_egadslib), Cint, (ego, Csize_t, Ptr{Cchar}, Ptr{ego}), context.ego, nbytes, stream, ptr)
    return Ego(ptr[] ; ctxt =  context)
end

end
