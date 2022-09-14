push!(LOAD_PATH, ".")

using Libdl
# SET GLOBAL PATHS TO DEPS. CHECK ESP_ROOT and include files exist !!
const ESP_ROOT        = ENV["ESP_ROOT"]
const ESP_ROOT_INC    = ESP_ROOT*"/include/"
const ESP_ROOT_LIB    = ESP_ROOT*"/lib/"
(!isdir(ESP_ROOT))     && error(" ESP_PATH ",ESP_ROOT," not fund. Have you installed ESP ??\n")
(!isdir(ESP_ROOT_INC)) && error(" Wrong ESP INCLUDE PATH ", ESP_ROOT_INC," not found !!\n")
(!isdir(ESP_ROOT_INC)) && error(" Wrong ESP LIBRARY PATH ",ESP_ROOT_LIB ," not found !!\n")

const FILE_EGADS_FUNS = "egads_Cwrap.jl"
const FILE_EGADS_CTT  = "egads_types_ctts.jl"
const EGADS_WRAP_PATH = normpath(dirname(@__FILE__) )

# Set egads library name. WARNING ! If you change C_egadslib name, all calls inside egads will fail.
const LIBEGADS        = "C_egadslib"

#Set library path
const egadsLIB = Base.Sys.isapple() ? "libegads.dylib" : Base.Sys.islinux() ?
                                      "libegads.so"    : Base.Sys.iswindows() ?
                                      "egads.dll"      : error("Unsupported operating system") 
C_egadslib = Libdl.find_library(egadsLIB, [ESP_ROOT_LIB])
isempty(C_egadslib) && error("EGADS library not found, update the location in egads.jl or update commonlocesp")
