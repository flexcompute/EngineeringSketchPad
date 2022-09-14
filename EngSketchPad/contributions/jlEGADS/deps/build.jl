using Clang.Generators, Clang.LibClang.Clang_jll

# File containing shared information with egads Module about library name and constants
@info " Setting global variables "
include("sharedFile.jl")

# STORE ALL .h FILES TO CONVERT C -> Julia
aux = [joinpath(ESP_ROOT_INC, header) for header in readdir(ESP_ROOT_INC)
       if endswith(header, ".h") && startswith(header,"egads")]
excluded_h = joinpath.(ESP_ROOT_INC,["egadsSplineFit.h", "egadsSplineVels.h"])
piv = findall(x -> x in excluded_h ,aux)
const cHeaders = deleteat!(aux, piv)


# Clang : set output files and create wrap
include_dir = joinpath(Clang_jll.artifact_dir, "include") |> normpath
clang_dir   = joinpath(include_dir, "clang-c")
options     = Dict("general" =>Dict("module_name" => "LibEgads",
               "library_name" => LIBEGADS,
               "jll_pkg_name" => "Clang_jll",
               "export_symbol_prefixes" => ["CX", "clang_"],
               "output_api_file_path" => "./"*FILE_EGADS_FUNS,
                "output_common_file_path" => "./"*FILE_EGADS_CTT))

# add compiler flags, e.g. "-DXXXXXXXXX"
args = get_default_args()
push!(args, "-I$include_dir")

# create context
ctx = create_context(cHeaders, args, options)

# run generator
build!(ctx)


# Now convert all C const : Int64 -> Int32
(tmppath, tmpio) = mktemp()
open(FILE_EGADS_CTT) do io
    for line in eachline(io, keep=true) # keep so the new line isn't chomped
        if occursin("EGADSPROP", line)
   	    s_beg = strip(split(line,"EGADSprop")[1])
            s_end = line[length(s_beg)+1:end-1]
            line = s_beg*"\" "*String(s_end)*"\" "
        elseif occursin("const",line)
            s_end = strip(split(line," ")[end])
            s_beg = line[1:end-length(s_end)-1]
            if occursin("e",s_end)
                write(tmpio, line)
                continue
            end
            isa(parse(Int, s_end), Int) && (line = s_beg*"Int32("*s_end*")")
        end
        write(tmpio, line)
    end
end
close(tmpio)
mv(tmppath, FILE_EGADS_CTT, force=true)
