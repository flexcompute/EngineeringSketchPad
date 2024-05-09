# Check Julia version
@static if VERSION < v"1.7" error(" jlEGADS requires Julia >= v1.7! Current Julia is ", string(VERSION), "...\n") end

# EGADS paths 
const ESP_ROOT        = ENV["ESP_ROOT"]               |> normpath
const ESP_ROOT_INC    = joinpath(ESP_ROOT, "include") |> normpath
(!isdir(ESP_ROOT))     && error(" ESP_ROOT ", ESP_ROOT, " not fund. Please fix the environment...\n")
(!isdir(ESP_ROOT_INC)) && error(" ESP INCLUDE PATH ", ESP_ROOT_INC ," not found !!\n")


# SET GLOBAL PATHS TO DEPS. CHECK ESP_ROOT and include files exist !!
const EGADS_COMMON_PATH = normpath(dirname(@__FILE__))
const FILE_EGADS_CTT    = joinpath(EGADS_COMMON_PATH,"egads_types_ctts.jl")

# Set egads library name. WARNING ! If you change C_egadslib name, all calls inside egads will fail.
const LIBEGADS        = "C_egadslib"

if !isfile(FILE_EGADS_CTT)
    using Pkg
    Pkg.add("Clang")    
    using Clang.Generators
    # STORE ALL .h FILES TO CONVERT C -> Julia
    cHeaders = [joinpath(ESP_ROOT_INC, "egadsTypes.h"),
                joinpath(ESP_ROOT_INC, "egadsErrors.h")]
               # joinpath(ESP_ROOT_INC, "egadsTris.h")]

    
    # add compiler flags, e.g. "-DXXXXXXXXX"
    args = get_default_args()
    push!(args, "-I$(ESP_ROOT_INC)")

    # dumb file for egads wrap
    EGADS_FUNS = "dumb.jl"
    # set output files and create wrap
    options = Dict("general" => Dict("module_name" => "egads",
               "library_name" => LIBEGADS,
               "export_symbol_prefixes" => ["EG_"],
               "output_api_file_path" => EGADS_FUNS,
               "output_common_file_path" => FILE_EGADS_CTT),
               "codegen.macro" => Dict("macro_mode" => "basic"))

    # create context

    ctx = create_context(cHeaders, args, options)

    # run generator
    build!(ctx)
    rm(EGADS_FUNS)

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
end