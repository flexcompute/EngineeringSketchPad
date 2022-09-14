 	jlEGADS: A Julia package wrapping EGADS functions
	_________________________________________________

		Author:  Julia Docampo Sanchez
		email:   julia.docampo@bsc.es
		Version: 1.0
	_________________________________________________


0. Note that this currently does not work with Windows


1. Install Julia 1.6.4 or newer

   OPT 1) pre-built: https://julialang.org/downloads/

   OPT 2) cloning directly
	$git clone git://github.com/JuliaLang/julia.git
	$cd julia
	$git checkout v1.6.x
	$make

    !! Add alias path_to_julia/julia/julia to ~/.bash_aliases (or somewhere)


2. Set up an environment variable jlEGADS in ESPenv.sh

   export jlEGADS=$ESP_ROOT/contributions/jlEGADS


3. Load & test the egads Package

   (3.1) open julia session inside your jlEGADS folder
	$cd jlEGADS
  	$julia
	NOTE: if you are in an editor, make sure the REPL is indeed in jlEGADS:
		julia> cd(ENV["jlEGADS"]

   (3.2) Enter Pkg environment

	julia> ]
	# you should now see something like @(v1.6) pkg>

   (3.3) Activate egads package

	(@v1.6) pkg> activate .   # The dot (.) tells Julia to look in local folder
	# you should see that (@v1.6) pkg> became (egads) pkg>

   	(3.3.1)	Julia version < 1.6
        	(egads) pkg> update     # Compile egads assuming header files prebuilt

	(3.3.2)	Julia version 1.6.x or higher
 		(egads) pkg> build	# Run the C-wrap and update the header files

	*** ERRORS ***
	1) If you had Clang v0.13 already installed. Skip build and run pkg> activate
	2) If errors persist, run pkg> resolve ; pkg> update

   (3.5) Test the package
	(egads) pkg>  test
	#If all went well, the screen should show 32 test PASSED


--------------------------------------------------------------------------
EXTRA: if you are new to Julia, here is a quick guide to generate scripts.
--------------------------------------------------------------------------

(1) Create a file: vi myfile.jl

   push!(LOAD_PATH, ".")
   push!(LOAD_PATH, ENV["jlEGADS"])	# tell Julia where to find egads
   using egads
   minR, maxR, occ = egads.revision()
   @info "********* EGADS version $maxR.$minR OCC $occ *********"
   context = egads.Ccontext()
   node = egads.makeTopology!(context, egads.NODE ; reals = [0.0,0.0,0.0])
   info = egads.getInfo(node)
   @info " NODE IS OCLASS ", info.oclass, " MTYPE ", info.mtype
   @info " OK. Bye"

(2) save & quit!

(3) Run the script:

   OPT 1) open Julia session
	$julia
  	julia> include("myfile.jl")

   OPT 2) directly on terminal. NOTE: Julia will exit after running the script !
	$ julia myfile.jl

  ** The advantage of leaving the terminal open (Opt 1) is that successive runs on the script
     are faster as it doesn't recompile anything. You can check this looking at what happens
     if you do

     $julia
     julia> include("myfile.jl")   #1st time: slower
     julia> include("myfile.jl")   #2nd time: faster

  ** The Revise package is useful for developing code. If you re-run the script in an open
     session, it won't recompile untouched modules. You just need to add

       using Revise

     at the top of the script. The 1st time you use it, you need to add it to Julia.

      julia> ]
      (@v1.6) pkg> add Revise
