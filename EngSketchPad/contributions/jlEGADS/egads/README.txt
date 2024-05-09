 	jlEGADS: A Julia package wrapping EGADS functions
	_________________________________________________

		Author: Julia Docampo Sanchez
		email: julia.docampo@bsc.es
		Version: 1.0
	_________________________________________________


1. Install Julia 1.7.X or newer

   OPT 1) pre-built: https://julialang.org/downloads/

   OPT 2) cloning directly
	$git clone git://github.com/JuliaLang/julia.git
	$cd julia
	$git checkout v1.7.x
	$make

    !! Add alias path_to_julia/julia/julia to ~/.bash_aliases (or somewhere)
    OR!! Modify your PATH variable 

3. Load & test the egads module

	3.1 Using Makefile
	  
	   (1) run make all 
		This will run the setup + the test files
	
	   (2) In this case, we are adding Julia packages directly 

	3.2 Using Package manager

	   (1) open julia session inside your jlEGADS folder
		$cd jlEGADS
	  	$julia

	   (2) Enter Pkg environment

		julia> ]
		# you should now see something like @(v1.7) pkg>

	   (3) Activate egads package

		(@v1.7) pkg> activate .   # The dot (.) tells Julia to look in local folder
		# you should see that (@v1.6) pkg> became (egads) pkg>

	   	(3.1)	Julia version < 1.7
			(egads) pkg> instantiate     # Compile egads assuming header files prebuilt

		(3.2)	Julia version 1.7.x or higher
	 		(egads) pkg> build	# Run the C-wrap and update the header files

		*** ERRORS ***
		1) If you had Clang v0.13 already installed. Skip build and run pkg> activate
		2) If errors persist, run pkg> resolve ; pkg> update

	   (4) Test the package
		(egads) pkg>  test
		#If all went well, the screen should show 209 test PASSED


--------------------------------------------------------------------------
EXTRA: if you are new to Julia, here is a quick guide to generate scripts.
--------------------------------------------------------------------------

(1) Create a file: vi myfile.jl

   push!(LOAD_PATH, ".")
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
	$julia --project=/path_to_jlEGADS
  	julia> include("myfile.jl") 

   OPT 2) directly on terminal. NOTE: Julia will exit after running the script !
	$ julia --project=/path_to_jlEGADS myfile.jl

  ** The advantage of leaving the terminal open (Opt 1) is that successive runs on the script
     are faster as it doesn't recompile anything. You can check this looking at what happens
     if you do

     $julia --project=/path_to_jlEGADS
     julia> include("myfile.jl")   #1st time: slower
     julia> include("myfile.jl")   #2nd time: faster

  ** The Revise package is useful for developing code. If you re-run the script in an open
     session, it won't recompile untouched modules. You just need to add

       using Revise

     at the top of the script. The 1st time you use it, you need to add it to Julia.

      julia> ]
      (@v1.6) pkg> add Revise
