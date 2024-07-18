using VrpSolver, JuMP, ArgParse

include("data.jl")
include("model.jl")
include("solution.jl")

function parse_commandline(args_array::Array{String,1}, appfolder::String)
   s = ArgParseSettings(usage="##### VRPSolver #####\n\n"*
	   "  On interactive mode, call main([\"arg1\", ..., \"argn\"])", exit_after_help=false)
   @add_arg_table s begin
      "instance"
         help = "Instance file path"
      "--cfg", "-c"
         help = "Configuration file path"
         default = "$appfolder/../config/HVRP.cfg"
      "--ub","-u"
         help = "Upper bound (primal bound)"
         arg_type = Float64
         default = 10000000.0
      "--out","-o"
         help = "Path to write the solution found"
      "--brandao","-B"
         help = "Brandao instances."
         action = :store_true
      "--batch", "-b" 
         help = "batch file path" 
   end
   return parse_args(args_array, s)
end
 
function run_hvrp(app::Dict{String,Any})

   println("Application parameters:");
   for (arg,val) in app
      println("  $arg  =>  $(repr(val))")
   end
   flush(stdout)

   instance_name = split(basename(app["instance"]), ".")[1]

   if app["brandao"]
      data = readHVRPBrandaoData(app["instance"])
   else
      data = readHVRPClassicData(app["instance"])
   end
   (model, x) = build_model(data)
   optimizer = VrpOptimizer(model, app["cfg"], instance_name)
   set_cutoff!(optimizer, app["ub"])

   (status, solution_found) = optimize!(optimizer)

   println("########################################################")
   if solution_found
      sol = getsolution(data, optimizer, x, get_objective_value(optimizer)) 
      print_routes(sol)
      println("Cost: $(sol.cost)")
      if app["out"] != nothing
         writesolution(app["out"], sol)
      end
   else
      if status == :Optimal
         println("Problem infeasible")
      else
         println("Solution not found")
      end
   end
   println("########################################################")

end

function main(args)
   appfolder = dirname(@__FILE__)
   app = parse_commandline(args, appfolder)
   isnothing(app) && return
   if app["batch"] != nothing
      for line in readlines(app["batch"])
         if isempty(strip(line)) || strip(line)[1] == '#'
            continue
         end
         args_array = [String(s) for s in split(line)]
         app_line = parse_commandline(args_array, appfolder)
         run_hvrp(app_line)
      end
   else
      run_hvrp(app)
   end
end

if isempty(ARGS)
   main(["--help"])
else
   main(ARGS)
end
