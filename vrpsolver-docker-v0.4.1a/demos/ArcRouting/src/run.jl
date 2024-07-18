using VrpSolver, JuMP, ArgParse

include("data.jl")
include("model.jl")
include("odcsep.jl")

function parse_commandline(args_array::Array{String,1}, appfolder::String)
   s = ArgParseSettings(usage="##### VRPSolver #####\n\n"*
	   "  On interactive mode, call main([\"arg1\", ..., \"argn\"])", exit_after_help=false)
   @add_arg_table s begin
      "instance"
         help = "Instance file path"
      "--cfg", "-c"
         help = "Configuration file path"
         default = "$appfolder/../config/CARP.cfg"
      "--ub","-u"
         help = "Upper bound (primal bound)"
         arg_type = Float64
         default = 10000000.0
      "--minr","-m"
         help = "Minimum number of routes in the solution"
         arg_type = Int
         default = 1 
      "--maxr","-M"
         help = "Maximum number of routes in the solution"
         arg_type = Int
         default = 999
      "--out","-o"
         help = "Path to write the solution found"
      "--batch","-b" 
         help = "batch file path" 
   end
   return parse_args(args_array, s)
end
 
function run_carp(app::Dict{String,Any})
   println("Application parameters:")
   for (arg,val) in app
      println("  $arg  =>  $(repr(val))")
   end
   flush(stdout)

   instance_name = split(basename(app["instance"]), ".")[1]

   data = readCARPData(app)

   (model, x, A, v⁻¹, source, sink) = build_model(data, app)
   optimizer = VrpOptimizer(model, app["cfg"], instance_name)
   set_cutoff!(optimizer, app["ub"])

   (status, solution_found) = optimize!(optimizer)
   if solution_found
      println("########################################################")
      n = data.dimension-1
      S = [a for a in A if get_value(optimizer, x[a]) > 0.5] # arcs of the solution

      # print the routes over original graph
      id_route = 1
      f = (app["out"] != nothing) ? open(app["out"], "w") : [] 
      while !isempty(S)
         print("Route #$(id_route): ")
         (app["out"] != nothing) && write(f, "Route #$(id_route): ")
         route, req_edges, aux_req = [], [], 0
         next, firstIter = source, true
         while next != sink || firstIter
            for i in 1:length(S)
               a = S[i]
               if a[1] == next
                  next = a[2]
                  r1,w1 = v⁻¹[a[1]]
                  r2,w2 = v⁻¹[a[2]]
                  w1 = (w1 == n+1) ? 0 : w1
                  w2 = (w2 == n+1) ? 0 : w2
                  # get subpath on original graph from arc a (network arc)
                  ori_path = copy(data.G′.D[ w1, (o(r2,w2)==n+1) ? 0 : o(r2,w2) ])
                  append!(route, ori_path)
                  push!(req_edges,aux_req + length(ori_path))
                  aux_req += length(ori_path)
                  deleteat!(S, i)
                  break
               end
            end
            firstIter = false
         end
         for i in 1:(length(route))
            if i in req_edges && route[i] != 0
               print("$(route[i]+1)=")
               (app["out"] != nothing) && write(f, "$(route[i]+1)=")
            elseif i > 1 && route[i] == 0
               print("$(route[i]+1)")
               (app["out"] != nothing) && write(f, "$(route[i]+1)")
            else
               print("$(route[i]+1)-")
               (app["out"] != nothing) && write(f, "$(route[i]+1)-")
            end
         end
         println()
         (app["out"] != nothing) && write(f, "\n")
         id_route += 1
      end
      println("Cost: $(get_objective_value(optimizer))")
      (app["out"] != nothing) && write(f, "Cost: $(get_objective_value(optimizer))\n")
      (app["out"] != nothing) && close(f)
      println("########################################################")
   elseif status == :Optimal
      println("Problem infeasible")
   else
      println("Solution not found")
   end
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
         run_carp(app_line)
      end
   else
      run_carp(app)
   end
end

if isempty(ARGS)
   main(["--help"])
else
   main(ARGS)
end
