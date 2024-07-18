function build_model(data::DataHVRP)

   E = edges(data) # set of edges of the input graph G′
   n = nb_customers(data)
   V = [i for i in 0:n] # set of vertices of the graphs G′ and G
   V⁺ = [i for i in 1:n] # set of customers of the input graph G′
   K = veh_types(data)

   # Formulation
   hvrp = VrpModel()
   @variable(hvrp.formulation, x[e in E, k in K], Int)
   @objective(hvrp.formulation, Min, sum(c(data, e, k) * x[e, k] for e in E, k in K))
   @constraint(hvrp.formulation, deg[i in V⁺], sum(x[e, k] for e in δ(data, i), k in K) == 2.0)

   #println(hvrp.formulation)

   # Build the model directed graph G=(V,A)
   function build_graph(k)

      v_source = v_sink = 0
      L, U = lowerBoundNbVehicles(data, k), upperBoundNbVehicles(data, k)

      # node ids of G from 0 to n
      G = VrpGraph(hvrp, V, v_source, v_sink, (L, U))
      cap_res_id = add_resource!(G, main = true) # R = Rᴹ = {cap_res_id}
      for i in V
         l_i, u_i = 0.0, Float64(veh_capacity(data, k)) # accumulated resource consumption interval [l_i, u_i] for the vertex i
         set_resource_bounds!(G, i, cap_res_id, l_i, u_i)
      end

      # Build set of arcs A from E (two arcs for each edge (i,j) assume that i < j)
      for (i,j) in E
         arc_id = add_arc!(G, i, j)
         add_arc_var_mapping!(G, arc_id, x[(i,j), k])
         set_arc_consumption!(G, arc_id, cap_res_id, (d(data, i) + d(data, j))/2)

         arc_id = add_arc!(G, j, i)
         add_arc_var_mapping!(G, arc_id, x[(i,j), k])
         set_arc_consumption!(G, arc_id, cap_res_id, (d(data, i) + d(data, j))/2)
      end

      return G
   end

   graphs = []
   for k in K
      G = build_graph(k)
      add_graph!(hvrp, G)
      #println(G)
      push!(graphs, G)
   end

   set_vertex_packing_sets!(hvrp, [[(G,i) for G in graphs] for i in V⁺])
   for G in graphs
      define_elementarity_sets_distance_matrix!(hvrp, G, [[distance(data, i, j) for j in V⁺] for i in V⁺])
   end
 
   add_capacity_cut_separator!(hvrp, [ ( [(G,i) for G in graphs], d(data, i)) for i in V⁺], maximum([veh_capacity(data, k) for k in K]))
 
   @expression(hvrp.formulation, edge[(i,j) in E], sum(x[(i,j), k] for k in K))
   set_branching_priority!(hvrp, edge, "edge", 1)

   @expression(hvrp.formulation, vehNum[k in K], sum(0.5 * x[(0,i), k] for i in V⁺))
   set_branching_priority!(hvrp, vehNum, "vehNum", 3)

   @expression(hvrp.formulation, assign[k in K, i in V⁺], sum(0.5 * x[e, k] for e in δ(data, i)))
   set_branching_priority!(hvrp, assign, "assign", 1)

   return (hvrp, x)
end
