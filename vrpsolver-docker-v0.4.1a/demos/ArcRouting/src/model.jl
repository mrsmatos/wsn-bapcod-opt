# cost for an arc of the VrpGraph
function c(data, a, v⁻¹)
   n = data.dimension-1
   w = v⁻¹[a[1]][2]
   r2,z = v⁻¹[a[2]]
   i,j = w < o(r2,z) ? (w, o(r2,z)) : (o(r2,z), w) # C[i,j] works when i < j
   i = (i == n+1) ? 0 : i # costs from n+1 (0') are equal those from 0 (depot)
   j = (j == n+1) ? 0 : j # costs to n+1 (0') are equal those to 0 (depot)
   if r2 == (0,n+1) # dummy required edge?
      return data.G′.C[i,j] # c_r2 = 0
   else
      return data.G′.C[i,j] + data.G′.cost[r2] # c_a = C(w,o(r2,z)) + c_r2
   end
end

# return incoming arcs at {i,j} (i,j ∈ V)
function δ⁻(A, i, j)
   in_arcs = Vector{Tuple}()
   for a in A
      if a[2] == i || a[2] == j
         push!(in_arcs, a)
      end
   end
   return in_arcs
end

# distance between required edges
function distance_req_edges(data, r1, r2)
    w1,w2 = r1
    z1,z2 = r2
    return minimum([data.G′.C[w1,z1], data.G′.C[w1,z2], data.G′.C[w2,z1], data.G′.C[w2,z2]])
end

function agg_ori_vertex_exp(data, n, v⁻¹, A, x)
   V′ = vertices(data)
   exps = [AffExpr() for i in V′]
   for a in A
      r1,w1 = v⁻¹[a[1]]
      r2,w2 = v⁻¹[a[2]]
      w1 = (w1 == n+1) ? 0 : w1
      w2 = (w2 == n+1) ? 0 : w2
      ori_path = copy(data.G′.D[ w1, (o(r2,w2)==n+1) ? 0 : o(r2,w2) ])
      push!(ori_path, w2)
      coeff = 0.0
      for (k,j) in enumerate(ori_path)
         if (k == 1 || k == length(ori_path)) #&& j == i
            coeff = 0.5
         else#if j == i
            coeff = 1.0
         end
         push!(exps[j+1], coeff, x[a])
      end
   end
   return exps
end

function agg_cvrp_edge_exp(data, S_0, x, r1, r2, v)
   exp = AffExpr()
   id1, id2 = findall(i->i==r1, S_0)[1], findall(i->i==r2, S_0)[1]
   if id1 >= id2 # pair already been treated (or skipped when equal)
      return exp
   end
   w1,w2 = r1
   z1,z2 = r2
   A′ = [] # set of 8 arcs between r1 and r2
   append!(A′,[(v[r1,w1],v[r2,z1]), (v[r1,w1],v[r2,z2]), (v[r1,w2],v[r2,z1]), (v[r1,w2],v[r2,z2])])
   append!(A′,[(v[r2,z1],v[r1,w1]), (v[r2,z2],v[r1,w1]), (v[r2,z1],v[r1,w2]), (v[r2,z2],v[r1,w2])])
   for a in A′
      push!(exp, 1.0, x[a])
   end
   return exp
end

function build_model(data::DataCARP, app::Dict{String,Any})

   V′ = vertices(data)
   S = required_edges(data)
   n = data.dimension-1
   r_0 = (0,n+1) # (0,0') (dummy required edge), 0'=n+1
   S_0 = vcat(S, r_0) # S ∪ r_0
   Q = veh_capacity(data)

   v = Dict() # access v[r,w] (vertex id in the VrpGraph)
   v⁻¹ = Dict() # reverse mapping from id to (r,w)
   vId = 0 # current id available for a new vertex

   V = Array{Int,1}() # set of vertices of the VrpGraph
   for r in S_0
      for w in r
         push!(V,vId)
         v[r,w] = vId
         v⁻¹[vId] = (r,w)
         vId += 1
      end
   end

   A = [] # set of arcs of the VrpGraph
   for r1 in S_0, r2 in S_0
      if r1 != r2
         w1,w2 = r1
         z1,z2 = r2
         # create 4 arcs from r1 nodes to r2 nodes
         append!(A, [(v[r1,w1],v[r2,z1]), (v[r1,w1],v[r2,z2]), (v[r1,w2],v[r2,z1]), (v[r1,w2],v[r2,z2])])
      end
   end

   # Formulation
   carp = VrpModel()
   @variable(carp.formulation, x[a in A], Int)
   @objective(carp.formulation, Min, sum(c(data, a, v⁻¹) * x[a] for a in A))
   for r in S
      w1, w2 = r
      @constraint(carp.formulation, sum(x[a] for a in δ⁻(A, v[r,w1], v[r,w2])) == 1.0)
   end

   #println(carp.formulation)

   # Build the model directed graph G=(V,A)
   function buildgraph(V, A)

      v_source, v_sink = v[r_0,0], v[r_0,n+1]

      L = max(app["minr"], lowerBoundNbVehicles(data))
      U = min(app["maxr"], length(S))

      G = VrpGraph(carp, V, v_source, v_sink, (L, U))
      cap_res_id = add_resource!(G, main = true) # R = Rᴹ = {cap_res_id}
      for i in V
         l_i, u_i = 0.0, Float64(Q) # accumulated resource consumption interval [l_i, u_i] for the vertex i
         set_resource_bounds!(G, i, cap_res_id, l_i, u_i)
      end

      for (i,j) in A
         arc_id = add_arc!(G, i, j)
         add_arc_var_mapping!(G, arc_id, x[(i,j)])
         r2 = v⁻¹[j][1]
         if r2 == (0,n+1) # dummy required edge
            set_arc_consumption!(G, arc_id, cap_res_id, 0.0)
         else
            set_arc_consumption!(G, arc_id, cap_res_id, d(data, r2))
         end
      end

      return G
   end

   G = buildgraph(V,A)
   add_graph!(carp, G)
   #println(G)

   set_vertex_packing_sets!(carp, [ [ (G,v[r,r[1]]), (G,v[r,r[2]]) ] for r in S])
   add_capacity_cut_separator!(carp, [ ( [ (G,v[r,r[1]]), (G,v[r,r[2]]) ], d(data,r) ) for r in S], Q)

   define_elementarity_sets_distance_matrix!(carp, G, [ [distance_req_edges(data, r1, r2) for r1 in S] for r2 in S])

   exps = agg_ori_vertex_exp(data, n, v⁻¹, A, x)
   @expression(carp.formulation, agg_ori_vertex[i in V′], exps[i+1])
   @expression(carp.formulation, agg_cvrp_edge[r1 in S_0, r2 in S_0], agg_cvrp_edge_exp(data, S_0, x, r1, r2, v))

   set_branching_priority!(carp, agg_ori_vertex, "agg_ori_vertex", 3) # higher priority
   set_branching_priority!(carp, agg_cvrp_edge, "agg_cvrp_edge", 2)
   set_branching_priority!(carp, "x", 1) # lower priority

   function odcut_callback()
      level(1) && println("Called user-defined odcut callback")

      # get the values of the path variables and aggregate by unordered pair
      # of source-sink
      n = length(V′)
      x_ = zeros(Float64, n+1, n+1)
      for a in A
         x_val = get_value(carp.optimizer, x[a])
         if x_val >= 0.000001
            r1,w1 = v⁻¹[a[1]]
            r2,w2 = v⁻¹[a[2]]
            i = (w1 == 0) ? n+1 : w1
            j = (o(r2,w2) == 0) ? n+1 : o(r2,w2)
            if i > j
               k = i
               i = j
               j = k
            end
            x_[i, j] += x_val
         end
      end

      # build the vectors that represent the path variable values for the
      # separation routine
      m = 0
      _ps = Vector{Int}()
      _pt = Vector{Int}()
      _pv = Vector{Float64}()
      for i = 1:n+1
         for j = i+1:n+1
            if x_[i, j] >= 0.000001
               push!(_ps, i)
               push!(_pt, j)
               push!(_pv, x_[i, j])
               m += 1
            end
         end
      end

      # build the vectors that represent the required edges for the
      # separation routine
      _ri = Vector{Int}()
      _rj = Vector{Int}()
      for re in S
         push!(_ri, (re[1] == 0) ? n+1 : re[1])
         push!(_rj, (re[2] == 0) ? n+1 : re[2])
      end
      _r = length(S)

      # call for the separation routine
      ncuts, cutSets = odcs_separate(n+1, m, _r, _ps, _pt, _pv, _ri, _rj, 50, 0.001)

      for c = 1:ncuts
         bits = zeros(Bool, n+1)
         for i in cutSets[c]
            bits[i] = true
         end
         arcs = Vector{Any}()
         coefs = Vector{Float64}()
         for a in A
            r1,w1 = v⁻¹[a[1]]
            r2,w2 = v⁻¹[a[2]]
            i = (w1 == 0) ? n+1 : w1
            j = (o(r2,w2) == 0) ? n+1 : o(r2,w2)
            if (bits[i] != bits[j])
               push!(arcs, a)
               push!(coefs, 1.0)
            end
         end
         add_dynamic_constr!(carp.optimizer, [x[a] for a in arcs], coefs, >=, 1.0, "odcut")
      end
   end
   add_cut_callback!(carp, odcut_callback, "odcut")

   return (carp, x, A, v⁻¹, v[r_0,0], v[r_0,n+1])
end
