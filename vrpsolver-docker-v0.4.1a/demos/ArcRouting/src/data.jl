import Unicode

mutable struct Vertex
   id_vertex::Int
   pos_x::Int
   pos_y::Int
   demand::Int
end

# Undirected graph
mutable struct InputGraph
   V′::Array{Int} # set of vertices 
   E::Array{Tuple{Int64,Int64}} # set of edges
   S::Array{Tuple{Int64,Int64}} # set of required edges
   cost::Dict{Tuple{Int64,Int64},Float64} # cost for each edge
   demand::Dict{Tuple{Int64,Int64},Float64} # demand for each edge (required when demand > 0)
   C::Dict{Tuple{Int64,Int64},Float64} # C[i,j]: shortest path cost from i to j
   D::Dict{Tuple{Int64,Int64},Array{Int}} # D[i,j]: shortest path from i to j
end

mutable struct DataCARP
   G′::InputGraph
   dimension::Int
   Q::Float64 # vehicle capacity
end

function getpath(next,u,v)
   if !haskey(next,(u,v))
      return []
   end
   path = [u]
   while u != v
      u = next[u,v]
      push!(path,u)
   end
   return path
end

# Floyd-Warshall Algorithm to find the shortest paths
function cheapest_path_costs(G)
   dist = Dict()
   next = Dict()
   for u in G.V′
      for v in G.V′
         if u == v
            dist[u,v] = 0.0
            next[u,v] = v
         else 
            dist[u,v] = Inf
         end
      end
   end
   for (i,j) in G.E
      dist[i,j] = G.cost[i,j]
      dist[j,i] = G.cost[i,j]
      next[i,j] = j
      next[j,i] = i
   end
   for k in G.V′
      for i in G.V′
         for j in G.V′
            if dist[i,j] > dist[i,k] + dist[k,j]
               dist[i,j] = dist[i,k] + dist[k,j]
               next[i,j] = next[i,k]
            end
         end
      end
   end
   return dist, next
end

contains(p, s) = findnext(s, p, 1) != nothing

# Read Eglese (1992) instances
function readCARPData(app::Dict{String,Any})

   str = Unicode.normalize(read(app["instance"], String); stripcc=true)
   breaks_in = [' '; ':'; '\n']
   aux = split(str, breaks_in; limit=0, keepempty=false)

   G′ = InputGraph([],[],[],Dict(),Dict(),Dict(),Dict())
   data = DataCARP(G′, 0, 0)

   dim = 0 # instance dimension
   nb_req_edges = 0
   nb_noreq_edges = 0

   for i in 1:length(aux)
      if contains(aux[i], "VERTICES")
         dim = parse(Int, aux[i+1])  # the method parse() convert the string to Int64
         data.dimension = dim
      end
      if contains(aux[i], "ARISTAS_REQ")
         nb_req_edges = parse(Float64, aux[i+1])  # the method parse() convert the string to Int64
      end
      if contains(aux[i], "ARISTAS_NOREQ")
         nb_noreq_edges = parse(Float64, aux[i+1])  # the method parse() convert the string to Int64
      end
      if contains(aux[i], "CAPACIDAD")
         data.Q = parse(Float64, aux[i+1])  # the method parse() convert the string to Int64
         break
      end
   end

   V′ = Set()
   ## getting edge/depot information
   for i in 1:length(aux)
      if contains(aux[i], "LISTA_ARISTAS_REQ")
         j, k = i, 0 
         while(k < nb_req_edges)
            v1 = parse(Int, split(aux[j+2],",")[1])-1
            v2 = parse(Int, split(aux[j+3],")")[1])-1
            cost = parse(Float64, aux[j+5])
            demand = parse(Float64, aux[j+7])
            push!(V′, v1)
            push!(V′, v2)
            e = (v1,v2)
            push!(G′.E, e) # add edge e
            data.G′.cost[e] = cost
            data.G′.demand[e] = demand
            if demand > 0
               push!(G′.S, e) # this is a required edge
            end
            j += 7 
            k += 1
         end
      end
      if contains(aux[i], "LISTA_ARISTAS_NOREQ")
         j, k = i, 0 
         while(k < nb_noreq_edges)
            v1 = parse(Int, split(aux[j+2],",")[1])-1
            v2 = parse(Int, split(aux[j+3],")")[1])-1
            cost = parse(Float64, aux[j+5])
            push!(V′, v1)
            push!(V′, v2)
            e = (v1,v2)
            push!(G′.E, e) # add edge e
            data.G′.cost[e] = cost
            data.G′.demand[e] = 0
            j += 5 
            k += 1
         end
         break
      end
   end

   data.G′.V′ = sort!(collect(V′)) # convert Set to Array
   # calculates C[i,j] (cheapest path from i to j)
   data.G′.C, next = cheapest_path_costs(data.G′)
   
   for i in data.G′.V′
      for j in data.G′.V′
         data.G′.D[i,j] = getpath(next, i, j)
      end
   end

   return data
end

required_edges(data::DataCARP) = data.G′.S # return set of edges
o(r, w) = (r[1] == w) ? r[2] : r[1]
d(data::DataCARP, e) = data.G′.demand[e] # return demand of i
veh_capacity(data::DataCARP) = data.Q
vertices(data::DataCARP) = data.G′.V′

function lowerBoundNbVehicles(data::DataCARP) 
   sum_demand = 0
   for e in data.G′.S
      sum_demand += data.G′.demand[e]
   end
   return trunc(Int,ceil(sum_demand/data.Q))
end

