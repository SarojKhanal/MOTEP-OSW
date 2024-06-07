mutable struct ThermalGenerator
  index::Int
  type::String
  location_node::Int
  Pg::Float64
  Qg::Float64
  Qmax::Float64
  Qmin::Float64
  Vg::Float64
  mBase::Float64
  status::Int
  Pmax::Float64
  Pmin::Float64
  c2::Float64
  c1::Float64
  c0::Float64
  SUcost::Float64
  SDcost::Float64
  RUrate::Float64
  RDrate::Float64
  UPtime::Float64
  DNtime::Float64
  Zone::String
  ave_marg_damages_isrm_LePeule::Float64
  CO2mtonperMWh::Float64
  SO2mtonperMWh::Float64
  NOxmtonperMWh::Float64
  PM25mtonperMWh::Float64
  PlannedRetirementYear::Int
  function ThermalGenerator(index, type, location_node, Pg, Qg, Qmax, Qmin, Vg, mBase, status, Pmax, Pmin, c2, c1, c0, SUcost, SDcost, RUrate, RDrate, UPtime, DNtime, Zone, ave_marg_damages_isrm_LePeule, CO2mtonperMWh, SO2mtonperMWh, NOxmtonperMWh, PM25mtonperMWh, PlannedRetirementYear)
     g = new()
     g.index = index
     g.type = type
     g.location_node = location_node
     g.Pg = Pg
     g.Qg = Qg
     g.Qmax = Qmax
     g.Qmin = Qmin
     g.Vg = Vg
     g.mBase = mBase
     g.status = status
     g.Pmax = Pmax
     g.Pmin = Pmin
     g.c2 = c2
     g.c1 = c1
     g.c0 = c0
     g.SUcost = SUcost
     g.SDcost = SDcost
     g.RUrate = RUrate
     g.RDrate = RDrate
     g.UPtime = UPtime
     g.DNtime = DNtime
     g.Zone = Zone
     g.ave_marg_damages_isrm_LePeule = ave_marg_damages_isrm_LePeule
     g.CO2mtonperMWh = CO2mtonperMWh
     g.SO2mtonperMWh = SO2mtonperMWh
     g.NOxmtonperMWh = NOxmtonperMWh
     g.PM25mtonperMWh = PM25mtonperMWh
     g.PlannedRetirementYear = PlannedRetirementYear
     return g
  end
end

mutable struct Bus
  index::Int
  node::Int
  Zone::String
  is_slack::Bool
  Pd::Float64
  Qd::Float64
  Vmax::Float64
  Vmin::Float64
  Gs::Float64
  Bs::Float64
  baseKV::Float64  
  inlist::Vector{Int}
  outlist::Vector{Int}
  generator::Vector{Int}
  function Bus(index, node, Zone, is_slack, Pd, Qd, Vmax, Vmin, Gs, Bs, baseKV)
     b = new()
     b.index = index
     b.node = node
     b.Zone = Zone
     b.is_slack = is_slack
     b.Pd = Pd
     b.Qd = Qd
     b.Vmax = Vmax
     b.Vmin = Vmin
     b.Gs = Gs
     b.Bs = Bs
     b.baseKV = baseKV
     b.inlist = Int[]
     b.outlist = Int[]
     b.generator = Int[]
     return b
  end
end

mutable struct Line
  index::Int
  from_node::Int
  to_node::Int
  r::Float64
  x::Float64
  b::Float64
  s_max::Float64
  function Line(index, from_node, to_node, r, x, b, s_max)
     l = new()
     l.index = index
     l.from_node = from_node
     l.to_node = to_node
     l.r = r
     l.x = x
     l.b = b
     l.s_max = s_max
     return l
  end
end

mutable struct Topology
 buses::Array{Bus}
 lines::Array{Line}
 generators::Array{ThermalGenerator}
 n_buses::Int
 n_lines::Int
 n_generators::Int
 slack_bus::Int
 function Topology(buses, lines, generators)
   f  = new()
   f.buses = buses
   f.lines = lines
   f.generators = generators
   f.n_buses = length(buses)
   f.n_lines = length(lines)
   f.n_generators = length(generators)
   for (i,b) in enumerate(buses)
     if b.is_slack
       f.slack_bus = i
       break
     end
   end
   return f
 end
end

function load_network(datadir)
# READ RAW DATA
println("Reading raw data from $(datadir)")

 nodes_raw = CSV.read(joinpath(datadir, "nodes.csv"), DataFrame)
 sum(nonunique(nodes_raw, :index)) != 0 ? @warn("Ambiguous Node Indices") : nothing

 buses = []
 for n in 1:nrow(nodes_raw)
     index = nodes_raw[n, :index]
     node = nodes_raw[n, :node]
     Zone = nodes_raw[n, :Zone]
     is_slack = nodes_raw[n, :is_slack]
     Pd = nodes_raw[n, :Pd]
     Qd = nodes_raw[n, :Qd]
     Vmax = nodes_raw[n, :Vmax]
     Vmin = nodes_raw[n, :Vmin]
     Gs = nodes_raw[n, :Gs]
     Bs = nodes_raw[n, :Bs]
     baseKV = nodes_raw[n, :baseKV]
     newb = Bus(index, node, Zone, is_slack, Pd, Qd, Vmax, Vmin, Gs, Bs, baseKV)
     push!(buses, newb)
 end

 generators_raw = CSV.read(joinpath(datadir,"generators_loc_emission_damages.csv"), DataFrame) # EIA-derived generator data by Christoph: ../generators_clustered_isone.csv # The ISO-NE Test System generator data: generators.csv
 sum(nonunique(generators_raw, :index)) != 0 ? @warn("Ambiguous Generator Indices") : nothing

 generators_raw[:,:Pg] .= 0
 generators_raw[:,:Qg] .= 0
 generators_raw[:,:Qmax] .= 0
 generators_raw[:,:Qmin] .= 0
 generators_raw[:,:Vg] .= 0
 generators_raw[:,:mBase] .= 1
 generators_raw[:,:status] .= 0

 generators = []
 for g in 1:nrow(generators_raw)
     index = generators_raw[g, :index]
     type = generators_raw[g, :type]
     location_node = generators_raw[g, :location_node]
     Pg = generators_raw[g, :Pg]
     Qg = generators_raw[g, :Qg]
     Qmax = generators_raw[g, :Qmax]
     Qmin = generators_raw[g, :Qmin]
     Vg = generators_raw[g, :Vg]
     mBase = generators_raw[g, :mBase]
     status = generators_raw[g, :status]
     Pmax = generators_raw[g, :Pmax]
     Pmin = generators_raw[g, :Pmin]
     c2 = generators_raw[g, :c2]
     c1 = generators_raw[g, :c1]
     c0 = generators_raw[g, :c0]
     SUcost = generators_raw[g, :SUcost]
     SDcost = generators_raw[g, :SDcost]
     RUrate = generators_raw[g, :RUrate]
     RDrate = generators_raw[g, :RDrate]
     UPtime = generators_raw[g, :UPtime]
     DNtime = generators_raw[g, :DNtime]
     Zone = generators_raw[g, :Zone]
     ave_marg_damages_isrm_LePeule = generators_raw[g, :ave_marg_damages_isrm_LePeule]
     CO2mtonperMWh = generators_raw[g, :CO2mtonperMWh]
     SO2mtonperMWh = generators_raw[g, :SO2mtonperMWh]
     NOxmtonperMWh = generators_raw[g, :NOxmtonperMWh]
     PM25mtonperMWh = generators_raw[g, :PM25mtonperMWh]
     PlannedRetirementYear = generators_raw[g, :PlannedRetirementYear]
     newg = ThermalGenerator(index, type, location_node, Pg, Qg, Qmax, Qmin, Vg, mBase, status, Pmax, Pmin, c2, c1, c0, SUcost, SDcost, RUrate, RDrate, UPtime, DNtime, Zone, ave_marg_damages_isrm_LePeule, CO2mtonperMWh, SO2mtonperMWh, NOxmtonperMWh, PM25mtonperMWh, PlannedRetirementYear)
     for n in 1:nrow(nodes_raw)
         if buses[n].node==newg.location_node
             push!(buses[n].generator, newg.index)
         end
     end
     push!(generators, newg)
 end

 lines_raw = CSV.read(joinpath(datadir,"lines.csv"), DataFrame)
 sum(nonunique(lines_raw, :index)) != 0  ? @warn("Ambiguous Line Indices") : nothing

 lines = []
 for l in 1:nrow(lines_raw)
     index = lines_raw[l, :index]
     from_node = lines_raw[l, :from_node]
     to_node = lines_raw[l, :to_node]
     r = lines_raw[l, :r]
     x = lines_raw[l, :x]
     b = lines_raw[l, :b]
     s_max = lines_raw[l, :s_max]
     newl = Line(index, from_node, to_node, r, x, b, s_max)
     for n in 1:nrow(nodes_raw)
         if buses[n].node == newl.from_node
             push!(buses[n].outlist, newl.index)
         elseif buses[n].node == newl.to_node
             push!(buses[n].inlist, newl.index)
         end
     end
     push!(lines, newl)
 end

 network = Topology(buses, lines, generators)

 return network
end

