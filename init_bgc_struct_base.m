function [surface_base, interior_base] = init_bgc_struct_base(sim)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

interior_base.domain.dzt    = sim.domain.dzt;
interior_base.tendency      = zeros(sim.bgc_struct_base.size.tendency)';

interior_base.tracer_name   = sim.bgc_struct_base.name.tracer;
interior_base.diag_name     = sim.bgc_struct_base.name.diag;
interior_base.forcing_name  = sim.bgc_struct_base.name.forcing;
interior_base.tracer_unit   = sim.bgc_struct_base.unit.tracer;

surface_base.tracer_name    = sim.bgc_struct_base.name.tracer;
surface_base.diag_name      = sim.bgc_struct_base.name.surf_diag;
surface_base.forcing_name   = sim.bgc_struct_base.name.surf_forcing;
surface_base.tracer_unit    = sim.bgc_struct_base.unit.tracer;

end
