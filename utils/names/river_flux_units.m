function [my_cell] = river_flux_units(lciso_on)
%tracer_names MARBL "units" of all tracers

% if MARBL is shutdown then mex_marbl_driver('tracer_sname', i) crashes...
%    just do a lot of typing...

    my_cell = strcat(tracer_units (lciso_on), '(cm/s)');
end