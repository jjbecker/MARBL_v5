function [bgc_struct] = init_bgc_struct(sim)


% size of one water column. Not efficient to transfer data to MARBL but
% claifies data struct code

% Create a struct that has meta data like size, names, and units.
% FIXME: possibly use this struct to hold values for a single water column
% to make it easy to create global arrays we actually use...

% Global grid size and number of water levels for all water columns.

bgc_struct.size.grd = size(sim.domain.M3d, [1 2]);
bgc_struct.size.lvl = size(sim.domain.M3d, 3);


% -- ASSUMED -- Size of MARBL data structs, cheched later at MARBL init!
%
% Inputs to MARBL
bgc_struct.size.kmt           = [1 1];
bgc_struct.size.tracer        = [bgc_struct.size.lvl, size( tracer_units( sim.lciso_on ), 2)];
% bgc_struct.size.old_tracer    = [bgc_struct.size.lvl, size( tracer_units( sim.lciso_on ), 2)];
bgc_struct.size.forcing       = [bgc_struct.size.lvl, size( interior_forcing_names(), 2)];
bgc_struct.size.surf_forcing  = size( surface_forcing_names( sim.lciso_on ));
bgc_struct.size.river_flux    = size( river_flux_names     ( sim.lciso_on ));

% FIXME: hardcoded
% Initial condition for pH, -and Output- diag for a single MARBL water column
bgc_struct.size.state         = [bgc_struct.size.lvl, 2];

% Output from MARBL for a single MARBL water column
bgc_struct.size.tendency      = bgc_struct.size.tracer;
bgc_struct.size.sfo           = size( surface_flux_output_names() );
bgc_struct.size.surf_flux     = [1, bgc_struct.size.tracer(2)];

% FIXME: hardcoded
% Diagnostics from MARBL for a single MARBL water column, no CISO
if (sim.lciso_on == 1)
    bgc_struct.size.diag      = [bgc_struct.size.lvl, 390]; % FIXME: hardcoded
    bgc_struct.size.surf_diag = [1, 43];                    % FIXME: hardcoded
else
    bgc_struct.size.diag          = [bgc_struct.size.lvl, 312]; % FIXME: hardcoded
    bgc_struct.size.surf_diag     = [1, 27];                    % FIXME: hardcoded
end

% bgc_struct.name.surf_diag     =
% bgc_struct.unit.surf_diag     =
%

% Values (aka actual data) for a single MARBL water column
% For global array, just repmat() this
bgc_struct.value.kmt          = zeros(bgc_struct.size.kmt);
bgc_struct.value.tracer       = zeros(bgc_struct.size.tracer);
% bgc_struct.value.old_tracer   = zeros(bgc_struct.size.tracer);
bgc_struct.value.forcing      = zeros(bgc_struct.size.forcing);
bgc_struct.value.surf_forcing = zeros(bgc_struct.size.surf_forcing);
bgc_struct.value.river_flux   = zeros(bgc_struct.size.river_flux);
bgc_struct.value.state        = zeros(bgc_struct.size.state);
bgc_struct.value.tendency     = zeros(bgc_struct.size.tendency);
bgc_struct.value.sfo          = zeros(bgc_struct.size.sfo);
bgc_struct.value.diag         = zeros(bgc_struct.size.diag);
bgc_struct.value.surf_diag    = zeros(bgc_struct.size.surf_diag);


% Amount of data for a single MARBL water column with all diags and outputs
if (or (sim.logDiags, sim.logTracers))
    disp(['single MARBL water column with all diags, forcing, and outputs: ', num2str(getMemSize(bgc_struct.value,1E+30)/1024,'%1.3f'), ' (KB)']);
end