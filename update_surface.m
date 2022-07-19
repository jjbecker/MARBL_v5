function [surface] = update_surface (dt, surface, interior, skip )

% use midpoint of time step to update surface?
% 
% FIXME: Midpoint update? How do we account for mixing and tendency?
% FIXME: further more update_interior DOES update the state, but NOT tracer.
% 
% surface.state  = (interior.state_old (:,1)  +interior.state (:,1) ) ./2;
% 
% FIXME: Use half of dt*tendency? Not perfect integration, and not doing 
% any not mixing...
% 
% delta_tracer   = dt *interior.tendency(:,1);
% new_tracer     = interior.tracer(:,1) +delta_tracer;
% surface.tracer = (interior.tracer_old(:,1) +new_tracer)' ./2;

surface.state  = (interior.state_old (:,1));
surface.tracer = (interior.tracer_old(:,1));

mex_marbl_driver ( 'restore_surface_flux_forcings',         surface.forcing );

% FIXME: FLUX limit?
mex_marbl_driver ( 'restore_tracers_at_surface',            surface.tracer  );
% mex_marbl_driver ( 'restore_tracers_at_surface',   max(eps, surface.tracer(:)));
% mex_marbl_driver ( 'restore_tracers_at_surface',     max(0, surface.tracer(:)));

% intial condition for calculation, not critical, just use previous value

mex_marbl_driver ( 'restore_surface_flux_saved_state', surface.state   ); 

% sign of flux is such that negative values are OUT water layer.

% marbl_log = mex_marbl_driver('surface_flux_compute');
mex_marbl_driver('surface_flux_compute');

surface.flux  = mex_marbl_driver ( 'surface_fluxes' );
% FIXME: comment next line out if not using all diags at all locations. 
% This copies a huge amount of data
% surface.diag  = mex_marbl_driver ( 'surface_flux_diags' );
% surface.sfo   = mex_marbl_driver ( 'sfo' );

% Mike Levi says using this output on next input speeds code up, but to me
% this output state is just a diagnostic, it's the pH. 
% save it for use as IC next iteration

surface.state = mex_marbl_driver ( 'surface_flux_saved_state' );


% FIXME: check flux for NaN and other errors...
% FIXME:    ...record (bad) tracer and surface_flux_forcing?

% Some updates are just startup spikes, like when CISO=1 and n=1.
% To avoid lots off by one bugs, just zero flux 0, but keepeverything else
% so we have same number of samples of flux, diags, etc.

if skip
    surface.flux = 0 * surface.flux;
end

% Need to use depth to convert surface flux to volume rate of change.
%
% Units of -METERS- are used in matlab domain.
%   convert depth here to units (cm?) used in MARBL

% FIXME: Hard coded units for depth for the moment. MARBL uses cm, we use m
MARBL_depth_per_m = 100;

surface.tendency = surface.flux' /( interior.domain.dzt(1) *MARBL_depth_per_m);

end % update_surface