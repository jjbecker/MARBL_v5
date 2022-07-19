function interior = update_interior ( interior, skip )

mex_marbl_driver ( 'restore_interior_tendency_forcings',    interior.forcing );

% state is just an intial condition for calculation, not critical, just use previous value

mex_marbl_driver ( 'set_depth',                             interior.domain.kmt );

% FIXME: FLUX limit? This should have been done in MARBL_loop_iteration()
mex_marbl_driver ( 'restore_tracers',                       interior.tracer     );
% mex_marbl_driver ( 'restore_tracers',              max(eps, interior.tracer(:)) );
% mex_marbl_driver ( 'restore_tracers',                max(0, interior.tracer)    );

mex_marbl_driver ( 'restore_interior_tendency_saved_state', interior.state      );

% compute tendency using MARBL

mex_marbl_driver('interior_tendency_compute');

% FIXME: check tendency for NaN and other errors...
% FIXME:    ...record record (bad) tracer and surface_flux_forcing

% read results

interior.tendency = mex_marbl_driver ( 'interior_tendencies' );
% FIXME: comment next line out if not using all diags at all locations. 
% This copies a huge amount of data
% interior.diag     = mex_marbl_driver ( 'interior_tendency_diags');
interior.state    = mex_marbl_driver ( 'interior_tendency_saved_state' );

% Some updates are just startup spikes, like when CISO=1 and n=1.
% To avoid lots off by one bugs, just zero flux 0, but keepeverything else
% so we have same number of samples of flux, diags, etc.

if skip
    interior.tendency = 0 * interior.tendency;
end

end % update_interior