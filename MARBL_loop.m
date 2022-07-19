function [bgc, time_series] = MARBL_loop (n, sim, bgc, time_series)

% A little note to me about the dimensions and shapes of all these variable..
%
% interior.forcing = 6x24                 bgc.forcing     (loc,1:24,6  )'
% interior.domain.kmt = 1x1               bgc.kmt         (loc)
% interior.tracer = 46x24                 bgc.tracer      (loc,1:24,46 )'
% interior.state = 2x24                   bgc.state       (loc,1:24,2  )'
% interior.tendency = 46x24 out           bgc.tendency    (loc,1:24,46 )'
% interior.diag = 386x24 out              bgc.diag        (loc,1:24,386)'
% interior.state = 2x24 out               bgc.state       (loc,1:24,2  )'
%
% restore_surface_flux_saved_state 2x1    bgc.state       (loc,1,:)'
% restore_tracers_at_surface 1x46         bgc.tracer      (loc,1,:)
% restore_surface_flux_forcings 1x13      bgc.surf_forcing(loc,1,:)
% surface_fluxes  1x46                    bgc.surf_flux   (loc,1,:)
% sfo 1x4                                 bgc.sfo         (loc,1,:)
% surface_flux_diags 43x1                 bgc.surf_diag   (loc,1,:)'

[surface_base, interior_base] = init_bgc_struct_base(sim);
surface  = surface_base;
interior = interior_base;

% sz = size(sim.domain.M3d);

bgc.tendency = nan(size(bgc.tendency));

for row = 1: sim.domain.num_wet_loc

    % Restore depth, state, and tracer at this location in structs

    interior.domain.kmt = bgc.kmt(row);
    interior.state      = squeeze(bgc.state (row,:,:))';
    interior.tracer     = squeeze(bgc.tracer(row,:,:))';

    % update forcing of interior and surface in structs

    surface.forcing     = squeeze(bgc.surf_forcing (row,:,:))';
    %     surface.river_flux  = squeeze(bgc.river_flux   (row,:,:))';
    interior.forcing    = squeeze(bgc.forcing      (row,:,:))';

    % Feed MARBL the structs and get updated structs with the result

    [surface, interior] = MARBL_loop_iteration(sim.dt, n, row, surface, interior);

    % record tendency at all levels for time step

    bgc.state    (row,:,:)  = interior.state';
    bgc.tendency (row,:,:)  = interior.tendency';

    % mex_marbl_driver('print_sfo');


    % log a time series at just one location. Gets very large very quickly

    if ((row == sim.time_series_loc) && sim.logTracers)
        surface.sfo   = mex_marbl_driver ( 'sfo' );

        if (sim.logDiags)
            interior.diag = mex_marbl_driver ( 'interior_tendency_diags');
            surface.diag  = mex_marbl_driver ( 'surface_flux_diags' );
        end
        update_log ();

    end
end % scan over all water columns in grid

    function update_log

        if (sim.logTracers)
            time_series.tracer        (:, :, n) = interior.tracer';
            time_series.sfo           (:, n)    = surface.sfo;      % -this- is sea->air gas flux
        end
        if (sim.logDiags)
            time_series.diag          (:, :, n) = interior.diag';
            time_series.surf_diag     (:, n)    = surface.diag';
        end
    end % update_log
end % MARBL_loop
