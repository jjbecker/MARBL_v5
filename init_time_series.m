function time_series = init_time_series(sim, bgc_struct)

time_series.dt                  = sim.dt;
time_series.nstep               = sim.num_time_steps;

time_series.tracer_name         = bgc_struct.name.tracer;
time_series.tracer_unit         = bgc_struct.unit.tracer;

if (sim.logTracers)
    time_series.moles               = squeeze( zeros([bgc_struct.size.surf_flux, sim.num_time_steps   ]));
    time_series.Dmoles              = squeeze( zeros([bgc_struct.size.surf_flux, sim.num_time_steps   ]));

    time_series.tracer              = zeros([bgc_struct.size.tracer, sim.num_time_steps               ] );
    time_series.sfo                 = squeeze( zeros([bgc_struct.size.sfo, sim.num_time_steps         ])); % e.g. air-sea-flux
    % time_series.accumulate          = zeros([bgc_struct.size.tracer, sim.num_time_steps               ] );
end

if (sim.logDiags)
    time_series.forcing_name        = bgc_struct.name.forcing;
    time_series.surf_forcing_name   = bgc_struct.name.surf_forcing;
    time_series.diag_name           = bgc_struct.name.diag;
    time_series.surf_diag_name      = bgc_struct.name.surf_diag;

    time_series.diag                = zeros([bgc_struct.size.diag, sim.num_time_steps                 ] );
    time_series.surf_diag           = squeeze( zeros([bgc_struct.size.surf_diag, sim.num_time_steps   ]));
    %     time_series.tendency            = zeros([bgc_struct.size.tracer, sim.num_time_steps               ] );
end

if (sim.verbose_debug)
    disp(['Log file of one location for all time steps uses ',num2str(getMemSize(time_series)/1024/1024, '%1.1f'),' MB'])
    disp(['Log file of one location uses ',num2str(getMemSize(time_series)/1024/sim.num_time_steps, '%1.1f'),' KB per time step'])
end
end
