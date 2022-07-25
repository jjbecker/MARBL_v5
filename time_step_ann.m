function [sim, bgc, time_series, n] = time_step_ann (sim, bgc, time_series, n, forcing, TR, month)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%%

persistent tracerCapture tendencyCapture
if (sim.captureAllSelectedTracers == 1)
    if isempty(tendencyCapture)
        tmp = packMarbl( bgc.tracer, sim.domain.iwet_JJ );
        tmp = tmp(:,sim.selection);
        myCaptureSz = [size(tmp(:),1), sim.num_time_steps];
        tracerCapture = zeros(myCaptureSz);
        fprintf('\n%s.m: HUGE tendency log created: %g bytes\n\n',mfilename, getMemSize(tracerCapture));
        tendencyCapture = zeros(myCaptureSz);
        fprintf('\n%s.m: HUGE tracer log created: %g bytes\n\n',mfilename, getMemSize(tendencyCapture));
    end
end

%%
dt               = sim.dt;
tsperDay         = round(sim.const.sec_d/dt);
days_in_month    = [31 28 31 30 31 30 31 31 30 31 30 31];
steps_per_period = tsperDay*days_in_month(month);
%%

% DEBUG: n = n+steps_per_period; return; % quick return to debug G

% forcing, river, and transport created monthly. We get only this month's
%  in this routine

bgc.forcing      = forcing.interior;
bgc.surf_forcing = forcing.surf_forcing;

% river flux is not a tendency until divided by surface layer thickness
% *** IN CM ***
% Yes! Units are mixed m and cm because MARBL depths are cm. 

surfaceLayerThickness_MARBL = (sim.domain.dzt(1) *sim.domain.MARBL_depth_per_m);

river_tendency = forcing.river_flux /surfaceLayerThickness_MARBL;

%%

C_1 = packMarbl( bgc.tracer, sim.domain.iwet_JJ );
[C_2, C_3, C_4] = deal( zeros( size(C_1)));

%%
% Advection Adams-Bashforth method
%
% source/sink term
% S = bgc.tendency + bgc.river_flux/thickness;
%
% LHS is the implicit diffusion
% RHS is the explitc advection
%
%     LHS  = d0(1+xi(:,2)) + dt.*(-TR.D);   % -CONSTANT- same LHS on all tracers
%     FLHS = mfactor(LHS);
%     RHS = C(:,k-1) + dt.*(A*w1(:,1) + H*C(:,k-1)+ S);
%
% Solve for C: (I - dt.*D)*C(:,k) = RHS
%
%     C(:,k) = mfactor(FLHS,RHS);           % aka time stepped tracers

LHS  = d0(1+dt*TR.dxidt) + dt.*(-TR.D);   % LHS on all tracers
FLHS = mfactor(LHS);    % run time = 0.44 sec for 3x3 any amount of tracers
A = TR.A;
H = TR.H;

for it = 1:steps_per_period

    n = n+1;
    [~,bgc] = calculate_forcing(sim, bgc, n);

    % FIXME     k = 4; if (it <= 2), k = it+1; end
    k = 4; if (it <= 2), k = it+1; end

    % calculate tendency
    tic  % comment this tic out to get cumulative time

    if (sim.runInParallel)
        [bgc, time_series] = MARBL_loop_parallel (n, sim, bgc,time_series);
    else
        [bgc, time_series] = MARBL_loop          (n, sim, bgc,time_series);
    end

    if (sim.captureAllSelectedTracers == 1)

        tmp = packMarbl( bgc.tracer, sim.domain.iwet_JJ );
        tracerCapture  (:,n) = tmp(:,sim.selection);

        tmp = packMarbl( bgc.tendency, sim.domain.iwet_JJ );
        tendencyCapture(:,n) = tmp(:,sim.selection);

    end
    % add river flux, as a tendency, to MARBL tendency

    bgc.tendency(:,1,:) = bgc.tendency(:,1,:) +river_tendency; % S = MARBL +river
    % bgc.tendency = 0 *bgc.tendency;                          % S = 0
    % bgc.tendency(:,1,:) = river_tendency;                    % S = 0 +river

    S  = packMarbl( bgc.tendency, sim.domain.iwet_JJ );

    %     toc
    %%
    % Capture time series of global volume integral of tracers and tendency
    % as part of the unpack
    switch (it)
        case 1
            w1 = C_1;
            RHS = C_1 + dt.*(A*w1 + H*C_1+ S);
            C_2 = mfactor(FLHS,RHS);
            bgc.tracer = unpackMarbl( C_2, sim.domain.iwet_JJ, size(bgc.tracer));
            %             bgc.accumulate = bgc.accumulate + ( S ./C_1 ); % FIXME: time stepped tracer or input tracer
            if (sim.logTracers)
                time_series.moles(:, n) = 1E-3 *sum( C_2 .* sim.domain.dVt_FP);
            end

        case 2
            w1 = (1/2).*(3*C_2 - C_1 );
            RHS = C_2 + dt.*(A*w1 + H*C_2+ S);
            C_3 = mfactor(FLHS,RHS);
            bgc.tracer = unpackMarbl( C_3, sim.domain.iwet_JJ, size(bgc.tracer));
            %             bgc.accumulate = bgc.accumulate + ( S ./C_2 ); % FIXME: time stepped tracer or input tracer
            if (sim.logTracers)
                time_series.moles(:, n) = 1E-3 *sum( C_3 .* sim.domain.dVt_FP);
            end

        otherwise
            w1 = ( (1/12).*( 23*C_3 - 16*C_2 + 5*C_1));
            RHS = C_3 + dt.*(A*w1 + H*C_3+ S);
            C_4 = mfactor(FLHS,RHS);
            bgc.tracer = unpackMarbl( C_4, sim.domain.iwet_JJ,size(bgc.tracer));
            %             bgc.accumulate = bgc.accumulate + ( S ./C_3 ); % FIXME: time stepped tracer or input tracer
            if (sim.logTracers)
                time_series.moles(:, n) = 1E-3 *sum( C_4 .* sim.domain.dVt_FP);
            end
    end

    if (sim.logTracers)
        time_series.moles (:, n) = global_moles(bgc.tracer          , sim);
        time_series.Dmoles(:, n) = 1E-3 *sum( S .* sim.domain.dVt_FP);
    end

    % Update temps storing previous time steps
    if (it > 2)
        C_1 = C_2; C_2 = C_3; C_3 = C_4;
    end

    %     disp(['Volume Integral: ',mfilename,'/test() end of Period #',num2str(n),' = ',num2str(time_series.moles(:, n)',12)])

    %     toc

    %%
    % Make occasional plots (slow) to monitor progress of sim and debug

    steps_between_plots = 1e+10;
    %    steps_between_plots = 8*7;
    %    steps_between_plots = 8*1;
    if mod(n, steps_between_plots) == 0
        if (sim.logTracers)
            small_plots(sim, time_series, n, sim.time_series_loc, sim.time_series_lvl);
            %         autoArrangeFigures(0,0);
        end
    end

    % fld=2; lvl=1;
    % Convert MARBL coordinates to FP
    % packed_data = packMarbl(bgc.forcing,sim.domain.iwet_JJ);%
    % Convert FP to xyz, but just the tracer/forcing of interest
    % data = nan(size(sim.domain.M3d)); data(sim.domain.iwet_FP) = packed_data(:,fld);
    % figure (100+fld); surf (data(:,:, 1), 'EdgeColor', 'none', 'FaceColor', 'interp'); view(2); colorbar;
    % % plot_global_interior_forcing(sim.domain.iwet_FP, sim.domain.M3d, sim.domain.zt, packed_data, 1)

    % [negativesFound] = negative_tracer_catcher(sim,bgc);
    % c0 = bgc2nsoli(sim, bgc.tracer);
    % fprintf('Fraction of tracers(:) <0 = %s\n', num2str(sum(c0<0)/numel(c0)));

end % steps_per_period
% mpiprofile viewer

if k <4 % error
    keyboard
end

% if (sim.logTracers)
% %     small_plots(sim, time_series, n, sim.time_series_loc, sim.time_series_lvl);
% %     disp(['Volume Integral: ',mfilename,' end of Month #',num2str(month),' = ',num2str(time_series.moles(:, n)',12)])
% end

if (sim.captureAllSelectedTracers == 1 && month == 12)

    myFileName = sprintf('%s/all_tracers_all_times.mat', sim.outputRestartDir);
    fprintf('%s.m: Saving "%s"...\n', mfilename, myFileName);
    tracer = tracerCapture;
% tmp_time_series = squeeze(time_series.tracer(sim.time_series_lvl,7,:));
% tmp = unpackMarbl(tracerCapture,sim.domain.iwet_JJ,[7881,60,730]);
% tmp_log = squeeze(tmp(sim.time_series_loc, sim.time_series_lvl,:));
% tst = [tmp_time_series tmp_log];
    save( myFileName, 'tracer' ,'-v7.3');

    myFileName = sprintf('%s/all_tendency_all_times.mat', sim.outputRestartDir);
    fprintf('%s.m: Saving "%s"...\n', mfilename, myFileName);
    tendency = tendencyCapture;
    save( myFileName, 'tendency' ,'-v7.3');

end

end % integration function
