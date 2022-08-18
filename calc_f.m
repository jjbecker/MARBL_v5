function tendency = calc_f( sim, bgc, time_series, forcing)

% calc_f: Calls MARBL with given tracers "C" in "fp" format [8500,32] etc.
%
% Input:
%
%   C : "fp" format tracers eg. sz=[8500,32], with MARBL unit eg (mmol/kg)
%
%   sim, and bgc: used for dimensions, M3d, and such
%
%   time_series: not changed but needs to have proper dimensions
%
% Does -NOT- change bgc, or time series when we call MARBL

% t1 = tic;

% make a "bgc" with "C" tracers passed in...
%
% MARBL wants sz=[545,20,32]

% Do -NOT- change bgc, or time series when we call MARBL

% tmp = bgc;
% x = replaceSelectedTracers(sim, ...
%     packMarbl(bgc.tracer,sim.domain.iwet_JJ), ...
%     C, sim.selection);
% tmp.tracer = nsoli2bgc(sim, bgc, x);   % marbl format

% make sure call to MARBL initializes it; aka n=1

tmp = bgc;
tmp.forcing      = forcing.interior;
tmp.surf_forcing = forcing.surf_forcing;
n = 1;
[~,tmp] = calculate_forcing(sim, tmp, n);

if (sim.runInParallel)
    [tmp, ~] = MARBL_loop_parallel ( n, sim, tmp, time_series);
else
    [tmp, ~] = MARBL_loop          ( n, sim, tmp, time_series);
end

tendency = tmp.tendency;

% tendency = packMarbl( tmp.tendency, sim.domain.iwet_JJ );

clear tmp
clear x x0

% elapsedTime = toc(t1);
% disp(' ');disp(['f: ', num2str(elapsedTime*1, '%1.3f'),' (s) for 1 tendency of all columns'])

end
