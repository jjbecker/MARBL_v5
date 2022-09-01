function tendency = calc_f( x0, sim, bgc, time_series, forcing, MTM, month)

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

c = packMarbl(bgc.tracer,sim.domain.iwet_JJ);
c = replaceSelectedTracers(sim, c(:), x0, sim.selection);
numWaterParcels = numel(sim.domain.iwet_JJ);
numTracers = sim.bgc_struct_base.size.tracer(2);
c = reshape(c, [numWaterParcels, numTracers]);
bgc.tracer = unpackMarbl(c, sim.domain.iwet_JJ,size(bgc.tracer));  

[~, tmp, ~, ~] = time_step_ann (sim, bgc, time_series, -1, forcing(month), MTM(month), month);

tendency = tmp.tendency;

clear c tmp

return

end
