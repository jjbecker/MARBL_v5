function J = calc_J_O2(sim, bgc, time_series, forcing, MTM)

fprintf('%s.m: Start at %s\n', mfilename, datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
%calc_J Summary of this function goes here
%   Detailed explanation goes here
%
% To avoid mysterous memory leaks, need to use "time_step_ann()" which
% requires a accurate struct for time_series, and obviously needs correct
% forcing and treansport.
%
%
%   x0 Input is O2 in "fp" format eg. sz=[8500,1] MARBL wants sz=[545,20,32]
%
% calculate partial of MARBL with respect to O2 in its own col
%
% f = MARBL() = tendency = time rate of change of all tracers at all loc
%
% J = simple numberical approximation of Jacobian(f)
%
% 3 index array:
%
%   first is location (aka water column),
%   second is index of water level
%   third is index of water lever of same tracer tracer being changed
%
%   J(:,i,j) = d(MARBL(O2(i))) / d(tracer(O2(j)) @ all water col, lvl i,
%   wrt O2 on lvl j

bgc_0 = bgc;

tName = tracer_names(0);    % no CISO tracers
tStr = string(tName(sim.selection));
fprintf('%s.m: Simplified Jacobian of "%s"; takes ~300 (s)\n',mfilename, tStr);

tStart = tic;

sz = size(bgc.tracer);
J = zeros([sz(1) sz(2) sz(2)]);

% calculate d(Tracer)/dt with tracer values x0
% need f0; size(f0)= [7881,60,32]
% f0 = calc_f(x0, sim, bgc, time_series, forcing, MTM,1);
month = 1;

[sim, bgc, time_series, n] = time_step_ann (sim, bgc, time_series, -1, forcing(month), MTM(month), month);

f0 = bgc.tendency;
f0 = squeeze(f0(:,:,sim.selection));   % size(f0) = [7881,60]

% FIXME: always tricky to guess "delta"
% h = sqrt(eps) *x0 but need to convert x0 from FP to MARBL coordinate
% size(h) =  [7881, 60, 1]
% h = unpackMarbl(x0 *sqrt(eps), sim.domain.iwet_JJ,[sz(1),sz(2),1]);

h = squeeze(bgc.tracer(:,:,sim.selection));   % size(f0) = [7881,60]
dh = eps + sqrt(eps) *h;

for lvl = 1:sz(2)       % loop over all lvl

    %     dx = h*0;
    %     dx(:,lvl) = h(:,lvl);
    %     dx = dx(sim.domain.iwet_JJ);
    %
    %     x1 = x0 +dx;
    %
    %     % calculate d(Tracer)/dt with tracer values x1
    %
    %     f1 = calc_f(x1, sim, bgc, time_series, forcing, MTM,1);

    bgc.tracer = bgc_0.tracer;
    bgc.tracer (:, lvl, sim.selection) = bgc.tracer (:, lvl, sim.selection) +dh(:,lvl);

    [sim, bgc, time_series, n] = time_step_ann (sim, bgc, time_series, -1, forcing(month), MTM(month), month);

    f1 = bgc.tendency;
    f1 = squeeze(f1(:,:,sim.selection));   % size(f1) = [7881,60]
    df = f1 -f0;
    df_dh = df ./ dh;  % size(df_dh) = [60 32] = d(tend(:,:))/dh(lvl,tr)

%     J(:,:,lvl) = df(:,:,sim.selection)./h(:,lvl);
    J(:,:,lvl) = df_dh;
end

elapsedTime = toc(tStart);
fprintf('%s.m: %1.3f (s) for partial of MARBL tendency(%s) on all levels, w.r.t. (%s) of all levels of same column, for all columns\n', mfilename, elapsedTime, tStr, tStr);

logJ=log10(abs(nonzeros(J(:))));
figure(400+sim.selection); histogram(logJ); title(sprintf("hist( log10( abs( J( %s ))))",tStr), 'Interpreter', 'none');

logJT=log10(sim.T *abs(nonzeros(J(:))));
figure(500+sim.selection); histogram(logJT); title(sprintf("hist( log10( abs( sim.T *J( %s ))))",tStr), 'Interpreter', 'none');

fprintf('End of %s.m: %s\n', mfilename, datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
end
