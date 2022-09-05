function J = calc_J_full(sim, bgc, time_series, forcing, MTM)

fprintf('%s.m: Start at %s\n', mfilename, datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
%calc_J Caclulate complete J, all tracers all levels, but at one location
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


tStart = tic;

% Code is confusing enough. Just make a copy of the imputs and then mess
% with them; keep copy of original just in case. But Matlab is call by
% vlaue so if we do not return sim or bgc then we don't need to restore
% them.

% sim_0 = sim;
bgc_0 = bgc;
sz = size(bgc.tracer);              % sz = [7881, 60, 32]

% want J at one water column, for all water levels and tracers wrt to all
% same column tracers at all water levels. Furthermore in the end we want
% a 2d array with all waterlevels of a single tracer, then same for next
% tracers, final size = [60*32, 60*32]

% FIXME: want to deal with slow speed of loading a sparse or use full?
J = zeros([sz(2) sz(3) sz(2) sz(3)]);   % size(J) = [60, 32, 60, 32]

% pick a location that is full depth, this finds all of them...
% full_depth_col = find(sim.domain.bottom_lvl(sim.domain.wet_loc) == 60);
% 
% This maps all of them in figure 666
% [iLat, iLon, kmt, lat, lon, water_depth_km] = col2latlon(sim, full_depth_col, 666, 'full_depth_col');
%
% Pick a water column east of Guam, near Mariana's Trench iCol = 3943
full_depth_col = 3943;
% [iLat, iLon, kmt, lat, lon, water_depth_km] = col2latlon(sim, full_depth_col, 666, 'full_depth_col');
col2latlon(sim, full_depth_col, 666, 'full_depth_col');

tName = tracer_names(0);    % no CISO tracers

month = 1;

[sim, bgc, time_series, n] = time_step_ann (sim, bgc, time_series, -1, forcing(month), MTM(month), month);

f0 = bgc.tendency;
f0 = squeeze(f0(full_depth_col,:,:));   % size(f0) = [60,32]

for iTr = 1:sz(3)

    sim.selection = iTr;

    tStr = string(tName(iTr));
    fprintf('%s.m: Full Jacobian of MARBL tendency for "%s"; takes ~600 (s)\n',mfilename, tStr);

    for lvl = 1:sz(2)       % loop over all lvl

        % h is a scaler = delta of this tracer on this lvl in this col
        h = squeeze(bgc.tracer(full_depth_col, lvl, iTr));
        dh = eps + sqrt(eps) *h;
        % dh = max(dh, sqrt(eps));

        bgc.tracer = bgc_0.tracer;
        bgc.tracer (full_depth_col, lvl, iTr) = bgc.tracer (full_depth_col, lvl, iTr) +dh;

        % calculate global tend = with x0+h = tracer values x1

        [sim, bgc, time_series, n] = time_step_ann (sim, bgc, time_series, -1, forcing(month), MTM(month), month);
        f1 = bgc.tendency;

        % keep only tendency for the choosen water col
        f1 = squeeze(f1(full_depth_col,:,:));
        df = f1 -f0;
        df_dh = df ./ dh;  % size(df_dh) = [60 32] = d(tend(:,:))/dh(lvl,tr)

        % first 2 cols specify location of --> h <--, col 3 and 4 are
        % change in tend for all tracers and lvl
        J(lvl, iTr, :, :) = df_dh;

        elapsedTime = toc(tStart);
        fprintf('%s.m: %1.3f (s) for partial of all MARBL tendency all levels, w.r.t. (%s) thru lvl %d of same column, for all columns. h %g, dh %g, norm(df_dh(:)) %g\n', mfilename, elapsedTime, tStr, lvl, h, dh, norm(df_dh(:)));
    end
end

elapsedTime = toc(tStart);
fprintf('%s.m: %1.3f (s) for partial of MARBL tendency(%s) on all levels, w.r.t. (%s) of all levels of same column, for all columns. \n', mfilename, elapsedTime, tStr, tStr);

% J_perm = permute(J, [2 1  3 4]);
J_1920x1920 = reshape(J, [], 32*60);
save('J_1920x1920', 'J_1920x1920');
spy(J_1920x1920);

fprintf('End of %s.m: %s\n', mfilename, datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
end
