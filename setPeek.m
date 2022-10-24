function [sim] = setPeek(sim)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% I think it's impossible to debug anything with out looking at time step
% data. If you can debug without intermediate results just turn off
% logging.
% %
% However MARBL diagnostics, while occasionally invaluable are HUGE and
% slow the sim down by at least a factor 2x, require huge memry etc; so I
% very rarely capture them...

if (or (sim.logDiags, sim.logTracers))
    disp('Setting up to log tracers and possibly MARBL diags...')
    % need grid dimensions for many things. Starting with possible log files...
    [~, sim.domain.wet_loc, sim.domain.iwet_FP, ~, sim.domain.bottom_lvl, ~]...
        = createIwetEtc(sim.domain.M3d);

    % Possibly Record everything at a single location, water level. But which?
    % FIXME: need a clean way to convert (lat,lon) to (iLat,iLon), but in the
    % mean time; iterate thise by hand...
    %
    % iLat = 50; iLon = 81; iLvl = 10;    % Galapagos(0.3N, 108.7W)   iFp = 6496
    % iLat = 77; iLon = 60; iLvl = 10;    % NorPac  (45.7N, 176.8E)   iFp = 4629
    % iLat = 50; iLon = 61; iLvl = 1;     % Dateline  ( 0.3N,  179.3E)   iFp = 4702
    % iLat = 49; iLon = 11; iLvl = 10;    % Zulu =  ( 0.3N, 0.7E)     iFp = 1049
    % iLat = 58; iLon = 50; iLvl = 10;    % Palau = ( 5.6N, 139.7E)   iFp = 3704
    % iLat = 50; iLon = 28; iLvl = 10;    % IO      ( 0.3N,  50.5E)   iFp = 2080
    %     iLat = 57; iLon =  3; iLvl = 10;    % AF 447 =  ( 4.7N, -29.5E)     iFp = 1049
    iLat = 20; iLon =  95; iLvl = 4;    % "-48" =  ( -45.695N, -58.3E)     iFp = 31045 iCol 7462
    %     iLat = 2; iLon =  95; iLvl = 1;    % 7445
    % Check that! Make a map!
    % first get iFp on level 1, Simpy put: on level 1, iFp = iCol...

    iFp = coordTransform_xyz2fp(iLat, iLon, 1, sim);
    [~, ~, ~, ~, ~, ~] = coordTransform_fp2xyz(iFp, sim, 999); title('Time Series "Peek" Localtion')
    %     [~, ~, ~, ~, ~, ~] = coordTransform_fp2xyz(iFp, sim);
    sim.time_series_loc = iFp ;
    % % % ... now set the level
    sim.time_series_lvl = iLvl;
    disp(['Time series(loc,lvl) = (', num2str(sim.time_series_loc), ', ', num2str(sim.time_series_lvl),')']);
    [~,~,~, lat, lon, ~] = col2latlon(sim, sim.time_series_loc);
    disp(['Time series(lat,lon, depth) = (', num2str(lat,'%.1f'), ' N, ', num2str(lon,'%.1f'),' E, ',num2str(sim.grd.zt(sim.time_series_lvl),'%.1f'),' m))'])
    disp(' ')
    clear iLat iLon iLvl lat lon iFp
else
    % avoid messy code in parallel;
    % just set a default legal array idx
    sim.time_series_loc = 1 ;
    sim.time_series_lvl = 1;
end
end