function [sim, bgc_struct] = calculate_depth_map_and_volumes(sim)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Define gird dimensions in meters

sim.domain.zt   = sim.grd.zt;
sim.domain.zw   = sim.grd.zw+sim.grd.dzt; % in MARBL zw is bottom, not top...
sim.domain.dzt  = sim.grd.dzt;
sim.domain.dzw  = sim.grd.dzw;

% ---> MARBL uses units of - CM - for depths, Matlab code uses meters!
sim.domain.MARBL_depth_per_m = 100;     % 1 meter here is 100 cm MARBL uses

% Create a struct that has meta data like size, names, and units.

bgc_struct = init_bgc_struct(sim);


[sim.domain.iwet_JJ, sim.domain.wet_loc, sim.domain.iwet_FP, sim.domain.num_wet_loc, ...
    sim.domain.bottom_lvl, sim.domain.iColLvlFromFp] = createIwetEtc(sim.domain.M3d);

% 2d matrix of "depth", e.g. 0:10900 meter <-- needed by MARBL
% FIXME: is this actually used anywhere???
% Only thing tricky is that Matlab array cant have index of 0...

zw = cumsum(sim.grd.dzt')';     % bottom of grid, not depth at given local

mapping = [1:length(zw); zw]';
mapping = [0 0; mapping];           % add row for "dry" with indx == 0
mapping(:, 1) = mapping(:, 1) +1;   % +1 for row with "dry" value
LUT(mapping(:, 1)) = mapping(:, 2);
sim.domain.bottom_depth = LUT(sim.domain.bottom_lvl+1);   % +1 for row with "dry" value

sim.domain.dVt    = sim.grd.dV;
sim.domain.dVt    = sim.domain.M3d .* sim.domain.dVt;
sim.domain.V      =         sum(sim.domain.dVt, 'all') ;
% with depth: layers get thicker, but area decreases. Volume per level
% almost constant but sort of random
sim.domain.V_lvl  = squeeze(sum(sim.domain.dVt, [1 2])); % vol at given iLvl
sim.domain.V_lat  = squeeze(sum(sim.domain.dVt, [3 2])); % vol at given iLat
sim.domain.V_lon  = squeeze(sum(sim.domain.dVt, [3 1])); % vol at given iLon

sim.domain.dVt_FP = sim.domain.dVt(sim.domain.iwet_FP);

end