function [iLat, iLon, kmt, lat, lon, water_depth_km] = col2latlon(sim, colNum, debugFigNum, debugStr)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% if nargin >2
%     colNum = colNum - (tracer-1)*sim.domain.num_wet_loc*size(sim.domain.dVt,3);
% end

if and(numel(colNum) >500, exist('debugFigNum','var'))
    disp('That is a whole lot of points to plot with this slow code...')
    keyboard
end
[iLat, iLon, ~] = ind2sub(size(sim.domain.M3d), sim.domain.wet_loc(colNum));
% tmp = sim.domain.bottom_lvl(:);
% % % % skip drama and use a for loop. Ohhhh no!!!

for i=numel(colNum):-1:1
    kmt(i)            = sim.domain.bottom_lvl(iLat(i), iLon(i));   % kmt is level of bottom
    lat(i)            = sim.grd.YT(iLat(i), iLon(i),1);
    lon(i)            = mod(180+sim.grd.XT(iLat(i), iLon(i),1),360)-180;
    water_depth_km(i) = sim.grd.ZT3d(iLat(i), iLon(i),kmt(i))/1e+3;
end
% kmt = tmp(colNum);


for idx = 1:numel(iLat)
    disp(['loc #',num2str(colNum(idx)), ...
        ' (lat,lon,bottom_depth): (', ...
        num2str(lat(idx), '%1.1f'),'N, ', ...
        num2str(lon(idx), '%1.1f'),'E, ', ...
        num2str(water_depth_km(idx), '%1.3f'),'km)', ...
        ' (iLat,iLon,kmt): (', ...
        num2str(iLat(idx)),',', ...
        num2str(iLon(idx)),',', ...
        num2str(kmt(idx)), ')' ])
end

if nargin<3
    return
else

    % Debug with a plot

    if (debugFigNum <1)
        return
    elseif numel(colNum) <= 0
        disp('col2latlon: No wet locations')
        return
    end

%     close(debugFigNum); % error if figure is not open
    figure(debugFigNum);

    myColors = 1:numel(iLat);
    geoscatter(lat, lon, 100, myColors,'^','filled');
    colorbar
%     b = num2str(myColors'); 
    b = num2str(colNum); 
    c = cellstr(b); % strings to label
    dx = 0.; dy = 0.; % displacement so text does not overlay data points
    text(lat+dx, lon+dy, c);
    title(debugStr, 'Interpreter', 'none');
    geobasemap topographic
end


end
