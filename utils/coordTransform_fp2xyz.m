function [iLat, iLon, iLvl, latitude, longitude, depth] = coordTransform_fp2xyz(idx_fp, sim, debugFigNum)
%coordTransform_fp2xyz Convert single tracer "water parcel" index to "iLat, iLon, iLvl".
%
% Optional: output geographical latitude, longitude, depth
% Optional: map location in figure(debugFigNum)
%
% If "debugFigNum" is input, geoscatter plot of locations is opened.
%
% If you want to convert "fp2xyz(idx,iTr)" just call fp2xyz(idx) because fp
% is the same for all tracers in a given water parcel.
%
% ==================
% Basically removing bottom levels made it at least tricky, not impossible,
% to convert "fp" idx back to (iLat, iLon, iLvl). The "find()" operation is
% culprit.
%
%   wet_loc         = find(M3d(:,:,1));           % "water colum locations"
%
% Convert both land and ocean to 2 subscripts: 3d->2d index = (sub2ind(lat,lon),lvl)
%   M3d_linear      = reshape(M3d, [sz(1)*sz(2), sz(3)]);
%
% Keep only ocean locations, with all levels: aka a water colum with all levels, even the bottom
%   M3d_linear_wet  = M3d_linear(wet_loc(:),:);
%   iwet_linear_wet = find(M3d_linear_wet(:));    % "" no bottom levels
%   iwet            = iwet_linear_wet;            % single index
% ==================
%
%   inputs:
%       idx_fp = array of single tracer "fp" indices
%       sim    = struct with all dimension and times for a simulation
%
%   outputs:
%
%       iLat, iLon, iLvl,       indices into M3d
%       latitude, longitude     deg (-90<lat<=90, -180<=lat<180) FIXME
%       depth                   meters
%
%   example:
%
%       [iLat, iLon] = coordTransform_fp2xyz([100 300 1], sim, 1)
%
%       [~,~,~,lat,long     ] = coordTransform_fp2xyz([100 300 1],sim,1)
%       [~,~,~,lat,lon,depth] = coordTransform_fp2xyz([100 300 1],sim,1)
%
%       coordTransform_fp2xyz([100 300 1],sim,1);  % returns iLat
%
%   and sorted uniquely...
%
%       [iLat, iLon] = coordTransform_fp2xyz(unique([100 101 1]), sim, 1)
%
%   See also: coordTransform_fp2nsoli, coordTransform_fp2bgc, coordTransform_fp2xyz,
%   coordTransform_nsoli2fp, coordTransform_nsoli2bgc, coordTransform_nsoli2xyz,
%   coordTransform_bgc2nsoli, coordTransform_bgc2fp, coordTransform_bgc2xyz,
%   coordTransform_xyz2nsoli, coordTransform_xyz2fp, coordTransform_xyz2bgc
%

% First step is the worst.
%
% iwet_FP = find(M3d(:));   % linear index of wet water parcel
%
%   --> "xyz" version of "iwet"

idx_M3d = sim.domain.iwet_FP(idx_fp); % linear index of M3d

% ...rest is easy. From here on, -IGNORE- wet dry stuff

sz = size(sim.domain.M3d,[1 2]);   % dimension of a single layer in M3d
num_location = prod(sz);    % or numel in a M3d lvl

iLvl = 1+ floor((idx_M3d-1)/num_location);  % layer of M3d
idx  = 1+ rem  ((idx_M3d-1),num_location);  % linear idx in a M3d layer
[iLat, iLon] = ind2sub(sz, idx );           % linear idx --> [row,col]

% latitude  = sim.grd.yt(iLat);               % +/- deg
% longitude = sim.grd.xt(iLon);               % 0 < lon <= 360
% depth     = sim.grd.zt(iLvl);               % meter

% latitude  = sim.grd.YT(iLat,iLon,iLvl);           % +/- deg
% longitude = sim.grd.XT(iLat,iLon,iLvl);
% longitude = mod(180+longitude,360)-180; % +/- deg
% depth = sim.grd.zt(iLvl);               % m

for i=numel(idx):-1:1
    latitude(i)  = sim.grd.YT(iLat(i),iLon(i),iLvl(i));     % +/- deg
    longitude(i) = sim.grd.XT(iLat(i),iLon(i),iLvl(i));     % 0 <= lon <360
    longitude(i) = mod(180+longitude(i),360)-180;           % -180 =< lon < 180
end
depth = sim.grd.zt(iLvl);               % m

if nargin<3
    return
else
    if (debugFigNum <1)
        return
    elseif numel(latitude) <= 0
        disp('coordTransform_fp2xyz: No wet locations')
        return
    end
    figure(debugFigNum);

    geoscatter(latitude, longitude, 100,'black','^','filled');
    geobasemap topographic

end

end
