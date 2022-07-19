function [iCol, iLvl, iLat, iLon, lat, lon, depth] = coordTransform_fp2bgc(idx_fp, sim, debugFigNum,debugStr)
% coordTransform_fp2bgc Convert single tracer "water parcel" index to "bgc" indices of "water column" and "water
% level".
%
% Optional: output iLat, iLon, lat, lon, depth
% Optional: map location in figure(debugFigNum)
%
% First two steps is the very tricky part, at least to understand fully.
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
%       idx_fp      array of single tracer "fp" indices
%       sim         struct with all dimension and times for a simulation
%
%   outputs:
%       iCol                linear index into a M3d layer
%       iLvl                layer of M3d
%       optional:
%           iLat, iLon                  indices into a M3d layer
%           lat, lon, dept   geographic (deg,ged, meter)
%
%   example:
%
%   [iCol, iLvl, iLat, iLon, lat, lon, depth] = coordTransform_fp2bgc(1,sim,345);
%
%
%   See also: coordTransform_fp2nsoli, coordTransform_fp2bgc, coordTransform_fp2xyz,
%   coordTransform_nsoli2fp, coordTransform_nsoli2bgc, coordTransform_nsoli2xyz,
%   coordTransform_bgc2nsoli, coordTransform_bgc2fp, coordTransform_bgc2xyz,
%   coordTransform_xyz2nsoli, coordTransform_xyz2fp, coordTransform_xyz2bgc


% number of wet surface locations,
% which is also number of water cols
% which is the cnt we want...

num_location = numel(sim.domain.wet_loc);
num_parcels  = numel(sim.domain.iwet_FP);

% idx_fp--> linear idx of bgc with all water levels in each col

idx = sim.domain.iwet_JJ( mod(idx_fp-1,num_parcels) +1 )';

if (idx_fp <=0)
    keyboard        % FIXME: return 0 ???  Error???
else


    % ...rest is easy.

    iLvl = 1+ floor((idx   -1)/num_location);
    iCol = 1+ rem  ((idx   -1),num_location); % linear idx in bgc

    iTr  = 1+ floor((idx_fp-1)/num_parcels);  % tracer idx

    % Need to convert from water column number to linear index into M3d
    %   -then- convert to lat/lon

    [iLat, iLon] = ind2sub(size(sim.domain.M3d,[1 2]), sim.domain.wet_loc(iCol)');

    for i=numel(idx):-1:1 % trick to effectively pre allocate the arrays
        lat(i) = sim.grd.YT(iLat(i),iLon(i),iLvl(i));   % +/- deg
        lon(i) = sim.grd.XT(iLat(i),iLon(i),iLvl(i));   % +/- deg
        lon(i) = mod(180+lon(i),360)-180;               % +/- deg
    end
    depth  = sim.domain.zt(iLvl);                       % m

    if nargin<3
        return
    else

        % Debug with a plot

        if (debugFigNum <1)
            return
        elseif numel(iCol) <= 0
            disp('coordTransform_fp2xyz: No wet locations')
            return
        end

        %     close(debugFigNum); % error if figure is not open
        figure(debugFigNum);

        %         myColors = 1:numel(iLat);
        myColors = iTr;
        geoscatter(lat, lon, 100, myColors,'^','filled');

        %         a = myColors;                  c = num2cell(a); % strings to label
        %         dx = 0.; dy = 0.; % displacement so the text does not overlay the data points
        %         text(lat+dx, lon+dy, c);

        colorbar;
        title(debugStr);
        geobasemap topographic
    end
end

