function [idx_fp] = coordTransform_bgc2fp(iCol, iLvl, sim)
%coordTransform_xyz2fp Convert "xyz" coordinate index vectors of "iLat,
% iLon, iLvl" to --SORTED-- single tracer "fp" index.
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
%       iLat, iLon, iLvl,       indices into M3d
%       sim    = struct with all dimension and times for a simulation
%
%   outputs:
%
%       idx_fp = --SORTED-- array of single tracer "fp" indices
%
%   example:
%
%       ifp = randperm(length(1:8500))';
%       [iLat, iLon, iLvl] = coordTransform_fp2xyz(ifp, sim, 345);
% %  Note: sort is needed for equality because ismember() sorts its input    
%       isequal(sort(ifp), coordTransform_xyz2fp(iLat, iLon, iLvl, sim))
%
%   See also: coordTransform_fp2nsoli, coordTransform_fp2bgc, coordTransform_fp2xyz,
%   coordTransform_nsoli2fp, coordTransform_nsoli2bgc, coordTransform_nsoli2xyz,
%   coordTransform_bgc2nsoli, coordTransform_bgc2fp, coordTransform_bgc2xyz,
%   coordTransform_xyz2nsoli, coordTransform_xyz2fp, coordTransform_xyz2bgc


[iLat, iLon] = ind2sub(size(sim.domain.M3d,[1 2]), sim.domain.wet_loc(iCol)');
[idx_fp] = coordTransform_xyz2fp(iLat, iLon, iLvl, sim);


end
