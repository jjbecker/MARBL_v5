function [idx_fp] = coordTransform_xyz2fp(iLat, iLon, iLvl, sim)
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


% First step is the worst.
%
% iwet_FP = find(M3d(:));   % linear index of wet water parcel
%
%   --> "xyz" version of "iwet"

sz = size(sim.domain.M3d,[1 2]);           % dimension of a single layer in M3d
iCol = sub2ind(sz,iLat,iLon);       % linear index of (ilat,ilon) in level
idx_xyz = iCol + prod(sz)*(iLvl-1); % linear index of xyz in all of M3d

% This is the tricky bit use sim.domain.iwet_FP,idx_xyz) to convert to 
% index in wet part of M3d. 
% 
% There are several clever ways to do this. "ismemeber" is fastest

xyz_is_wet = ismember(sim.domain.iwet_FP,idx_xyz); % This sorts idx_xyz
idx_fp = find(xyz_is_wet);                              
if (isempty(idx_fp))
    keyboard
end

% [~,idx_fp] = intersect(sim.domain.iwet_FP, idx_xyz);  % 3x slower. obscure. sorts...
% idx_fp = sort(dsearchn(sim.domain.iwet_FP, idx_xyz)); % clean. REALLY slow. SCRAMBLES idx_fp

end
