function [iwet_JJ, wet_loc, iwet, num_wet_loc, bottom_lvl, iColLvlFromFp] = createIwetEtc(M3d)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% ================= Mess to convert from e.g. (91 x180 x24) grid to various
% packed coordinates

% iwet is list of idx into M3d(:) that are 1; IOW iwet picks wet data from 
% x(iLat,iLon,iLvl)
% 
% x(idx) = x(iwet) = x(iLat,iLon,iLvl) at (iLat,iLon,iLvl) that are wet
% 
% iwet is accidently sorted by level, since M3d is 
% 
% IOW iwet converts 3d grid e.g. (91,180,24) to wet only Francois (379913)

iwet = find(M3d(:));    % linear idx of (iLat,iLon,iLvl) for water parcels
sz   = size(M3d);

% convert x from (iLat, iLon) to x(iCol) simply with x(wet_loc) 

wet_loc     = find(M3d(:,:,1)); % sub2ind(sz,lat,lon) of water cols
num_wet_loc = size(wet_loc,1);  % just number water cols

% disp(['fraction of grid that is a water column: ', num2str(num_wet_loc/numel(M3d(:,:,1)))]);

% indx of deepest ocean level aka water level above bottom. IOW this is a 
% wet location, the land under ocean is +1 of this.
bottom_lvl= sum(M3d,3); % 2d(ilat, ilon), includes land

% kmt = bottom_lvl(wet_loc);    % deepest wet level for each water column
% [iLat,iLon] = ind2sub(sz,wet_loc); % (iLat,iLon) of 7881 water col


% Next is critical to invert fully packed to lat,lon,lvl
% 
% Convert M3d and iwet to other coordinate systems...
% M3d w/2d index of (sub2ind(sz,iLat,iLon),iLvl)
M3d_linear      = reshape(M3d, [sz(1)*sz(2), sz(3)]);   % (116*100, 60)

% iwet(iCol,iLvl) but full depth for MARBL
M3d_linear_wet 	= M3d_linear(wet_loc(:),:);             % (7881, 60)

% convert (iCol,1:60) to idx = (sub2ind(sz,iLat,iLon),iLvl)
iwet_linear_wet = find(M3d_linear_wet(:));  % converts from (iCol,iLvl) to iLat,iLon,iLvl)
iwet_JJ = iwet_linear_wet;                  % converts from (iCol,iLvl) to iLat,iLon,iLvl)

% Note '*' in the above comment; i.e (10441,24) -> (200160)

% get idx_FP into M3d from (iCol,iLvl) 
% iParcel = iColLvlFromFp(iCol, iLvl);
iColLvlFromFp          = zeros(size(M3d_linear_wet));
iColLvlFromFp(iwet_JJ) = 1:numel(iwet_JJ);

end % createIwetEtc