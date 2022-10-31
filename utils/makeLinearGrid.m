function [gridOfArray] = makeLinearGrid(wet_loc, gridSize, myArray)
% makeGridFromArray: Duplicate myArray into a grid.

%  grdSize is dimensions of a 2D array
%  myArray might be a scaler, 1D array, or matrix
%  size(myArray) = [1 1], or [1 N], or [M N]
% 
%  result is 2D or 3D global array

arraySize = size( myArray );

% FIXME: must be a clever "Matlab" way to do all this in one step...

% Make something with right basic idea, but wrong shape (lvl, tracer, lat, lon)...

gridOfArray = repmat( myArray, [1, 1, gridSize] );

% ...then make right shape (lat, lon, lvl, tracer)...

gridOfArray = permute( gridOfArray, [3 4 1 2] );      

% ...but that has subscripts for (iLat, iLon).
%
% Convert to linear index which is easy to divide among threads

gridOfArray = reshape( gridOfArray, [ prod( gridSize ), arraySize] );

% keep only wet locations

gridOfArray = gridOfArray( wet_loc, :, : );

end

