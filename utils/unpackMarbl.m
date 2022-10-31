function [unpacked] = unpackMarbl(packed,iwet_JJ,sz_unpacked)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% num wet node  = size(packed,1); % ~num_rows*num_lvls, -BUT- not seafloor
% num_tracers = size(packed,2);
% 
% num_rows    = sz_unpacked(1);
% num_lvls    = sz_unpacked(2);
% num_tracers = sz_unpacked(3);

unpacked = nan( [ prod(sz_unpacked(1:2)) , sz_unpacked(3) ] );
unpacked(iwet_JJ,:) = packed;
unpacked = reshape(unpacked,sz_unpacked);


% % % % Check it
% % % 
% % % % preallocate for speed
% % % 
% % % num_tracers = sz_unpacked(3);
% % % sz_unpacked(3) = 1;
% % % tst_unpacked = zeros(sz_unpacked);  % one tracer
% % % 
% % % for iTracer = 1: num_tracers
% % %     my_slice = zeros(sz_unpacked);      % restore bottom of water col
% % %     my_slice(iwet_JJ) = packed(:,iTracer);    
% % %     tst_unpacked(:,:,iTracer) = my_slice;   
% % % end
% % % 
% % % if ~isequaln(unpacked,tst_unpacked)
% % %     keyboard
% % % end

end

