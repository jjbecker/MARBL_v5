function [packed] = packMarbl(rectangular,iwet_JJ)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% num_rows    = size(rectangular,1);
% num_lvls    = size(rectangular,2);
% num_tracers = size(rectangular,3);

sz  = size(rectangular);
tmp = reshape(rectangular, [sz(1)*sz(2), sz(3)]);
packed = tmp(iwet_JJ,:);


% % % % Check it
% % % 
% % % % preallocate for speed
% % % 
% % % tst_packed  = zeros(size(iwet_JJ,1),num_tracers);
% % % 
% % % for iTracer = 1: num_tracers
% % %     my_slice   = rectangular(:,:,iTracer);
% % %     tst_packed (:,iTracer) = my_slice(iwet_JJ);
% % % end
% % % 
% % % if ~isequaln(packed,tst_packed)
% % %     keyboard
% % % end
% % % 
end
