function [sim, bgc] = saveRestartFiles(sim, bgc, tracer, newRestartFileName)
% function [sim, bgc] = saveRestartFiles(sim, bgc, tracer_0, years_gone_by)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if sim.debug_disable_phi
    fprintf('\n\n\t%s.m: ********* phi() is short circuited skip inputRestartFile save  *********\n\n',mfilename)
    return
end

% copy sim.inputRestartFile, and replace it's "tracer" with input tracer
%
% Surprisingly fast


fprintf('%s.m: Saving "%s"...\n', mfilename, newRestartFileName);

if ~exist(sim.inputRestartFile, 'file')
    keyboard
end
copyfile( sim.inputRestartFile, newRestartFileName);

save( newRestartFileName, 'tracer',  '-append' );    % overwrites tracer ONLY, keep forcing from init
if ~exist(newRestartFileName, 'file')
    keyboard
end

% save "x1" or final state file...
% tName = tracer_names(0);    % no CISO tracers
% myRestartFile = sprintf('%s/restart_%d_%s_x1.mat', sim.outputRestartDir, round(sim.start_yr+years_gone_by),strjoin(tName(sim.selection)));

%
%
% % save "x0" or initial state file...
%
% myRestartFile = sprintf('%s/restart_%d_%s_x0.mat', sim.outputRestartDir, round(sim.start_yr+years_gone_by),strjoin(tName(sim.selection)));
% fprintf('%s.m: Saving "%s"...\n', mfilename, newRestartFileName);
% % copy original restart file, then replace original "tracer" with
% % current bgc.tracer. Surprisingly fast!
% copyfile( sim.inputRestartFile, newRestartFileName);
%
% tracer = tracer_0;                            % --- these are x0 everywhere  ---
% tracer(:,:,sim.selection) = bgc.tracer(:,:,sim.selection);   % only selected x1
% save( myRestartFile, 'tracer',  '-append' );    % overwrites tracer ONLY, keep forcing from init

end
