function [sim, bgc] = saveRestartFiles(sim, bgc, tracer_0, years_gone_by)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%         toc(timer_total);
%         save the entire workspace. Surprisingly slow. Perhaps v7.3 compress of so much data is slow...
%         allFile = sprintf('%s/all_%d.mat', sim.outputRestartDir, round(1+years_gone_by));
%         fprintf('%s.m: Saving "%s"...\n', mfilename, allFile);
%         save(allFile,'-v7.3','-nocompression');
%
% Matlab load() has trouble with filenames that space and so on.
% KISS


tName = tracer_names(0);    % no CISO tracers

% save "x1" or final state file...

myRestartFile = sprintf('%s/restart_%d_%s_x1.mat', sim.outputRestartDir, round(sim.start_yr+years_gone_by),strjoin(tName(sim.selection)));
fprintf('%s.m: Saving "%s"...\n', mfilename,myRestartFile);
% copy original restart file, then replace original "tracer" with
% the current bgc.tracer. Surprisingly fast!
copyfile( sim.inputRestartFile, myRestartFile);

tracer = bgc.tracer;                            % --- these are x1 everywhere  ---
save( myRestartFile, 'tracer',  '-append' );    % overwrites tracer ONLY, keep forcing from init



% save "x0" or initial state file...

myRestartFile = sprintf('%s/restart_%d_%s_x0.mat', sim.outputRestartDir, round(sim.start_yr+years_gone_by),strjoin(tName(sim.selection)));
fprintf('%s.m: Saving "%s"...\n', mfilename,myRestartFile);
% copy original restart file, then replace original "tracer" with
% the current bgc.tracer. Surprisingly fast!
copyfile( sim.inputRestartFile, myRestartFile);

tracer = tracer_0;                            % --- these are x0 everywhere  ---
tracer(:,:,sim.selection) = bgc.tracer(:,:,sim.selection);   % only selected x1
save( myRestartFile, 'tracer',  '-append' );    % overwrites tracer ONLY, keep forcing from init

end
