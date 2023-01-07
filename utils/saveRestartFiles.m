function saveRestartFiles(sim, tracer, newRestartFileName)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fprintf('%s.m: Saving "%s"...\n', mfilename, newRestartFileName);
if sim.debug_disable_phi
    fprintf('\n\n\t%s.m: ********* phi() is short circuited. just touch file  *********\n\n',mfilename);
    system(sprintf('touch %s', newRestartFileName));
    return
end


% copy sim.inputRestartFile, and replace it's "tracer" with input tracer
%
% Surprisingly fast

if ~exist(sim.inputRestartFile, 'file')
    keyboard
end
copyfile( sim.inputRestartFile, newRestartFileName);

save( newRestartFileName, 'tracer',  '-append' );    % overwrites tracer ONLY, keep forcing from init
if ~exist(newRestartFileName, 'file')
    keyboard
end

end
