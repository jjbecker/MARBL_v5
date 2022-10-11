function [restartDir] = myRestartDir(ck_years)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if isunix && ~ismac
    % Code to run on Linux platform
    restartDir = sprintf('/DFS-L/DATA/primeau/jjbecker/restart_0_%s',int2str(ck_years));
    % restartDir = '/DFS-L/SCRATCH/primeau/jjbecker/';
elseif ismac
    % Code to run on Mac platform
    restartDir = sprintf('Data/restart_0_%s',int2str(ck_years));
elseif ispc
    % Code to run on Windows platform
    disp('Platform not supported')
    keyboard
else
    disp('Platform not supported')
    keyboard
end

restartDir = sprintf('%s%s',restartDir,'_output');

if ~exist(restartDir, 'dir')
    mkdir(restartDir)
end

end