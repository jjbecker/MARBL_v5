function [dataDir] = myDataDir()
% function [restartDir] = myDataDir(ck_years)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if isunix && ~ismac
    % Code to run on Linux platform
    dataDir = sprintf('/DFS-L/DATA/primeau/jjbecker/restart_0_%s',int2str(ck_years));
    % dataDir = '/DFS-L/SCRATCH/primeau/jjbecker/';
elseif ismac
    % Code to run on Mac platform
%     restartDir = sprintf('Data/restart_0_%s',int2str(ck_years));
    dataDir = strcat(pwd,'/Data/');
elseif ispc
    % Code to run on Windows platform
    disp('Platform not supported')
    keyboard
else
    disp('Platform not supported')
    keyboard
end

if ~exist(dataDir, 'dir')
    mkdir(dataDir)
end

end