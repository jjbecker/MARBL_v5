function killAndCleanThreads()

% kill any leftover threads still running

delete(gcp('nocreate'))

myCluster = parcluster;     % writes empty garbage file to trash ?#$%^&*
delete(findJob(myCluster,'State','failed'));
delete(myCluster.Jobs);
clear myCluster

% clean out any leftover thread dumps from mex crashes%
% get the unix directory where the dumps are located...
% myCluster.JobStorageLocation
unix('rm -rf /Users/jj/Library/Application\ Support/MathWorks/MATLAB//local_cluster_jobs/R2020a/Job*');

% Whew, done with clean up!!!

end

