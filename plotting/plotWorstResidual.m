function plotWorstResidual(sim, fname)

load(fname, 'r');
figure (700)
plot(r)

% rO2 = unpackMarbl(r,sim.domain.iwet_JJ,[7881,60,1]);
% rO2_tst = rO2(sim.time_series_loc,:);

myRng = 1:100;
[maxR,idxR] = sort(   (r),"descend",'MissingPlacement','last');
maxR(1)
[~, ~, ~, ~, ~, ~] = coordTransform_fp2xyz(idxR(myRng), sim, 701); title('Most postive')

[minR,idxR] = sort(   (r),"ascend",'MissingPlacement','last');
minR(1)
[~, ~, ~, ~, ~, ~] = coordTransform_fp2xyz(idxR(myRng), sim, 702);  title('Most Negative')

[maxAbsR,idxAbs] = sort(abs(r),"descend",'MissingPlacement','last');
[~, ~, ~, ~, ~, ~] = coordTransform_fp2xyz(idxAbs(myRng), sim, 703); title('Largest Abs')

% 
% [doubleSort,idxAbs] = sort(maxAbsR(:,1),"descend");
% doubleSort(1:30)
% idxAbs(1:30)
%
end