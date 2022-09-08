function J = calc_J_Single_Tracer(sim, bgc, time_series, forcing, MTM)

fprintf('%s.m: Start at %s\n', mfilename, datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
%calc_J_O2 Summary of this function goes here
%   Detailed explanation goes here
%
% To avoid mysterous memory leaks, need to use "time_step_ann()" which
% requires a accurate struct for time_series, and obviously needs correct
% forcing and treansport.
%
% calculate partial of single selected MARBL tracer single output with 
% respect to its tracer in its own col
%
% f = MARBL() = tendency = time rate of change of tracers at all loc
%
% J = simple numerical approximation of Jacobian(f)
%
% 3 index array:
%
%   first is location (aka water column),
%   second index is water level where input tracer was changed
%   third is water lever index of output tendecy
%
%   J(:,i,j) = d(MARBL(:,i,idxTr)) / d(tracer(:,j,idxTr) @ all water col

bgc_0 = bgc;

tName = tracer_names(0);    % no CISO tracers
tStr = string(tName(sim.selection));
fprintf('%s.m: Simplified Jacobian of "%s"; takes ~300 (s)\n',mfilename, tStr);

tStart = tic;

sz = size(bgc.tracer);
J = zeros([sz(1) sz(2) sz(2)]);                 % [7881,60, 60]

x0 = squeeze(bgc.tracer(:,:,sim.selection));    % size(x0) = [7881,60]
dx = eps + sqrt(eps) *x0;                       % x0 = NaN below bottom

% set x0, dx below bottom to avoid divide by Nan (not necessary)
% x0(~isfinite(x0(:))) = 0;                       % set x0=0 below bottom
% dx(~isfinite(dx(:))) = 0;                       % set dx=0 below bottom

% FIXME: always tricky to guess "delta x"
% dx = eps +sqrt(eps) *x0                       % size(dx) = [7881, 60]

% calculate d(Tracer)/dt with x0 tracer values
% need f0; size(f0)= [7881,60,32]

month = 1;
[sim, bgc, time_series, n] = time_step_ann (sim, bgc, time_series, -1, forcing(month), MTM(month), month);

f0 = bgc.tendency;
f0 = squeeze(f0(:,:,sim.selection));            % size(f0) = [7881,60]

for h_lvl = 1:sz(2)       % loop over all levels

    % add finite difference to tracer in all water col, but only this level
    %
    % change the tracer of every water column, but only on one level..
    % can be confusing, next line is correct but...
    % bgc.tracer(:, h_lvl, sim.selection) = bgc.tracer(:, h_lvl, sim.selection)  +dx(:, h_lvl);
    
    % remember that x0, dx and x1 are NaN below bottom
    
    x1 = x0;
    x1(:,h_lvl) = x0(:,h_lvl)  +dx(:, h_lvl);   % x1=x0+h; only on j_lvl

    % calculate tendency every where with x1 = x0+dx

    bgc.tracer = bgc_0.tracer;  % use saved initial value for tracer
    bgc.tracer(:,h_lvl,sim.selection) = x1(:,h_lvl);

    [sim, bgc, time_series, n] = time_step_ann (sim, bgc, time_series, -1, forcing(month), MTM(month), month);

    f1 = bgc.tendency;                      % f1(x1) = tendency all tracers
    f1 = squeeze(f1(:,:,sim.selection));    % tendency of selected tracer 

    df = f1 -f0;                            % delta(tend)

    % remember that f0, f1 and df are zero below the bottom

    % divide df, the change of tendency of every water column on every 
    % level by dh the change of of its own water columnm tracer but ONLY 
    % on the given level
    
    % only do calculation on columns where kmt >= h_lvl (aka wet at this
    % depth)

    valid_col = find(bgc.kmt >= h_lvl);

    dh = x1 - x0;
    dh = dh(:,h_lvl);                       % [7881,1]

    df_dx = 0 * f0;
    df_dx(valid_col,:) = df(valid_col,:) ./dh(valid_col);

    % check it!
    if sum(~isfinite(df_dx(:))) > 0
        keyboard
    end

    % To be clear!
    % 
    % Calculated how much tendency changed in every water column when O2
    % changed a small amount in every water column, but only on a given 
    % SINGLE level "h_lvl"
    %
    % J(10,20,30) is change in tendency of water column #10, at level
    % #20 for a small change in tracer on column #10 level #30

    J(:,:,h_lvl) = df_dx;                           % [7881, 60, 60] = [iCOl,jLvl, iLvl]

end

J_perm   = permute(J,[1 3 2]);                      % [7881, 60, 60] = [iCOl,hLvl, iLvl]
J_packed = packMarbl(J_perm,sim.domain.iwet_JJ);    % [379913, 60)
% [iCol, iLvl, iLat, iLon, lat, lon, depth] = coordTransform_fp2bgc(379893, sim, 666,'');

elapsedTime = toc(tStart);
fprintf('%s.m: %1.3f (s) for partial of MARBL tendency(%s) on all levels, w.r.t. (%s) of all levels of same column, for all columns\n', mfilename, elapsedTime, tStr, tStr);

logJ=log10(abs(nonzeros(J(:))));
figure(400+sim.selection); histogram(logJ); title(sprintf("hist( log10( abs( J( %s ))))",tStr), 'Interpreter', 'none');

logJT=log10(sim.T *abs(nonzeros(J(:))));
figure(500+sim.selection); histogram(logJT); title(sprintf("hist( log10( abs( sim.T *J( %s ))))",tStr), 'Interpreter', 'none');

J = J_perm;


fprintf('%s.m: nnz(J) = %f\n', mfilename, nnz(J));
fprintf('%s.m: nnz(J)/numel(x0(:) = %f\n', mfilename, nnz(J)/numel(x0(:)));

[iCol, iLvl] = coordTransform_fp2bgc(1:379913, sim);

numRows = numel(iCol);
J_2d = sparse(numRows,numRows);

tic;
for row = 1:numRows
    myLvl = find(J_packed(row,:) ~=0);
    if numel(myLvl) > 0
        rows = repelem(row, numel(myLvl));
        cols = coordTransform_bgc2fp(iCol(rows), myLvl, sim);
        vals = squeeze(J_packed(row,myLvl));
        J_2d(row, cols) = vals;
    end
    if mod(row,7881) == 0
        fprintf('row %d ', row)
        toc;
    end
end
toc;
save('J', 'J', 'J_packed', 'J_2d');
toc;

spy(J_2d)
title('Partial of d(O2)/dt wrt to O2 everywhere');
xlabel('water level FP index)')
ylabel('water level FP index)')

fprintf('End of %s.m: %s\n', mfilename, datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));

end
