function [x0_sol, c0, sim, bgc, time_series, forcing, MTM, PQ_inv, myRestartFile_x0] = marbl_solve(~, x0, c0, sim, bgc, time_series, forcing, MTM, PQ_inv)
fprintf('%s: Parameters nsoli()... \n', mfilename)
% The max POSSIBLE total calls to G (or phi) ~ maxit*maxarm
%   realistically: ~maxit <=  total calls to G <= ~maxit+maxarm
%   because line search rarely actually works
maxit  = 1+7;   % Max number of calls to phi() where norm(G) decreases (+1).
% >2 because need at least 1 for f0;
% maxmium number of nonlinear iterations *** NOT in line search ***
maxdim = 40;    % maximum number of x in history (aka emperical Jacobian)
maxarm = 10;     % >0 or always fails (ierr=2) or any one given itc, *** 1+max allowed calls to phi() *** with norm(G) not decreaasing
%   = maximum number of steplength reductions;
% POSSIBLE calls to phi COULD be ~maxit*maxit, if steplength reductions actually do any good,
% which could be a
%   VERY long time
% And == 1 is pointless because first one is
% evaluating new value of f0, so 1 means waste a
% call to phi and then quit anyway
% used only by nsoli()
% etamax = 0.9;           % maximum error tol for residual in inner iteration, default = 0.9
% lmeth  = 2;             % Nsoli() method 2 = GMRES(m), not used by brsola().
% restart_limit = 10;     % max number of restarts for GMRES if lmeth = 2, default = 20;

% first 3 of these parms are used by (modified) brsola, rest are specific to nsoli
parms  = [maxit,maxdim, maxarm];

% Get current drift of all tracers to pick a sensible rtol for selected tracer
% [r0,G0,x1] = calc_G(x0,c0,sim,bgc,time_series,forcing,MTM,PQ_inv);

atol = sqrt(eps);       % stop when norm(drift,2) < sqrt(eps) (numerical noise)
rtol = 1e-2;            % stop when norm(drift,2) <10% of of G(x0)
% atol = eps;     % DEBUG it stops when residual is < this, so perfect residual ==0 will not stop run :-)
% rtol = 0;       % DEBUG
% atol = 150
% % DEBUG it stops when residual is < this, so perfect residual ==0 will not stop run :-)

% remember that "sol" of nsoli() is an x0 value !!!

[x0_sol,it_hist,ierr,x_hist] = brsola(x0, @(x) calc_G(x,c0,sim,bgc,time_series,forcing,MTM,PQ_inv), [atol,rtol], parms);

% if ierr == 0, or 1, last column of x_hist is x0_sol
%
% if ierr == 2 then last column x_hist is x when iarm == maxarm,
% and next to last column x_hist is x when armijo iterations
% started, and second to last column, might be x if there was no
% restart, or only one line search; in short if ierr = 2 x_hist is
% hard to understand...

% save complete bgc with -sol- tracers
x0_bgc  = replaceSelectedTracers(sim, c0, x0_sol, sim.selection);
bgc.tracer = nsoli2bgc(sim, bgc, x0_bgc);   % marbl format x0
% bgc_sol_x0 = bgc;  % this has c0 for all tracer

tName = tracer_names(0);    % no CISO tracers
sol_fname = sprintf('%s/sol_%s_ierr_%d_x0', sim.outputRestartDir, string(tName(sim.selection)), ierr);
fprintf('%s.m: Saving just x0_sol, ierr, it_hist, x_hist in %s\n', mfilename, sol_fname);
save(sol_fname, 'x0_sol', 'ierr', 'it_hist', 'x_hist', '-v7.3','-nocompression');

myRestartFile_x0 = sprintf('%s/restart_%d_%s_sol_x0.mat', sim.outputRestartDir, round(sim.start_yr), strjoin(tName(sim.selection)));
[sim, bgc] = saveRestartFiles(sim, bgc, bgc.tracer, myRestartFile_x0);

% "sol" from brsola() will be last x that was smallest --not-- last
% one caclualated which is "restart_x0.mat"
%
% confusing? YES
% tracer = [7881, 60, 32]
% x, sol = [379913]
% c      = [379913, 32]
% c1     = 12157216]

%         % FIXME: this is NOT x1 of sol, it is last x1 calulated
%                 IOW it might be x1 from a line search and hence terrible
%
%         bgc_sol_x1.tracer = load(sprintf('%s/restart_x1.mat', sim.outputRestartDir), 'tracer').tracer;     % [78881, 60, 32]
%
%         c1 = bgc2nsoli(sim, bgc_sol_x1.tracer);    % nsoli format; unitless; aka scaled FP
%         c1 = reshape(c1, sz);                  % [393913,32]
%         x1_sol = c1(:,sim.selection);           % [393913] initial condition for Nsoli()
%         % x1_sol = x1_sol(:);                     % unitless
%         myRestartFile_x1 = sprintf('%s/restart_%d_%s_sol_x1.mat', sim.outputRestartDir, round(sim.start_yr), strjoin(tName(sim.selection)))
%         [sim, bgc] = saveRestartFiles(sim, bgc, bgc.tracer, myRestartFile_x1);

end