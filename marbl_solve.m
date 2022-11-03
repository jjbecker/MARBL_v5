function [ierr, myRestartFile_x0, x0_sol, c0, sim, bgc, time_series, forcing, MTM, PQ_inv] = marbl_solve(x0, c0, sim, bgc, time_series, forcing, MTM, PQ_inv)
fprintf('%s: Parameters nsoli()... \n', mfilename)

% Get initial G of all tracers to pick sensible rtol for selected tracer???
% [r0,G0,x1] = calc_G(x0,c0,sim,bgc,time_series,forcing,MTM,PQ_inv);


% DEBUG stops when residual < tol, perfect residual = 0? does not stop run
atol = sqrt(eps);         % stop if norm(drift,2) < sqrt(eps) (noise)
rtol = 1e-2;              % stop if norm(drift,2) < 10% of G(x0)
% atol = eps;
% rtol = 0;       % DEBUG
% atol = 150;     %

% maxfeval or maxit == 1 is pointless.
% 
% Nsoli() always evaluates value of f(x0), so 1 means wasted a call to phi,
% t always quits AFTER that call.
%
% maxfeval is the absolute max number of calls to G. must be >2; need at 
% least 1 for f0 and then 1 for a Newton;
%
% maxitl is NOT maxmium number of nonlinear iterations because of possible
% line searches; maybe many line searches. maxit = Max number of calls to G
% or phi() IF AND ONLY IF G always decreases. 
% Must be>2 because need at least 1 for a call to get f0 and 1 for a Newton
% 
% COULD BE maxitl = maxit*maxit, (if steplength reductions do work on last 
% try), and that could be a VERY long time. maxfeval bounds the total cnt 
% of calls to G to definite value.
%
%   realistically: maxit ~ total calls to G ~ ~maxit+maxarm
%   because line search rarely works, and search fails quickly...
%
% maxdim = maximum number of x in history (aka emperical Jacobian). 
% size(x_hist) = [393913, maxdim];
%
% Given x is a large array, do NOT want large maxdim, and it could never be
% case we need maxdim > maxfeval which is the max POSSIBLE total calls to G
%
% first 3 of these parms are used by (modified) brsola,
% rest are specific to nsoli

maxfeval = 1+4;

maxit    = maxfeval;

maxdim   = 40;            % default is 40 in brsola()
maxdim   = min(maxfeval, maxdim);

% used only by nsoli()
% etamax = 0.9;           % maximum error tol for residual in inner iteration, default = 0.9
% lmeth  = 2;             % Nsoli() method 2 = GMRES(m), not used by brsola().
% restart_limit = 10;     % max number of restarts for GMRES if lmeth = 2, default = 20;

parms  = [maxit, maxdim, maxfeval];

[x0_sol,it_hist,ierr,x_hist] = brsola(x0, @(x) calc_G(x,c0,sim,bgc,time_series,forcing,MTM,PQ_inv), [atol,rtol], parms);
x0_sol_norm=norm(x0_sol)
norm_x_hist = vecnorm(x_hist)
if ierr >1
    fprintf('\n\nn%s.m: #$@ ierr = %d\n\n', mfilename,ierr);
%     keyboard
end

% remember that "sol" of nsoli() is an x0 value !!!
%
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