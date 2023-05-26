function [ierr, fnrm, myRestartFile_x0, sol, c0, sim, bgc] = marbl_solve(x0, c0, sim, bgc, f)

% Get initial G of all tracers to pick sensible rtol for selected tracer???
% [r0,G0,x1] = calc_G(x0,c0,sim,bgc,time_series,forcing,MTM,PQ_inv);


% DEBUG stops when residual < tol, perfect residual = 0? does not stop run
atol = sim.epsilon;     % stop if norm(drift,2) < sqrt(eps) (noise)
rtol = sim.rtol;        % stop if norm(drift,2) < 10% of G(x0)

% maxfeval or maxit == 1 is pointless, if we do not input f(x0).
% 
% Nsoli() always evaluates value of f(x0), so 1 means wasted a call to phi,
% t always quits AFTER that call, unless we pass in f(x0).
%
% % However! When setting limits used to test convergence (atol and rtol)
% it is crucial to remember that these tests are done on (probably 
% preconditioned) "f(x)" not "x"!!! 
% 
% So it is extrememtly likely "f(x0)" is needed before we ever get here to 
% set useful values for atol and rtol. Since it takes an hour or 3 or 24 to
% calculate f(x0) it is very inefficient, but clearer code, to just pass in
% x0 to nsoli() even if f(x0) is already know because we might know that 
% nsoli() is going to evaluate f(x0 a 1,000 times so you cares if nsoli() 
% calculates f(x0) an extra time? 
% 
% But additional complication is that for some scheme that call nsoli(), 
% might only need a few iterations of f(x); probably to have practical run
% times. 
% 
% IN addition, it can be case that for a vector valued search, we 
% always start with same x0 even if we are solving only for certain 
% components of vector so calculating f(x0) repeatedly is a real waste 
% of time.
% 
% In all of these cases we want to set a small limit for numbers of calls 
% to f(x); like 5, and it is very valuable to hand in f(x0) along with x0.
%
% maxit is NOT maxmium number of nonlinear iterations because of possible
% line searches; maybe many line searches. maxit = Max number of G() or 
% phi() IF AND ONLY IF G always decreases. Must be>2 because need at least 
% 1 for a call to get f0 and 1 for a Newton, if we are calculating f(x0) in
% nsoli().
% 
% COULD BE maxit = maxit*maxit, (if steplength reductions do work on last 
% try), and that could be a VERY long time. maxfeval bounds total cnt of 
% calls to G to definite value.
%
%   realistically: maxit ~ total calls to G ~ ~maxit+maxarm
%   because line search rarely works, and search fails quickly...
%
% maxdim = maximum number of x in history (aka emperical Jacobian). 
% size(x0_hist) = [393913, maxdim];
%
% Given x is a large array, do NOT want large maxdim, and it could never be
% case we need maxdim > maxfeval which is max POSSIBLE total calls to G
%
% first 3 of these parms are used by (modified) brsola,
% rest are specific to nsoli

maxit    = sim.maxfeval;

maxdim   = 40;            % default is 40 in brsola()
maxdim   = min(sim.maxfeval, maxdim);

% used only by nsoli()
% etamax = 0.9;           % maximum error tol for residual in inner iteration, default = 0.9
% lmeth  = 2;             % Nsoli() method 2 = GMRES(m), not used by brsola().
% restart_limit = 10;     % max number of restarts for GMRES if lmeth = 2, default = 20;

tName = tracer_names(0);    % no CISO tracers
tracerStr   = strjoin(tName(sim.selection));

parms  = [maxit, maxdim, sim.maxfeval, sim.selection];

% NOTE: f has one input, x, and a bunch of parameters; e.g. bgc etc. If G()
% changes parameters subsequent cals to f do NOT get new parameter value 
% calculated in bgc;  e.g. if G() changes bgc internally bgc in calls to 
% f() are ones here, right now, not updated ones in G().


% if we or caller  already knows value of f(x0)=f0 vector, pass it in.
%
% However, if we do not include f0 on call of brsola(), brsola will
% calculate f0 itself. 
% 
% This is useful  depending of if this code is being called in a loop
% where each of 32 tracers has exactly same x0 and have exactly same f0.
% 
% But if we are not in a loop or debugging repeatedly using same intial 
% condition, let brsola do call.


% This call could take several hours or even days...
% f0=feval(f,x0);


% This is finally call to Newton Krylov solver!

fprintf('%s.m: (%s) ###### norm(x0) = %g\n', mfilename, tracerStr, norm(x0));

% [sol,it_hist,ierr,x0_hist] = brsola(x0, f, [atol,rtol],parms,tracerStr,f0);
[sol,it_hist,ierr,x0_hist] = brsola(x0,f,[atol,rtol],parms,tracerStr);

fnrm = it_hist(end,1);
fprintf('%s.m: (%s) ###### fnrm = %g\n', mfilename, tracerStr, fnrm)

% Save results, which are many and sort of tricky. Might want x0, the
% initial condition that solves f, or one forward integral of that (x1),
% or...

% remember that "sol" of nsoli() is an x0 value !!!
%
% if ierr == 0, or 1, last column of x0_hist is sol
%
% if ierr == 2 then last column x0_hist is x when iarm == maxarm,
% and next to last column x0_hist is x when armijo iterations
% started, and second to last column, might be x if there was no
% restart, or only one line search; in short if ierr = 2 x0_hist is
% hard to understand...

if ierr >1
    fprintf('%s.m: (%s) ###### ierr = %d replacing sol with x0_hist(:,end-1)\n', mfilename, tracerStr, ierr);
%     keyboard
    norm(sol)
    norm(x0_hist( :,end  ))
    norm(x0_hist( :,end-1))
    sol = x0_hist(:,end-1);
    [myMin, myIdx] = min(it_hist(:,1));
    sol = x0_hist(:,myIdx);
end

tmpStr=sprintf('%g ', vecnorm(x0_hist));
fprintf('%s.m: (%s) ###### norm(x0_hist) = %s\n', mfilename, tracerStr, tmpStr)

tmpStr=sprintf('%g ', it_hist(:,1));
fprintf('%s.m: (%s) ###### it_hist(:,1)  = iterations of norm(preconditioned(drift = x1-x0)\n', mfilename, tracerStr)
fprintf('%s.m: (%s) ###### it_hist(:,1)  = %s\n', mfilename, tracerStr, tmpStr)
% fprintf('%s.m: (%s) ###### it_hist(:,1) = norm(%s\n', mfilename, tracerStr, tmpStr)

fprintf('%s.m: (%s) ###### norm(sol)     = %g\n', mfilename, tracerStr, norm(sol))

% save complete bgc with -sol- tracers
x0_bgc  = replaceSelectedTracers(sim, c0, sol, sim.selection);
bgc.tracer = nsoli2bgc(sim, bgc, x0_bgc);   % marbl format x0

myRestartFile_x0 = sprintf('%s/restart_%d_%s_sol_x0.mat', sim.outputRestartDir, round(sim.start_yr), strjoin(tName(sim.selection)));
saveRestartFiles(sim, bgc.tracer, myRestartFile_x0);

sol_fname = sprintf('%s/sol_%s_ierr_%d_x0', sim.outputRestartDir, string(tName(sim.selection)), ierr);
fprintf('%s.m: Saving just sol, ierr, it_hist, x0_hist in %s\n', mfilename, sol_fname);
if sim.debug_disable_phi
    fprintf('%s.m: (%s) ###### phi() is short circuited skip sol_fname save\n',mfilename, tracerStr)
    return
end
save(sol_fname, 'sol', 'ierr', 'it_hist', 'x0_hist', '-v7.3','-nocompression');

% "sol" from brsola() will be last x0 that had smaller norm than previous 
% call of f(x), --not-- last one caclualated, which is "restart_x0.mat". In
% other words, last x might be a failed line search and sol is last value
% of x that was converging.
% 
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
%         saveRestartFiles(sim, bgc.tracer, myRestartFile_x1);

end