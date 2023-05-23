function [tracerError, tracerNorm, c_sol]  = unscrambleAndSaveParFor(newRestartFileName, sim, bgc, tmp_ierr, tmp_xsol, tmp_fnrm, parforIdxRange, tracer_cell )

% Unscramble results captured in random order by parfor loop. No need
% for a for loop, use array index on slices. Then replace initial values in
% restart file with single column solution.

tracerError = zeros([1, size(bgc.tracer,3)]);
tracerNorm  = zeros([1, size(bgc.tracer,3)]);
% singleColTracerSolution = 0 *c0;

tracerRange = sim.tracer_loop_idx (parforIdxRange);
tracerError (:, tracerRange) = tmp_ierr (:, parforIdxRange);
tracerNorm  (:, tracerRange) = tmp_fnrm (:, parforIdxRange);
% tmp_c_sol          (:, tracerRange) = tmp_xsol (:, parforIdxRange);

ierrLimit = 1;
ierrLimit = 2;

bad_par_idx = parforIdxRange( tmp_ierr >  ierrLimit )
badTracers  = sim.tracer_loop_idx ( bad_par_idx )

good_par_idx = parforIdxRange;
good_par_idx( bad_par_idx ) = []
goodTracers  = sim.tracer_loop_idx ( good_par_idx )

%%%
% Need moles for calc_global_moles_and_means() for bgc2nsoli()
% Need initial tracers for restart file

c0_nsoli = bgc2nsoli(sim, bgc.tracer);    % nsoli format; unitless; aka scaled FP
sz_bgc = [ numel(sim.domain.iwet_JJ) , size(bgc.tracer,3) ];
c_sol = reshape(c0_nsoli, sz_bgc);
c_sol (:, goodTracers) = tmp_xsol (:, good_par_idx);

tracer = nsoli2bgc(sim, bgc, c_sol);
% result is saved in a file, bgc and sim returned are bogus...
saveRestartFiles(sim, tracer, newRestartFileName);

% debug output
tmp = reshape(c0_nsoli, sz_bgc);
for idx = parforIdxRange
    col = sim.tracer_loop_idx (idx);
    fprintf('%s.m: (%s)\tierr = %d fnrm(r) = %-#13.7g norm(sol) = %-#10.7g norm(x0) = %-#10.7g norm(sol-x0) = %-#10.7g\n', ...
        mfilename, string(tracer_cell(idx)), tracerError(col), ...
        tracerNorm(col), norm(c_sol(:,col)), norm(tmp(:,col)), norm( tmp(:,col)-c_sol(:,col) ))
end

end % unscrambleAndSaveParFor function
