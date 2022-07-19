function moles = global_moles(tracer, sim)
% global_moles calculates total moles of all tracers
% MARBL fields or ---mmole---

c = packMarbl( tracer, sim.domain.iwet_JJ);

mmoles = sum( c .* sim.domain.dVt_FP);

% MARBL uses ---mmol--- we want MOLES

moles = 1e-3 *mmoles;  % mmol ->mol

end
