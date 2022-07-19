function [x] = bgc2nsoli(sim, bgc_tracer)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


x = packMarbl(bgc_tracer, sim.domain.iwet_JJ);
x = x(:)./sim.tracerScaleFactor;


end

