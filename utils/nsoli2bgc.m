function [tracer] = nsoli2bgc(sim, bgc, x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

sz  = [ numel(sim.domain.iwet_JJ) , size(bgc.tracer,3) ];                 

% restore units and reshape to [iwet, num_tr]

tmp = reshape(x.*sim.tracerScaleFactor, sz);                  

% reshape to [col,lvl,tr]

tracer = unpackMarbl(tmp, sim.domain.iwet_JJ, size(bgc.tracer)); 

end

