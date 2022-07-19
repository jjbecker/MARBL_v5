function c = selectedTracers(sim, c0, selection)
%selectedTracers select tracers from Nsoli() format
% c0 is ALL tracers, in 1d format.
% Need to reshape c0 to size = [ numWaterParcels, 32 ].
% Then select tracers 

numWaterParcels = numel(sim.domain.iwet_JJ);
numTracers = sim.bgc_struct_base.size.tracer(2);

c = reshape(c0, [numWaterParcels, numTracers]);

c = c(:,selection);

% back to 1d vector

c = c(:);

end