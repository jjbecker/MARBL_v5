function c = replaceSelectedTracers(sim, c0, x0, selection)
%replaceSelectedTracers replace selected with selected tracers
% c0 is ALL tracers, in 1d format.
% x0 is a subset of tracers; probably from from Nsoli()...
% Need to reshape c0 to size = [ numWaterParcels, 32 ].
% Then replace the selected tracers with the ones from Nsoli

numWaterParcels = numel(sim.domain.iwet_JJ);
numTracers = sim.bgc_struct_base.size.tracer(2);

c = reshape(c0, [numWaterParcels, numTracers]);

% x0 is a subset of tracers; probably from from Nsoli()...
% Need to reshape x0 from 1d to size = [ numWaterParcels, selection ]

x = reshape(x0, [numWaterParcels, numel(selection)] );

% Now replace selected tracers with tracers with Nsoli()

c(:,selection) = x;

% back to 1d vector

c = c(:);

end