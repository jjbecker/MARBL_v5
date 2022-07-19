function sim = calc_global_moles_and_means(bgc, sim)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

disp('Calculate global average mean of all tracers');
% MARBLuses ---mmol---
sim.globalMean     = 1e3 *global_moles          (bgc.tracer, sim)./ sim.domain.V;       % (millimole m^3)
sim.globalMean_lvl = 1e3 *global_moles_per_layer(bgc.tracer, sim)./ sim.domain.V_lvl;   % (millimole m^3)

% use all levels to calculate scale factor
sz  = [ numel(sim.domain.iwet_JJ) , size(bgc.tracer,3) ];
sim.tracerScaleFactor = ones(sz).*sim.globalMean;
sim.tracerScaleFactor = sim.tracerScaleFactor(:);

% use a per level mean to calculate scale factor
tmp = reshape(sim.globalMean_lvl, [1,size(sim.globalMean_lvl)]); % [1,20,32]
tmp = repmat( tmp, [sim.domain.num_wet_loc, 1] );   % [545,20,32]
tmp = packMarbl(tmp, sim.domain.iwet_JJ);           % [8500,32]
sim.tracerScaleFactor = tmp(:);                     % [8500*32]



sim.tracerScaleFactor = 1;




    function moles_lvl = global_moles_per_layer(tracer, sim)
        % global_moles calculates total moles of all tracers
        % MARBL fields or ---mmole---

        % convert fld  from e.g. (10441,24) to (ilat, ilon, lvl) = (91,180,24)

        c = nan+sim.domain.M3d;
        % mmoles(size(tracer,3)) = 0;

        for iTr = size(tracer,3):-1:1
            fld = tracer(:,:,iTr);
            % convert fld  from e.g. (iCol,iLvl) to (ilat, ilon, lvl) = (91,180,24)
            c(sim.domain.iwet_FP) = fld(sim.domain.iwet_JJ);    % c (iCol,iLvl) -> [iLat, iLon, iLvl]
            % c is now at [iLat, iLon, iLvl] like M3d

            % mmoles(iTr)       = sum( c .*sim.domain.dVt, 'all', 'omitnan');
            mmoles_lvl(:,iTr) = sum( c .*sim.domain.dVt, [1 2], 'omitnan');
        end

        % MARBLuses ---mmol---
        %
        % convert to moles
        %
        % moles     = mmoles     * 1e-3;  % mmol ->mol

        moles_lvl = 1e-3 *mmoles_lvl;  % mmol ->mol

    end

end

