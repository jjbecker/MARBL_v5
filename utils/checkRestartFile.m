function [] = checkRestartFile(sim, bgc, forcing)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fprintf('\n%s.m: Display a WHOLE LOT of debug info to get a sanity check of input grids and tracers...\n\n',mfilename)
% Check it!
% squeeze(bgc.tracer(:,1,:));
% squeeze(bgc.forcing(:,1,:));
% squeeze(bgc.surf_forcing);
% squeeze(bgc.river_flux);
% Check it!

% average them over year
bgc.forcing      = 0*forcing(1).interior;
bgc.surf_forcing = 0*forcing(1).surf_forcing;
bgc.river_flux   = 0*forcing(1).river_flux;

% FIXME: obviously we should check if this input struct has monthly, annual
% yearly, or whatever members for forcing...
for ii=1:12 
    bgc.forcing      = forcing(ii).interior     +bgc.forcing      /12;
    bgc.surf_forcing = forcing(ii).surf_forcing +bgc.surf_forcing /12;
    bgc.river_flux   = forcing(ii).river_flux   +bgc.river_flux   /12;
end


% Next statement works only because "Areat" is organized as lat,lon, depth. 
% So wet locations which only has enough elements for first page of "Areat"
% accidently is top layer.
dA_m  = sim.grd.Areat(:,:,1);      % Surface Area of land and water
dA_m  = dA_m(sim.domain.wet_loc);  % Surface Area of water cols in m^2
dA_cm = dA_m *(1E+2)^2;            % cm^2
disp(['Surface Area of Ocean = ', num2str(sum(dA_m),3), ' m^2'])

% dV_m  = sim.domain.dVt;           % Volume of parcels land and water  m^3
tmp = reshape( sim.domain.dVt, [numel(sim.domain.M3d(:,:,1)), numel(sim.domain.M3d(1,1,:))] );
dV_m  = tmp(sim.domain.wet_loc,:,:);% Volume of parcels in water cols m^3
disp(['Volume of Ocean = ', num2str(sum(dV_m,[1 2]),3), ' m^3'])
disp(['Depth of Ocean = ', num2str(sum(dV_m,[1 2])/sum(dA_m),3), ' m'])

dNOx  = bgc.surf_forcing(:,7);      % nM/cm^2/s
mySum = sum(dA_cm.*dNOx);           % nM/s
mySum = mySum *sim.const.sec_y;     % nM/y
mySum = mySum /1e+9 /1e+12;         % Tmol/y
gm_mole = 14+16*2;                  % molecular weight of NOx as NO2
mySum = mySum *gm_mole;             % Tg/y
%  https://doi.org/10.1038/s41467-018-03553-w
disp(['global(surface NOx Flux) Absorbed = ', num2str(mySum,3), ' (Tg/y), expected 16 ??? '])

dNHy  = bgc.surf_forcing(:,8);      % nM/cm^2/s
mySum = sum(dA_cm.*dNHy);           % nM/s
mySum = mySum *sim.const.sec_y;     % nM/y
mySum = mySum /1e+9 /1e+12;         % nmole -> mol -> Tmol/y
gm_mole = 14+1*4;                   % molecular weight of NHy as NH4
mySum = mySum *gm_mole;             % Tg/y
disp(['global(surface NHy Flux) Emission(????) = ', num2str(mySum,3), ' (Tg/y), expected 13 ??? '])

dDust = bgc.surf_forcing(:,5);      % g/cm^2/s
mySum = sum(dA_cm.*dDust);          % g/s
disp(['mean(surface Dust Flux) = ', num2str(mySum/sum(dA_cm),3), ' (gm/cm^2/s), expected 3.485E-12 '])
mySum = mySum *sim.const.sec_y;     % g/y
mySum = mySum /1e+15;               % Pg/y
disp(['global(surface Dust Flux) = ', num2str(mySum,3), ' (Pg/y), '])

% % "...Fe Inputs dFe atms 5.0, seds 20.5, vents 5.0, rivers 0.3 Gmol/yr

dFe   = bgc.surf_forcing(:,6);      % IRON_FLUX mmol/m^2/s
% NOTE: dA_m not cm
mySum = sum(dA_m.*dFe);             % mM/s
globalAvgFeFlux = mySum / sum(dA_m);     % mmole/m^2/s
disp(['mean(surface Iron Flux) = ', num2str(globalAvgFeFlux,3), ' (nmol/s/cm^2), expected 4.9e-10?? '])
mySum = mySum *sim.const.sec_y;     % mM/y
mySum = mySum /1e+3 /1e+9;          % mmol -> mol -> Gmol/y
disp(['global surface Fe Flux = ', num2str(mySum,3), ' (Gmol/y), expected 5.4?? or 13.62?? '])

% % "...Fe Inputs dFe atms 5.0, seds 20.5, vents 5.0, rivers 0.3 Gmol/yr

dFeSed = bgc.forcing(:,:,6);        % FESEDFLUX nmol/cm^2/s
globalAvgFeSedFlux = sum(dFeSed .* dA_cm, 'all','omitnan')/sum(dA_cm);
disp(['mean(FESEDFLUX) = ', num2str(globalAvgFeSedFlux,3), ' (nmol/cm^2/s), expected 2.24e-07 (nmol/cm^2/s)?? '])

% Convert flux at surface between each and every layers into volume
% tendency in each layer. e.g. just divide value of flux by layer thickness
% Make sure we get same answer for annual average!

dFlux = dFeSed * 1E4;               % nmol/m^2/s
dFlux = dFlux ./sim.domain.dzt;    % nmol/m^3/s
mySum = sum(dV_m.*dFlux,'all','omitnan');   % nmole/s 
globalAvgFESEDFLUX = mySum / sum(dA_cm);    % nmole/s/cm^2
disp(['mean(FESEDFLUX) = ', num2str(globalAvgFESEDFLUX,3), ' (nmol/cm^2/s), expected 2.24e-07 (nmol/cm^2/s)?? '])
mySum = mySum *sim.const.sec_y;     % nmol/y
mySum = mySum /1E9 /1E9;       % nm -> m -> Gm
disp(['global(interior FESEDFLUX) = ', num2str(mySum,3), ' (Gmol/y), expected 25.5 = 20.5+5.0'])

% % "...Fe Inputs dFe atms 5.0, seds 20.5, vents 5.0, rivers 0.3 Gmol/yr

dRiver = bgc.river_flux(:,:);       % Fe_RIV_FLUX mmol/m^3 cm/s
mySum = sum(dA_m.*dRiver) /100;     % m^2 *mmol/m^3 cm/s *m/cm = mmole/s
% disp(['global(River Flux(mmol/s)) = ', num2str(mySum,3), ' (mmol/s), expected ??? '])
mySum = mySum *sim.const.sec_y;     % mmol/y
mySum = mySum /1e+3;                % mol/y
mySum = mySum /1e+9;                % Gmol/y
disp(['global(River Fe Flux) = ', num2str(mySum(5),3), ' (Gmol/y), expected 0.345???'])

tmp=tracer_names(0); cstr = tmp(1:17)';
disp(strcat('global(River Flux = ', {' '}, string(strjoin(cstr))))
disp(['River Flux(Gmol/y) = ', num2str(mySum(1:17),3), ' (Gmol/y), expected ??? '])

fprintf('...end of sanity check\n\n')

end