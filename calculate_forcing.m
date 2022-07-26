function [sim,bgc] = calculate_forcing(sim,bgc, timeStep)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% All of interior forcing is constant, for simplicity, except for PAR.
%
%  Calulate PAR everywhere, then "nan" from bottom to last element of array

% kmt          = bgc.kmt;
% [iLat, iLon] = ind2sub(size(sim.domain.M3d), sim.domain.wet_loc);
% lat          = sim.grd.YT(iLat);
% lon          = mod(180+sim.grd.XT(iLon),360)-180;
% 
% % PAR is surface 400-700(nm) -AIR-
% % Surface Shortwave (W/m^2)
% %
% % MARBL attenuates PAR with depth including "blooms" in layers above.
% %
% % Should have a GCM tell us weather (clouds, aka air mass, etc)
% %
% % https://www.pveducation.org/pvcdrom/properties-of-sunlight/air-mass
% 
% solar_constant = 1367;      % TOA (W/m^2)
% visible_400_700 = 0.383;    % 400-700 nm wavelength == PAR?
% % FIXME: should air mass be something like 1/cos(lat)
% air_mass = 0;               % Fraction lost to atmo
% 
% K = solar_constant *visible_400_700 *(1-air_mass) ; % 523.6 W/m^2
% 
% % calculate cos of zenith angle of sun everywhere, multiply by intensity.
% %
% % irradiance = K *cos_zenith(sim.JD, lat, lon)' ;
% irradiance = K *cos_zenith(0, lat, lon)' ;
% 
% % need this on every level for some reason, "nan" below bottom
% 
% num_lvl = size(bgc.forcing(:,:,2),2);
% irradiance = repmat(irradiance',[1,num_lvl]);
% 
% % % skip drama and use a for loop. Ohhhh no!!!
% %
% % FIXME: need to add bottom CORRECTLY
% % FIXME: why do we need to add bottom????
% % keyboard
persistent solar_force
persistent num_lvl
if isempty(solar_force)
    solar_force = load(strcat(myDataDir(),'/SOLAR_3hr_forcing.mat'), '-mat').solar_force;
    num_lvl = size(bgc.forcing(:,:,2),2);
end

timeStepInYear = 1 +mod( timeStep-1, round(sim.const.sec_y /sim.dt));

% Check if "dt" and size of solar match...
if ( round(sim.const.sec_y /sim.dt) == size(solar_force,1))
    idx = timeStepInYear;
else
    myFrac = (timeStepInYear -1) / round(sim.const.sec_y /sim.dt);
%     myFrac = myFrac + 1/2/round(sim.const.sec_y /sim.dt);
    idx = 1 +round(myFrac *size(solar_force,1));
end

% irradiance = solar_force(timeStepInYear,:);
% irradiance = repmat(irradiance',[1,num_lvl]);
irradiance = repmat(solar_force(idx,:)',[1,num_lvl]);

% for i=1:numel(iLat)
%     irradiance(i,bgc.kmt(i)+1:num_lvl) = nan;
% end


bgc.forcing(:,:,2) = irradiance;

end

