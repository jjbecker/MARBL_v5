% %%
% 
% file_a = '/Users/jj/Desktop/relax/restart_790_relax.mat';
% file_b = '/Users/jj/Desktop/relax/restart_791_relax.mat';
% 
% tracer_a = load(file_a, '-mat').tracer;
% tracer_b = load(file_b, '-mat').tracer;
% 
% sim      = load(file_a, '-mat').sim;
% M3d      = sim.domain.M3d;
% 
% % tracer_names(0)
% fld  = 5;   % Fe
% lvl  = 1;
% 
% fld  = 1;   % PO4
% lvl  = 10;
% 
% %%
% % load tracer from files, regardless of the file extension..
% 
% %% On jj's mac...
% file_a = '/Users/jj/Desktop/restart_781_relax.mat';
% file_b = '/Users/jj/Desktop/restart_791_relax.mat';
% 
% %% On GP
% % file_a = '/DFS-L/DATA/primeau/jjbecker/Data/restart_0_1_output/restart_781_relax.mat';
% % file_b = '/DFS-L/DATA/primeau/jjbecker/Data/restart_0_1_output/restart_791_relax.mat';
% 
% 
% %%
% disp 'loading large files... (takes ~60 sec)';
% 
% tracer_a = load(file_a, '-mat').tracer;
% tracer_b = load(file_b, '-mat').tracer;
% 
% % DEBUG
% fprintf('size(''%s'').tracer = %s\n', file_a, sprintf('%d ', size(tracer_a)))
% fprintf('size(''%s'').tracer = %s\n', file_b, sprintf('%d ', size(tracer_b)))
% 
% 
% %%
% % Need M3d from either restart file..
% sim      = load(file_a, '-mat').sim;
% M3d      = sim.domain.M3d;
% 
% 
% %%
% % Plot the data
% 
% myDebugPlots(tracer_a, tracer_b, lvl, fld, M3d)
% 

%%
function [] = myDebugPlots(tracer_a, tracer_b, lvl, fld, M3d)

%% Contour plot of Tracer A

myContourPlot(tracer_a, lvl, fld, M3d, 201, "Tracer A");

%% Contour plot of Tracer B

myContourPlot(tracer_b, lvl, fld, M3d, 202, "Tracer B");
% mySurfacePlot(tracer_b, lvl, fld, M3d, 202, "Interpolated 'surf plot' Tracer B");

%% Contour plot of Tracer B - Tracer A

myContourPlot(tracer_b -tracer_a, lvl, fld, M3d, 203, "(Tracer B - Tracer A)");

%% Make a cool 3D surf plot of Tracer A - Tracer B we can spin around...

mySurfacePlot(tracer_b -tracer_a, lvl, fld, M3d, 204, "Interpolated (Tracer B - Tracer A)");
% mySurfacePlot(tracer_b , lvl, fld, M3d, 204, "3D Surface Tracer B-A");

end % myDebugPlots


%%
function [] = myContourPlot(tracer, lvl, fld, M3d, figNum, myTitle)

% Convert from MARBL coordinates to xyz

wet_loc  = find(M3d(:,:,1));

tmp = nan(size(M3d,1),size(M3d,2)); 
tmp(wet_loc) = tracer(:,lvl,fld);

figure(figNum)

contourf (tmp,'LineColor','none'); colorbar;
% y = 1:size(tmp,1);
% x = 1:size(tmp,2);
% Z = tmp(:);
% contourf (x,y,tmp); colorbar;
% surf (tmp, 'EdgeColor', 'none', 'FaceColor', 'interp'); view(2); colorbar;
surf (tmp, 'EdgeColor', 'none', 'FaceColor', 'flat'); view(2); colorbar;

% Title and labels...

name  = tracer_names(0);
units = tracer_units(0);
nameUnits = strcat(name," ",units);

% title(myTitle+" at level "+lvl+" of "+nameUnits(fld), 'Interpreter', 'none');
title(nameUnits(fld)+" at level "+lvl+" of "+myTitle, 'Interpreter', 'none');
xlabel ("ilon");
ylabel ("iLat");

end % myContourPlot


%%
function [] = mySurfacePlot(tracer, lvl, fld, M3d, figNum, myTitle)

% Convert from MARBL coordinates to xyz

wet_loc  = find(M3d(:,:,1));

tmp = nan(size(M3d,1),size(M3d,2)); 
tmp(wet_loc) = tracer(:,lvl,fld);

figure(figNum)

% surf (tmp, 'EdgeColor', 'none', 'FaceColor', 'interp'); view(15,55); colorbar;
surf (tmp, 'EdgeColor', 'none', 'FaceColor', 'interp'); view(-20,27); colorbar;

% Title and labels...

name  = tracer_names(0);
units = tracer_units(0);
nameUnits = strcat(name," ",units);

% title(myTitle+" at level "+lvl+" of "+nameUnits(fld), 'Interpreter', 'none');
title(nameUnits(fld)+" at level "+lvl+" of "+myTitle, 'Interpreter', 'none');
xlabel ("ilon");
ylabel ("iLat");

end % myContourPlot


%%
function [my_cell] = tracer_names(lciso_on)
%tracer_names MARBL "short name" of all tracers

% if MARBL is shutdown then mex_marbl_driver('tracer_sname', i) crashes...
%    just do a lot of typing...

my_cell{ 1} = 'PO4';
my_cell{ 2} = 'NO3';
my_cell{ 3} = 'SiO3';
my_cell{ 4} = 'NH4';
my_cell{ 5} = 'Fe';
my_cell{ 6} = 'Lig';
my_cell{ 7} = 'O2';
my_cell{ 8} = 'DIC';
my_cell{ 9} = 'DIC_ALT_CO2';
my_cell{10} = 'ALK';
my_cell{11} = 'ALK_ALT_CO2';
my_cell{12} = 'DOC';
my_cell{13} = 'DON';
my_cell{14} = 'DOP';
my_cell{15} = 'DOPr';
my_cell{16} = 'DONr';
my_cell{17} = 'DOCr';
my_cell{18} = 'zooC';
my_cell{19} = 'spChl';
my_cell{20} = 'spC';
my_cell{21} = 'spP';
my_cell{22} = 'spFe';
my_cell{23} = 'spCaCO3';
my_cell{24} = 'diatChl';
my_cell{25} = 'diatC';
my_cell{26} = 'diatP';
my_cell{27} = 'diatFe';
my_cell{28} = 'diatSi';
my_cell{29} = 'diazChl';
my_cell{30} = 'diazC';
my_cell{31} = 'diazP';
my_cell{32} = 'diazFe';
if lciso_on == 1
    my_cell{33} = 'DI13C';
    my_cell{34} = 'DO13Ctot';
    my_cell{35} = 'DI14C';
    my_cell{36} = 'DO14Ctot';
    my_cell{37} = 'zootot13C';
    my_cell{38} = 'zootot14C';
    my_cell{39} = 'sp13C';
    my_cell{40} = 'sp14C';
    my_cell{41} = 'spCa13CO3';
    my_cell{42} = 'spCa14CO3';
    my_cell{43} = 'diat13C';
    my_cell{44} = 'diat14C';
    my_cell{45} = 'diaz13C';
    my_cell{46} = 'diaz14C';
end
end % tracer_names


%%
function [my_cell] = tracer_units(lciso_on)
%tracer_names MARBL "units" of all tracers

% if MARBL is shutdown then mex_marbl_driver('tracer_sname', i) crashes...
%    just do a lot of typing...

my_cell{ 1} = '(mmol/m^3)';
my_cell{ 2} = '(mmol/m^3)';
my_cell{ 3} = '(mmol/m^3)';
my_cell{ 4} = '(mmol/m^3)';
my_cell{ 5} = '(mmol/m^3)';
my_cell{ 6} = '(mmol/m^3)';
my_cell{ 7} = '(mmol/m^3)';
my_cell{ 8} = '(mmol/m^3)';
my_cell{ 9} = '(mmol/m^3)';
my_cell{10} = '(meq/m^3)';
my_cell{11} = '(meq/m^3)';
my_cell{12} = '(mmol/m^3)';
my_cell{13} = '(mmol/m^3)';
my_cell{14} = '(mmol/m^3)';
my_cell{15} = '(mmol/m^3)';
my_cell{16} = '(mmol/m^3)';
my_cell{17} = '(mmol/m^3)';
my_cell{18} = '(mmol/m^3)';
my_cell{19} = '(mg/m^3)';
my_cell{20} = '(mmol/m^3)';
my_cell{21} = '(mmol/m^3)';
my_cell{22} = '(mmol/m^3)';
my_cell{23} = '(mmol/m^3)';
my_cell{24} = '(mg/m^3)';
my_cell{25} = '(mmol/m^3)';
my_cell{26} = '(mmol/m^3)';
my_cell{27} = '(mmol/m^3)';
my_cell{28} = '(mmol/m^3)';
my_cell{29} = '(mg/m^3)';
my_cell{30} = '(mmol/m^3)';
my_cell{31} = '(mmol/m^3)';
my_cell{32} = '(mmol/m^3)';
if lciso_on == 1
    my_cell{33} = '(mmol/m^3)';
    my_cell{34} = '(mmol/m^3)';
    my_cell{35} = '(mmol/m^3)';
    my_cell{36} = '(mmol/m^3)';
    my_cell{37} = '(mmol/m^3)';
    my_cell{38} = '(mmol/m^3)';
    my_cell{39} = '(mmol/m^3)';
    my_cell{40} = '(mmol/m^3)';
    my_cell{41} = '(mmol/m^3)';
    my_cell{42} = '(mmol/m^3)';
    my_cell{43} = '(mmol/m^3)';
    my_cell{44} = '(mmol/m^3)';
    my_cell{45} = '(mmol/m^3)';
    my_cell{46} = '(mmol/m^3)';
end
end % tracer_names

