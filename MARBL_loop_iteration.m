function [ surface, interior ] = MARBL_loop_iteration(dt, n, col_num, surface, interior)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% t1 = tic;

% FIXME: time step can generate tracers that are slightly negative, however
% FIXME: negative tracer is nonsense. MARBL might limit to 0. Force that.
% FIXME: those negative values get here so global mass is conserved.
% FIXME: limiting needs to be  done here before we call MARBL.
%
% Code should be
%
%      interior.tracer = max ( -sim.epsilon, interior.tracer);
%
% but no globals allowed in parallel code. So it's hard coded...
% 
% interior.tracer = max ( -sqrt(eps), interior.tracer);
% 
% neither of surface_update or interior_update changes tracers, 
% so don't have to check after updates

% Save tracers and state before update.
% Use to get mid point of time step for surface flux update

% FIXME: MARBL did NOT change tracers during update_interior, so old_tracer
% == new tracer, but it does update pH or state...

interior.state_old  = interior.state;
interior.tracer_old = interior.tracer;

% Integrate forcing to get tendency
% Calculate surface flux and update interior with it

if (n == 1  && col_num == 1)
    
    % FIXME: first call to update can be garbage.
    % FIXME: Uninitialized vars in mex ?
    
    ignore = 1;
    % FIXME: use these throw away states and diags????
    %     interior = update_interior ( interior, ignore );
    %     surface  = update_surface  ( surface, interior, ignore );
    update_interior ( interior, ignore );   % does NOT change tracers
    update_surface  ( surface,  interior, ignore ); % "" no change
end
ignore = 0;

% neither update NOT changes tracers
interior = update_interior(           interior, ignore );
surface  = update_surface ( surface,  interior, ignore );

interior.state   (:,1) = surface.state;
interior.tendency(:,1) = interior.tendency(:,1) +surface.tendency;

% elapsedTime = toc(t1);
% disp(' ');disp(['t1: ', num2str(elapsedTime*1000, '%1.1f'),' (ms) for 1 tendency'])

% FIXME: do NOT Flux limit output from MARBL. right???

end
