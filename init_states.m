function [surf_state, state]= init_states ()

% Note: surface and interior have "level" and tracer indices swapped...
% Read correct sized array of garbge from MARBL, initialize that, 
% write it back...

state       = mex_marbl_driver('interior_tendency_saved_state')';
surf_state  = mex_marbl_driver('surface_flux_saved_state')';

% surface_flux_saved_state      has 2 states but just one level; ph, alt_ph
% interior_tendency_saved_state has 2 states for each level

state(:,:)      = 8.1821049785262; % (pH)
surf_state(:)   = 7.9830634698986; % (pH)

% FIXME: hack in a wild IC for testing
% state      = 0*state;       % (pH)
% surf_state = 0*surf_state;  % (pH)

end % init_states.m
