function [PQ_inv, J_FP] = calc_PQ_inv(sim, bgc, time_series, forcing, MTM)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% 
% From NK-EXAMPLES: make a sparse operator that restores the surface to zero with a time scale of tau
% msk  = sim.domain.M3d;     % wet == 1, dry == 0
% iwet = sim.domain.iwet_FP;
% % tau = 24 * 60^2;                % (sec/d)
% % tau = 7*sim.const.sec_d;        % (sec/wk)
% tau = sim.const.sec_h;          % (sec/hr)
% temp = msk;
% temp(:,:,2:end) = 0;
% R =  d0( temp(iwet) / tau );   % (1/sec)
% T = 12*num_step_per_month*dt;

tName = tracer_names(0);    % no CISO tracers
fprintf('%s.m: calculating PQ_inv for tracers: [%s]\n', mfilename, strjoin(tName(sim.selection)) )

fprintf('\n%s.m: Takes about 20 mins to average transport, compute a single tracer J, and mfactor PQ...\n', mfilename)
fprintf('%s.m: Averaging transport...\n', mfilename)
tic
Q =  MTM(1).A + MTM(1).H + MTM(1).D;
for k = 2:12
    Q = Q + MTM(k).A + MTM(k).H + MTM(k).D;
end
Q = Q/12; % annually averaged transport and surface restoring (aka birth)
fprintf('%s.m: %1.3f (s) to average transport\n', mfilename, toc);


J_FP = calc_J_Single_Tracer(sim, bgc, time_series, forcing, MTM);

PQ = Q +J_FP;

tStart = tic;
fprintf('%s.m: Factoring 23 GB preconditioner PQ = Q +J = mean( MTM(k).A + MTM(k).H + MTM(k).D ) +J_FP...\n', mfilename)

%     keyboard
PQ_inv = mfactor(sim.T*PQ);     % 22 GB!  aka 2.4133e+10 bytes
elapsedTime = toc(tStart);
fprintf('%s.m: %1.3f (s) to mfactor of PQ\n', mfilename, toc(tStart));

fprintf('%s.m: Saving 23 GB preconditioner...\n', mfilename)
save (strcat('sol/',strjoin(tName(sim.selection)),'_QJ'), 'PQ_inv', 'PQ', 'Q', 'J_FP')
clear Q J PQ

elapsedTime = toc(tStart);
fprintf('%s.m: %1.0f (s) to factor and save PQinv \n',mfilename, toc(tStart));
end
