function [PQ_inv, J_FP] = calc_PQ_inv(sim, bgc, time_series, forcing, MTM)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
% From NK-EXAMPLES: make a sparse operator that restores surface to zero with a time scale of tau
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
if sim.recalculate_PQ_inv && ~(numel(sim.disabledPreconditoners)>0 && ismember(tName(sim.selection), sim.disabledPreconditoners))
    fprintf('%s.m: Start at %s\n', mfilename, datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
    tName = tracer_names(0);    % no CISO tracers
    fprintf('%s.m: Calculating PQ_inv for tracers: [%s]\n', mfilename, strjoin(tName(sim.selection)) )
    fprintf('%s.m: Takes about 20 mins to average transport, compute a single tracer J, and mfactor PQ...\n', mfilename)
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
    % elapsedTime = toc(tStart);
    fprintf('%s.m: %1.3f (s) to mfactor PQ\n', mfilename, toc(tStart));

    fprintf('%s.m: Saving 23 GB preconditioner...\n', mfilename)
    % save(strcat(myDataDir(),'sol/',strjoin(tName(sim.selection)),'_QJ'), 'PQ_inv', 'PQ', 'Q', 'J_FP','-v7.3','-nocompression')
    save(strcat(myDataDir(),'sol/',strjoin(tName(sim.selection)),'_QJ'), 'PQ_inv','-v7.3','-nocompression')
    % save(strcat(myDataDir(),'sol/',strjoin(tName(sim.selection)),'_QJ'), 'PQ_inv','-v7.3')
    clear Q J PQ

    % elapsedTime = toc(tStart);
    fprintf('%s.m: %1.0f (s) to factor and save PQinv \n',mfilename, toc(tStart));
    fprintf('End of %s.m: %s\n\n', mfilename, datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
else
    tStart = tic;
    if sim.disable_ALL_Preconditioner || (numel(sim.disabledPreconditoners)>0 && ismember(tName(sim.selection), sim.disabledPreconditoners))
        fprintf('\n\n\t%s.m: ********* Replace preconditioner with 1 *********\n\n',mfilename)
        PQ_inv = 1
    else
        fprintf('\n%s.m: Loading ~30 GB(!) mfactored preconditioner PQ_inv from %s solution...\n', mfilename, strcat(string(tName(sim.selection))))
        load (strcat(myDataDir(),'sol/',strjoin(tName(sim.selection)),'_QJ'), 'PQ_inv');
    end
    fprintf('%s.m: %1.0f (s) to init sim and load PQinv \n',mfilename, toc(tStart));
end % calculate or load PQ_inv

end