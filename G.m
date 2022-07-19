function [accept,r,x0] = G(x0,Npt,c0,sim,bgc,time_series,forcing,MTM,PQ_inv)
%UNTITLED Take an initial value of tracers, return change at end of a year
%   Detailed explanation goes here

persistent gFileCnt x0_prev
if isempty(gFileCnt)
    gFileCnt = 1;
%     fprintf('call #%d to G\n', gFileCnt);
    fprintf('first norm(x0   ) = %f\n', norm(x0         ));
% checkNegAndHisto(sim, x0, 100.0, 'x0', 900);
else
    gFileCnt = gFileCnt +1;
    fprintf('call #%d to G\n', gFileCnt);
    fprintf('norm(x0_prev    ) = %f\n', norm(x0_prev    ));
    fprintf('norm(x0         ) = %f\n', norm(x0         ));
    fprintf('norm(x0 -x0_prev) = %.6f\n', norm(x0 -x0_prev));
end
x0_prev = x0;


% Check for negative tracers

accept = 1;
% c0 = bgc2nsoli(sim, bgc.tracer);    % nsoli format; unitless; aka scaled FP
% [~,c0] = ck_constraints(c0,-2);     % corrected any negatives
% 
% tmp = x0;
% [accept,x0] = ck_constraints(x0,Npt);
% if (accept == false)
%     disp('accept = false; histo our input in fig(999)')
%     figure(999); histogram(tmp)
%     keyboard
%     r = NaN(size(x0));
%     return
% end % if not within constraints


% Combine current tracers being changed by Nsoli, aka "sim.selection",
% with ones not being fiddled with to create full set of tracer to input
% to phi() aka MARBL.

numWaterParcels = numel(sim.domain.iwet_JJ);
numTracers = sim.bgc_struct_base.size.tracer(2);
sz = [numWaterParcels, numTracers];

% "_bgc" in tracer name means it has all 32 tracers

x0_bgc = replaceSelectedTracers(sim, c0, x0, sim.selection);
bgc.tracer = nsoli2bgc(sim, bgc, x0_bgc);   % marbl format

    initial_moles = global_moles(bgc.tracer, sim);  % DEBUG
[sim, bgc, ~] = phi(sim, bgc, time_series, forcing, MTM);
    final_moles = global_moles(bgc.tracer, sim);    % DEBUG

x1_bgc = bgc2nsoli(sim, bgc.tracer); % unitless end of year values
% checkNegAndHisto(sim, selectedTracers(sim, x1_bgc, sim.selection), 100.0, 'x', 900+gFileCnt);


% depending on preconditioner used, might need all residuals in actual
% units, or just selected ones, or something else.

res = -reshape(x1_bgc -x0_bgc, sz);     % needed res size is sz; aka 32 col
res = res(:,sim.selection);                 % just selected cols
res = res(:);                           % nsoli format

% keep res, actual diff, and find the preconditioned(res), r
r = res(:);

% Precondition the residual

% % Normal case. Plot histo of difference of starting and ending tracers.
% %     figure(700); histogram(log10(abs(dx_selected(abs(dx_selected)>eps)))); xlabel('log10(all global mean normalized tracers'); title('G: Histogram in sim.selection x1-x0')
% % figure(700); histogram(log10(abs(dx_selected))); xlabel('log10(all global mean normalized tracers'); title('G: Histogram in sim.selection x1-x0')
% 
% % % Debug output
% % if (sim.logTracers)
% % %     small_plots(sim, time_series, sim.num_time_steps, sim.time_series_loc, sim.time_series_lvl);
% % end
% % disp(['Npt = ', int2str(Npt),': norm(x-x0)[selected] = ', num2str(norm(dx_selected), '%1.7g')])
% % figure(1123);
% % ecdf(log10(abs(dx_selected))); grid on
% % xlabel('log10(global mean normalized tracers');
% % title('G: Empirical CDF in sim.selection tracers (x-x0)')
% 
% 
% % Precondition residual...
% 
% % Apply precondition to residual. Don't waste time preconditioning tracers
% % not going back to nsoli. Code is a little tricky. Breaks full nsoli
% % vector into chunks by tracer.
% 
% % c3d = NaN(size(sim.domain.M3d));
% 
% % keyboard
% % for idx = 1:sz(2)
% %
% %     % convert tr from e.g. (10441,24) to (ilat, ilon, lvl) = (91,180,24)
% %
% %     if (~ismember(idx, sim.selection))
% %         % figure(1); histogram(log10(abs(tr))); title(int2str(idx));
% %     else
% %         % P   = inv(Q) -I
% %         % P*z = inv(Q)*resid -resid
% %         %     = mfactor(Q,resid) -resid
% %
% %         z = mfactor(factoredQ(idx),res(:,idx));
% %
% % figure(400+idx);scatter(res(:,idx),z); title('Scatter of residual and "z"');xlabel('Residual');ylabel('z = mfactor(Q,Residual)');
% %
% %         % P*z = z-resid?
% %
% %         w = P{idx} *z;
% % figure(500+idx);scatter(res(:,idx),w); title('Scatter of residual and "w=Pz"');xlabel('Residual');ylabel('z = mfactor(Q,Residual)');
% % figure(600+idx);scatter(res(:,idx),res(:,idx)+w); title('Scatter of residual and "w=Pz"');xlabel('Residual');ylabel('z = mfactor(Q,Residual)');
% %
% %         res(:,idx)  = w;            % FIXME: must be more to it....
% %
% %         %         figure(600+idx);
% %         %         d = tr;
% %         %         histogram(log10(abs(d(abs(d)>1e-6))));
% %         %         xlabel('log10(all global mean normalized tracers'); title(int2str(idx))
% %
% %     end
% %
% %     % c3d(sim.domain.iwetXyzFrancois) = tr;
% %     % c3d = (c3d.*sim.domain.dVt /sim.domain.V);
% %     % tr = c3d(sim.domain.iwetXyzFrancois);
% %
% % end
% %
% % % Keep only tracers we are manipulating and return to nsoli format.
% % res = res(:,sim.selection);
% % r = res(:);
% 
% % r = PQ_inv*(dx_selected .* sim.tracerScaleFactor);
% 
% 
% % dx_selected = dx_selected .* sim.tracerScaleFactor;
% % if (0)
% %     checkSum = 0;
% %     checkIdx = 0*dTr;
% %     r        = 0*dTr;
% %     max_r    = zeros(size(sim.domain.wet_loc'));
% %     for k = 1:sim.domain.num_wet_loc;
% %         % find elements in "nsoli() vector for this water column:
% %         % index into nsoli of everything in iCol "k"? there are numel()/32
% %         % water parcels in col, aka number of wet levels in this column
% %         %
% %         % Location of these elements in nsoli vector is hard to picture...
% % 
% %         % Have to be very careful. "Q_inv" is not full inverse, just the
% %         % elements of J for this column... So we can --NOT-- take whole row.
% % 
% %         tmp = bgc.tracer*0;
% %         tmp(k,:,sim.selection) = 1;
% %         idx_iCol = find(bgc2nsoli(sim,tmp));
% % 
% %         r(idx_iCol) = PQ_inv(idx_iCol,idx_iCol) *dTr(idx_iCol);
% % 
% %         % debug
% %         checkSum    = checkSum+numel(idx_iCol);
% %         checkIdx(idx_iCol) = 1;
% %         max_r(k)    = max(abs(norm(r(idx_iCol))));
% % 
% %     end
% %     % debug
% %     if (nnz(checkIdx) ~= numel(dx_selected))
% %         disp(['checkSum = ', num2str(checkSum)]);
% %         keyboard
% %     end
% %     figure(51);plot(log10(max_r+eps));title('log10(max_r)', 'Interpreter', 'none');xlabel('col');
% % 
% % end


% DEBUG
disp([mfilename,'.m: Moles  start of phi() = ',num2str(initial_moles,7)])
disp([mfilename,'.m: Moles  end of phi()   = ',num2str(final_moles,7)])
disp([mfilename,'.m: Moles  delta          = ',num2str(final_moles-initial_moles,7)])
ppm = ((final_moles-initial_moles)./ final_moles *1e6);
disp([mfilename,'.m: Moles  delta (ppm)    = ',num2str(ppm,7)])

fprintf('%s.m: Npt = %d P*dx norm    = %1.10g\n', mfilename, Npt, norm(r));
fprintf('%s.m: Npt = %d P*dx max     = %1.7g\n',  mfilename, Npt, max((r)));
fprintf('%s.m: Npt = %d P*dx min     = %1.7g\n',  mfilename, Npt, min((r)));
fprintf('%s.m: Npt = %d P*dx median  = %1.7g\n',  mfilename, Npt, median(r));
fprintf('%s.m: Npt = %d P*dx mean    = %1.7g\n',  mfilename, Npt, mean(r));
fprintf('%s.m: Npt = %d P*dx std     = %1.7g\n',  mfilename, Npt, std(r));
fprintf('%s.m: Npt = %d P*dx mad(avg)= %1.7g\n',  mfilename, Npt, mad(r,0));
fprintf('%s.m: Npt = %d P*dx mad(med)= %1.7g\n',  mfilename, Npt, mad(r,1));
tmp = replaceSelectedTracers(sim, c0, r, sim.selection);
res_moles = global_moles(nsoli2bgc(sim, bgc, tmp), sim);
res_moles = res_moles(sim.selection);
% res_moles ./ final_moles(sim.selection) *1e6
fprintf('%s.m: Npt = %d P*dx moles   = %1.7g\n',  mfilename, Npt, res_moles);

% if (0)
%     myGfile = sprintf('%s/G_%d.mat', sim.outputRestartDir, round(gFileCnt));
%     fprintf('%s.m: Saving "%s"...\n', mfilename,myGfile);
%     tracer = bgc.tracer;
%     copyfile( sim.inputRestartFile, myGfile);
%     save( myGfile, 'tracer', '-append' ); % overwrites tracer from input
%     save( myGfile, 'c', '-append' );      % overwrites tracer from input
%     save(myGfile,'-v7.3');
% end


end % G()




% % check for bogus output from MARBL, like sun is down everywhere in world.
%
% % [iParcel,iTr] = coordTransform_nsoli2fp(find(abs(dx)>0.010*x0_bgc), sim, 666);
% % myLimit = 0.001;
% % myTr=5;
% % [iParcel,iTr] = coordTransform_nsoli2fp(find(abs(dx)>myLimit), sim, 483);   % bad, all tracers
% % idx_selected  = ismember(iTr, sim.selection);                               % in sim.selection and bad
% % iParcel = iParcel(idx_selected);                                        % ""
% % iTr = iTr(idx_selected);                                                % tracers at bad place
% % tracer_names(0); disp(ans(unique(iTr)'))                                % names of bad and in sim.selection
% % figure(2233); histogram(iTr);title(['x-x0 is > ',num2str(myLimit)]);xlabel('Tracer number');
% % [iLat, iLon, iLvl, latitude, longitude, depth] = coordTransform_fp2xyz(iParcel(find(iTr==myTr))', sim, 800); title(['x-x0 is > ',num2str(myLimit),' and ITr ',int2str(myTr)])
% if (sum(abs(time_series.diag(:,:,151)),'all')<=0)
%     % MARBL crashed!
%     figure(999); histogram(x0)
%     % small_plots(sim, time_series, sim.num_time_steps, sim.time_series_loc, sim.time_series_lvl);
%     disp(' '); disp('MARBL appears to have crashed: interior diag(151) = "PAR()" is all zeroes');disp(' ');
%     keyboard
% end


% % % % create UNITLESS aka scaled with all 32 tracers initial condition sz=[num_col,32]
% % % sz = [ numel(sim.domain.iwet_JJ) , size(bgc.tracer,3) ];
% % % x0_bgc  = reshape(c0,sz);
% % % % replace IC values for selected tracers with UNITLESS aka SCALED value
% % % % from Nsoli()...
% % % x0_bgc(:,sim.selection) = reshape(x0, [size(sim.domain.iwet_JJ,1), numel(sim.selection)] );
% % % x0_bgc = x0_bgc(:);             % back to 1d vector
% if (Npt == 1)
%     %     recalculate Jacobian in Preconditioner...
%     sim = calc_global_moles_and_means(bgc, sim);
%     J = calc_J_full(@calc_f, packMarbl(bgc.tracer, sim.domain.iwet), sim, bgc, time_series);
%     Q_inv = calc_Q_inv(J, bgc, sim);
% end

% run sim for 1 year
% m = matfile(filename);
% [nrows,ncols] = size(m,'c.c1');
% listOfVariables = who('-file', 'census.mat');
% ismember('pop', listOfVariables) % returns true
% ismember('doesNotExist', listOfVariables) % returns false[sim, bgc, time_series] = phi(sim, bgc, time_series, forcing, MTM);


