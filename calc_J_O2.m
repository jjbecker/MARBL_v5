function J = calc_J_O2(x0, sim, bgc, time_series, forcing, MTM)
fprintf('%s.m: Start at %s\n', mfilename, datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
%calc_J Summary of this function goes here
%   Detailed explanation goes here
%
% To avoid mysterous memory leaks, need to use "time_step_ann()" which
% requires a accurate struct for time_series, and obviously needs correct
% forcing and treansport.
%
%
%   x0 Input is O2 in "fp" format eg. sz=[8500,1] MARBL wants sz=[545,20,32]
%
% calculate partial of MARBL with respect to O2 in its own col
%
% f = MARBL() = tendency = time rate of change of all tracers at all loc
%
% J = simple numberical approximation of Jacobian(f)
%
% 3 index array:
%
%   first is location (aka water column),
%   second is index of tendency
%   third is index of tracer being changed
%
%   J(:,i,j) = d(MARBL(O2(i))) / d(tracer(O2(j)) @ all water col, lvl i,
%   wrt O2 on lvl j

%
%   J(1,2,3) = d(MARBL(2)) / d(tracer(3) @ loc 1
%
%            = d(d(NO3)/dt)/ d(SiO3)     @ loc 1

fprintf('%s.m: Simplified Jacobian; takes ~600 seconds/tracer\n',mfilename);

tStart = tic;

sz = size(bgc.tracer);
J = zeros([sz(1) sz(2) sz(2)]);


% FIXME: always tricky to guess "delta"
% FIXME: need to have correct uinit for h; what is that? median x0(iTr)??? or...

% h = x0 /1e3;
% x0 = packMarbl(bgc.tracer,sim.domain.iwet_JJ);
% x0 = x0(:,7);
h = x0 *sqrt(eps);
h = unpackMarbl(h, sim.domain.iwet_JJ,[7881,60,1]);
% h = sqrt(eps) *ones(size(x0));

% calculate d(Tracer)/dt with tracer values x0

% tName = tracer_names(0);    % no CISO tracers
f0 = calc_f(x0, sim, bgc, time_series, forcing, MTM,1);
for lvl = 1:sz(2)       % loop over all lvl

    dx = h*0;
    dx(:,lvl) = h(:,lvl);
    dx = dx(sim.domain.iwet_JJ);

    x1 = x0 +dx;

    % calculate d(Tracer)/dt with tracer values x1

    f1 = calc_f(x1, sim, bgc, time_series, forcing, MTM,1);

    df = f1 -f0;
    idx = 7;

    %     its too late to be cute, just use a loop...
    J(:,:,lvl) = df(:,:,idx)./h(:,lvl);
    %     for j=1:sz(2)
    % %         J(:,j,lvl) = df(:,j) ./(dx+eps);
    %         J(:,j,lvl) = df(:,j,idx)./h(:,j);
    % %         if idx == 7
    % %             fprintf('%s:m: norm of d(%s)\t/d(%s) = %.2g\n', mfilename, string(tName(j)), string(tName(idx)),norm(df(:,j)./dx));
    % %         end
    %     end
    % https://duckduckgo.com/?q=181+This+sparse+indexing+expression+is+likely+to+be+slow.&t=osx&ia=web
end

elapsedTime = toc(tStart);
fprintf('%s.m: %1.3f (s) for partial of MARBL (O2 only) tendency WRT tracers (O2 only), on all levels of same column, and for all columns\n', mfilename, elapsedTime);

logJ=log10(abs(nonzeros(J(:))));
figure(455); histogram(logJ); title("hist(log10(J))");

fprintf('End of %s.m: %s\n', mfilename, datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
end
