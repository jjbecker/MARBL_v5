function J = calc_J_O2(x0,sim, bgc, time_series, forcing)
%calc_J Summary of this function goes here
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

disp([mfilename,'.m: Starting simplified Jacobian; takes 30 seconds/tracer (maybe?) ...']);

tic

sz = size(bgc.tracer);
J = zeros([sz(1) sz(2) sz(2)]);


% FIXME: always tricky to guess "delta"
% FIXME: need to have correct uinit for h; what is that? median x0(iTr)??? or...

% h = x0 /1e3;
x0 = packMarbl(bgc.tracer,sim.domain.iwet_JJ);
x0 = x0(:,7);
h = x0 *sqrt(eps);
% h = sqrt(eps) *ones(size(x0));

% calculate d(Tracer)/dt with tracer values x0

tmp = time_series;

% initial_moles = global_moles(bgc.tracer, sim);  % DEBUG
% f0 = f(x0, sim, bgc, time_series);                % sz=[379913,32]
f0 = calc_f(sim, bgc, time_series, forcing);    % sz=[379913,32]
% final_moles = global_moles(bgc.tracer, sim);    % DEBUG
% tend_moles = global_moles(unpackMarbl(f0, sim.domain.iwet_JJ,size(bgc.tracer)), sim);    

tName = tracer_names(0);    % no CISO tracers
for lvl = 1:sz(2)       % loop over all lvl

    dx = h(:,idx);

    x1 = x0;
    x1(:,idx) = x0(:,idx) +dx;

    % calculate d(Tracer)/dt with tracer values x1

%     f1 = f(x1,sim, bgc, time_series);
    f1 = calc_f(x1, sim, bgc, time_series, forcing);

    % deltaF = d(Tracer)
    df = f1 -f0;

%     fprintf('%s:m: norm of d(%s) = %.2g\n', mfilename, string(tName(idx)), norm(df(:,idx)));
    fprintf('%s:m: norm of d(%s)\t/d(%s) = %.2g\n', mfilename, string(tName(idx)), string(tName(idx)),norm(df(:,idx)./dx));
%     if idx == 7
%         keyboard
%     end

    %     its too late to be cute, just use a loop...
    for j=1:sz(2)
        J(:,j,idx) = df(:,j) ./dx;
        if idx == 7
            fprintf('%s:m: norm of d(%s)\t/d(%s) = %.2g\n', mfilename, string(tName(j)), string(tName(idx)),norm(df(:,j)./dx));
        end
    end
    % debug deltaF at loc 1 all tracers wrt to Fe
    %     loc = 152; iTr = 5;
    %     foo = deltaF./dx./x0;
    %     tmp = [f0(loc,:)', f1(loc,:)', deltaF(loc,:)', J(loc,:,idx)', foo(loc,:)', J(loc,:,iTr)'];

    % https://duckduckgo.com/?q=181+This+sparse+indexing+expression+is+likely+to+be+slow.&t=osx&ia=web
%     toc
end
% tmp = [f0(loc,:)', f1(loc,:)', deltaF(loc,:)', J(loc,:,idx)', foo(loc,:)', J(loc,:,iTr)'];

elapsedTime = toc;
disp(['J: ', num2str(elapsedTime*1, '%1.3f'),' (s) for partial of MARBL tendency (all tracers), WRT all tracers, at all locations']);

logJ=log10(abs(nonzeros(J(:))));
figure(504); histogram(logJ); title("hist(log10(J))");

% sz = size(J)tmp = zeros(prod(size()));
% tmp = zeros([sz(1)*sz(2),sz(3));
% for idx = 1:sz(2)       % loop over tracers
%     tmp(:,idx) = J(:,idx,:);
% end
end
