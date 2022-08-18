function J = calc_J_simplified(fsdfsfs,x0, sim, bgc, time_series, forcing)
%calc_J Summary of this function goes here
% calculate partial of MARBL with respect to tracers
%
%   Input:
%
% f = MARBL() = tendency = time rate of change of all tracers at all loc
% x0 = one or more tracers is in "nsoli" format eg. sz=[8500,32] 
%           MARBL wants sz=[545,20,32]!!
%
% sim, bgc, time_series, forcing are the usual inputs to G, phi, and
% time_step
%
%   Output:
%
% J = simple numberical approximation of Jacobian(f)
%
% 3 index array:
%
%   first is location (aka water column),
%   second is index of tendency
%   third is index of tracer being changed
%
%   J(:,i,j) = d(MARBL(i)) / d(tracer(j) @ all locations
%
%   J(1,2,3) = d(MARBL(2)) / d(tracer(3) @ loc 1
%
%            = d(d(NO3)/dt)/ d(SiO3)     @ loc 1

disp([mfilename,'.m: Starting simplified Jacobian; takes 30 seconds/tracer (maybe?) ...']);

tic

x0 = packMarbl(bgc.tracer,sim.domain.iwet_JJ);
x0 = x0(:,7);
sz = size(x0);              % e.g. [num_wet_loc, numel(selected)]
J = zeros([sz, sz(2)]);     % e.g. [num_wet_loc, numel(selected), numel(selected)]


% FIXME: always tricky to guess "delta"
% FIXME: need to have correct uinit for h; what is that? median x0(iTr)??? or...

% h = x0 /1e3;
% h = x0 *sqrt(eps);
h = sqrt(eps);

% calculate d(Tracer)/dt with tracer values x0

initial_moles = global_moles(bgc.tracer, sim);  % DEBUG
% f0 = f(x0, sim, bgc, time_series);                % sz=[379913,32]
f0 = calc_f(x0, sim, bgc, time_series, forcing);    % sz=[379913,32]
final_moles = global_moles(bgc.tracer, sim);    % DEBUG
tend_moles = global_moles(unpackMarbl(f0, sim.domain.iwet_JJ,size(bgc.tracer)), sim);    

tName = tracer_names(0);    % no CISO tracers
for idx = 1:numel(bgc.tracer)       % loop over all tracers

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
