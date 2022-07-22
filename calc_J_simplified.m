function [J] = calc_J_simplified(f,x0,sim, bgc, time_series)
%calc_J Summary of this function goes here
%
%   Input is in "fp" format eg. sz=[8500,32] MARBL wants sz=[545,20,32]
%
% calculate partial of MARBL with respect to tracers
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
%   J(:,i,j) = d(MARBL(i)) / d(tracer(j) @ all locations
%
%   J(1,2,3) = d(MARBL(2)) / d(tracer(3) @ loc 1
%
%            = d(d(NO3)/dt)/ d(SiO3)     @ loc 1

disp('Starting calc_J; takes 60 seconds...');

tic

[f0] = f(x0,sim, bgc, time_series);     % sz=[8500,32]

sz = size(x0);
J = zeros([sz, sz(2)]);


% FIXME: always tricky to guess "delta"
% FIXME: need to have correct uinit for h; what is that? median x0(iTr)??? or...

h = x0 /1e3;

for idx = 1:sz(2)       % loop over tracers
    
    dx = h(:,idx);
    
    x1 = x0;
    x1(:,idx) = x1(:,idx) +dx;
    
    f1 = f(x1,sim, bgc, time_series);
    
    df = f1 -f0;
    %     its too late to be cute...
    for j=1:sz(2)
        J(:,j,idx) = df(:,j) ./dx;
    end
    % debug df at loc 1 all tracers wrt to Fe
%     loc = 152; iTr = 5;
%     foo = df./dx./x0;
%     tmp = [f0(loc,:)', f1(loc,:)', df(loc,:)', J(loc,:,idx)', foo(loc,:)', J(loc,:,iTr)'];


% https://duckduckgo.com/?q=181+This+sparse+indexing+expression+is+likely+to+be+slow.&t=osx&ia=web                

end
% tmp = [f0(loc,:)', f1(loc,:)', df(loc,:)', J(loc,:,idx)', foo(loc,:)', J(loc,:,iTr)'];

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
