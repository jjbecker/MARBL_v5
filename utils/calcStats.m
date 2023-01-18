function [x0,tmpG_all] = calcStats(x0, tmpG_all, G_bgc, initial_moles, final_moles, selection, current_yr, callersName)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% Print annual stats for DEBUG and learning!

if nargin < 7
    callersName = mfilename;
end

tmpG = tmpG_all(:,selection);       % just selected cols
tmpG = tmpG(:);                         % nsoli format


tName = tracer_names(0);    % no CISO tracers
tendStr   = strjoin(tName(selection));
gStr = sprintf('G( %s )', tendStr);

fprintf(        '%s.m: Year %d                        %s\n',callersName,current_yr,strjoin(pad(tName,14)));
disp([callersName,  '.m: Year ',num2str(current_yr),' ',['sol(' tendStr ')'],' Moles start = ',num2str(initial_moles,'%-#15.7g')])
disp([callersName,  '.m: Year ',num2str(current_yr),' ',['sol(' tendStr ')'],' Moles end   = ',num2str(final_moles,'%-#15.7g')])
disp([callersName,  '.m: Year ',num2str(current_yr),' ',['sol(' tendStr ')'],' Moles delta = ',num2str(final_moles-initial_moles,'%-#15.7g')])
ppm = ((final_moles-initial_moles)./ final_moles *1e6);
disp([callersName,  '.m: Year ',num2str(current_yr),' ',['sol(' tendStr ')'],' Moles (ppm) = ',num2str(ppm,'%-#15.7g')])
normG = vecnorm (tmpG_all);
%         disp([callersName,  '.m: Year ',num2str(current_yr),' ',['sol(' tendStr ')'],' norm G      = ', num2str(max(abs(tmpG_all)),'%-#15.7g')])
disp([callersName,  '.m: Year ',num2str(current_yr),' ',['sol(' tendStr ')'],' norm(G,2)   = ', num2str(normG,'%-#15.7g')])
normG = vecnorm (tmpG_all,inf);
disp([callersName,  '.m: Year ',num2str(current_yr),' ',['sol(' tendStr ')'],' norm(G,inf) = ', num2str(normG,'%-#15.7g')])
fprintf(        '%s.m: Year %d                        %s\n',callersName,current_yr,strjoin(pad(tName,14)));

% selection = [ ...
%     find( strcmp(tName,'SiO3') ) ];     % #3
figure (700); scatter(x0(:,selection),tmpG); title(strjoin(["scatter(",gStr,", ",strjoin(tName(selection)),")"]));    xlabel(strjoin(tName(selection)));   ylabel(gStr); grid on
figure (701); plot(tmpG);       title(strjoin(["plot(",gStr,")"]));         xlabel('idx FP');                        ylabel(gStr); grid on
figure (702); qqplot(tmpG);     title(strjoin(["qqplot(",gStr,")"])); grid on
figure (601); histogram(tmpG);  title(strjoin(["histogram(",gStr,")"]));    xlabel(gStr);                            ylabel('Count'); grid on
figure (602); clf(602)
% whisker plot each water column
% cross at median
% https://stackoverflow.com/questions/20224615/matlab-boxplots
maxData = max((G_bgc),[],2,'omitnan')';
minData = min((G_bgc),[],2,'omitnan')';
medData = median((G_bgc),2,'omitnan')';
stdData = std(G_bgc(:),'omitnan');
x = 1:numel(maxData);

hold on

line([x; x], [minData; maxData])    % vertical whiskers at each col#
plot(x, medData, '+')               % ticks on whiskers at median

yline(+stdData)                     % reference line at std
yline(-stdData)

hold off

xlabel('Water Column');
tmpStr = strjoin(["Extremes of",gStr,"in each water column"]);
ylabel(tmpStr);
title(tmpStr);


% figure(211); plot(colMax /std(colMax));
% title(strjoin(["max ( abs (",gStr,")) / std()"]));         xlabel('water col#');                        ylabel(strjoin(["abs ( ",gStr,") / std(abs)"])); grid on
% figure(212); histogram( colMax /std(colMax));   set(gca,'YScale','log')
% title(strjoin(["max ( abs (",gStr,")) / std()"]));  xlabel(strjoin(["max ( abs (",gStr,")) / std()"]));    ylabel("log(cnt)"); grid on
% figure(213); qqplot( colMax );
% colMax = max((G_bgc),[],2);
% figure(111); plot(colMax /std(colMax)); title(strjoin(["max ( ",gStr,")) / std()"]));         xlabel('water col#');                        ylabel(strjoin([" ",gStr,") / std(abs)"])); grid on
% figure(112); histogram( colMax /std(colMax));   set(gca,'YScale','log')
% title(strjoin(["max ( ",gStr,")) / std()"]));   xlabel(strjoin(["max ( ",gStr,")) / std()"]));   ylabel("log(cnt)");     grid on
% figure(113); qqplot( colMax );
%
% % Find the abs of each water column, and max histogram
% colMax = max(abs(G_bgc),[],2);
% figure(211); plot(colMax /std(colMax));
% title(strjoin(["max ( abs (",gStr,")) / std()"]));         xlabel('water col#');                        ylabel(strjoin(["abs ( ",gStr,") / std(abs)"])); grid on
% figure(212); histogram( colMax /std(colMax));   set(gca,'YScale','log')
% title(strjoin(["max ( abs (",gStr,")) / std()"]));  xlabel(strjoin(["max ( abs (",gStr,")) / std()"]));    ylabel("log(cnt)"); grid on
% figure(213); qqplot( colMax );
