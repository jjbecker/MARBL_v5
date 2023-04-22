% function plotWorstResidual(sim, fname)
function plotWorstResidual(x0, G, G_bgc_all, current_yr, callersName)

% load(fname, 'r');
myFig = 800;
figure (myFig)
figure (myFig); 
clf(myFig)

tl = tiledlayout('flow','TileSpacing','compact','Padding','compact');

myInc = 4;
for i=1:myInc:size(G_bgc_all,3)

    G_bgc = G_bgc_all(:,:,i);       % just selected cols
%     tmpG = tmpG(:);              % nsoli format


    tName = tracer_names(0);    % no CISO tracers
    tendStr   = strjoin(tName(i));
    gStr = sprintf('G( %s )', tendStr, 'Interpreter', 'none');


    ax(i) = nexttile(tl);

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


    grid on
end
% linkaxes(ax,'x')

% title(tl, myTitle, 'Interpreter', 'none');
% if (depthNotTime)
%     ylabel (tl,"Depth (m)");
% else
%     xlabel (tl,"Time (day)");
% end

end