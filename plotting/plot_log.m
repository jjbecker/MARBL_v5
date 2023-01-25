function myFig = plot_log(myFig, myTitle, t_or_z, y_array, name, idx, depthNotTime)

% plot_time_series plot time evolution of array of tracers

figure(myFig)
myFig = myFig +1;

% tiledlayout requires R2019b or later

% a = 1;  b = 1;  c = -size(y_array,2); m = (-b+sqrt(b^2-4*a*c))/2/a; m = ceil(m); m*(m+b)
% tl = tiledlayout(m, m+b, 'TileSpacing','compact', 'Padding','compact');
tl = tiledlayout('flow','TileSpacing','compact','Padding','compact');

for i=1:size(y_array,2)
    ax(i) = nexttile(tl);
    plot(t_or_z, y_array(:,i)); % FIXME: Surprisingly this takes first of a multiple column array!
    %     ylabel (idx(i)+". "+name(i), 'Interpreter', 'none')
    ylabel (idx(i)+". "+name(i), 'Interpreter', 'none')
    if (depthNotTime)
        view([90 90]);
    else
        hold on
        plot([t_or_z(1) t_or_z(end)],[y_array(1,i) y_array(end,i)])
        if floor(t_or_z(end)/365) >0
            xline([1:floor(t_or_z(end)/365)]*365);
        end
        hold off
    end
    grid on
end
linkaxes(ax,'x')

title(tl, myTitle, 'Interpreter', 'none');
if (depthNotTime)
    ylabel (tl,"Depth (m)");
else
    xlabel (tl,"Time (day)");
end

