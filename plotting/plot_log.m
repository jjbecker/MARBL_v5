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
    else                % annote plot with lines at each year
        hold on
%         plot([t_or_z(1) t_or_z(end)],[y_array(1,i) y_array(end,i)])
        if floor(t_or_z(end)/365) >0
            
            % xline([1:floor(t_or_z(end)/365)]*365, '-.k'); % doesn't work on pre 2022 Matlab
            for x_yrs = [1:floor(t_or_z(end)/365)]*365
                xline(x_yrs, '-.k')
            end

            dt = t_or_z(2)-t_or_z(1);
            steps_yr = 365/dt;
            for my_yr = 1:floor(t_or_z(end)/365)
                x0 = 1+(my_yr-1)*steps_yr;    x1 = 1+(my_yr-0)*steps_yr;
                plot(t_or_z([x0 x1]), y_array([x0 x1],i))
            end
            x0 = x1; x1 = numel(t_or_z);
            plot(t_or_z([x0 x1]), y_array([x0 x1],i))

        else
            plot([t_or_z(1) t_or_z(end)],[y_array(1,i) y_array(end,i)])
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

