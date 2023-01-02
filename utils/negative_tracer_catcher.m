function [negativesFound] = negative_tracer_catcher(sim,bgc)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% ad hoc checking for crazy values
nameUnits = tracer_names(sim.lciso_on);

tracers_fp = packMarbl(bgc.tracer,sim.domain.iwet_JJ);

negativesFound = 0;

for fld=1:size(nameUnits')

    % scale tracer by global mean of that tracer
    % bottom NaN already removed by pack()

    %     scaledTr = tracers_fp(:,fld) ./ sim.globalMean(fld);
    scaledTr = tracers_fp(:,fld) ;

    % FIXME: find "bad points" but what does "bad" mean?

    negLimit = 0;
    negLimit = -sim.epsilon;
    negLimit = -10 *sqrt(eps);
    negLimit = -1e-6;

    idx = find(scaledTr < negLimit);
    cnt = numel(idx);

    % look for obvious errors...

    if max(scaledTr) > +1e+6
        disp(' ');
        disp '>>>>>>>>>> CRAZY positive tracer!'
        keyboard
    end
    if min(scaledTr) < -1e+6
        disp(' ');
        disp '>>>>>>>>>> CRAZY negative tracer!'
        %         keyboard
    end

    % find all bad, get worst of bad
    if cnt>0

        negativesFound = negativesFound +cnt;

        % index into whole array, not negatives.

        [~,sorted_idx] = sort(scaledTr(idx));
        idx = idx(sorted_idx);

        % Show some of worst of bad locations.
        % Possibly want first bad, rather than worst.

        [iCol, iLvl, iLat, iLon, lat, lon, depth] = coordTransform_fp2bgc( idx, sim);

        badBottomLvl = [];badBottomDepth= [];
        for i=1:cnt
            badBottomLvl(i)   = sim.domain.bottom_lvl(iLat(i),iLon(i));
            badBottomDepth(i) = sim.domain.bottom_depth(iLat(i),iLon(i));
        end

        disp(['All ', nameUnits{fld}, ': ', ...
            '  scaled value min: ', num2str(min(scaledTr)), ...
            ' med: ', num2str(median(scaledTr)), ...
            ' max: ', num2str(max(scaledTr))]);

        disp(['Bad ', nameUnits{fld}, ': cnt ', num2str(cnt), ...
            ' (', num2str(100*cnt/numel(scaledTr),2), '%) ', ...
            ', scaled value min: ', num2str(min(scaledTr(idx))), ...
            ' med: ', num2str(median(scaledTr(idx))), ...
            ' max: ', num2str(max(scaledTr(idx)))]);

        formatSpec = '%.0f';
        disp(['Bad ',  nameUnits{fld}, ': depth (m) ', ...
            ' min: ', num2str(min(depth),formatSpec), ...
            ' med: ', num2str(median(depth),formatSpec), ...
            ' max: ', num2str(max(depth),formatSpec)]);

        scaledZ = depth./badBottomDepth;
        formatSpec = '%.2f';
        disp(['Bad ',  nameUnits{fld}, ': depth as fraction of water depth', ...
            ' min: ', num2str(min(scaledZ),formatSpec), ...
            ' med: ', num2str(median(scaledZ),formatSpec), ...
            ' max: ', num2str(max(scaledZ),formatSpec)]);

        scaledZ = iLvl./badBottomLvl;
        disp(['Bad ',  nameUnits{fld}, ': depth as fraction of bottom lvl ', ...
            ' min: ', num2str(min(scaledZ),formatSpec), ...
            ' med: ', num2str(median(scaledZ),formatSpec), ...
            ' max: ', num2str(max(scaledZ),formatSpec)]);

        disp(['Most negative at loc #',num2str(iCol(1)), ...
            ' (iLat,iLon,iLvl): (', ...
            num2str(iLat(1)),',', ...
            num2str(iLon(1)),',', ...
            num2str(iLvl(1)), ')', ...
            ' (lat,lon,depth): (', ...
            num2str(lat(1), '%1.1f'),'N, ', ...
            num2str(lon(1), '%1.1f'),'E, ', ...
            num2str(depth(1), '%1.0f'),'m)'])

        numPrint = 5;
        rng = 1:min(numPrint, cnt);
        fprintf('Most negative %d locations...\n',numPrint);
        fprintf('(loc #,lvl #): ')
        fprintf('(%g, %g) ', [iCol(rng);iLvl(rng);]); disp ' ';

        fprintf('(lat,lon,depth): ')
        fprintf('(%1.1f, %1.1f, %1.0fm) ', [lat(rng);lon(rng);depth(rng)]); disp ' '

        % optionally plot location of bad points, but not a zillion

        if (cnt > 1e6)
            fprintf('It is going to take a few minutes to plot %d locations...\n', cnt);
            disp 'skipping plot'
        else
            debugFigMapNum = 1300+fld;
figure(debugFigMapNum);
tl = tiledlayout('flow','TileSpacing','compact','Padding','compact');
ax(1) = nexttile(tl);
            debugStr = sprintf('%s: Negative %s', mfilename, nameUnits{fld});
            [iCol, iLvl, iLat, iLon, lat, lon, depth] = coordTransform_fp2bgc( idx, sim, debugFigMapNum,'');
            title(debugStr, 'Interpreter', 'none'); % get ride of subscripts

%             debugFigHistNum = 1400+fld;
%             figure(debugFigHistNum)
ax(1) = nexttile(tl);
            histType = 'cumcount';
            debugStr = sprintf('%s<%s "%s" vs Depth (m)', ...
                nameUnits{fld}, num2str(negLimit,2),histType);
            histogram(depth,'Normalization',histType); grid on
            title(debugStr, 'Interpreter', 'none'); % get ride of subscripts

            %             pause(1)
        end
        disp ' '
    end % of cnt >0
end % of fld

fprintf('Found %d (%.0f%%) negatives in %.2e tracers\n', ...
    negativesFound, 100*negativesFound/numel(bgc.tracer(:)), numel(bgc.tracer(:)))
