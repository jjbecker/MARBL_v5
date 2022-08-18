function J = calc_J_full(f,x0,sim, bgc, time_series, forcing)
%calc_J Summary of this function goes here

numWetLoc = numel(sim.domain.wet_loc);  % = 8500 in 8x8 grid, 379913 in 3x3
nmWetLvl  = numel(sim.domain.zt);       % = 20 in 8x8 grid, 60 in 3x3
numTracer = size(bgc.tracer,3);         % = 32 in MARBL circa 2020

%   Input is in "fp" format eg. 8x8 grid sz=[8500,32] MARBL wants
%   sz=[545,20,32] = [col,lvl,tr]
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
%   J(:,i,j,idx_lvl) = d(MARBL(i)) / d(tracer(j on level k(i)) @ all locations
%
%   J(1,2,3,4) = d(MARBL(2)) / d(tracer(3 on level 4(1)) @ loc(1)
%
%            = d(d(NO3)/dt)/ d(SiO3 on level 4(1))     @ loc(1)

dbstop if error

disp([mfilename,'.m: Starting full Jacobian; takes long time...]']);
ticAll = tic;

% Do one easy part of this index festival...

f0 = f(x0,sim, bgc, time_series, forcing);     % = size(M3d)

% Simple so far, use ppm of each tracer as step size, then forward diff...

% dh = 1e-6;    % ???

dh = sqrt(eps);

% h = dh *abs(x0);
%
% But! Ultimately we need to make a "stencil" that picks out part of h on
% certain levels and for certain tracers. Avoid much confusion and don't
% even try to do this in "packed" aka "iwet" coordinates. for each level we
% need an "h" that only has that part of "h"
%
% Convert h to coordinates of [waterColumnIdx, waterLevel, tracerIdx]

% FIXME: need to address x0() elements == 0; supposed to be 1e-15
h_3d  = unpackMarbl( dh*abs(x0), sim.domain.iwet, size(bgc.tracer)); % [545,20,32]

% x0_3d = unpackMarbl( x0,         sim.domain.iwet, size(bgc.tracer)); % [545,20,32]
% x0_3d(isnan(x0_3d)) = 0;

% Tricky! We need PER LEVEL -and- per tracer. Need to pack/unpack, etc.

% J_3d = zeros( [545 20 32 20 32] );  % size = [8500 32 [20 32]] = 272000*20*32 = 175M for 8x8 grid!!!
% delta_tendency = zeros([545 20 32]);

J = sparse(8500*32, 8500*32);
rows_visted = zeros(8500*32,1);
myName = tracer_names(sim.lciso_on);
for idx_tr = 1: size(x0,2)  % loop over input tracers, but not all !!
% for idx_tr = 1:23 % FIXME: 23? why 23?
    ticLevel = tic;
    
    parcelTouched = zeros(8500,1);
    for idx_lvl = 1:20
        disp(['dMARBL()/d(',char(myName(idx_tr)),'(lvl ',int2str(idx_lvl),'))']);
        
        % find fp coordinates of each parcel in each col on this level
        iParcel = coordTransform_bgc2fp(1:545, idx_lvl, idx_tr, sim);
        % find iParcel ~= 0, then pick them out of array
        wet_iParcel = iParcel(iParcel);
        
        % double check that we end up with all elements after all levels
        parcelTouched(wet_iParcel) = 1;
        numWetParcelThisLevelAndTracer = sum(sim.domain.M3d(:,:,idx_lvl),'all');
        if (numel(wet_iParcel) ~= numWetParcelThisLevelAndTracer)
            disp(['numel, min, max of wet_iParcel: ',int2str(numel(wet_iParcel)),', ', int2str(min(wet_iParcel)),', ',  int2str(max(wet_iParcel))]);
            keyboard
        end
        
        
        % Cut out "stencil" for this level and tracer and replace "nan"
        % with 0 as we are in effect going to divide entire MARBL output
        % with a fraction of tracers in other parts so nans will not line
        % up. Trust me. See below in "dx_lvl_tr" caclulation
        %
        % Use dx which is zero except for "slice" in call to f(x)
        
        dx_slice = h_3d(:,idx_lvl,idx_tr);              % [545,1]
        dx = zeros(size(h_3d));                         % [545,20,32]
        dx(:,idx_lvl,idx_tr) = dx_slice;                % 3d zero+patch of h
        
        dx = packMarbl(dx, sim.domain.iwet);            % [8500,32]
        
        % Now comes tricky part! Need to divide all "tiles" in result of
        % "f(x)" with slice. Many ways to do that but replicating tiles
        % everywhere is easy to debug and fast "enough
        %
        % This is tricky. One way to think about it is that we are doing
        % "block" matrix math. Need to divide blocks in f(x+dx) by block
        % dx. We do that be making a matrix with block of dx everywhere.
        % However matrix "f(x)" is not actually in any special format, just
        % plain old double. There are many ways, this is easier to debug.
        
        dx_full   = repmat( dx_slice, 1, 20, 32);   % madness! [545,20,32]
        
        dx_lvl_tr = packMarbl( dx_full, sim.domain.iwet);   % [8500,32]
        
        
        % The rest is easy (er)...
        %
        % THe small slice is used in call to f, and its divided by full
        %
        % Mostly dividing 0 by h from a small slice. Other choices exist.
        
        % Forward difference
        
        fp = f( x0 +dx, sim, bgc, time_series);         % f(x+dx)
        df = fp -f0;
        
        % Central difference seems to give exactly same result, skip cost?
        % fm = f( x0 -dx, sim, bgc, time_series);       % f(x-dx)
        % df = (fp -fm)/2;
        % df_dx = df ./dx_lvl_tr;
        
        
        % Finally we get Jacobian of all tracers on all levels, wrt this
        % tracer on this level, by dividing d(all tracers) on all  levels
        % by that input delta on one level and tracer.
        
        df_dx = df ./dx_lvl_tr;         % [8500, 32] / [8500 32]
        % FIXME: remove nan without checking where they are??? dangerous.
        df_dx(isnan(df_dx)) = 0;
        
        % Unpack blocks we so carefully made. We can store them in a way
        % that is efficient to use inside nsoli(). Maybe. Need [545,20,32]
        
        df_dx_3d  = unpackMarbl( df_dx, sim.domain.iwet, size(bgc.tracer));
        % FIXME: remove nan without checking where they are??? dangerous.
        df_dx_3d(isnan(df_dx_3d)) = 0;
        
        % What is this intermediate result? We added h to one tracer on one
        % level and captured change on all levels for all tracers.
        %
        % So we turned a number in to 20x32 numbers at all 8500 locations.
        %
        % To get result at one location for one tracer we would
        % multiple residual for all tracers at a single location and
        % all levels tracers
        
        
        % FIXME: need to pack J_3d into [8500,32] format; somehow...
%         disp(['old nnz(J_3d,,,,) = ', int2str(nnz(J_3d))]);
%         disp(['nnz(df_dx) = ', int2str(nnz(df_dx))]);
%         disp(['nnz(df_dx_3d) = ', int2str(nnz(df_dx_3d))]);
%         J_3d(:,:,:,idx_lvl,idx_tr) = df_dx_3d;
        my_rows = coordTransform_bgc2fp(1:545, idx_lvl, idx_tr, sim);
        for k = 1:545
            if (my_rows(k)>0)
                
                foo = 0*df_dx_3d;
                foo(k,:,:) = df_dx_3d(k,:,:);
                this_row = packMarbl(foo, sim.domain.iwet);            % [8500,32]
                
                rowInJ = my_rows(k)+(idx_tr-1)*8500;
                if (rows_visted(rowInJ) ~=0)
                    keyboard
                end
                if (nnz(df_dx_3d(k,:,:)) ~= nnz(this_row))
                    keyboard
                end
                
                rows_visted(rowInJ)=1;
                tmp = reshape(this_row,1,8500*32);
                % super slow to write zeros to a sparse matrix which then
                % deletes them...
%                 J(rowInJ,:) = tmp;
                % Find nonzero in this row of J, then write them
                myNnz = find(tmp);
% https://www.mathworks.com/matlabcentral/answers/363341-sparse-indexing-expression-is-likely-to-be-slow                J(rowInJ,myNnz) = tmp(myNnz);
                J(rowInJ,myNnz) = tmp(myNnz);
            end
            
        end
%         nnz_J = nnz(J_3d);
%         disp(['nnz(J_3d,,,,) = ', int2str(nnz_J)]);
        nnz_J = nnz(J);
        disp(['nnz(J) = ', int2str(nnz_J)]);
        disp(['nnz(rows_visted) = ', int2str(nnz(rows_visted))]);
%         if (nnz_J ~= nnz_J)
%             keyboard
%         end
        figure(666); spy(J); title('Jacobian of MARBL Tendency w.r.t. Tracers')
        
        disp(['nnz dMARBL() /dTracer(idx_lvl ',int2str(idx_lvl),', idx_tr ',int2str(idx_tr),')              : ', int2str(nnz(df_dx_3d(:))), ' of possible ', int2str(numel(df_dx_3d(:))/1e6),' k 3d (NOT all wet)' ])
        
        disp(['avg over col: nnz dMARBL() /dTracer(idx_lvl ',num2str(idx_lvl),', idx_tr ',int2str(idx_tr),'): ', num2str(nnz(df_dx_3d(:))/numWetParcelThisLevelAndTracer,'%.2f'), ' on average in each of ',int2str(numWetParcelThisLevelAndTracer),' water parcels' ])
        disp(['total       : nnz J(:)                                  : ', num2str(nnz(J(:))/1e6),' M' ])
        
        max_df_dx = max(abs(df_dx),[],'all');
        if (max_df_dx>0)
            debug_J_full();
        else
            disp([...
                '"h" on level ', int2str(idx_lvl), ...
                ' and tracer ',   int2str(idx_tr), ...
                ' does NOT create terms in J on other levels'])
            
        end % debug
        disp(' ')
        
        logJ = log10(abs(nonzeros(df_dx_3d(:))));
        figure(400); histogram(logJ); title(['log10 abs(df_dx_3d) Tracer ',int2str(idx_tr),' all loc, thru Level ',int2str(idx_lvl),'...'], 'Interpreter', 'none');xlabel('log10 abs(J)');ylabel('Cnt');
        %     autoArrangeFigures(0,0);
        
        %         pause(1)
        
    end % of level loop
    
    if (sum(parcelTouched) ~= numel(sim.domain.iwet))
        keyboard % FIXME: num of iParcel ~= num of iwet is a big bummer
    end
    elapsedTime = toc(ticLevel);
    
    disp(['J: ', num2str(elapsedTime*1, '%1.3f'),' (s), for partial of MARBL tendency, all loc and levels for tracer ',int2str(idx_tr)]);
    disp(' ')
    
    logJ = log10(abs(nonzeros(J(:))));
    figure(500); histogram(logJ); title(['log10 abs(J) all loc and levels, thru tracer ',int2str(idx_tr),' ...'], 'Interpreter', 'none');xlabel('log10 abs(J)');ylabel('Cnt');
    autoArrangeFigures(0,0);
    
    % keyboard
    
end % of tracer loop

elapsedTime = toc(ticAll);

disp(['J: ', num2str(elapsedTime*1, '%1.3f'),' (s) for partial of MARBL tendency, all tracers on all levels...']);

logJ = log10(abs(nonzeros(J(:))));
figure(500); histogram(logJ); title("log10 abs(J) all levels", 'Interpreter', 'none');xlabel('log10 abs(J)');ylabel('Cnt');
autoArrangeFigures(0,0);

dbclear if error

    function debug_J_full()
        
        % Find value and location of biggest element in J.
        idx_nsoli = min(unique(find(abs(df_dx) == max_df_dx)));
        
        % Unpack...
        iTr = 1 +floor(idx_nsoli /8500);
        iFp = 1 +rem(fix(idx_nsoli),8500);
        
        % Convert to human readable
        debugFigNum = 0;    % do NOT plot location, for speed
        [iCol, iLvl, iLat, iLon, latitude, longitude, depth] = coordTransform_fp2bgc(iFp, sim, debugFigNum);
        
        disp(['max: J = d(MARBL(packed iFp ',int2str(iFp),' = (iCol, iLvl, idx_tr) = (',int2str(iCol),', ',int2str(iLvl),', ',int2str(iTr),...
            ')) wrt (idx_lvl ',int2str(idx_lvl),', idx_tr ',int2str(idx_tr),')'])
        disp(['max: J = ',num2str(max_df_dx), ...
            ', at iFp ', int2str(iFp),' = (iLat: ',num2str(iLat), ', iLon: ',num2str(iLon), ', iLvl: ',num2str(iLvl),')'])
        disp(['max: J = ',num2str(max_df_dx), ...
            ', at iFp ', int2str(iFp),' = (latitude ',num2str(latitude,'%.1f'), 'N, longitude ',num2str(longitude,'%.1f'), 'E, depth ',num2str(depth,'%.0f'), ' m)'])
        
        debugCrossLevelTerms();
        
        function debugCrossLevelTerms()
            
            other_lvl = setdiff(1:20, idx_lvl);
            tmp2 = squeeze(df_dx_3d(:,other_lvl,:));                            % all on OTHER levels
            tmp3 = squeeze(max(abs(squeeze(df_dx_3d)),[],'omitnan'));           % sum over all levels max(abs(
%             tmp4 = tmp3(other_lvl,setdiff(1:32, idx_tr));                          % max(abs( of OTHER tr on OTHER idx_lvl
%             tmp5 = sum(squeeze(sum(abs(tmp2),'omitnan')),2)';                   % sum on OTHER levels over locs and tracers
            % tmp5>0 indicates levels with some change in OTHER tracer
            
            % FIXME: need to count number of wet locations on this level and divide
            % into number of non zero J to see if it makes sense
            
            if (sum(abs(tmp2),'all','omitnan'))  % there are cross level terms
                
                disp([...
                    '"h" on level ', int2str(idx_lvl), ...
                    ' and tracer ',   int2str(idx_tr), ...
                    ' makes ',       int2str(nnz(tmp3(other_lvl,:))),' terms of J -on OTHER levels- non zero, for ',      ...
                    int2str(nnz(sum(tmp3(other_lvl,setdiff(1:32, idx_tr)),1))), ' OTHER tracers.'])
                
                
                % tmp6(idx_tr)>0 indicates that tracer changed on some other level
                tmp6 = sum(squeeze(sum(abs(tmp2),'omitnan')),1)';               % sun of all OTHER level of sum of tracer
                
                tmp7 = tmp3; tmp7(idx_lvl,:) = 0;   % find big "off level" elements
                tmp8 = squeeze(sum(abs(tmp7),'omitnan'));   % find tracer
                tmp9 = squeeze(sum(abs(tmp7),2,'omitnan')); % find idx_lvl
                
                disp(['h(idx_lvl ',int2str(idx_lvl),', idx_tr ',int2str(idx_tr),'): idx of tracers that change on other level: ', int2str(setdiff(find(tmp6), idx_lvl)')])
                disp(['h(idx_lvl ',int2str(idx_lvl),', idx_tr ',int2str(idx_tr),'): idx of levels with tracers that changed  : ', int2str(setdiff(find(tmp9), idx_lvl)')])
                
                disp([...
                    'Max abs(J) on OTHER level: J([', ...
                    int2str(find(tmp9 == max(tmp9))'),'],[', ...
                    int2str(find(tmp8 == max(tmp8)) ),']) = ', ...
                    num2str(max(tmp8))])
                disp([...
                    'Max(other lvl)/max(all) = : J([', ...
                    int2str(find(tmp9 == max(tmp9))'),'],[', ...
                    int2str(find(tmp8 == max(tmp8)) ),']) ' ...
                    'wrt (idx_lvl ',int2str(idx_lvl),', idx_tr ',int2str(idx_tr),') = ', ...
                    num2str(max(tmp8)/max_df_dx, 3)])
                
                %             disp('J([idx_lvl,tracer])) != 0')
                %             disp(sparse(tmp3))
                %
                %                 [iParcel,~] = coordTransform_bgc2fp(iCol, idx_lvl, -987777, sim);
            else
                disp([...
                    '"h" on level ', int2str(idx_lvl), ...
                    ' of tracer ',   int2str(idx_tr), ...
                    ' does NOT create terms in J on other levels'])
                
            end % of if there are cross level terms
            
        end % of debugCrossLecavelTerms()
        
    end % end of debug_J_all

end % of calc_J_full


%         for i=1:545
%             for j=1:20
%                 for k=1:32
%                     delta_tendency(i,j,k)=0;
%                     for l=1:32
%                         for m=1:20
%                             delta_tendency(i,j,k)=delta_tendency(i,j,k) + J(i,j,k,l,m)*x0_3d(i,m,l);
%                         end
%                     end
%                 end
%             end
%         end
%         foo = sum(squeeze(df_dx_3d(1,:,:)) * squeeze(x0_3d(1,:,:))',2);
%         delta_tracer = delta_tendency *sim.dt;
%         dTracerFrac = delta_tracer./x0_3d;
%         logDelta = log10(abs(nonzeros(dTracerFrac(:))));
%         figure(200); histogram(logDelta); title("log10 abs(dT *delta_tendency /x0) all levels", 'Interpreter', 'none');
%         autoArrangeFigures(0,0);
%    keyboard
% done but for endless debug and test

