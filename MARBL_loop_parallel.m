function [bgc, time_series] = MARBL_loop_parallel ( n, sim, bgc, time_series )

% these variables need to be OUTSIDE of spmd loops to a normal struct of
% normal doubles, otherwise we end up with "composite", "distributed" or
% "codistrubted" depending exactly where we declare them, and worse a
% struct of entire grid is opposite of what we want.
%
% We want a struct of data from exactly one location.

[surface_base, interior_base] = init_bgc_struct_base(sim);

interior_base_client = interior_base;
surface_base_client  = surface_base;

interior_time_series_loc = interior_base;
surface_time_series_loc  = surface_base;

bgc.tendency = nan(size(bgc.tendency));

% Likewise we want to send a normal matrix, not a structure of matrix

Tracers_client      = bgc.tracer;
kmt_client          = sim.domain.bottom_lvl  (sim.domain.wet_loc);
Forcing_client      = bgc.forcing;
Surf_Forcing_client = bgc.surf_forcing;
States_client       = bgc.state ;

S_in_sz     = size ( States_client );
S_out_sz    = S_in_sz;
tend_sz     = size ( bgc.tendency );
% do NOT read diags in spmd, those huge arrays get copied by each worker
% diag_sz = size ( bgc.diag );
% surf_sz = size ( bgc.surf_diag );

% do NOT copy in a 100 MB of data for 3 numbers...

dt_client               = sim.dt;
num_threads             = sim.number_of_threads;
time_series_loc_client  = sim.time_series_loc;
logDiags_client         = sim.logDiags;
time_series_lab         = 0;


% distribute local copy to all workers

% FIXME: It would be a little complex, but "state" is unused anywhere but
% in MARBL, so S_out could be copied over S_in, and skip scatter/gather...

% FIXME: kmt and almost all of forcing are constant in time, should not
% need to scatter ever time step, not that this is particularily slow.
% However it does NOT work to just scatter on first pass for some reason

% constant inputs; -SHOULD- only need to be initialized and sent once.

% if (n == 1); ticBytes(gcp); end
spmd(num_threads)
    kmt_in = codistributed(kmt_client          , codistributor1d(1)); % 1 is by row
    F_in   = codistributed(Forcing_client      , codistributor1d(1));
    SF_in  = codistributed(Surf_Forcing_client , codistributor1d(1));

    % deterime which of "lab" or worker has data for location being logged
    % determine row on this worker of client's row for time_series_loc.

    global_rows_this_worker_stores = globalIndices(kmt_in, 1);
    row_time_series_loc = find( global_rows_this_worker_stores == time_series_loc_client);
    if (row_time_series_loc>0)
        time_series_lab = spmdIndex;
    end
    % end
    % if (n == 1)
    %     whos States_client Tracers_client kmt_client
    %     whos kmt_in F_in SF_in
    %     disp('bytes sent by scatter of kmt_in, F_in, SF_in: '); tocBytes(gcp)
    %     toc
    % end
    %
    % % input every time step: State and current Tracers
    % % Using 'distributed" is odd because a) takes as long as codistributed even
    % % if it sends way less (1/numWorks) data and b) does not work because it
    % % distributes over last index instead of first.
    % % if (n == 1); ticBytes(gcp); end;
    % % S_in2 = distributed(States_client );
    % % T_in2 = distributed(Tracers_client);
    % % if (n == 1)
    % %     whos S_in2 T_in2
    % %     disp('bytes sent by scatter of S_in2 and T_in2 : '); tocBytes(gcp)
    % %     toc
    % % end
    % if (n == 1); ticBytes(gcp); end
    % spmd(num_threads)
    S_in   = codistributed(States_client       , codistributor1d(1));
    T_in   = codistributed(Tracers_client      , codistributor1d(1));
    % end
    % if (n == 1)
    %     whos S_in T_in
    %     disp('bytes sent by scatter of S_in and T_in: '); tocBytes(gcp)
    %     toc
    % end
    %
    % % output every time step: State and Tendency
    % if (n == 1); ticBytes(gcp); end
    % spmd(num_threads)
    S_out_globalSize = S_out_sz ;
    S_out_codistr    = codistributor1d ( 1, codistributor1d.unsetPartition, S_out_globalSize );

    tend_globalSize = tend_sz ;
    tend_codistr    = codistributor1d ( 1, codistributor1d.unsetPartition, tend_globalSize );

    % Currently NOT using D, do NOT waste time: mex read, D_out writes.
    % diag_globalSize = diag_sz ;
    % diag_codistr    = codistributor1d ( 1, codistributor1d.unsetPartition, diag_globalSize );
    % surf_globalSize = surf_sz ;
    % surf_codistr    = codistributor1d ( 1, codistributor1d.unsetPartition, surf_globalSize );
    % % % end
    % % % spmd(num_threads)
    % inputs need to grab local part so that we do NOT get codistributed
    % which can NOT be processed (sensibly) as part of a struct like
    % "interior" and "surface"
    %
    % e.g 7881 water columns, and e.g. 6 workers, then each worker has 1314

    T_local   = getLocalPart( T_in   ); % size(T_local)     % 1314x60x32
    kmt_local = getLocalPart( kmt_in );                     % 1314x1
    F_local   = getLocalPart( F_in   ); % size(F_local)     % 1314x60x6
    SF_local  = getLocalPart( SF_in  ); % size(SF_local)    % 1314x1x11
    S_local   = getLocalPart( S_in   ); % size(S_local)     % 1314x60x2

    num_rows_for_this_worker = size(kmt_local,1);

    interior = interior_base_client;
    surface  = surface_base_client;

    % build local storage on each thread for outputs. Do NOT write to
    % client in for loop!

    tend_localSize = tend_globalSize;
    tend_localSize ( tend_codistr.Dimension ) = tend_codistr.Partition ( spmdIndex );
    tend_L = zeros ( tend_localSize );          % What we write to in loop

    S_out_localSize = S_out_globalSize;
    S_out_localSize ( S_out_codistr.Dimension ) = S_out_codistr.Partition ( spmdIndex );
    S_out_L = zeros ( S_out_localSize );          % What we write to in loop

    % FIXME: currently we do NOT use these, probably should NOT waste time
    % FIXME: collecting them in mex read and/or D_out and SD_out writes
    %     diag_localSize = diag_globalSize;
    %     diag_localSize ( diag_codistr.Dimension ) = diag_codistr.Partition ( labindex );
    %     diag_L = zeros ( diag_localSize );          % What we write to in loop

    %     surf_localSize = surf_globalSize;
    %     surf_localSize ( surf_codistr.Dimension ) = surf_codistr.Partition ( labindex );
    %     surf_L = zeros ( surf_localSize );          % What we write to in loop
    % end
    % if (n == 1)
    %     whos S_out_globalSize S_out_codistr tend_globalSize tend_codistr
    %     disp('bytes sent by creation of S_out and T_end: '); tocBytes(gcp)
    %     toc
    % end
    %
    %
    % % Finally! Lets calculate something in parallel!
    %
    % spmd(num_threads)

    for row = 1: num_rows_for_this_worker

        % disp(['water column #', num2str(row),' runs on thread #', num2str(labindex)])

        % Restore depth, state, and tracer at this location in structs

        interior.domain.kmt = kmt_local      (row);         % size(kmt_local
        interior.state      = squeeze(S_local(row,:,:))';   % size(S_local(row,:,:))    % 1x60x2
        interior.tracer     = squeeze(T_local(row,:,:))';   % size(T_local(row,:,:))    % 1x60x32

        % update forcing of interior and surface in structs

        surface.forcing     = squeeze(SF_local      (row,:,:))';  % size(SF_local(row,:,:))   % 1x1x11
        %         surface.river_flux  = squeeze(bgc.river_flux(row,:,:))';
        interior.forcing    = squeeze(F_local       (row,:,:))';  % size(F_local(row,:,:))    % 1x60x6

        % Run actual MARBL calulations

        [surface, interior] = MARBL_loop_iteration ( dt_client, n, row, surface, interior );

        % record result at all levels for time step

        tend_L(row,:,:) = interior.tendency';
        S_out_L(row,:,:) = interior.state';
        % FIXME: currently we do NOT use these, probably should NOT waste time
        % FIXME: collecting them in mex read and/or D_out and SD_out writes
        %         diag_L(row,:,:) = mex_marbl_driver ( 'interior_tendency_diags')'; % large and slow to store
        % surf_L(row,:,:) = mex_marbl_driver ( 'surface_flux_diags'             )'; % neither large or slow
        % bgc.sfo       (row,:,:) = mex_marbl_driver ( 'sfo'                    )'; % neither large or slow

        % log time series at just 1 location. Gets very large, very quickly

        global_row = global_rows_this_worker_stores(row);
        if ( (global_row == time_series_loc_client)  && sim.logTracers)
            %             disp(['logging col: ', num2str(row)])
            if (logDiags_client)
                surface.diag    = mex_marbl_driver ( 'surface_flux_diags' );
                interior.diag   = mex_marbl_driver ( 'interior_tendency_diags');
            end
            surface.sfo     = mex_marbl_driver ( 'sfo' );
            %             disp(['PAR forcing and diag: (', num2str(interior.forcing(2,1), '%1.1f'), ', ',num2str(interior.diag(151,1), '%1.1f'),') (W/m^2)']);
            surface_time_series_loc  = surface ;
            interior_time_series_loc = interior;
        end

    end % scan over all water columns in grid
end
% spmd
%     % efficient gather(), done in parallel
%     % t_spmd_inner = tic;
%     T_out  = codistributed.build( tend_L, tend_codistr,'noCommunication');
%     S_out  = codistributed.build( S_out_L, S_out_codistr,'noCommunication');
%
%     % FIXME: currently we do NOT use these, probably should NOT waste time
%     % FIXME: collecting them in mex read and/or D_out and SD_out writes
%     %     D_out  = codistributed.build( diag_L, diag_codistr);
%     %     SD_out = codistributed.build( surf_L, surf_codistr);
%
%     %  disp(['worker: MARBL build: ', num2str(toc(t_spmd_inner)*1000, '%1.1f'),' (ms) for ', num2str(num_rows_for_this_worker),' tendency'])
%
% end % spmd


if (sim.logTracers)
    % FIXME: deal with VERY frustrating and undesireable creation of composite
    % Note: curly braces of time_series_lab{} -not- parentheses ()
    time_series_lab = max                      ([time_series_lab{:}]);
    surface         = surface_time_series_loc  { time_series_lab };
    interior        = interior_time_series_loc { time_series_lab };
    update_log ();
end

% if (n == 1); ticBytes(gcp); end
% tendency must be gathered because mfactor() does NOT work on distributed
% bgc.tendency  = gather ( T_out);
% bgc.state     = gather ( S_out);

if (num_threads == 4)
    bgc.state     = [   (S_out_L{1}); (S_out_L{2}); (S_out_L{3}); (S_out_L{4}) ];
    bgc.tendency  = [   (tend_L{1}); (tend_L{2}); (tend_L{3}); (tend_L{4}) ];
elseif (num_threads == 12)
    bgc.state     = [   (S_out_L{1}); (S_out_L{2}); (S_out_L{3}); (S_out_L{4}); ...
        (S_out_L{5}); (S_out_L{6}); ...
        (S_out_L{7}); (S_out_L{8}); (S_out_L{9}); ...
        (S_out_L{10});(S_out_L{11});(S_out_L{12}); ...
        ];
    bgc.tendency  = [   (tend_L{1}); (tend_L{2}); (tend_L{3}); (tend_L{4}); ...
        (tend_L{5}); (tend_L{6}); ...
        (tend_L{7}); (tend_L{8}); (tend_L{9}); ...
        (tend_L{10});(tend_L{11});(tend_L{12}); ...
        ];
else
    keyboard;
end

% FIXME: currently we do NOT use these, probably should NOT waste time
% FIXME: collecting them in mex read and/or D_out and SD_out writes
% bgc.diag      = gather ( D_out);
% bgc.surf_diag = gather (SD_out);

% if (n == 1)
%     whos T_out S_out
%     disp('bytes received by gather: '); tocBytes(gcp)
%     toc
% end

    function update_log
        if (sim.logTracers)
            time_series.tracer        (:, :, n) = interior.tracer';
            time_series.sfo           (:, n)    = surface.sfo;
        end
        if (sim.logDiags)
            time_series.diag      (:, :, n) = interior.diag';
            time_series.surf_diag (:, n)    = surface.diag';
        end
    end % update_log

end % MARBL_loop_parallel
