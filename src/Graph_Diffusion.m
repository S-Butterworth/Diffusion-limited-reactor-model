%% Importing inital conditions and constants 
[C,~] = Initial_Conditions_and_Constants();

%% File directory for figures and .mat files
iptDir = fileparts(mfilename('fullpath'));
outDirfig = fullfile(iptDir, 'figures');
if ~isfolder(outDirfig)
    mkdir(outDirfig);
end

iptDir = fileparts(mfilename('fullpath'));
outDirdata = fullfile(iptDir, 'data');
if ~isfolder(outDirdata)
    mkdir(outDirdata);
end

%% Funciton to produce the full radius and concentration values when there is a dead core region
function [r_full, CA_full, CB_full] = physicalProfile(sol, rc)
    if ~isnan(rc)
        r_active = rc + (1-rc)*sol.x;
        r_dead = linspace(0, rc, 50);
        CA_full = [zeros(1,50), sol.y(1,:)];
        CB_full = [sol.y(3,1)*ones(1,50), sol.y(3,:)];
        r_full = [r_dead, r_active];
    else
        r_full = sol.x;
        CA_full = sol.y(1,:);
        CB_full = sol.y(3,:);
    end
end

%% Single Temperature diffusion
resp = '';
fprintf('\n--- Single Temperature Diffusion: normalised concentraions and pellet radius at fixed temperature ---\n');
while ~any(strcmpi(resp, {'R','S'}))
    fprintf('\n')
    resp = input('Run (R) or Skip (S)?: ', 's');
    if resp == 'R'
        temp = NaN;
        while isnan(temp)
            temp = str2double(input('Enter the temperature to evaluate at: ', 's'));
            C.t = str2double(input('Enter the time on stream to evaluate at: ', 's'));
        end
    end
    fprintf('\n')
end

switch upper(resp)
    case 'R'
        T0_steps = linspace(0, temp, 200);
        [~, ~, rc,sol] = Diffusion(T0_steps, C);
        [r_full, CA_full, CB_full] = physicalProfile(sol, rc(end));
    case 'S'
        disp('Skipped.'); 
end
if resp == 'R'
    f = figure('Visible','on');
    theme(f,"light");
    ax = axes(f);

    p = plot(ax, r_full, CA_full, r_full, CB_full);
    if ~isnan(rc(end))
        xline(ax, rc(end), '--', 'r_c', ...
            'LabelVerticalAlignment','middle', 'LabelHorizontalAlignment','center');
    end

    legend(ax, p, {'A','B'});
    xlabel(ax,'Normalised Radius'); ylabel(ax,'Normalised Concentration')
    title(ax,sprintf('Normalised Concentration Diffusion in Catalyst at: %g °C after %g days on stream', temp,C.t))

    grid(ax,'on'); grid(ax,'minor')
    ax.MinorGridLineStyle = '-';
    ax.MinorGridAlpha = 0.05;
    ax.XDir = 'reverse';

    fname = sprintf('Diffusion_T_%gC_Day_%g.png', T0_steps(end), C.t);
    exportgraphics(f, fullfile(outDirfig, fname), ...
        'Resolution',600, 'BackgroundColor','white', 'Padding',60);
end

%% Temperature Sweep at input mole fractions
T0_steps  = linspace(0, 1200, 200);
sweepFile = fullfile(outDirdata, 'temperature_sweep.mat');

if isfile(sweepFile)
    resp = '';
    fprintf('\n--- Temperature Sweep: eta and r_c vs T at fixed input composition ---\n');
    showTableInfo(sweepFile)
    while ~any(strcmpi(resp, {'L', 'R', 'S'}))
        resp = input('Load (L), Rerun (R) or Skip (S)?: ', 's');
        if isempty(resp), resp = 'S'; end
    end
else
    resp = 'R';
end

switch upper(resp)
    case 'L'
        S = load(sweepFile);
        T0_steps    = S.T0_steps;
        eta_history = S.eta_history;
        rc_history  = S.rc_history;
        tableInfo   = S.tableInfo;
        fprintf('Loaded temperature sweep from %s (generated %s)\n', sweepFile, tableInfo.dateGenerated);
    case 'R'
        [T0_steps, eta_history, rc_history, ~] = Diffusion(T0_steps, C);

        tableInfo = struct();
        tableInfo.dateGenerated = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');
        tableInfo.N_T = length(T0_steps);
        tableInfo.T_range = [min(T0_steps), max(T0_steps)];
        tableInfo.t_use = C.t;
        tableInfo.k = C.k; tableInfo.EA = C.EA;
        tableInfo.kd = C.kd; tableInfo.n_deactivation = C.n;
        tableInfo.dp = C.dp; tableInfo.pc = C.pc; tableInfo.D = C.D;
        tableInfo.PT = C.PT; tableInfo.x = C.x;
        tableInfo.notes = 'Single-composition temperature sweep at fixed input mole fractions.';

        save(sweepFile, 'T0_steps', 'eta_history', 'rc_history', 'tableInfo');
        fprintf('Sweep complete, saved to %s\n', sweepFile);
    case 'S'
        disp('Skipped.');     
end

fprintf('\n')

if resp == 'L' || resp == 'R'
    f = figure('Visible','on');
    theme(f,"light");
    tl = tiledlayout(f,2,1,'TileSpacing','compact');

    data   = {eta_history, rc_history};
    ylabs  = {'\eta', 'r_c'};
    titles = {'Effectiveness factor vs temperature', 'Dead-core radius vs temperature'};

    for k = 1:2
        ax = nexttile(tl);
        plot(ax, T0_steps, data{k})
        xlabel(ax,'T_0 (°C)'); ylabel(ax, ylabs{k}); title(ax, titles{k})
        grid(ax,'on'); grid(ax,'minor')
        ax.MinorGridLineStyle = '-';
        ax.MinorGridAlpha = 0.05;
        xlim(ax, [T0_steps(1) T0_steps(end)])
    end

    exportgraphics(f, fullfile(outDirfig,"Effectiveness_factor_and_dead_core_radius_vs_temperature.png"), ...
        'Resolution',600,'BackgroundColor','white','Padding',60);
end

%% Changing conversion values, independent for stochiometic feeds where CAs/CBs ratio is constant
pellet_eta = fullfile(outDirdata, 'pellet_eta_table.mat');

if isfile(pellet_eta)
    resp = '';
    fprintf('\n--- Composition Table: eta and r_c vs (T, X) along the reactor conversion trajectory ---\n');
    showTableInfo(pellet_eta)
    while ~any(strcmpi(resp, {'L', 'R', 'S'}))
        resp = input('Load (L), Rerun (R) or Skip (S)?: ', 's');
        if isempty(resp), resp = 'S'; end
    end
else
    resp = 'R';
end

switch upper(resp)
    case 'L'
        S = load(pellet_eta);
        eta_table = S.eta_table;
        rc_table  = S.rc_table;
        T_grid    = S.T_grid;
        X_grid    = S.X_grid;
        xA_traj   = S.xA_traj;
        xB_traj   = S.xB_traj;
        tableInfo = S.tableInfo;

        fprintf('Loaded eta table %s (generated %s)\n', pellet_eta, tableInfo.dateGenerated);
    case 'R'
        tic;
        N_T = 500;
        N_X = 200;

        pool = gcp('nocreate');
        if isempty(pool)
            pool = parpool;
        end

        X_grid = linspace(0, 0.95, N_X);
        T_grid = linspace(0,1200,N_T);

        FA = C.F(1)*(1-X_grid);
        FB = C.F(2) - 0.5*C.F(1)*X_grid;
        FC = C.F(3) + C.F(1)*X_grid;
        Finert = C.F(4)*ones(size(X_grid));
        Ftotal = FA + FB + FC + Finert;

        xA_traj = FA./Ftotal;
        xB_traj = FB./Ftotal;

        eta_flat = nan(N_T, N_X);
        rc_flat  = nan(N_T, N_X);

        futures = parallel.FevalFuture.empty;
        for idx = 1:N_X
            Clocal = C;
            Clocal.x(1) = xA_traj(idx);
            Clocal.x(2) = xB_traj(idx);
            futures(idx) = parfeval(pool, @Diffusion, 4, T_grid, Clocal);
        end
        cleanup = onCleanup(@() cancel(futures));

        barWidth = 30;
        msgLen = 0;
        tStart = tic;
        for n = 1:N_X
            [idx, ~, eta_hist, rc_hist, ~] = fetchNext(futures);
            eta_flat(:,idx) = eta_hist(:);
            rc_flat(:,idx)  = rc_hist(:);

            frac    = n / N_X;
            filled  = round(frac * barWidth);
            elapsed = toc(tStart);
            eta     = elapsed / frac - elapsed;
            msg = sprintf('[%s%s] %3.0f%%  %d/%d  ETA %s', ...
                repmat('=', 1, filled), repmat(' ', 1, barWidth - filled), ...
                100*frac, n, N_X, formatTime(eta));
            out = [msg repmat(' ', 1, max(0, msgLen - length(msg)))];
            fprintf([repmat('\b', 1, msgLen) '%s'], out);
            msgLen = length(out);
        end
        fprintf('\n');

        eta_table = eta_flat;
        rc_table  = rc_flat;
        nanMask = isnan(eta_table);

        fprintf('%d of %d table entries failed to converge (%.1f%%)\n', ...
                sum(nanMask(:)), numel(eta_table), 100*sum(nanMask(:))/numel(eta_table));
        outFile = 'pellet_eta_table.mat';

        tableInfo = struct();
        tableInfo.dateGenerated = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');
        tableInfo.elapsedMinutes = toc/60;

        % Grid definition
        tableInfo.N_T = N_T; tableInfo.N_X = N_X;
        tableInfo.T_range = [min(T_grid), max(T_grid)]; tableInfo.X_range = [min(X_grid), max(X_grid)];
        tableInfo.t_use = C.t;

        % Fixed physical/kinetic parameters used to build this table
        tableInfo.k = C.k;   tableInfo.EA = C.EA; tableInfo.kd = C.kd; tableInfo.n_deactivation = C.n;
        tableInfo.dp = C.dp; tableInfo.pc = C.pc; tableInfo.D = C.D; tableInfo.PT = C.PT; tableInfo.x = C.x;

        % Solver settings
        tableInfo.threshold = C.threshold;   tableInfo.minRcToSwitch = C.minRcToSwitch;
        tableInfo.RelTol_stage1 = C.RelTol1; tableInfo.RelTol_stage2 = C.RelTol2;

        % Convergence summary
        tableInfo.nanCount = sum(nanMask(:));
        tableInfo.totalPoints = numel(eta_table);
        tableInfo.nanPercent = 100*tableInfo.nanCount/tableInfo.totalPoints;

        % Assumptions
        tableInfo.notes = ['Composition varies via conversion X, tied to stoichiometric feed ratios. ' ...
                            'If FA0=2*FB0 exactly, eta is independent of X (see derivation notes) — ' ...
                            'this table remains useful for non-stoichiometric feed variants.'];

        save(pellet_eta, 'eta_table', 'rc_table', 'T_grid', 'X_grid', 'xA_traj', 'xB_traj', 'tableInfo');
        fprintf('Eta table saved as %s\n', outFile);
    case 'S'
        disp('Skipped.');
end

if resp == 'L' || resp == 'R'
    f = figure('Visible','on');
    theme(f,"light");
    tl = tiledlayout(f,2,1,'TileSpacing','compact','Padding','compact');

    Z      = {eta_table, rc_table};
    zlabs  = {'\eta', 'r_c'};
    titles = {'Effectiveness factor vs conversion and temperature', ...
            'Dead-core radius vs conversion and temperature'};

    for k = 1:2
        ax = nexttile(tl);
        surf(ax, X_grid, T_grid, Z{k}, ...
            'FaceColor','interp', 'EdgeColor','k', 'EdgeAlpha',0.15);

        title(ax, titles{k})
        xlabel(ax,'Conversion X'); ylabel(ax,'T_0 (°C)'); zlabel(ax, zlabs{k})

        grid(ax,'on'); grid(ax,'minor')
        ax.MinorGridLineStyle = '-';
        ax.MinorGridAlpha = 0.05;
        ax.Box = 'on';
        ax.FontSize = 10;
        ax.LineWidth = 0.8;
        view(ax, -35, 30)
        colormap(ax, parula)
        colorbar(ax);
    end
    exportgraphics(f, fullfile(outDirfig,"Effectiveness_factor_and_dead_core_radius_vs_temperature_and_conversion.png"), ...
        'Resolution',600,'BackgroundColor','white','Padding',60);
end

%% Deactivation impact on diffusion
pellet_eta_act = fullfile(outDirdata, 'pellet_eta_table_activity.mat');

if isfile(pellet_eta_act)
    resp = '';
    fprintf('\n--- Activity/Time-in-Use Table: eta and r_c vs (T, t) at fixed catalyst age ---\n');
    showTableInfo(pellet_eta_act)
    while ~any(strcmpi(resp, {'L', 'R', 'S'}))
        resp = input('Load (L), Rerun (R) or Skip (S)?: ', 's');
        if isempty(resp), resp = 'S'; end
    end
else
    resp = 'R';
end

switch upper(resp)
    case 'L'
        S = load(pellet_eta_act);
        eta_table = S.eta_table;
        rc_table  = S.rc_table;
        T_grid    = S.T_grid;
        t_grid    = S.t_grid;
        tableInfo = S.tableInfo;

        fprintf('Loaded eta table %s (generated %s)\n', pellet_eta_act, tableInfo.dateGenerated);
    case 'R'
        tic;
        pool = gcp('nocreate');
        if isempty(pool)
            pool = parpool;
        end

        N_T = 500;
        N_t = 200;

        t_grid = linspace(0, 100, N_t);
        T_grid = linspace(0, 1200, N_T);

        eta_flat = nan(N_T, N_t);
        rc_flat  = nan(N_T, N_t);

        futures = parallel.FevalFuture.empty;
        for idx_t = 1:N_t
            Clocal = C;
            Clocal.t  = t_grid(idx_t);
            futures(idx_t) = parfeval(pool, @Diffusion, 4, T_grid, Clocal);
        end
        cleanup = onCleanup(@() cancel(futures));

        barWidth = 30;
        msgLen = 0;
        tStart = tic;
        for n = 1:N_t
            [idx_t, ~, eta_hist, rc_hist, ~] = fetchNext(futures);
            eta_flat(:,idx_t) = eta_hist(:);
            rc_flat(:,idx_t)  = rc_hist(:);

            frac    = n / N_t;
            filled  = round(frac * barWidth);
            elapsed = toc(tStart);
            eta     = elapsed / frac - elapsed;
            msg = sprintf('[%s%s] %3.0f%%  %d/%d  ETA %s', ...
                repmat('=', 1, filled), repmat(' ', 1, barWidth - filled), ...
                100*frac, n, N_t, formatTime(eta));

            out = [msg repmat(' ', 1, max(0, msgLen - length(msg)))];
            fprintf([repmat('\b', 1, msgLen) '%s'], out);
            msgLen = length(out);
        end
        fprintf('\n');

        eta_table = eta_flat;
        rc_table  = rc_flat;

        nanMask = isnan(eta_table);
        fprintf('%d of %d table entries failed to converge (%.1f%%)\n', ...
                sum(nanMask(:)), numel(eta_table), 100*sum(nanMask(:))/numel(eta_table));
                
        tableInfo = struct();
        tableInfo.dateGenerated = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');
        tableInfo.elapsedMinutes = toc/60;

        % Grid definition
        tableInfo.N_T = N_T; tableInfo.N_t = N_t;
        tableInfo.T_range = [min(T_grid), max(T_grid)]; tableInfo.t_range = [min(t_grid), max(t_grid)];

        % Fixed physical/kinetic parameters used to build this table
        tableInfo.k = C.k;   tableInfo.EA = C.EA; tableInfo.kd = C.kd; tableInfo.n_deactivation = C.n;
        tableInfo.a0 = C.a0; tableInfo.dp = C.dp; tableInfo.pc = C.pc; tableInfo.D = C.D;
        tableInfo.PT = C.PT; tableInfo.x = C.x;

        % Solver settings
        tableInfo.threshold = C.threshold;   tableInfo.minRcToSwitch = C.minRcToSwitch;
        tableInfo.RelTol_stage1 = C.RelTol1; tableInfo.RelTol_stage2 = C.RelTol2;

        % Convergence summary
        tableInfo.nanCount = sum(nanMask(:));
        tableInfo.totalPoints = numel(eta_table);
        tableInfo.nanPercent = 100*tableInfo.nanCount/tableInfo.totalPoints;
        
        % Assumptions
        tableInfo.notes = ['2D (T,t) table. X dimension omitted: feed is exactly stoichiometric ' ...
                    '(FA0=2*FB0), so CAs/CBs ratio is invariant with conversion and eta is ' ...
                    'independent of X for this feed. See derivation notes if feed changes.'];

        save(pellet_eta_act, 'eta_table', 'rc_table', 'T_grid', 't_grid', 'tableInfo');
        fprintf('Eta table saved as %s\n', pellet_eta_act);
    case 'S'
        disp('Skipped.');
end

if resp == 'L' || resp == 'R'
    f = figure('Visible','on');
    theme(f,"light");
    tl = tiledlayout(f,2,1,'TileSpacing','compact','Padding','compact');

    Z      = {eta_table, rc_table};
    zlabs  = {'\eta', 'r_c'};
    titles = {'Effectiveness factor vs time in use and temperature', ...
            'Dead-core radius vs time in use and temperature'};

    [t_grid, T_grid] = meshgrid(t_grid, T_grid); 

    for k = 1:2
        ax = nexttile(tl);
        surf(ax, t_grid, T_grid, Z{k},t_grid, ...
            'FaceColor','interp', 'EdgeColor','none');

        title(ax, titles{k})
        xlabel(ax,'Time in use (days)'); ylabel(ax,'T_0 (°C)'); zlabel(ax, zlabs{k})
        ylim(ax, [min(T_grid(:)) max(T_grid(:))])
        grid(ax,'on'); grid(ax,'minor')
        ax.MinorGridLineStyle = '-';
        ax.MinorGridAlpha = 0.05;
        ax.Box = 'on';
        ax.FontSize = 10;
        ax.LineWidth = 0.8;
        view(ax, 90, 0)
        colormap(ax, parula)
    end
    clim(ax, [min(t_grid(:)) max(t_grid(:))])
    cb = colorbar(ax);
    cb.Layout.Tile = 'east';
    cb.Label.String = 'Time on stream (days)';

    exportgraphics(f, fullfile(outDirfig,"Effectiveness_factor_and_dead_core_radius_vs_temperature_and_time_on_stream.png"), ...
        'Resolution',600,'BackgroundColor','white','Padding',60);
end

%% Time formatting
function s = formatTime(t)
    t = round(t);
    if t < 60
        s = sprintf('%ds', t);
    elseif t < 3600
        s = sprintf('%dm %02ds', floor(t/60), mod(t, 60));
    else
        s = sprintf('%dh %02dm', floor(t/3600), floor(mod(t, 3600)/60));
    end
end

%% Table information function
function showTableInfo(filename)
    info = load(filename, 'tableInfo');
    ti = info.tableInfo;

    fprintf('=== Table: %s ===\n', filename);

    if isfield(ti, 'dateGenerated')
        fprintf('Generated:        %s\n', ti.dateGenerated);
    end
    if isfield(ti, 'elapsedMinutes')
        fprintf('Build time:       %.1f min\n', ti.elapsedMinutes);
    end

    dims = {};
    if isfield(ti,'N_T'), dims{end+1} = sprintf('%d (T)', ti.N_T); end
    if isfield(ti,'N_X'), dims{end+1} = sprintf('%d (X)', ti.N_X); end
    if isfield(ti,'N_t'), dims{end+1} = sprintf('%d (t)', ti.N_t); end
    if ~isempty(dims)
        fprintf('Grid:             %s\n', strjoin(dims, ' x '));
    end
    if isfield(ti,'T_range')
        fprintf('T range:          %.0f to %.0f C\n', ti.T_range(1), ti.T_range(2));
    end
    if isfield(ti,'X_range')
        fprintf('X range:          %.3f to %.3f\n', ti.X_range(1), ti.X_range(2));
    end
    if isfield(ti,'t_range')
        fprintf('t range:          %.0f to %.0f days\n', ti.t_range(1), ti.t_range(2));
    end
    if isfield(ti,'t_use')
        fprintf('Fixed t (days):   %.1f\n', ti.t_use);
    end
    if isfield(ti,'nanPercent')
        fprintf('Convergence:      %.1f%% failed (%d of %d points)\n', ...
                 ti.nanPercent, ti.nanCount, ti.totalPoints);
    end

    paramFields = {'k','EA','kd','n_deactivation','dp','pc','PT'};
    printed = false;
    for i = 1:length(paramFields)
        if isfield(ti, paramFields{i})
            if ~printed
                fprintf('Parameters:\n');
                printed = true;
            end
            fprintf('  %-16s %g\n', paramFields{i}, ti.(paramFields{i}));
        end
    end

    if isfield(ti,'notes')
        fprintf('Notes: %s\n', ti.notes);
    end
    fprintf('\n')
end
