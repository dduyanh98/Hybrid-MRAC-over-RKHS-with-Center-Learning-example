%% ======================================================================
%  Animation and snapshot export.
%  1) Video: phase portrait, active library, and used centers
%     - Current-state circles match trajectory colors (real vs ref)
%     - Active library split into original vs refined/new centers
%  2) Paper figure: 2x3 snapshots at t = [0.04 0.10 0.50 3.50 6.35 9.79]
%     - Same color conventions as video
% ======================================================================

clearvars -except LOG Xi_hist center_times Xi_lib_hist active_idx_hist lib_times
close all; clc;

% ===================== GLOBAL LATEX STYLE =====================
set(groot, ...
    'defaultTextInterpreter','latex', ...
    'defaultLegendInterpreter','latex', ...
    'defaultAxesTickLabelInterpreter','latex', ...
    'defaultAxesFontSize',18);

font_size = 18;
font_size_title = 20;

%% ===================== REQUIRED VARIABLES =====================
assert(exist('LOG','var')==1 && isfield(LOG,'T') && isfield(LOG,'X'), 'LOG.T / LOG.X not found.');
assert(exist('Xi_hist','var')==1 && exist('center_times','var')==1, 'Xi_hist / center_times not found.');
assert(exist('Xi_lib_hist','var')==1 && exist('active_idx_hist','var')==1 && exist('lib_times','var')==1, ...
    'Xi_lib_hist / active_idx_hist / lib_times not found (add snapshot logging).');
assert(size(LOG.X,2) >= 4, ...
    'LOG.X must have at least 4 columns: [theta, theta_dot, theta_ref, theta_dot_ref].');

%% ===================== EXTRACT TRAJECTORIES =====================
t_all     = LOG.T(:);
theta_all = LOG.X(:,1);
thd_all   = LOG.X(:,2);

theta_ref_all = LOG.X(:,3);
thd_ref_all   = LOG.X(:,4);

th_deg      = rad2deg(theta_all);
thd_deg     = rad2deg(thd_all);
th_ref_deg  = rad2deg(theta_ref_all);
thd_ref_deg = rad2deg(thd_ref_all);

U = numel(Xi_hist);               % number of "used-centers" updates
S = numel(Xi_lib_hist);           % number of library snapshots

% used centers in degrees
Xi_used_deg = cell(U,1);
for u = 1:U
    Xi_used_deg{u} = rad2deg(Xi_hist{u});
end

%% ===================== COLORS + ORIGINAL/REFINED SPLIT =====================
cReal = [0 0 0];          % real trajectory + marker
cRef  = [1 0 0];          % reference trajectory + marker

% Use the first snapshot to identify original library centers.
if S >= 1 && ~isempty(Xi_lib_hist{1})
    N0 = size(Xi_lib_hist{1}, 1);   % original centers count
else
    N0 = 0;
end

% Library colors: original vs refined/new
cLibOld = [0 0.6 0];   % original active library centers
cLibNew = [0.95 0.50 0.10];   % refined/new active library centers

%% ===================== AXIS BOUNDS (GLOBAL) =====================
xmin = min([th_deg;  th_ref_deg]);
xmax = max([th_deg;  th_ref_deg]);
ymin = min([thd_deg; thd_ref_deg]);
ymax = max([thd_deg; thd_ref_deg]);

for s = 1:S
    Xi_s = Xi_lib_hist{s};
    idx  = active_idx_hist{s};
    if isempty(Xi_s) || isempty(idx), continue; end
    idx  = idx(idx>=1 & idx<=size(Xi_s,1));
    if isempty(idx), continue; end
    Xact = rad2deg(Xi_s(idx,:));
    xmin = min(xmin, min(Xact(:,1)));
    xmax = max(xmax, max(Xact(:,1)));
    ymin = min(ymin, min(Xact(:,2)));
    ymax = max(ymax, max(Xact(:,2)));
end

pad = 0.05;
dx = xmax-xmin; if dx==0, dx=1; end
dy = ymax-ymin; if dy==0, dy=1; end

XL = [xmin-pad*dx, xmax+pad*dx];
YL = [ymin-pad*dy, ymax+pad*dy];

%% ======================================================================
%  PART 1: VIDEO ANIMATION
% ======================================================================

%% --- Video settings
outFile = 'pendulum_library_used_centers_dynamic.mp4';
fps     = 25;
stride  = 3;
doSave  = true;

%% --- Figure setup
fig = figure('Color','w','Renderer','opengl');
ax  = axes(fig);
hold(ax,'on'); grid(ax,'on');

xlabel(ax,'$\theta$ [deg]','Interpreter','latex');
ylabel(ax,'$\dot{\theta}$ [deg/s]','Interpreter','latex');
set(ax,'TickLabelInterpreter','latex');

xlim(ax, XL);
ylim(ax, YL);
axis(ax,'normal');

try
    fig.WindowState = 'maximized';
catch
    set(fig,'Units','pixels','Position',get(0,'ScreenSize'));
end

%% --- Plot handles
hTraj = plot(ax, nan, nan, 'LineWidth', 1.5, 'Color', cReal);
hNow  = scatter(ax, nan, nan, 70, 'filled', ...
    'MarkerFaceColor', cReal, 'MarkerEdgeColor', cReal);

hRefTraj = plot(ax, nan, nan, '--', 'LineWidth', 1.5, 'Color', cRef);
hRefNow  = scatter(ax, nan, nan, 70, 'filled', ...
    'MarkerFaceColor', cRef, 'MarkerEdgeColor', cRef);

% Active library split into original vs refined/new
hLibOld = scatter(ax, nan, nan, 14, 'filled', ...
    'MarkerFaceColor', cLibOld, 'MarkerEdgeColor', cLibOld);
hLibNew = scatter(ax, nan, nan, 14, 'filled', ...
    'MarkerFaceColor', cLibNew, 'MarkerEdgeColor', cLibNew);

% Used centers + box
hUsed = scatter(ax, nan, nan, 120, 'filled');
hPoly = plot(ax, nan, nan, '--', 'LineWidth', 1.4);

legend(ax, {'Trajectory','Current state','Reference traj','Reference state', ...
            'Active library (orig)','Active library (refined)', ...
            'Used centers','Used box'}, ...
       'Location','best');

uistack(hRefTraj,'bottom');

%% --- Video writer
if doSave
    v = VideoWriter(outFile, 'MPEG-4');
    v.FrameRate = fps;
    open(v);
end

%% --- Animation loop
ct = center_times(:);
lt = lib_times(:);

u = 0;  % used-centers index
s = 0;  % library snapshot index

for i = 1:stride:numel(t_all)
    t = t_all(i);

    % advance used-centers index
    while (u < numel(ct)) && (t >= ct(u+1))
        u = u + 1;
    end

    % advance library snapshot index
    while (s < numel(lt)) && (t >= lt(s+1))
        s = s + 1;

        Xi_s = Xi_lib_hist{s};
        idx  = active_idx_hist{s};

        if isempty(Xi_s) || isempty(idx)
            set(hLibOld,'XData',nan,'YData',nan);
            set(hLibNew,'XData',nan,'YData',nan);
        else
            idx = idx(idx>=1 & idx<=size(Xi_s,1));

            if isempty(idx)
                set(hLibOld,'XData',nan,'YData',nan);
                set(hLibNew,'XData',nan,'YData',nan);
            else
                idx_old = idx(idx <= N0);
                idx_new = idx(idx >  N0);

                if isempty(idx_old)
                    set(hLibOld,'XData',nan,'YData',nan);
                else
                    Xold = rad2deg(Xi_s(idx_old,:));
                    set(hLibOld,'XData',Xold(:,1),'YData',Xold(:,2));
                end

                if isempty(idx_new)
                    set(hLibNew,'XData',nan,'YData',nan);
                else
                    Xnew = rad2deg(Xi_s(idx_new,:));
                    set(hLibNew,'XData',Xnew(:,1),'YData',Xnew(:,2));
                end
            end
        end
    end

    % update trajectory + current state
    set(hTraj, 'XData', th_deg(1:i), 'YData', thd_deg(1:i));
    set(hNow,  'XData', th_deg(i),   'YData', thd_deg(i));

    % update reference trajectory + current state
    set(hRefTraj, 'XData', th_ref_deg(1:i), 'YData', thd_ref_deg(1:i));
    set(hRefNow,  'XData', th_ref_deg(i),   'YData', thd_ref_deg(i));

    % update used centers + axis-aligned rectangle (horizontal/vertical)
    if u >= 1
        XU = Xi_used_deg{u};
        set(hUsed,'XData',XU(:,1),'YData',XU(:,2));

        xmin_u = min(XU(:,1)); xmax_u = max(XU(:,1));
        ymin_u = min(XU(:,2)); ymax_u = max(XU(:,2));

        xbox = [xmin_u xmax_u xmax_u xmin_u xmin_u];
        ybox = [ymin_u ymin_u ymax_u ymax_u ymin_u];
        set(hPoly,'XData',xbox,'YData',ybox);
    else
        set(hUsed,'XData',nan,'YData',nan);
        set(hPoly,'XData',nan,'YData',nan);
    end

    title(ax, sprintf('Pendulum  |  t = %.2f s', t), 'Interpreter','none');

    drawnow limitrate;

    if doSave
        writeVideo(v, getframe(fig));
    end
end

if doSave
    close(v);
    fprintf('Saved video: %s\n', outFile);
end

%% ======================================================================
%  PAPER FIGURE: 2x3 snapshots
% ======================================================================

t_snap = [0.04, 0.10, 0.50, 3.50, 6.35, 9.79];

ct = center_times(:);
lt = lib_times(:);

fig2 = figure('Color','w','Renderer','opengl');
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

for k = 1:numel(t_snap)

    t = t_snap(k);

    i = find(t_all <= t, 1, 'last');
    if isempty(i), i = 1; end

    u = find(ct <= t, 1, 'last');
    if isempty(u), u = 0; end

    s = find(lt <= t, 1, 'last');
    if isempty(s), s = 0; end

    nexttile; hold on; grid on;
    xlim(XL); ylim(YL);

    xlabel('$\theta$ [deg]','Interpreter','latex');
    ylabel('$\dot{\theta}$ [deg/s]','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    %% ===== Trajectories =====
    hTraj_k    = plot(th_deg(1:i),     thd_deg(1:i),     ...
        'LineWidth',1.5,'Color',cReal);

    hRefTraj_k = plot(th_ref_deg(1:i), thd_ref_deg(1:i), ...
        '--','LineWidth',1.5,'Color',cRef);

    %% ===== Current States =====
    hNow_k = scatter(th_deg(i), thd_deg(i), 60, 'filled', ...
        'MarkerFaceColor',cReal,'MarkerEdgeColor',cReal);

    hRefNow_k = scatter(th_ref_deg(i), thd_ref_deg(i), 60, 'filled', ...
        'MarkerFaceColor',cRef,'MarkerEdgeColor',cRef);

    %% ===== Active Library (orig = green, refined = orange) =====
    % Initialize legend handles for both library groups.
    hLibOld_k = scatter(nan,nan,12,'filled', ...
        'MarkerFaceColor',cLibOld,'MarkerEdgeColor',cLibOld);

    hLibNew_k = scatter(nan,nan,12,'filled', ...
        'MarkerFaceColor',cLibNew,'MarkerEdgeColor',cLibNew);

    if s >= 1
        Xi_s = Xi_lib_hist{s};
        idx  = active_idx_hist{s};

        if ~isempty(Xi_s) && ~isempty(idx)
            idx = idx(idx>=1 & idx<=size(Xi_s,1));

            idx_old = idx(idx <= N0);
            idx_new = idx(idx >  N0);

            if ~isempty(idx_old)
                Xold = rad2deg(Xi_s(idx_old,:));
                set(hLibOld_k,'XData',Xold(:,1),'YData',Xold(:,2));
            end

            if ~isempty(idx_new)
                Xnew = rad2deg(Xi_s(idx_new,:));
                set(hLibNew_k,'XData',Xnew(:,1),'YData',Xnew(:,2));
            end
        end
    end

    %% ===== Used Centers + Used Box =====
    if u >= 1
        XU = Xi_used_deg{u};

        hUsed_k = scatter(XU(:,1), XU(:,2), 90, 'filled');

        xmin_u = min(XU(:,1)); xmax_u = max(XU(:,1));
        ymin_u = min(XU(:,2)); ymax_u = max(XU(:,2));

        xbox = [xmin_u xmax_u xmax_u xmin_u xmin_u];
        ybox = [ymin_u ymin_u ymax_u ymax_u ymin_u];

        hBox_k = plot(xbox,ybox,'--','LineWidth',1.4);
    else
        hUsed_k = scatter(nan,nan,90,'filled');
        hBox_k  = plot(nan,nan,'--','LineWidth',1.4);
    end

    title(sprintf('Pendulum  |  t = %.2f s', t),'Interpreter','none');

    %% ===== LEGEND (FIRST TILE ONLY) =====
    if k == 1
        legend([hTraj_k, hNow_k, hRefTraj_k, hRefNow_k, ...
                hLibOld_k, hLibNew_k, hUsed_k, hBox_k], ...
            {'Trajectory','Current state','Reference traj','Reference state', ...
             'Active library (orig)','Active library (refined)', ...
             'Used centers','Used box'}, ...
            'Location','best');
    end
end

set(fig2,'PaperPositionMode','auto');
print(fig2,'pendulum_rkhs_snapshots_2x3.png','-dpng','-r300');
print(fig2,'pendulum_rkhs_snapshots_2x3.pdf','-dpdf','-r300');

fprintf('Saved:\n  pendulum_rkhs_snapshots_2x3.png\n  pendulum_rkhs_snapshots_2x3.pdf\n');