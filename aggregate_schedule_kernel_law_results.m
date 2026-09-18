function TBL = aggregate_schedule_kernel_law_results()
%AGGREGATE_SCHEDULE_KERNEL_LAW_RESULTS Rank saved cases and plot best runs.

addpath(genpath(pwd));

in_dir = fullfile('Results', 'schedule_kernel_law_cases');
files = dir(fullfile(in_dir, '*.mat'));
if isempty(files)
    error('No case result files found in %s', in_dir);
end

results = repmat(empty_result(), numel(files), 1);
for i = 1:numel(files)
    S = load(fullfile(files(i).folder, files(i).name), 'result');
    r = S.result;
    fields = fieldnames(results(i));
    for f = 1:numel(fields)
        field = fields{f};
        if isfield(r, field)
            results(i).(field) = r.(field);
        end
    end
end

TBL = struct2table(rmfield(results, {'T', 'e'}));
TBL = sortrows(TBL, {'schedule', 'l2_5_10'}, {'ascend', 'ascend'});

best_runs = struct();
schedules = {'aggressive', 'current'};
for s = 1:numel(schedules)
    schedule = schedules{s};
    valid = find(strcmp({results.schedule}, schedule) & ...
                 cellfun(@isempty, {results.error}) & ...
                 ~isnan([results.l2_5_10]));
    if isempty(valid)
        warning('No valid completed cases for schedule %s', schedule);
        continue;
    end
    [~, best_local] = min([results(valid).l2_5_10]);
    best_runs.(schedule) = results(valid(best_local));
end

if isfield(best_runs, 'aggressive') && isfield(best_runs, 'current')
    plot_best_error_comparison(best_runs);
end

save(fullfile('Results', 'schedule_kernel_laws_compare.mat'), ...
    'results', 'TBL', 'best_runs');

disp(TBL(:, {'schedule','case_name','kernel_mode','kernel_l2_high','kernel_l2_low', ...
    'l2_5_10','theta_rmse_5_10','theta_dot_rmse_5_10','max_abs_theta_err_5_10', ...
    'kernel_updates','final_centers','mean_centers','max_centers','error'}));
end

function plot_best_error_comparison(best_runs)
figure('Color', 'w', 'Visible', 'off');
hold on;
plot(best_runs.aggressive.T, rad2deg(best_runs.aggressive.e(:,1)), ...
    'r-', 'LineWidth', 1.6, 'DisplayName', best_label(best_runs.aggressive));
plot(best_runs.current.T, rad2deg(best_runs.current.e(:,1)), ...
    'b--', 'LineWidth', 1.6, 'DisplayName', best_label(best_runs.current));
xlabel('Time [s]');
ylabel('\theta error [deg]');
title('\theta Tracking Error: Best Aggressive vs Best Mild/Current');
legend('show', 'Location', 'best');
grid on;
saveas(gcf, fullfile('Results', 'schedule_kernel_laws_theta_error.png'));
savefig(gcf, fullfile('Results', 'schedule_kernel_laws_theta_error.fig'));

figure('Color', 'w', 'Visible', 'off');
hold on;
plot(best_runs.aggressive.T, rad2deg(best_runs.aggressive.e(:,2)), ...
    'r-', 'LineWidth', 1.6, 'DisplayName', best_label(best_runs.aggressive));
plot(best_runs.current.T, rad2deg(best_runs.current.e(:,2)), ...
    'b--', 'LineWidth', 1.6, 'DisplayName', best_label(best_runs.current));
xlabel('Time [s]');
ylabel('\theta dot error [deg/s]');
title('\theta dot Tracking Error: Best Aggressive vs Best Mild/Current');
legend('show', 'Location', 'best');
grid on;
saveas(gcf, fullfile('Results', 'schedule_kernel_laws_theta_dot_error.png'));
savefig(gcf, fullfile('Results', 'schedule_kernel_laws_theta_dot_error.fig'));
end

function label = best_label(run_data)
label = sprintf('%s best: %s, L2_{5-10}=%.4g', ...
    run_data.schedule, run_data.case_name, run_data.l2_5_10);
end

function result = empty_result()
result = struct( ...
    'schedule', '', ...
    'case_name', '', ...
    'kernel_mode', NaN, ...
    'kernel_l2_high', NaN, ...
    'kernel_l2_low', NaN, ...
    'l2_5_10', NaN, ...
    'l2_full', NaN, ...
    'theta_rmse_5_10', NaN, ...
    'theta_dot_rmse_5_10', NaN, ...
    'max_abs_theta_err_5_10', NaN, ...
    'control_rms', NaN, ...
    'control_max_abs', NaN, ...
    'kernel_updates', NaN, ...
    'final_centers', NaN, ...
    'mean_centers', NaN, ...
    'max_centers', NaN, ...
    'T', [], ...
    'e', [], ...
    'error', '');
end
