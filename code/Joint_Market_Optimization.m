%  Joint_Market_Optimization
clear; 
%clc; 
close all;
warning('off', 'all');

data_dir = fullfile(pwd, 'mat_data');
agg_file = fullfile(data_dir, 'Aggregation_Results.mat'); 
load(agg_file, 'H_avg', 'h_avg_store', 'timesave' , 'Mean_beta' , ...
    'beta_CVaR_agg', 't_CVaR_agg' , 'beta_CVaR_save' ,  't_CVaR_save' , ...
    'beta_c_agg' , 't_c_agg' , 'beta_c_save' ,  't_c_save' , ...
    'beta_n_agg', 't_n_agg' ,'beta_Bonf_agg' , 't_Bonf_agg' ,...
    'K', 'cfg', 'cluster_indices', 'batts');
for k = 1:K
    t_Bonf_agg{k}(t_Bonf_agg{k} < 1e-04) = 0;
end

price_file = fullfile(data_dir, 'PJM_Market_Data_Processed.mat');
load(price_file , 'X_LMP', 'X_Res');

T = 24;
dt = 1; 
Samples = [X_LMP; X_Res]'; 

rng(1);
S = 8; 
[idx, C] = kmeans(Samples, S);
LMP_Scenarios = C(:, 1:24)';
Pr_Res_Scenarios = C(:, 25:48)';
LMP_Scenarios = max(0, LMP_Scenarios);
Pr_Res_Scenarios = max(0, Pr_Res_Scenarios);
Scen_Probs = zeros(S, 1);
N_total_days = size(Samples, 1);
for s = 1:S
    Scen_Probs(s) = sum(idx == s) / N_total_days;
end
LMP_Base = sum(LMP_Scenarios .* Scen_Probs', 2);
Price_Res_Base = sum(Pr_Res_Scenarios .* Scen_Probs', 2);
C_rep = 600000;   
Prob_Res_Call = 0.1; 
psi = 0;       
alpha = 0.95; 
  
%%  Case 1: Uncoordinated Charging
P_dumb_total = zeros(T, 1);
eta_ch_assumed = 0.95; 

for i = 1:cfg.EVsNum
    E_cap = batts.Ecap(i);
    E_init = batts.E0_exp_ratio(i) * E_cap;
    E_target = batts.Etarget_ratio(i) * E_cap;
    P_max = batts.Pmax_ch_base(i, 1);
    
    T_arr = round(batts.Tarrive_exp(i)); 
    T_dep = round(batts.Tdepart_exp(i));
    E_need_batt = max(0, E_cap - E_init);
    
    curr_t = T_arr;
    
    while E_need_batt > 0 && curr_t < T_dep && curr_t <= T        
        P_grid_req = E_need_batt / eta_ch_assumed / dt;
        p_act_grid = min(P_max, P_grid_req);
        
        idx = curr_t + 1; 
        if idx > T, break; end
        
        P_dumb_total(idx) = P_dumb_total(idx) + p_act_grid;
        E_charged_batt = p_act_grid * dt * eta_ch_assumed;
        E_need_batt = E_need_batt - E_charged_batt;
        
        curr_t = curr_t + 1;
    end
end
P_dumb_MW = P_dumb_total / 1000;
Cost_Energy_C1 = sum(LMP_Base .* P_dumb_MW) * dt; 
Stats_C1.Revenue_Energy = -Cost_Energy_C1; 
Stats_C1.Revenue_Res = 0;
Stats_C1.Cost_Deg = 0; 
Stats_C1.Total_Profit = Stats_C1.Revenue_Energy;

fprintf('  Case 1  $%.2f\n', Stats_C1.Total_Profit);

%% 3. Case 2: Energy
[Stats_C2, Res_C2] = run_market_opt_cvar_simple(true, LMP_Scenarios, Pr_Res_Scenarios, Scen_Probs, ...
                                         C_rep, Prob_Res_Call, H_avg, h_avg_store, ...
                                         beta_c_agg, t_c_agg, K, T, dt, batts, cluster_indices, psi, alpha);
fprintf('  Case 2 : $%.2f\n', Stats_C2.Total_Profit);

%% 4. Case 3: Joint Market With SOC uncertainty psi = 0
[Stats_C3, Res_C3] = run_market_opt_cvar_simple(false, LMP_Scenarios, Pr_Res_Scenarios, Scen_Probs, ...
                                         C_rep, Prob_Res_Call, H_avg, h_avg_store, ...
                                         beta_CVaR_agg, t_CVaR_agg, K, T, dt, batts, cluster_indices, psi, alpha);
fprintf('  Case 3 : $%.2f\n', Stats_C3.Total_Profit);

%% 5. Case 4: (Homothetic) psi = 0
% -------------------------------------------------------------------------

[Stats_C4, Res_C4] = run_market_opt_cvar_Nocluster(false, LMP_Scenarios, Pr_Res_Scenarios, Scen_Probs, ...
                                         C_rep, Prob_Res_Call, H_avg, h_avg_store, ...
                                         beta_n_agg, t_n_agg, T, dt, batts, psi, alpha);

fprintf('  Case 4 : $%.2f\n', Stats_C4.Total_Profit);

%% 6. Case 5: Joint Market Without SOC uncertainty psi = 0
[Stats_C5, Res_C5] = run_market_opt_cvar_simple(false, LMP_Scenarios, Pr_Res_Scenarios, Scen_Probs, ...
                                         C_rep, Prob_Res_Call, H_avg, h_avg_store, ...
                                         beta_c_agg, t_c_agg, K, T, dt, batts, cluster_indices, psi, alpha);
fprintf('  Case 5 : $%.2f\n', Stats_C5.Total_Profit);

%% 7. Case 6: Joint Market Without SOC uncertainty psi = 1
[Stats_C6, Res_C6] = run_market_opt_cvar_simple(false, LMP_Scenarios, Pr_Res_Scenarios, Scen_Probs, ...
                                         C_rep, Prob_Res_Call, H_avg, h_avg_store, ...
                                         beta_CVaR_agg, t_CVaR_agg, K, T, dt, batts, cluster_indices, 1, alpha);
fprintf('  Case 6 : $%.2f\n', Stats_C6.Total_Profit);

%% 8. Case 7: Joint Market-Bonf
[Stats_C7, Res_C7] = run_market_opt_cvar_simple(false, LMP_Scenarios, Pr_Res_Scenarios, Scen_Probs, ...
                                         C_rep, Prob_Res_Call, H_avg, h_avg_store, ...
                                         beta_Bonf_agg, t_Bonf_agg, K, T, dt, batts, cluster_indices, psi, alpha);
fprintf('  Case 7 : $%.2f\n', Stats_C7.Total_Profit);




%%  Visual

color_scheme = {
    '#9bbf8a', '#82afda', '#f79059', '#e7dbd3', '#c2bdde', ...
    '#8dcec8', '#add3e2', '#3480b8', '#ffbe7a', '#fa8878', '#c82423'
};

colors = color_scheme(1:S);
line_width = 2;

figure('Position', [100, 100, 1200, 500]);

subplot(1, 2, 1);
hold on;
for s = 1:S
    plot(1:24, LMP_Scenarios(:, s), ...
        'LineWidth', line_width, ...
        'Color', colors{s}, ...
        'DisplayName', sprintf('Scenario %d (%.1f%%)', s, Scen_Probs(s)*100));
end
hold off;

title('LMP Scenarios (Day-Ahead Market)', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Hour of Day', 'FontSize', 12);
ylabel('Price ($/MWh)', 'FontSize', 12);
grid off;
xlim([1 24]);
xticks(1:24);
set(gca, 'FontSize', 11);

lgd1 = legend('Location', 'best', 'FontSize', 10);
lgd1.Title.String = 'Scenarios with Probabilities';
lgd1.Title.FontWeight = 'bold';

subplot(1, 2, 2);
hold on;
for s = 1:S
    plot(1:24, Pr_Res_Scenarios(:, s), ...
        'LineWidth', line_width, ...
        'Color', colors{s}, ...
        'DisplayName', sprintf('Scenario %d (%.1f%%)', s, Scen_Probs(s)*100));
end
hold off;

title('Reserve MCP Scenarios (Ancillary Services)', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Hour of Day', 'FontSize', 12);
ylabel('Price ($/MW)', 'FontSize', 12);
grid off;
xlim([1 24]);
xticks(1:24);
set(gca, 'FontSize', 11);

lgd2 = legend('Location', 'best', 'FontSize', 10);
lgd2.Title.String = 'Scenarios with Probabilities';
lgd2.Title.FontWeight = 'bold';

sgtitle('Market Price Scenarios with Probabilities (CVaR Clustering)', ...
    'FontSize', 16, 'FontWeight', 'bold', 'Color', '#2E4053');

fprintf('\n=================================================================================\n');
fprintf('%-15s %-15s %-15s %-15s %-15s\n', 'Case', 'Energy Rev($)', 'Reserve Rev($)', 'Aging Cost($)', 'Net Profit($)');
fprintf('%-15s %-15.2f %-15.2f %-15.2f %-15.2f\n', '1(Uncoord)', Stats_C1.Revenue_Energy, Stats_C1.Revenue_Res, Stats_C1.Cost_Deg, Stats_C1.Total_Profit);
fprintf('%-15s %-15.2f %-15.2f %-15.2f %-15.2f\n', '2(Energy)', Stats_C2.Revenue_Energy, Stats_C2.Revenue_Res, Stats_C2.Cost_Deg, Stats_C2.Total_Profit);
fprintf('%-15s %-15.2f %-15.2f %-15.2f %-15.2f\n', '3(Joint_cluster_Uncertainty)psi = 0', Stats_C3.Revenue_Energy, Stats_C3.Revenue_Res, Stats_C3.Cost_Deg, Stats_C3.Total_Profit);
fprintf('%-15s %-15.2f %-15.2f %-15.2f %-15.2f\n', '4(Joint_Nocluster)', Stats_C4.Revenue_Energy, Stats_C4.Revenue_Res, Stats_C4.Cost_Deg, Stats_C4.Total_Profit);
fprintf('%-15s %-15.2f %-15.2f %-15.2f %-15.2f\n', '5(Joint_cluster_No_Uncertainty)', Stats_C5.Revenue_Energy, Stats_C5.Revenue_Res, Stats_C5.Cost_Deg, Stats_C5.Total_Profit);
fprintf('%-15s %-15.2f %-15.2f %-15.2f %-15.2f\n', '6(Joint_cluster_Uncertainty)psi = 1', Stats_C6.Revenue_Energy, Stats_C6.Revenue_Res, Stats_C6.Cost_Deg, Stats_C6.Total_Profit);
fprintf('%-15s %-15.2f %-15.2f %-15.2f %-15.2f\n', '6(Joint_Bonf)', Stats_C7.Revenue_Energy, Stats_C7.Revenue_Res, Stats_C7.Cost_Deg, Stats_C7.Total_Profit);

fprintf('=================================================================================\n');

figure('Color', 'w', 'Position', [100, 50, 1000, 600]);
x_ax = 1:T;

yyaxis left;
pos_energy = max(0, Res_C3.P_net);
reserve_cap = Res_C3.R; 
neg_energy = min(0, Res_C3.P_net);

b = bar(x_ax, [neg_energy, pos_energy, reserve_cap], 'stacked');


b(1).FaceColor = [0.4 0.7 0.9]; 
b(1).DisplayName = 'Charging Power (MW)';

b(2).FaceColor = [0.2 0.6 0.8]; 
b(2).DisplayName = 'Discharging Power (MW)';

b(3).FaceColor = [0.2 0.8 0.2]; 
b(3).DisplayName = 'Reserve Capacity (MW)';
b(3).FaceAlpha = 0.7; 

ylabel('Power (MW)', 'FontSize', 12, 'FontWeight', 'bold');
yline(0, 'k-', 'HandleVisibility', 'off');
ylim([-15, 15]); 


yyaxis right;
p1 = plot(x_ax, LMP_Base, 'r-o', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'LMP ($/MWh)');
hold on;
p2 = plot(x_ax, Price_Res_Base, 'm--s', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Reserve Price ($/MW)');

ylabel('Price ($)', 'FontSize', 12, 'FontWeight', 'bold');
ax = gca;
ax.YColor = 'r'; 
ylim([0, 100]); 

xlabel('Hour of Day', 'FontSize', 12);
title('Optimal Dispatch vs. Market Prices (Zhang et al. 2022 Data)', 'FontSize', 14);
grid on;
xlim([0.5, 24.5]);
set(gca, 'XTick', 1:24);

legend([b(1), b(2), b(3), p1, p2], 'Location', 'northwest', 'NumColumns', 2);

E_cap_total = sum(batts.Ecap);        
E_init_total = sum(batts.E0_exp_ratio .* batts.Ecap);  
E_target_total = sum(batts.Etarget_ratio .* batts.Ecap); 
E_min_total = sum(batts.Emin_ratio .* batts.Ecap);  

P_max_ch_total = sum(batts.Pmax_ch_base(:, 1));  
P_max_dis_total = sum(batts.Pmax_dis_base(:, 1)); 

T_arrive = round(batts.Tarrive_exp(1));
T_depart = round(batts.Tdepart_exp(1))+1;

time_hours_fine = T_arrive:0.1:T_depart;

E_upper = zeros(size(time_hours_fine));
E_lower = zeros(size(time_hours_fine));

for i = 1:length(time_hours_fine)
    t = time_hours_fine(i);

    dt_from_arr = max(0, t - T_arrive);
    E_upper(i) = min(E_cap_total, E_init_total + P_max_ch_total * dt_from_arr);
    
    E_discharge_path = max(E_min_total, E_init_total - P_max_dis_total * dt_from_arr);
    time_remaining = max(0, T_depart - t);
    E_required_for_target = E_target_total - P_max_ch_total * time_remaining;
    
    E_lower(i) = max(E_discharge_path, E_required_for_target);
    E_lower(i) = min(E_lower(i), E_cap_total);
end

E_upper_MWh = E_upper / 1000;
E_lower_MWh = E_lower / 1000;

eta_ch_avg = mean(batts.eta_ch);
eta_dis_avg = mean(batts.eta_dis);

E_traj_C1 = zeros(T, 1);
E_curr = E_init_total;
for t = 1:T
    if t >= T_arrive && t <= T_depart
        P_charge_kW = P_dumb_MW(t) * 1000; 
        E_curr = E_curr + (P_charge_kW * eta_ch_avg) * dt;
        E_curr = min(E_curr, E_cap_total);
    end
    E_traj_C1(t) = E_curr;
end
E_traj_C1_MWh = E_traj_C1 / 1000;

E_traj_C2 = zeros(T, 1);
E_curr = E_init_total;
for t = 1:T
    if t >= T_arrive && t <= T_depart
        P_net_kW = Res_C2.P_net(t) * 1000;
        if P_net_kW >= 0 
            E_curr = E_curr - (P_net_kW / eta_dis_avg) * dt;
        else 
            E_curr = E_curr - (P_net_kW * eta_ch_avg) * dt;
        end
        E_curr = min(max(E_curr, E_min_total), E_cap_total);
    end
    E_traj_C2(t) = E_curr;
end
E_traj_C2_MWh = E_traj_C2 / 1000;

E_traj_C3 = zeros(T, 1);
E_curr = E_init_total;
for t = 1:T
    if t >= T_arrive && t <= T_depart
        P_net_kW = Res_C3.P_net(t) * 1000;
        P_res_exp_kW = Res_C3.R(t) * 1000 * Prob_Res_Call; 
        P_total_kW = P_net_kW + P_res_exp_kW;
        
        if P_total_kW >= 0 
            E_curr = E_curr - (P_total_kW / eta_dis_avg) * dt;
        else 
            E_curr = E_curr - (P_total_kW * eta_ch_avg) * dt;
        end
        E_curr = min(max(E_curr, E_min_total), E_cap_total);
    end
    E_traj_C3(t) = E_curr;
end
E_traj_C3_MWh = E_traj_C3 / 1000;

E_traj_C4 = zeros(T, 1);
E_curr = E_init_total;
for t = 1:T
    if t >= T_arrive && t <= T_depart
        P_net_kW = Res_C4.P_net(t) * 1000;
        P_res_exp_kW = Res_C4.R(t) * 1000 * Prob_Res_Call; 
        P_total_kW = P_net_kW + P_res_exp_kW;
        
        if P_total_kW >= 0 
            E_curr = E_curr - (P_total_kW / eta_dis_avg) * dt;
        else 
            E_curr = E_curr - (P_total_kW * eta_ch_avg) * dt;
        end
        E_curr = min(max(E_curr, E_min_total), E_cap_total);
    end
    E_traj_C4(t) = E_curr;
end
E_traj_C4_MWh = E_traj_C4 / 1000;

E_traj_C5 = zeros(T, 1);
E_curr = E_init_total;
for t = 1:T
    if t >= T_arrive && t <= T_depart
        P_net_kW = Res_C5.P_net(t) * 1000;
        P_res_exp_kW = Res_C5.R(t) * 1000 * Prob_Res_Call; 
        P_total_kW = P_net_kW + P_res_exp_kW;
        
        if P_total_kW >= 0 
            E_curr = E_curr - (P_total_kW / eta_dis_avg) * dt;
        else 
            E_curr = E_curr - (P_total_kW * eta_ch_avg) * dt;
        end
        E_curr = min(max(E_curr, E_min_total), E_cap_total);
    end
    E_traj_C5(t) = E_curr;
end
E_traj_C5_MWh = E_traj_C5 / 1000;

E_traj_C6 = zeros(T, 1);
E_curr = E_init_total;
for t = 1:T
    if t >= T_arrive && t <= T_depart
        P_net_kW = Res_C6.P_net(t) * 1000;
        P_res_exp_kW = Res_C6.R(t) * 1000 * Prob_Res_Call;
        P_total_kW = P_net_kW + P_res_exp_kW;
        
        if P_total_kW >= 0 
            E_curr = E_curr - (P_total_kW / eta_dis_avg) * dt;
        else 
            E_curr = E_curr - (P_total_kW * eta_ch_avg) * dt;
        end
        E_curr = min(max(E_curr, E_min_total), E_cap_total);
    end
    E_traj_C6(t) = E_curr;
end
E_traj_C6_MWh = E_traj_C6 / 1000;


figure('Color', 'w', 'Position', [100, 100, 1000, 500]);
hold on;
custom_colors_rgb = [
    155, 191, 138;  % #9bbf8a (Green)
    130, 175, 218;  % #82afda (Blue)
    247, 144, 89;   % #f79059 (Orange)
    231, 219, 211;  % #e7dbd3 (Beige)
    194, 189, 222;  % #c2bdde (Light Purple)
    141, 206, 200;  % #8dcec8 (Teal)
    173, 211, 226;  % #add3e2 (Light Blue)
    52,  128, 184;  % #3480b8 (Dark Blue)
    255, 190, 122;  % #ffbe7a (Peach)
    250, 136, 120;  % #fa8878 (Salmon)
    200, 36,  35;   % #c82423 (Dark Red)
] / 255;
fill([time_hours_fine, fliplr(time_hours_fine)], ...
     [E_upper_MWh, fliplr(E_lower_MWh)], ...
     [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.4, ... 
     'DisplayName', 'Feasible Energy Envelope');
plot(time_hours_fine, E_upper_MWh, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
plot(time_hours_fine, E_lower_MWh, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
yline(E_init_total/1000, 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
yline(E_target_total/1000, 'g-.', 'LineWidth', 1.5, 'DisplayName', 'Target Energy');
valid_idx = T_arrive:T_depart;
% % Case 1: Uncoordinated 
% plot(valid_idx, E_traj_C1_MWh(valid_idx), '-^', 'LineWidth', 2, 'MarkerSize', 6, ...
%      'Color', custom_colors_rgb(1,:),'DisplayName', 'Case 1: Uncoordinated');
% 
% % Case 2: Energy Only 
% plot(valid_idx, E_traj_C2_MWh(valid_idx), '-o', 'LineWidth', 2, 'MarkerSize', 6, ...
%      'Color', custom_colors_rgb(2,:) , 'DisplayName', 'Case 2: Energy Only');

% Case 3: Joint Market-uncertainty 
plot(valid_idx, E_traj_C3_MWh(valid_idx), '-s', 'LineWidth', 2, 'MarkerSize', 6, ...
     'Color', custom_colors_rgb(3,:) ,'DisplayName', 'Case 3: Uncertainty pis=0 (Expected)');

% % Case 4: Joint Market-nocluster 
% plot(valid_idx, E_traj_C4_MWh(valid_idx), '-s', 'LineWidth', 2, 'MarkerSize', 6, ...
%      'Color', custom_colors_rgb(4,:) , 'DisplayName', 'Case 4: Homothetic (Expected)');
% 
% % Case 5: Joint Market -cluster
% plot(valid_idx, E_traj_C5_MWh(valid_idx), '-s', 'LineWidth', 2, 'MarkerSize', 6, ...
%      'Color', custom_colors_rgb(5,:) ,'DisplayName', 'Case 5: Cluster (Expected)');
% 
% % Case 6: Joint Market -uncertainty psi=1
% plot(valid_idx, E_traj_C6_MWh(valid_idx), '-s', 'LineWidth', 2, 'MarkerSize', 6, ...
%      'Color', custom_colors_rgb(6,:) ,'DisplayName', 'Case 6: Uncertainty pis=1 (Expected)');
title('Aggregated Energy Trajectory vs. Flexibility Envelope', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Time (Hour)', 'FontSize', 12);
ylabel('Total Energy (MWh)', 'FontSize', 12);
legend('Location', 'best', 'NumColumns', 2);
grid off;
xlim([T_arrive-0.5, T_depart+0.5]);
xticks(T_arrive:T_depart);
text(T_arrive, E_init_total/1000, '  Start', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
text(T_depart, E_target_total/1000, '  Target', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
hold off;



figure('Color', 'w', 'Position', [100, 50, 1000, 600]);
x_ax = 1:T;
if ~isfield(Res_C3, 'P_net_k')
    error('请先运行修改后的优化函数以获取簇级功率数据 (Results.P_net_k)');
end
P_cluster_matrix = Res_C7.P_net_k; 
P_cluster_pos = max(0, P_cluster_matrix); 
P_cluster_neg = min(0, P_cluster_matrix); 

custom_colors_rgb = [
    155, 191, 138;  % #9bbf8a (Green)
    130, 175, 218;  % #82afda (Blue)
    247, 144, 89;   % #f79059 (Orange)
    231, 219, 211;  % #e7dbd3 (Beige)
    194, 189, 222;  % #c2bdde (Light Purple)
    141, 206, 200;  % #8dcec8 (Teal)
    173, 211, 226;  % #add3e2 (Light Blue)
    52,  128, 184;  % #3480b8 (Dark Blue)
    255, 190, 122;  % #ffbe7a (Peach)
    250, 136, 120;  % #fa8878 (Salmon)
    200, 36,  35;   % #c82423 (Dark Red)
] / 255; 

num_colors = size(custom_colors_rgb, 1);
hold on;
b_pos = bar(x_ax, P_cluster_pos, 'stacked', 'EdgeColor', 'none', 'BarWidth', 0.7);
b_neg = bar(x_ax, P_cluster_neg, 'stacked', 'EdgeColor', 'none', 'BarWidth', 0.7);
for k = 1:K

    color_idx = mod(k-1, num_colors) + 1;
    this_color = custom_colors_rgb(color_idx, :);
    
    b_pos(k).FaceColor = this_color;
    b_neg(k).FaceColor = this_color;
    

    b_pos(k).FaceAlpha = 0.9;
    b_neg(k).FaceAlpha = 0.9;
    

    b_pos(k).DisplayName = sprintf('Cluster %d', k);
    b_neg(k).HandleVisibility = 'off'; 
end
ylabel('Cluster Power Dispatch (MW)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.2]);
yline(0, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', [0.4 0.4 0.4]); 
ax = gca;
ax.YColor = [0.2 0.2 0.2]; 
title('Optimal Dispatch by Cluster vs. Market Prices', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Hour of Day', 'FontSize', 12);
grid off;
ax.GridAlpha = 0.3; 
xlim([0.5, 24.5]);
set(gca, 'XTick', 1:24);
lgd = legend('Location', 'eastoutside');
title(lgd, 'Clusters');
hold off;




counts = histcounts(cluster_indices, 1:K+1);
[max_num, k_target] = max(counts);
P_net_agg_MW = Res_C3.P_net_k(:, k_target);
P_net_agg_kW = P_net_agg_MW * 1000; 
P_dis_agg = max(0, P_net_agg_kW);
P_ch_agg  = max(0, -P_net_agg_kW);
x_agg = [P_ch_agg; P_dis_agg]; 
Beta_sum_k = beta_c_agg(k_target);
Chi_sum_k  = t_c_agg{k_target}; 
idx_cars = find(cluster_indices == k_target);
num_cars_in_cluster = length(idx_cars);
P_individual_profiles = zeros(T, num_cars_in_cluster);

E_caps = batts.Ecap(idx_cars);
P_maxs = batts.Pmax_dis_base(idx_cars, 1);
betas_i = beta_c_save{k_target};
chis_i = t_c_save{k_target};
for n = 1:num_cars_in_cluster

    b_i = betas_i(n);
    c_i = chis_i(:, n);

    x_i = (b_i / Beta_sum_k) * (x_agg - Chi_sum_k) + c_i;

    p_ch_i = x_i(1:T);
    p_dis_i = x_i(T+1:end);

    p_ch_i = max(0, p_ch_i);
    p_dis_i = max(0, p_dis_i);

    P_individual_profiles(:, n) = p_dis_i - p_ch_i;
end



% -------------------------------------------------------------------------
figure('Color', 'w', 'Position', [100, 100, 1200, 600]);
hold on;

custom_colors = [
    0.61 0.75 0.54;  
    0.51 0.69 0.85;  
    0.97 0.56 0.35;  
    0.91 0.86 0.83;  
    0.76 0.74 0.87;  
    0.55 0.81 0.78;  
    0.68 0.83 0.89;  
    0.20 0.50 0.72;  
    1.00 0.75 0.48;  
    0.98 0.53 0.47;  
    0.78 0.14 0.14;  
];

num_plot = min(num_cars_in_cluster, 300); 
plot_indices = 1:num_plot; 

time_vec = 1:T;

time_stairs = [time_vec, time_vec(end)];

for i = 1:num_plot
    idx = plot_indices(i);

    power_data = P_individual_profiles(:, idx);
    

    power_stairs = [power_data; power_data(end)];
    

    color_idx = mod(i-1, size(custom_colors, 1)) + 1;
    line_color = custom_colors(color_idx, :);
    

    stairs(time_stairs, power_stairs, ...
        'Color', line_color, ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('EV %d', idx));
end


P_avg_profile = mean(P_individual_profiles, 2);
P_avg_stairs = [P_avg_profile; P_avg_profile(end)];
stairs(time_stairs, P_avg_stairs, ...
    'k-', 'LineWidth', 3, ...
    'DisplayName', 'Average Profile');

yline(0, 'k-', 'LineWidth', 1.5, 'Alpha', 0.3, 'HandleVisibility', 'off');


grid off;
box off;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3);

title(sprintf('Disaggregated EV Dispatch Profiles - Step Plot (Cluster #%d, %d EVs Shown)', ...
      k_target, num_plot), 'FontSize', 16, 'FontWeight', 'bold');
xlabel('Time (Hour)', 'FontSize', 14);
ylabel('Individual Power (kW)', 'FontSize', 14);
xlim([1, T]);
set(gca, 'XTick', 1:T);

if num_plot <= 15
    legend('Location', 'best', 'FontSize', 10, 'NumColumns', 2);
else
    legend('Average Profile', 'Location', 'best', 'FontSize', 12);
end

y_lim = get(gca, 'YLim');

if y_lim(1) < 0
    area_x = [1, T, T, 1];
    area_y = [0, 0, y_lim(1), y_lim(1)];
    fill(area_x, area_y, [0.93, 0.96, 1.0], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.2, ...
        'HandleVisibility', 'off', 'DisplayName', 'Charging Region');
    

    area_y2 = [0, 0, y_lim(2), y_lim(2)];
    fill(area_x, area_y2, [1.0, 0.96, 0.93], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.2, ...
        'HandleVisibility', 'off', 'DisplayName', 'Discharging Region');
    

    uistack(findobj(gca, 'Type', 'patch'), 'bottom');
    

    text(T-2, y_lim(1)*0.7, 'Charging', ...
        'Color', [0.2, 0.4, 0.8], 'FontSize', 12, 'FontWeight', 'bold');
    text(T-2, y_lim(2)*0.7, 'Discharging (V2G)', ...
        'Color', [0.8, 0.4, 0.2], 'FontSize', 12, 'FontWeight', 'bold');
end

hold off;


annotation('textbox', [0.02, 0.02, 0.3, 0.05], ...
    'String', sprintf('Total EVs in cluster: %d | Time slots: %d', num_cars_in_cluster, T), ...
    'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'none');



%% 
%  This function forces the optimization of all EVs as a whole (Cluster k = 1).
%  The input parameters beta_agg should be beta_n_agg, and t_agg should be t_n_agg

function [Stats, Results] = run_market_opt_cvar_Nocluster(is_energy_only, LMP_Scens, Pr_Res_Scens, Scen_Probs, ...
                                                C_rep, Prob_Res_Call, H_avg, h_avg_store, ...
                                                beta_agg, t_agg, T, dt, batts, psi, alpha)

    P_net_total = sdpvar(T, 1);
    R_total     = sdpvar(T, 1);

    K_sim = 1; 
    
    P_ch_k  = sdpvar(T, 1);
    P_dis_k = sdpvar(T, 1);
    R_k     = sdpvar(T, 1); 

    theta = sdpvar(1, 1);
    phi   = sdpvar(size(LMP_Scens, 2), 1);
    
    N_seg = 64;
    Constraints = [];
    
    Cost_Aging_Total_Expression = 0;
    

    Constraints = [Constraints, P_ch_k >= 0, P_dis_k >= 0, R_k >= 0];
    
    if is_energy_only
        Constraints = [Constraints, R_k == 0];
    end

    E_cap_total = sum(batts.Ecap); 
    P_max_total = sum(batts.Pmax_dis_base(:, 1));
    

    omega_total = calculate_aging_coefficients(E_cap_total, C_rep, 0.95, N_seg);
    

    P_dis_seg = sdpvar(T, N_seg);
    Constraints = [Constraints, P_dis_k == sum(P_dis_seg, 2)];
    limit_per_seg = P_max_total / N_seg;
    Constraints = [Constraints, 0 <= P_dis_seg <= limit_per_seg];
    
    for m = 1:N_seg
        Cost_Aging_Total_Expression = Cost_Aging_Total_Expression + ...
            sum(omega_total(m) * (P_dis_seg(:, m) / 1000)) * dt;
    end
    

    
    H = H_avg;          
    if iscell(h_avg_store)
        h0 = h_avg_store{1}; 
    else
        h0 = h_avg_store;
    end

    if length(beta_agg) > 1
        warning('If the length of the input beta_agg is greater than 1, only the first element will be taken as the parameter for a single cluster.');
        beta = beta_agg(1);
    else
        beta = beta_agg;
    end
    
    if iscell(t_agg)
        chi = t_agg{1};
    else
        chi = t_agg;
    end
    
    x_state = [P_ch_k; P_dis_k]; 
    x_res_vec = [zeros(T, 1); R_k];   
    

    Constraints = [Constraints, H * (x_state - chi) <= beta * h0];

    Constraints = [Constraints, H * (x_state + Prob_Res_Call .* x_res_vec - chi) <= beta * h0];
    

    Constraints = [Constraints, P_net_total == P_dis_k - P_ch_k];
    Constraints = [Constraints, R_total == R_k];
    

    P_net_MW = P_net_total / 1000;
    R_MW     = R_total / 1000;
    
    S_scen = length(Scen_Probs);
    Profit_Scenarios = [];
    
    for s = 1:S_scen
        lambda_en  = LMP_Scens(:, s);
        lambda_res = Pr_Res_Scens(:, s);
        
        Rev_En = sum(lambda_en .* P_net_MW) * dt;
        Rev_Res = sum(lambda_res .* R_MW);
        Rev_Res_Call = sum(Prob_Res_Call * (R_MW .* lambda_en)) * dt; 
        
        Profit_s = Rev_En + Rev_Res + Rev_Res_Call - Cost_Aging_Total_Expression;
        
        Profit_Scenarios = [Profit_Scenarios; Profit_s];
        Constraints = [Constraints, phi(s) >= theta - Profit_s];
    end
    Constraints = [Constraints, phi >= 0];
    
    Expected_Profit = sum(Scen_Probs .* Profit_Scenarios);
    CVaR_Term = theta - (1 / (1 - alpha)) * sum(Scen_Probs .* phi);
    Objective = (1 - psi) * Expected_Profit + psi * CVaR_Term;

    ops = sdpsettings('solver', 'gurobi', 'verbose', 0);
    sol = optimize(Constraints, -Objective, ops);

    if sol.problem == 0
        lambda_en_avg = LMP_Scens * Scen_Probs;
        lambda_res_avg = Pr_Res_Scens * Scen_Probs;
        
        Stats.Revenue_Energy = value(sum(lambda_en_avg .* P_net_MW) * dt);
        Stats.Revenue_Res    = value(sum(lambda_res_avg .* R_MW));
        Stats.Cost_Deg       = value(Cost_Aging_Total_Expression);
        Stats.Total_Profit   = value(Expected_Profit);
        Stats.CVaR           = value(CVaR_Term);
        
        Results.P_net = value(P_net_MW);
        Results.R     = value(R_MW);
    else
        warning('求解失败 (No Cluster Case)');
        Stats = struct('Revenue_Energy',0,'Revenue_Res',0,'Cost_Deg',0,'Total_Profit',-Inf,'CVaR',0);
        Results = struct('P_net',zeros(T,1),'R',zeros(T,1));
    end
end