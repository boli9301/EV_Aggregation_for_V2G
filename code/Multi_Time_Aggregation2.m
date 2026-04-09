%% Section 1
clc;
clear;
ops = sdpsettings('verbose', 0, 'solver', 'gurobi' , 'gurobi.TimeLimit', 120);
cfg.T = 24;         
cfg.deltaT = 1;     
cfg.EVsNum = 1000;  
cfg.w1 = 0.5;       
cfg.w2 = 1 - cfg.w1;
cfg.epsilon_soc = 0.1; 
cfg.numerical_tolerance = 0;

rng(1); 

batts.Tarrive_exp = 8*ones(cfg.EVsNum , 1);  
batts.Tdepart_exp = 20*ones(cfg.EVsNum , 1); 
trunc_norm = @(mu, sigma, low, high, sz) ...
    max(low, min(high, normrnd(mu, sigma, sz)));
Type1 = [0.50, 40, 3.0, 15 , 60    11, 2.5 ,7 , 15]; 
Type2 = [0.50, 40, 3.0, 15 , 60 ,  90.0, 20 , 60 , 120]; 
Type3 = [0, 350, 50, 280 , 570 , 500.0, 50.0 , 300 , 600];
rand_types = rand(cfg.EVsNum, 1);
veh_types = ones(cfg.EVsNum, 1);
veh_types(rand_types > Type1(1)) = 2;
Ecap_vec = zeros(cfg.EVsNum, 1);
Pmax_vec = zeros(cfg.EVsNum, 1);
for i = 1:cfg.EVsNum
    if veh_types(i) == 1
        param = Type1;
    elseif veh_types(i) == 2
        param = Type2;
    else
        param = Type3;
    end
    Ecap_vec(i) = trunc_norm(param(2), param(3), param(4), param(5), [1,1]);
    Pmax_vec(i) = trunc_norm(param(6), param(7), param(8), param(9), [1,1]);
end
batts.Ecap = Ecap_vec;
batts.Pmax_ch_base  = repmat(Pmax_vec, 1, cfg.T);
batts.Pmax_dis_base = repmat(Pmax_vec, 1, cfg.T);
batts.E0_exp_ratio = trunc_norm(0.4, 0.15, 0.2, 0.6, [cfg.EVsNum, 1]);
batts.eta_ch  = 0.95 * ones(cfg.EVsNum, 1);
batts.eta_dis = 0.95 * ones(cfg.EVsNum, 1); 
batts.Emin_ratio = 0.2 * ones(cfg.EVsNum, 1); 
batts.Emax_ratio = 1.0 * ones(cfg.EVsNum, 1);
batts.Etarget_ratio = 0.95 * ones(cfg.EVsNum, 1); 


%% Section 2
H_i = cell(cfg.EVsNum, 1);
h_i = cell(cfg.EVsNum, 1);
valid_ev_count = 0;

for i = 1:cfg.EVsNum
    E0_exp_val = batts.E0_exp_ratio(i) * batts.Ecap(i);
    try
        [H_i{i}, h_i{i}] = local_build_polytope(cfg, batts, i, E0_exp_val, batts.Tarrive_exp(i), batts.Tdepart_exp(i));
        if ~isempty(H_i{i}) && ~any(isinf(h_i{i})) && ~any(isnan(h_i{i}))
            valid_ev_count = valid_ev_count + 1;
        else
            fprintf('Error \n');
        end
    catch ME
        fprintf('Error: EV %d : %s\n', i, ME.message);
        H_i{i} = [];
        h_i{i} = [];
    end
end

eta_ch_avg = mean(batts.eta_ch);
eta_dis_avg = mean(batts.eta_dis);
A_EV_avg = tril(ones(cfg.T)).*(cfg.deltaT*eta_ch_avg);
B_EV_avg = tril(ones(cfg.T)).*(cfg.deltaT/eta_dis_avg);
C_EV_avg = -eta_ch_avg*cfg.deltaT.*ones(1,cfg.T);
D_EV_avg = -cfg.deltaT/eta_dis_avg.*ones(1,cfg.T);
H_avg = [eye(cfg.T)  , zeros(cfg.T) ;   
    -eye(cfg.T) , zeros(cfg.T) ;   
    zeros(cfg.T) , eye(cfg.T)   ;   
    zeros(cfg.T) , -eye(cfg.T)  ;   
    A_EV_avg    ,  -B_EV_avg   ;   
    -A_EV_avg   ,  B_EV_avg    ;   
    C_EV_avg    ,  -D_EV_avg   ];  
fprintf('Success! ValidNum: %d / %d\n', valid_ev_count, cfg.EVsNum);

%%  Section 3 -- homothetic 
tic;
h_avg = mean(cat(2, h_i{:}), 2);
s_n_val = zeros(cfg.EVsNum, 1);
r_n_val = zeros(size(H_avg,2), cfg.EVsNum);
valid_n_flags = false(cfg.EVsNum, 1); 

parfor n = 1:cfg.EVsNum
    H_current = H_i{n};
    h_current = h_i{n};
    s_n_var = sdpvar(1, 1, 'full');
    r_n_var = sdpvar(size(H_avg, 2), 1, 'full');
    G_n_var = sdpvar(size(H_current, 1), size(H_current, 1), 'full');
   
    constraints_n = [];
    constraints_n = [constraints_n, G_n_var >= 0];
    constraints_n = [constraints_n, s_n_var >= 0];
    constraints_n = [constraints_n, G_n_var * H_avg == H_current];
    constraints_n = [constraints_n, G_n_var * h_avg <= s_n_var * h_current + H_current * r_n_var];

    obj_n = s_n_var;
    sol = optimize(constraints_n, obj_n, ops);
    
    if sol.problem == 0
        s_n_val(n) = value(s_n_var);
        r_n_val(:, n) = value(r_n_var);
        valid_n_flags(n) = true;
    else
        s_n_val(n) = NaN;
        disp('Error!')
    end
end

success_count = sum(valid_n_flags);
fprintf('Success: %d / %d\n', success_count, cfg.EVsNum);


valid_indices = find(valid_n_flags);
if isempty(valid_indices)
    error('All vehicle solutions have failed, rendering aggregation unfeasible ! ');
end

s_valid = s_n_val(valid_indices);
r_valid = r_n_val(:, valid_indices);

h_n_homothetic = cell(success_count,1);
beta_n_calc = 1 ./ s_valid;
t_n_calc = zeros(size(r_valid));
for k = 1:length(s_valid)
    t_n_calc(:, k) = -r_valid(:, k) ./ s_valid(k);
end
beta_n_agg = sum(beta_n_calc);
t_n_agg = sum(t_n_calc, 2);
fprintf('Method 1 (Homothetic) Beta: %.4f\n', beta_n_agg);
Mean_beta.Homothetic = beta_n_agg/cfg.EVsNum;
disp('Homothetic Time:');
timesave.Homothetic = toc  
yalmip('clear');

%% Section 4 -- The proposed Method 


%% Section 5 -- The proposed Method With SOC Uncertainty(Proposed Approximation DRCCP)


%% Section 6 -- The proposed Method With SOC Uncertainty(Bonf-Approximation DRJCC)
tic;

cfg.epsilon = 0.05; 

uncertain_rows = [(4*cfg.T + 1) : (6*cfg.T), size(H_avg,1)];
num_uncertain_constraints = length(uncertain_rows);

E0_ratio_std = 0.02;
E0_ratio_scenarios = zeros(cfg.EVsNum, N_scenarios);
rng(2);
for i = 1:cfg.EVsNum
    E0_ratio_scenarios(i,:) = normrnd(batts.E0_exp_ratio(i), E0_ratio_std, [1, N_scenarios]);
    E0_ratio_scenarios(i, E0_ratio_scenarios(i,:) < 0.1) = 0.1;
    E0_ratio_scenarios(i, E0_ratio_scenarios(i,:) > 1.0) = 1.0;
end

Hrd_EV_Scenarios = zeros(size(H_avg,1), cfg.EVsNum, N_scenarios);
for i = 1:cfg.EVsNum
    Emax = batts.Emax_ratio(i) * batts.Ecap(i);
    Emin = batts.Emin_ratio(i) * batts.Ecap(i);
    Etarget = batts.Etarget_ratio(i) * batts.Ecap(i);
    for s = 1:N_scenarios
        E0_s = E0_ratio_scenarios(i,s) * batts.Ecap(i);
        h_energy_up_s = ones(cfg.T,1) * (Emax - E0_s);
        h_energy_down_s = ones(cfg.T,1) * (E0_s - Emin);
        h_energy_depart_s = E0_s - Etarget;
        h_i_s = [batts.Pmax_ch_base(i,:)';        
                 zeros(cfg.T,1);            
                 batts.Pmax_dis_base(i,:)';         
                 zeros(cfg.T,1);             
                 h_energy_up_s;               
                 h_energy_down_s;             
                 h_energy_depart_s];           
        Hrd_EV_Scenarios(:, i, s) = h_i_s;
    end
end

gamma_EV = sqrt(var(Hrd_EV_Scenarios, 0, 3)); 
beta_Bonf_agg = zeros(K, 1);
t_Bonf_agg = cell(K, 1);

for k = 1:K
    ev_indices_in_cluster = find(cluster_indices == k);
    num_ev_in_cluster = length(ev_indices_in_cluster);
    if num_ev_in_cluster == 0, continue; end

    h_i_cluster = h_i(ev_indices_in_cluster);
    h_avg_cluster = mean(cat(2, h_i_cluster{:}), 2);

    H_i_cluster = H_i(ev_indices_in_cluster);
    H_avg_cluster = zeros(size(H_i_cluster{1},1), size(H_i_cluster{1},2));
    for n = 1:num_ev_in_cluster
        H_avg_cluster = H_avg_cluster + H_i_cluster{n};
    end
    H_avg_cluster = H_avg_cluster ./ num_ev_in_cluster;

    s_u_val = zeros(num_ev_in_cluster, 1);
    r_u_val = zeros(size(H_avg_cluster,2), num_ev_in_cluster);

    parfor n = 1:num_ev_in_cluster
        global_idx = ev_indices_in_cluster(n);
        s_u_local = sdpvar(1, 1, 'full');
        G_u_local = sdpvar(size(H_i_cluster{n},1), size(H_i_cluster{n},1), 'full');
        r_u_local = sdpvar(size(H_avg_cluster,2), 1, 'full');
        zs_local = sdpvar(num_uncertain_constraints, 1, 'full');
        vs_local = sdpvar(num_uncertain_constraints, 1, 'full');
        epsilons_local = sdpvar(num_uncertain_constraints, 1, 'full');

        constraints_local = [];
        constraints_local = [constraints_local, s_u_local >= 0];
        constraints_local = [constraints_local, G_u_local >= 0];
        constraints_local = [constraints_local, G_u_local * H_avg_cluster == H_i_cluster{n}];
        for power_idx = 1:4*cfg.T
            constraints_local = [constraints_local, ...
                G_u_local(power_idx,:) * h_avg_cluster <= s_u_local * h_i_cluster{n}(power_idx) + H_i_cluster{n}(power_idx,:) * r_u_local];
        end

        for uncert_idx = 1:num_uncertain_constraints
            constraint_row = uncertain_rows(uncert_idx);
            current_gamma = gamma_EV(constraint_row, global_idx);

            constraints_local = [constraints_local, ...
                zs_local(uncert_idx) * current_gamma + G_u_local(constraint_row,:) * h_avg_cluster <= ...
                s_u_local * h_i_cluster{n}(constraint_row) + H_i_cluster{n}(constraint_row,:) * r_u_local];

            constraints_local = [constraints_local, ...
                norm([2*vs_local(uncert_idx); 1]) <= 2*epsilons_local(uncert_idx) + 1];
            constraints_local = [constraints_local, ...
                norm([2; zs_local(uncert_idx) - vs_local(uncert_idx)]) <= zs_local(uncert_idx) + vs_local(uncert_idx)];
        end

        constraints_local = [constraints_local, sum(epsilons_local) <= cfg.epsilon];
        constraints_local = [constraints_local, zs_local >= 0];
        constraints_local = [constraints_local, vs_local >= 0];
        constraints_local = [constraints_local, epsilons_local >= 0];

        obj = s_u_local;
        sol_local = optimize(constraints_local, obj, ops);

        if sol_local.problem == 0
            s_u_val(n) = value(s_u_local);
            r_u_val(:, n) = value(r_u_local);
        elseif sol_local.problem == 11
            warning('EV %d Time Out!', global_idx);
            s_u_val(n) = Inf; 
        else
            warning('EV %d Error', global_idx);
            s_u_val(n) = Inf; 
        end
    end

    failed_indices = find(isinf(s_u_val));
    vaild_Bonf_flags = cfg.EVsNum - length(failed_indices);
    if ~isempty(failed_indices)
        fprintf('  Cluster %d contains %d EV Error\n', k, length(failed_indices));
        for idx = failed_indices'
            s_u_val(idx) = 1; 
            r_u_val(:, idx) = zeros(size(H_avg_cluster,2), 1); 
        end
    end

    beta_u = 1 ./ s_u_val;
    beta_Bonf_agg(k) = sum(beta_u);

    t_u = zeros(size(H_avg_cluster,2), num_ev_in_cluster);
    for j = 1:num_ev_in_cluster
        t_u(:,j) = -r_u_val(:,j) ./ s_u_val(j);
    end
    t_Bonf_agg{k} = sum(t_u, 2);
end
Mean_beta.Bonf = sum(beta_Bonf_agg)/ vaild_Bonf_flags;
fprintf('Time:');
timesave.Bonf = toc  
yalmip('clear');
%% Section 7. Save
%save_path = fullfile(pwd, 'mat_data', 'Aggregation_Results.mat');
if ~exist(fullfile(pwd, 'mat_data'), 'dir')
    mkdir(fullfile(pwd, 'mat_data'));
end
%save(save_path, 'H_avg', 'h_avg_store', 'timesave' , 'Mean_beta' , 'beta_CVaR_agg', 't_CVaR_agg' , 'beta_CVaR_save' ,  't_CVaR_save' , 'beta_c_agg' , 't_c_agg' , 'beta_c_save' ,  't_c_save' , 'beta_n_agg', 't_n_agg' ,'beta_Bonf_agg' , 't_Bonf_agg' , 'K', 'cfg', 'cluster_indices', 'batts');
yalmip('clear');

