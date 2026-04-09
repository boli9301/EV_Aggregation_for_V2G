function [Stats, Results] = run_market_opt_cvar_simple(is_energy_only, LMP_Scens, Pr_Res_Scens, Scen_Probs, ...
                                                C_rep, Prob_Res_Call, H_avg, h_avg_store, ...
                                                beta_agg, t_agg, K, T, dt, batts, cluster_indices, psi, alpha)
    
    % --- 1. 变量定义 ---
    P_net_total = sdpvar(T, 1);
    R_total     = sdpvar(T, 1);
    
    % 簇级变量
    P_ch_k  = cell(K, 1);
    P_dis_k = cell(K, 1);
    P_dis_seg_k = cell(K, 1);
    R_k     = cell(K, 1); 
    
    % CVaR 变量
    theta = sdpvar(1, 1);
    phi   = sdpvar(size(LMP_Scens, 2), 1);
    
    N_seg = 64;
    Constraints = [];
    
    P_dis_sum_all = 0; P_ch_sum_all = 0; R_sum_all = 0;
    Cost_Aging_Total_Expression = 0;
    
    % --- 2. 簇级约束 ---
    for k = 1:K
        if beta_agg(k) <= 1e-6, continue; end
        
        P_ch_k{k}  = sdpvar(T, 1);
        P_dis_k{k} = sdpvar(T, 1);
        R_k{k}     = sdpvar(T, 1); 
        P_dis_seg_k{k} = sdpvar(T, N_seg);
        
        % 非负约束
        Constraints = [Constraints, P_ch_k{k} >= 0, P_dis_k{k} >= 0, R_k{k} >= 0];
        
        if is_energy_only
            Constraints = [Constraints, R_k{k} == 0];
        end
        
        % 老化模型
        idx_in_cluster = find(cluster_indices == k);
        E_cap_agg_k = sum(batts.Ecap(idx_in_cluster)); 
        P_max_agg_k = sum(batts.Pmax_dis_base(idx_in_cluster, 1));
        omega_k = calculate_aging_coefficients(E_cap_agg_k, C_rep, 0.95, N_seg);
        
        Constraints = [Constraints, P_dis_k{k} == sum(P_dis_seg_k{k}, 2)];
        limit_per_seg = P_max_agg_k / N_seg;
        Constraints = [Constraints, 0 <= P_dis_seg_k{k} <= limit_per_seg];
        
        for m = 1:N_seg
            Cost_Aging_Total_Expression = Cost_Aging_Total_Expression + ...
                sum(omega_k(m) * (P_dis_seg_k{k}(:, m) / 1000)) * dt;
        end
        
        % --- 几何灵活性约束 ---
        H = H_avg;          
        h0 = h_avg_store{k}; 
        beta = beta_agg(k);
        chi = t_agg{k};     
        
        x_state_k = [P_ch_k{k}; P_dis_k{k}]; 
        x_res_vec = [zeros(T, 1); R_k{k}];   
        
        % (a) 基准功率可行
        Constraints = [Constraints, H * (x_state_k - chi) <= beta * h0];
        
        % (b) 备用调用可行
        Constraints = [Constraints, H * (x_state_k + Prob_Res_Call .* x_res_vec - chi) <= beta * h0];
        
        % 累加
        P_dis_sum_all = P_dis_sum_all + P_dis_k{k};
        P_ch_sum_all  = P_ch_sum_all + P_ch_k{k};
        R_sum_all     = R_sum_all + R_k{k};
    end
    
    % --- 3. 市场层约束 ---
    Constraints = [Constraints, P_net_total == P_dis_sum_all - P_ch_sum_all];
    Constraints = [Constraints, R_total == R_sum_all];
    
    % --- 4. 目标函数 ---
    P_net_MW = P_net_total / 1000;
    R_MW     = R_total / 1000;
    
    S = length(Scen_Probs);
    Profit_Scenarios = [];
    
    for s = 1:S
        lambda_en  = LMP_Scens(:, s);
        lambda_res = Pr_Res_Scens(:, s);
        
        % (a) 能量收益
        Rev_En = sum(lambda_en .* P_net_MW) * dt;
        
        % (b) 备用容量收益
        Rev_Res = sum(lambda_res .* R_MW);
        
        % (c) 备用调用预期收益 (保留收益，忽略成本)
        Rev_Res_Call = sum(Prob_Res_Call * (R_MW .* lambda_en)) * dt;
        
        % (d) [已删除] 备用调用预期老化成本
        
        % 总利润 = 能量 + 容量 + 调用预期收益 - 仅基准老化
        Profit_s = Rev_En + Rev_Res + Rev_Res_Call - Cost_Aging_Total_Expression;
        
        Profit_Scenarios = [Profit_Scenarios; Profit_s];
        Constraints = [Constraints, phi(s) >= theta - Profit_s];
    end
    Constraints = [Constraints, phi >= 0];
    
    Expected_Profit = sum(Scen_Probs .* Profit_Scenarios);
    CVaR_Term = theta - (1 / (1 - alpha)) * sum(Scen_Probs .* phi);
    Objective = (1 - psi) * Expected_Profit + psi * CVaR_Term;
    
    % --- 5. 求解 ---
    ops = sdpsettings('solver', 'gurobi', 'verbose', 0);
    sol = optimize(Constraints, -Objective, ops);
    
    % --- 6. 结果 ---
    if sol.problem == 0
        lambda_en_avg = mean(LMP_Scens, 2);
        lambda_res_avg = mean(Pr_Res_Scens, 2);
        
        Stats.Revenue_Energy = value(sum(lambda_en_avg .* P_net_MW) * dt);
        Stats.Revenue_Res    = value(sum(lambda_res_avg .* R_MW));
        Stats.Cost_Deg       = value(Cost_Aging_Total_Expression);
        Stats.Total_Profit   = value(Expected_Profit);
        Stats.CVaR           = value(CVaR_Term);
        
        Results.P_net = value(P_net_MW);
        Results.R     = value(R_MW);
        P_net_k_val = zeros(T, K);
        for k = 1:K
            if beta_agg(k) > 1e-6
                % 计算该簇净功率: (放电 - 充电) / 1000 转 MW
                p_k_mw = (value(P_dis_k{k}) - value(P_ch_k{k})) / 1000;
                P_net_k_val(:, k) = p_k_mw;
            end
        end
        Results.P_net_k = P_net_k_val; % 将矩阵存入 Results 结构体
    else
        warning('求解失败');
        Stats = struct('Revenue_Energy',0,'Revenue_Res',0,'Cost_Deg',0,'Total_Profit',-Inf,'CVaR',0);
        Results = struct('P_net',zeros(T,1),'R',zeros(T,1));
    end
end