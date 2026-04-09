function [H, h] = local_build_polytope(cfg, batts, ev_idx, E0, Tarrive, Tdepart)
%  构建单个EV多面体的函数
    try
        T = cfg.T;
        deltaT = cfg.deltaT;
        eta_ch = batts.eta_ch(ev_idx);
        eta_dis = batts.eta_dis(ev_idx);
        Ecap = batts.Ecap(ev_idx);
        Emin = batts.Emin_ratio(ev_idx) * Ecap;
        Emax = batts.Emax_ratio(ev_idx) * Ecap;
        Etarget = batts.Etarget_ratio(ev_idx) * Ecap;

        
        % 确保时间参数的有效性
        Tarrive = max(0, min(T-4, round(Tarrive)));
        Tdepart = max(Tarrive+3, min(T-1, round(Tdepart)));% 至少停留3小时
        
        % 构建功率约束矩阵
        Pmax_ch_v = zeros(T, 1);
        Pmax_dis_v = zeros(T, 1);
        active_indices = (Tarrive+1):(Tdepart+1);
        active_indices = active_indices(active_indices >= 1 & active_indices <= T);
        
        if ~isempty(active_indices)
            Pmax_ch_v(active_indices) = batts.Pmax_ch_base(ev_idx , active_indices);
            Pmax_dis_v(active_indices) = batts.Pmax_dis_base(ev_idx , active_indices);
        end
        
        % 充放电功率约束 (4T)
        % P_ch <= Pmax_ch, P_ch >= 0, P_dis <= Pmax_dis, P_dis >= 0
        H_power = [eye(T), zeros(T); -eye(T), zeros(T); zeros(T), eye(T); zeros(T), -eye(T)];
        h_power = [Pmax_ch_v; zeros(T,1); Pmax_dis_v; zeros(T,1)];
        
        % SOC 约束 (2T)
        A_ch = tril(ones(T)) * (deltaT * eta_ch);
        A_dis = tril(ones(T)) * (-deltaT / eta_dis);
        
        % SOC上限: E0 + sum(eta_ch*P_ch*dt - P_dis/eta_dis*dt) <= Emax
        % => sum(eta_ch*P_ch*dt - P_dis/eta_dis*dt) <= Emax - E0
        H_soc_max = [A_ch, A_dis];
        h_soc_max = ones(T, 1) * (Emax - E0);
        
        % SOC下限: E0 + sum(eta_ch*P_ch*dt - P_dis/eta_dis*dt) >= Emin  
        % => -sum(eta_ch*P_ch*dt - P_dis/eta_dis*dt) <= E0 - Emin
        H_soc_min = [-A_ch, -A_dis];
        h_soc_min = ones(T, 1) * (E0 - Emin);
        
        % 离网SOC目标约束 (1)
        H_target = zeros(1, 2*T);
        h_target = E0 - Etarget;
        
        if Tdepart < T-1
            % 如果在调度期内离网，需要满足目标SOC
            % E(Tdepart) >= Etarget
            % => E0 + sum_{t=1}^{Tdepart}(eta_ch*P_ch(t) - P_dis(t)/eta_dis)*dt >= Etarget
            % => sum_{t=1}^{Tdepart}(eta_ch*P_ch(t) - P_dis(t)/eta_dis)*dt >= Etarget - E0
            % => -sum_{t=1}^{Tdepart}(eta_ch*P_ch(t) - P_dis(t)/eta_dis)*dt <= E0 - Etarget
            active_target_indices = 1:(Tdepart+1);
            active_target_indices = active_target_indices(active_target_indices >= 1 & active_target_indices <= T);
            if ~isempty(active_target_indices)
                H_target(1, active_target_indices) = -deltaT * eta_ch;
                H_target(1, T + active_target_indices) = deltaT / eta_dis;
            end
        end
        
        % 组合矩阵 (6T+1 rows)
        H = [H_power; H_soc_max; H_soc_min; H_target];
        h = [h_power; h_soc_max; h_soc_min; h_target];

        
        
        % 检查基本可行性
        if any(h < -1e10) || any(isinf(h)) || any(isnan(h))
            warning('EV %d 的多面体参数异常', ev_idx);
        end
        
    catch ME
        fprintf('构建EV %d 多面体时发生错误: %s\n', ev_idx, ME.message);
        H = [];
        h = [];
    end
end