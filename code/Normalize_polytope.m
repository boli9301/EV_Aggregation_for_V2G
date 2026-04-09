function [H_norm, h_norm, success] = Normalize_polytope(H, h, H_box, h_box)
    s_inv = sdpvar(1, 1, 'full');
    r_inv = sdpvar(size(H,2), 1, 'full');
    G_inv = sdpvar(size(H,1), size(H,1), 'full');
    
    % 约束条件
    constraints = [];
    constraints = [constraints ,s_inv >= 0];
    constraints = [constraints, G_inv >= 0];
    constraints = [constraints, G_inv * H == H_box];
    constraints = [constraints, G_inv * h <= s_inv * h_box + H_box * r_inv];
    
    % 求解
    ops = sdpsettings('verbose', 0, 'solver', 'gurobi');
    result = optimize(constraints, s_inv, ops);
    
    if result.problem == 0
        s_val = value(s_inv);
        r_val = value(r_inv);
    
        if s_val > 0 && ~any(isinf(r_val)) && ~any(isnan(r_val))
            beta_inv = 1 / s_val;
            ts_inv = -r_val / s_val;
    
            % 保存归一化后的多面体参数
            H_norm = H;
            h_norm = beta_inv * h + H * ts_inv;
            success = true;
        else
            H_norm = H;
            h_norm = h;
            success = false;
        end
    else
        disp('归一化失败\n');
        H_norm = H;
        h_norm = h;
        success = false;
    end
end