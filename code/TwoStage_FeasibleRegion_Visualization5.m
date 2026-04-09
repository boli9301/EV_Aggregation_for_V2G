function [] = TwoStage_FeasibleRegion_Visualization5(P_compare, P_homothetic, P_gen_affine, P_cluster, P_box)
    %% 绘图 - 五种方法对比 (Minkowski, Homothetic, Gen. Affine, Cluster, Box)
    disp('--- 两阶段灵活性可视化 (5种方法对比) ---');
    
    % --- 投影到 P_ch(1) vs P_ch(2) 平面 ---
    dims_to_project = [1, 2];  
    
    % 辅助函数：安全投影和计算面积
    function [P_2d, area] = safe_proj_area(P, name)
        if isempty(P) || P.isEmptySet
            P_2d = Polyhedron(); area = 0;
            fprintf('%s 为空\n', name);
        else
            P_2d = P.projection(dims_to_project);
            area = P_2d.volume;
            fprintf('%s 面积: %.2f\n', name, area);
        end
    end

    [P_comp_2d, area_comp] = safe_proj_area(P_compare, '精确 Minkowski');
    [P_homo_2d, area_homo] = safe_proj_area(P_homothetic, 'Homothetic');
    [P_ga_2d, area_ga] = safe_proj_area(P_gen_affine, 'Gen. Affine');
    [P_clu_2d, area_clu] = safe_proj_area(P_cluster, 'Cluster-based');
    [P_box_2d, area_box] = safe_proj_area(P_box, 'Box Approx');

    fprintf('=============================\n');

    % --- 绘图 ---
    figure('Name', '聚合方法全面对比', 'Position', [100, 100, 1000, 800], 'Color', 'w');
    hold on;
    
    % 定义样式
    style_comp = {'Color', 'k', 'LineWidth', 2, 'LineStyle', '-'};      % 黑实线
    style_homo = {'Color', [0, 0.45, 0.74], 'LineWidth', 2, 'LineStyle', '-'}; % 蓝实线
    style_ga   = {'Color', [0.85, 0.33, 0.10], 'LineWidth', 2, 'LineStyle', '-'}; % 橙实线
    style_clu  = {'Color', [0.47, 0.67, 0.19], 'LineWidth', 2, 'LineStyle', '-'}; % 绿实线
    style_box  = {'Color', 'k', 'LineWidth', 2, 'LineStyle', '--'};    % 黑虚线 (Box)

    legend_handles = [];
    legend_str = {};

    % 辅助绘图函数
    function h = plot_poly(P_2d, style, label, area_val, ref_area)
        h = [];
        if ~P_2d.isEmptySet
            V = P_2d.V;
            if ~isempty(V) && size(V, 2) >= 2
                k = convhull(V(:,1), V(:,2));
                h = plot(V(k,1), V(k,2), style{:});
                
                pct = 0;
                if ref_area > 1e-6, pct = (area_val/ref_area)*100; end
                
                if contains(label, 'Minkowski')
                    legend_str{end+1} = sprintf('%s (Area=%.0f)', label, area_val);
                else
                    legend_str{end+1} = sprintf('%s (Area=%.0f, %.1f%%)', label, area_val, pct);
                end
            end
        end
    end

    % 按顺序绘制 (大的在下，小的在上，或者线型区分)
    % 1. 精确 (基准)
    h1 = plot_poly(P_comp_2d, style_comp, 'Exact Minkowski', area_comp, area_comp);
    if ~isempty(h1), legend_handles = [legend_handles, h1]; end
    
    % 2. Gen Affine
    h2 = plot_poly(P_ga_2d, style_ga, 'Gen. Affine', area_ga, area_comp);
    if ~isempty(h2), legend_handles = [legend_handles, h2]; end

    % 3. Cluster
    h3 = plot_poly(P_clu_2d, style_clu, 'Cluster-based', area_clu, area_comp);
    if ~isempty(h3), legend_handles = [legend_handles, h3]; end

    % 4. Homothetic
    h4 = plot_poly(P_homo_2d, style_homo, 'Homothetic', area_homo, area_comp);
    if ~isempty(h4), legend_handles = [legend_handles, h4]; end
    
    % 5. Box (最上层，虚线)
    h5 = plot_poly(P_box_2d, style_box, 'Box Approx', area_box, area_comp);
    if ~isempty(h5), legend_handles = [legend_handles, h5]; end

    % 美化
    title('Flexibility Aggregation Methods Comparison (T=2)', 'FontSize', 16, 'FontWeight', 'bold');
    xlabel('P_{ch}(t=1) [kW]', 'FontSize', 14);
    ylabel('P_{ch}(t=2) [kW]', 'FontSize', 14);
    
    if ~isempty(legend_handles)
        legend(legend_handles, legend_str, 'Location', 'best', 'FontSize', 10);
    end
    
    grid on; box on; axis equal;
    set(gca, 'FontSize', 12);
    hold off;
end
