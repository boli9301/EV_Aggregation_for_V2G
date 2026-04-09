function [] = TwoStage_FeasibleRegion_Visualization2(P_input , P_compare)
    %% 绘图
    disp('--- 两阶段灵活性可视化 ---');
    
    % --- 投影到 P_ch(1) vs P_ch(2) 平面 ---
    dims_to_project = [1, 2];  % 对应 P_ch(t=1) 和 P_ch(t=2)
    P_input_2d = P_input.projection(dims_to_project);
    P_compare_2d = P_compare.projection(dims_to_project);
    % --- 计算面积 ---
    area_input = P_input_2d.volume;
    area_compare = P_compare_2d.volume;
    % --- 打印结果 ---
    fprintf('灵活性面积为: %.4f\n', area_input);
    if area_compare > 1e-6
        fprintf('聚合灵活性面积:       %.4f (占精确值的 %.2f%%)\n', area_input, (area_input/area_compare)*100);
    else
        fprintf('精确面积出错，聚合灵活性面积:   %.4f\n', area_input);
    end
    fprintf('=============================\n');

    % --- 使用顶点方法绘制轮廓线 ---
    figure('Name', '聚合方法精度对比', 'Position', [100, 100, 900, 700]);
    hold on;
    % 定义颜色和线型
    color_compare = [0, 0, 0];              % 黑色
    color_input = [0, 0.45, 0.74];    % 蓝色
    line_width = 3;

    % 获取顶点并绘制轮廓
    legend_handles = [];
    legend_str = {};

    % 绘制精确闵可夫斯基和 (黑色实线)
    if ~P_compare_2d.isEmptySet
        V_compare = P_compare_2d.V;
        if ~isempty(V_compare) && size(V_compare, 2) >= 2
            k = convhull(V_compare(:,1), V_compare(:,2));
            h1 = plot(V_compare(k,1), V_compare(k,2), '-', ...
                'Color', color_compare, 'LineWidth', line_width);
            legend_handles = [legend_handles, h1];
            legend_str{end+1} = sprintf('精确闵可夫斯基和 (面积=%.0f)', area_compare);
        end
    end

    % 绘制平移+缩放结果 (蓝色实线)
    if ~P_input_2d.isEmptySet
        V_input = P_input_2d.V;
        if ~isempty(V_input) && size(V_input, 2) >= 2
            k = convhull(V_input(:,1), V_input(:,2));
            h2 = plot(V_input(k,1), V_input(k,2), '-', ...
                'Color', color_input, 'LineWidth', line_width);
            legend_handles = [legend_handles, h2];
            if area_input > 1e-6
                legend_str{end+1} = sprintf('聚合灵活性 (面积=%.0f, %.1f%%)', area_input, (area_input/area_compare)*100);
            else
                legend_str{end+1} = sprintf('精确面积出错，聚合灵活性 (面积=%.0f)', area_input);
            end
        end
    end

    % 图表美化
    title('内逼近聚合方法精度对比 (T=2)', 'FontSize', 16, 'FontWeight', 'bold');
    xlabel('充电功率 P_ch(t=1) [kW]', 'FontSize', 14);
    ylabel('充电功率 P_ch(t=2) [kW]', 'FontSize', 14);

    % 添加图例
    if ~isempty(legend_handles)
        legend(legend_handles, legend_str, 'Location', 'northeast', 'FontSize', 12);
    end

    % 网格和坐标轴美化
    grid on;
    box on;
    axis equal;
    set(gca, 'FontSize', 12);
    hold off;
end