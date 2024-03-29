function Part1_driver()
    mat_size = [10^1, 10^2, 10^3, 10^4, 10^5, 10^6];
    tol = -2; % We are only interested in the tendency in part1
    max_iter = 100;
    norm_type = 2;
    weight = [1.0, 1.3, 1.6]; % Weights for SOR method

    for i = 1:length(mat_size)
        fprintf('Processing Matrix Size: %d\n', mat_size(i));
        [table_data, plotHandle, matSizeLabel] = main_part1(mat_size(i), tol, max_iter, norm_type, weight);
        
        % Make the figure visible if not already
        display(table_data)
        figure(plotHandle);
        
        % Customization of figure properties
        set(plotHandle, 'Name', sprintf('Results for Matrix Size %d', matSizeLabel), 'NumberTitle', 'off');
    
    end
end

function [table_data, plotHandle, matSizeLabel] = main_part1(mat_size, tol, max_iter, norm_type, weight)
    [A, D, L, U] = mat_creation(mat_size);
    b = ones(mat_size, 1);
    x0 = zeros(mat_size, 1);
    
    % Solver instances
    Jacobi_solver = Jacobi_method(A, D, L, U, b, x0, tol, max_iter, norm_type);
    SOR_solver = SOR_method(A, D, L, U, b, x0, tol, max_iter, norm_type);
    CG_solver = Conjugate_Gradient(A, b, x0, tol, max_iter, norm_type);
    
    % Solve and collect results
    [ ~ , ~, Jacobi_result] = Jacobi_solver.main();
    [ ~ , ~, SOR_1_0_result] = SOR_solver.main(weight(1));
    [ ~ , ~, SOR_1_3_result] = SOR_solver.main(weight(2));
    [ ~ , ~, SOR_1_6_result] = SOR_solver.main(weight(3));
    [ ~ , ~, CG_result] = CG_solver.main();

    
    % Plot results
    table_data = tabulize_results(Jacobi_result, SOR_1_0_result, SOR_1_3_result, SOR_1_6_result, CG_result, max_iter);
    plotHandle = plot_result(Jacobi_result, SOR_1_0_result, SOR_1_3_result, SOR_1_6_result, CG_result, max_iter);
    matSizeLabel = mat_size;
end

function table_data = tabulize_results(result1, result2, result3, result4, result5, max_iter)
    iteration_indices = [1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 100];
    method_names = {'Jacobi', 'SOR(1.0)', 'SOR(1.3)', 'SOR(1.6)', 'CG'};
    table_data = zeros(length(iteration_indices), length(method_names));
    for i = 1:length(iteration_indices)
        iteration = iteration_indices(i);
        table_data(i, 1) = result1(iteration);
        table_data(i, 2) = result2(iteration);
        table_data(i, 3) = result3(iteration);
        table_data(i, 4) = result4(iteration);
        table_data(i, 5) = result5(iteration);
    end

    % Display the table with row and column labels
    row_labels = arrayfun(@(x) sprintf('Iteration %d', x), iteration_indices, 'UniformOutput', false);
    column_labels = method_names;
    table_data = array2table(table_data, 'RowNames', row_labels, 'VariableNames', column_labels);
end


function plot_handle = plot_result(result1, result2, result3, result4, result5, max_iter)
    
    plot_handle = figure;
    
    hold on;
    
    methodColors = {'k-', 'b--', 'r-.', 'g:', 'm-'}; 
    
    x_axis = 1:1:max_iter;
    semilogy(x_axis, result1, methodColors{1}, 'LineWidth', 2, 'DisplayName', 'Jacobi');
    semilogy(x_axis, result2, methodColors{2}, 'LineWidth', 2, 'DisplayName', 'SOR(1.0)');
    semilogy(x_axis, result3, methodColors{3}, 'LineWidth', 2, 'DisplayName', 'SOR(1.3)');
    semilogy(x_axis, result4, methodColors{4}, 'LineWidth', 2, 'DisplayName', 'SOR(1.6)');
    semilogy(x_axis, result5, methodColors{5}, 'LineWidth', 2, 'DisplayName', 'CG');
    
    hold off;
    
    title('Relative Residual Error vs Iteration with Interval: 50');
    xlabel('Iteration');
    ylabel('Relative Residual Error');
    legend('show', 'Location', 'best'); % Show legend with 'DisplayName' of each plot
end

function [A, D, L, U] = mat_creation(n)
    D = sparse(1:n, 1:n, repmat(2.1, 1, n), n, n); % diagonal
    L = sparse(2:n,1:n-1,ones(1,n-1),n,n); % lower triangle
    U = L'; % upper triangle
    
    A = D - L - U ;
end



