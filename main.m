
function Part1_driverFunctionForPlots()
    
    mat_size = [10^1, 10^2, 10^3, 10^4, 10^5, 10^6];
    tol = -1; % Tolerance for convergence
    max_iter = 1000; % Maximum number of iterations
    norm_type = 2; % Norm type for convergence check
    weight = [1.0, 1.3, 1.6]; % Weights for SOR method




        
    
    end


    % Call the main_part1 function with the initialized parameters
    plotHandles = main_part1(tol, max_iter, norm_type, weight);

    % Display each plot generated for the different matrix sizes
    for i = 1:length(plotHandles)
        figureHandle = plotHandles{i}; % Retrieve the figure handle
        
        % Make the figure visible if not already
        figure(figureHandle);
        
        % Customization of figure properties can be done here if needed
        set(figureHandle, 'Name', sprintf('Results for Matrix Size %d', 10^i), 'NumberTitle', 'off');
        
        % Ensure the figure is correctly drawn and updated
        drawnow;
        
        % Optionally, you can add pauses or user prompts to control the display flow
        % pause; % Uncomment to require pressing a key before showing the next plot
    end
end


function Part1 = main_part1(mat_size, tol, max_iter, norm_type, weight) % weight for SOR(possibly array; multiple omegas)
    
    %structure: plot_each_matsize = {plot_handle1, plot_handle2, ..., plot_handle6}
    %where plot_handle1: 10^1, plot_handle2: 10^2, ..., plot_handle6: 10^6
    %each plot_handle is a figure consists of 5 plots for each method(in a single grid)

    %mat_size = [10^1, 10^2, 10^3, 10^4, 10^5, 10^6];
    % plot_each_matsize = cell(1, length(mat_size))

    for i = 1:length(mat_size)

        [A, D, L, U] = mat_creation(mat_size(i));
        b = ones(mat_size(i), 1);
        x0 = zeros(mat_size(i), 1);
        
        Jacobi_solver = Jacobi_method(A, D, L, U, b, x0, tol, max_iter, norm_type);
        SOR_solver = SOR_method(A, D, L , U, b, x0, tol, max_iter, norm_type); % weight for SOR(possibly array; multiple omegas)
        CG_solver = Conjugate_Gradient(A, b, x0, tol, max_iter, norm_type);
        
        [ ~ , flag, Jacobi_result] = Jacobi_solver.main(); % ~ for x_min, but we are not interested in the value itself for this purpose
        [ ~ , flag, SOR_1_0_result] = SOR_solver.main(weight(1));
        [ ~ , flag, SOR_1_3_result] = SOR_solver.main(weight(2));
        [ ~ , flag, SOR_1_6_result] = SOR_solver.main(weight(3));
        [ ~ , flag, CG_result] = CG_solver.main();

        result_each_method = {Jacobi_result, SOR_1_0_result, SOR_1_3_result, SOR_1_6_result, CG_result};
        Part1 =  plot_result(Jacobi_result, SOR_1_0_result, SOR_1_3_result, SOR_1_6_result, CG_result, max_iter);
        
    end
end

function plot_handle = plot_result(result1, result2, result3, result4, result5, max_iter)
    
    plot_handle = figure;

    hold on;

    methodColors = {'k-', 'b--', 'r-.', 'g:', 'm-'}; 
    
    x_axis = 1:50:max_iter;
    semilogy(x_axis, cell2mat(result1), methodColors{1}, 'LineWidth', 2, 'DisplayName', 'Jacobi');
    semilogy(x_axis, cell2mat(result2), methodColors{2}, 'LineWidth', 2, 'DisplayName', 'SOR(1.0)');
    semilogy(x_axis, cell2mat(result3), methodColors{3}, 'LineWidth', 2, 'DisplayName', 'SOR(1.3)');
    semilogy(x_axis, cell2mat(result4), methodColors{4}, 'LineWidth', 2, 'DisplayName', 'SOR(1.6)');
    semilogy(x_axis, cell2mat(result5), methodColors{5}, 'LineWidth', 2, 'DisplayName', 'CG');

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

% function plot_handle = plot_result(result, max_iter)
    
%     plot_handle = figure;

%     hold on;

%     x_axis = 1: 50 : max_iter;
%     semilogy(x_axis, cell2mat(result{1}));
%     semilogy(x_axis, cell2mat(result{2}));
%     semilogy(x_axis, cell2mat(result{3}));
%     semilogy(x_axis, cell2mat(result{4}));
%     semilogy(x_axis, cell2mat(result{5}));

%     hold off;
    
%     % Setting the title, labels, and legend outside the loop
%     title('Relative Residual Error vs Iteration with Interval: 50');
%     xlabel('Iteration');
%     ylabel('Relative Residual Error');
%     legend({'Jacobi', 'SOR(1.0)', 'SOR(1.3)', 'SOR(1.6)', 'CG'}, 'Location', 'best');
% end