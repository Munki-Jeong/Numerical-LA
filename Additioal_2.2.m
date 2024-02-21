function Part2_driver()
    mat_size = [10^1, 10^2, 10^3, 10^4, 10^5, 10^6]; % Adjusted for practical computation
    tol = 1e-12;
    max_iter = 700;
    norm_type = 2;
    weights = [1.0, 1.3, 1.6]; % Weights for SOR method

    % Initialize array to store iteration counts for each method and matrix size
    iter_num = cell(length(weights) + 2, length(mat_size)); % +2 for Jacobi and CG


    for i = 1:length(mat_size)
        fprintf('Processing Matrix Size: %d\n', mat_size(i));
        iter_num(:, i) = main_part2(mat_size(i), tol, max_iter, norm_type, weights); %each iter_num(:, i): consist of Jacobi, SOR1.0, SOR1.3, SOR1.6, and CG
    end
    
    plot_part2(mat_size, iter_num, weights);
    
    % Plot the results
end

function iter_num = main_part2(mat_size, tol, max_iter, norm_type, weights)
    [A, D, L, U] = mat_creation(mat_size);
    b = ones(mat_size, 1);
    x0 = zeros(mat_size, 1);

    % Initialize array to store iterations for Jacobi, 3xSOR, and CG methods
    iterations = cell(1, length(weights) + 2); % +2 for Jacobi and CG
    
    % Jacobi solver
    Jacobi_solver = Jacobi_method(A, D, L, U, b, x0, tol, max_iter, norm_type);
    [~, flag, ~, Jacobi_iter] = Jacobi_solver.main();
    iterations{1} = Jacobi_iter;
    
    
    % SOR solver for each weight
    for w = 1:length(weights)
        SOR_solver = SOR_method(A, D, L, U, b, x0, tol, max_iter, norm_type);
        [~, flag, ~, SOR_iter] = SOR_solver.main(weights(w));
        iterations{w + 1} = SOR_iter;
    end
    
    % Conjugate Gradient solver
    CG_solver = Conjugate_Gradient(A, b, x0, tol, max_iter, norm_type);
    [~, flag, ~, CG_iter] = CG_solver.main();
    iterations{5} = CG_iter;
    
    iter_num = iterations;
end

function plot_part2(mat_size, iter_num, weights)
    plot_handle = figure;
    
    % Define method markers for the plot
    methodMarkers = {'ko', 'b^', 'rs', 'gd', 'm*'};
    methodNames = {'Jacobi', sprintf('SOR(%g)', weights(1)), sprintf('SOR(%g)', weights(2)), sprintf('SOR(%g)', weights(3)), 'CG'};
    
    hold on;
    
    % Determine the x-coordinates for each method's plot on the log scale
    x_values = log10(mat_size);
    
    % Get the maximum iteration number for y-axis limits
    max_iter_num = max(cellfun(@(c) max(c), iter_num(:)));
    
    % Draw vertical lines and plot iteration numbers for each method
    for idx = 1:length(x_values)
        % Plot iteration numbers for each method
        for methodIdx = 1:length(methodNames)
            y_value = iter_num{methodIdx, idx};
            plot(x_values(idx) * ones(size(y_value)), y_value, methodMarkers{methodIdx}, 'LineWidth', 2, 'DisplayName', methodNames{methodIdx});
        end
    end
    
    hold off;
    
    title('Iterations Needed for Convergence at Different Matrix Sizes');
    xlabel('log_{10}(Matrix Size)');
    ylabel('Iterations');
    
    % Set the x-ticks to be the log of matrix sizes and label them with the actual matrix sizes
    set(gca, 'XTick', x_values, 'XTickLabel', arrayfun(@num2str, mat_size, 'UniformOutput', false));
    ylim([0, 700]);
    
    % Set y-axis limits if needed
% Set the y-axis limits to zoom out

% Add the remaining code here

legend(methodNames, 'Location', 'northeastoutside');
grid on;
    
    legend(methodNames, 'Location', 'northeastoutside');
    grid on;
end




function [A, D, L, U] = mat_creation(n)
    D = sparse(1:n, 1:n, repmat(2.1, 1, n), n, n); % diagonal
    L = sparse(2:n,1:n-1,ones(1,n-1),n,n); % lower triangle
    U = L'; % upper triangle
    
    A = D - L - U ;
end
