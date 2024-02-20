% Test Script for SOR_Method Class with Different Weights

% Define a symmetric positive definite matrix A
n = 5; % Size of the matrix
A = gallery('lehmer', n);
A = A'*A; % Ensure A is symmetric positive definite to make SOR applicable

% Decompose A into D, L, and U components
D = diag(diag(A));
L = -1 * tril(A, -1);
U = -1 * triu(A, 1);

% Define b
b = rand(n, 1);

% Initial guess x0
x0 = zeros(n, 1);

% Tolerance and maximum number of iterations
tol = 1e-6;
max_iter = 1000;

% Norm type
norm_type = 2;

% Weights for SOR
weights = [1.0, 1.3, 1.6];

% Create an instance of the SOR_method class
for weight = weights
    sor_solver = SOR_method(A, D, L, U, b, x0, tol, max_iter, norm_type);
    
    % Call the main method to solve the system with the current weight
    [x_min, flag, ~] = sor_solver.main(weight);
    
    % Compare the solution with MATLAB's built-in solver
    x_actual = A\b;
    
    % Display results
    disp(['SOR Method with weight = ', num2str(weight), ' Converged: ', num2str(flag)]);
    disp(['SOR Solution: ', num2str(x_min')]);
    disp(['Actual Solution: ', num2str(x_actual')]);
    
    % Check the accuracy of the solution
    error_norm = norm(x_min - x_actual, norm_type);
    disp(['Error Norm: ', num2str(error_norm)]);
    
    % Assess convergence based on the specified tolerance
    if error_norm < tol && flag == 1
        disp(['Test Passed: The SOR solver with weight = ', num2str(weight), ' converged to the correct solution.']);
    else
        disp(['Test Failed: The SOR solver with weight = ', num2str(weight), ' did not converge to the correct solution.']);
    end
    disp(' '); % Add a newline for readability
end
