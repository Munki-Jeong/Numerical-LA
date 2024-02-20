% Test Script for Conjugate_Gradient Class

% Define a symmetric positive definite matrix A
n = 5; % Size of the matrix
A = gallery('lehmer', n);
A = A'*A; % Ensure A is symmetric positive definite

% Define b
b = rand(n, 1);

% Initial guess x0
x0 = zeros(n, 1);

% Tolerance and maximum number of iterations
tol = 1e-6;
max_iter = 1000;

% Norm type
norm_type = 2;

% Create an instance of the Conjugate_Gradient class
cg_solver = Conjugate_Gradient(A, b, x0, tol, max_iter, norm_type);

% Call the main method to solve the system
[x_min, flag, ~] = cg_solver.main();

% Compare the solution with MATLAB's built-in solver
x_actual = A\b;

% Display results
disp(['Conjugate Gradient Method Converged: ', num2str(flag)]);
disp(['CG Solution: ', num2str(x_min')]);
disp(['Actual Solution: ', num2str(x_actual')]);

% Check the accuracy of the solution
error_norm = norm(x_min - x_actual, norm_type);
disp(['Error Norm: ', num2str(error_norm)]);

% Assess convergence based on the specified tolerance
if error_norm < tol && flag == 1
    disp('Test Passed: The Conjugate Gradient solver converged to the correct solution.');
else
    disp('Test Failed: The Conjugate Gradient solver did not converge to the correct solution.');
end
