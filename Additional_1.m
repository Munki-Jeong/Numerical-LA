A = [4 1 -1 1 ; 1 4 -1 -1 ; -1 -1 5 1; 1 -1 1 3];
b = [-2; -1; 0 ;1];
x0 = [0; 0; 0; 0];
tol = 1e-6;
max_iter = 2;
norm_type = 2;

CG_solver = Conjugate_Gradient(A, b, x0, tol, max_iter, norm_type);
x_exact = A \ b;
[x_min, flag, result, ~] = CG_solver.main();


disp("exact solution: ");
disp(x_exact);
disp("Solution from Conjugate Gradient after two iterations: ");
disp(x_min);

if flag == 0
    disp('Not Converged with in 1e-6');
else
    disp('Converged');
end

