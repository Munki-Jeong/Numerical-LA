%test whether each method class is properly implemented


A = [4, -1, 0; -1, 4, -1; 0, -1, 4];

D = diag(diag(A));
L = tril(A) - D;
U = triu(A) - D;


b = [1; 2; 5];
x0 = [0; 0; 0];
tol = 1e-6 ;
max_iter = 1000;
norm_type = 2;

Jacobi_solver = Jacobi_method(A, D, L, U, b, x0, tol, max_iter, norm_type);
CG_solver = Conjugate_Gradient(A, b, x0, tol, max_iter, norm_type);
SOR_solver = SOR_method(A, D, L , U, b, x0, tol, max_iter, norm_type); 

[x_min_Jacobi, flag, ~] = Jacobi_solver.main();
[x_min_SOR_1_3, flag, ~ ] = SOR_solver.main(1.3);
[x_min_CG, flag, ~] = CG_solver.main();


disp(x_min_Jacobi);
disp(x_min_SOR_1_3);
disp(x_min_CG);
