%Jacobi_solver = Jacobi(D, L, U, b, x0, tol, max_iter, norm_type);
classdef Jacobi_method
    properties
        A;
        D;
        L;
        U;
        b;
        x0;
        tol;
        max_iter;
        norm_type;
    end
    
    methods
        function obj = Jacobi_method(A, D, L, U, b, x0, tol, max_iter, norm_type)
            obj.A = A;
            obj.D = D;
            obj.L = L;
            obj.U = U;
            obj.b = b;
            obj.x0 = x0;
            obj.tol = tol;
            obj.max_iter = max_iter;
            obj.norm_type = norm_type;
        end

        function [x_min, flag, result, conv_iter] = main(obj)

            flag = 0;
            result = zeros(1, max(1, ceil(obj.max_iter / 1)));
            resultIndex = 1;

            x = obj.x0;
            for iter = 1:obj.max_iter
                x = obj.D \ (obj.b + (obj.L + obj.U)* x);

                r = obj.b - obj.A * x;
                criteria = norm(r, obj.norm_type)/norm(obj.b, obj.norm_type);
                if mod(iter, 1) == 0
                    result(resultIndex) = criteria;
                    resultIndex = resultIndex + 1;
                end

                if criteria < obj.tol
                    flag = 1;
                    conv_iter = iter;
                    break;
                end


            end
            x_min = x;

            if flag == 0
                conv_iter = obj.max_iter+1;
            end
        end
    end

    
    % Define additional class elements here
    
end





% function result = jacobiMethod(obj, A, b, x0, maxIterations, tolerance)
%     % Jacobi method implementation
    
%     % Initialize variables
%     n = size(A, 1);
%     x = x0;
%     xPrev = x0;
%     iteration = 0;
    
%     % Perform Jacobi iterations
%     while iteration < maxIterations
%         for i = 1:n
%             sigma = 0;
%             for j = 1:n
%                 if j ~= i
%                     sigma = sigma + A(i, j) * xPrev(j);
%                 end
%             end
%             x(i) = (b(i) - sigma) / A(i, i);
%         end
        
%         % Check convergence
%         if norm(x - xPrev) < tolerance
%             break;
%         end
        
%         % Update variables for next iteration
%         xPrev = x;
%         iteration = iteration + 1;
%     end
    
%     % Return the result
%     result = x;
% end