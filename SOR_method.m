%SOR_solver = SOR(A, D, L , U, b, x0, tol, max_iter, norm_type); % weight for SOR(possibly array; multiple omegas)

classdef SOR_method
    properties
        A
        D
        L
        U
        b
        x0
        tol
        max_iter
        norm_type
    end
    
    methods
        function obj = SOR_method(A, D, L, U, b, x0, tol, max_iter, norm_type)
            
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
        
        function [x_min, flag, result, conv_iter] = main(obj, weight)
            flag = 0 ;
            result = zeros(1, max(1, ceil(obj.max_iter / 1)));
            resultIndex = 1;
            
            LHS = (obj.D - weight * obj.L);
            RHS = (1-weight)*obj.D + weight*obj.U;
            b_w = weight*obj.b;
            
            x = obj.x0;
            for iter = 1: obj.max_iter
                x = LHS \ (RHS * x + b_w);

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
end
