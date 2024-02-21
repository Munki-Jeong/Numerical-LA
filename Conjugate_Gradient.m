classdef Conjugate_Gradient %flag 1: converged, 0: not converged
    %CONJUGATE_GRADIENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        A
        b
        x0
        tol
        max_iter
        norm_type
    end

    methods
        function obj = Conjugate_Gradient(A, b, x0, tol, max_iter, norm_type) %constructor
            obj.A = A;
            obj.b = b;
            obj.x0 = x0;
            obj.tol = tol;
            obj.max_iter = max_iter;
            obj.norm_type = norm_type;
        end
        
        function t = stepsize(obj,v, u, r)
            t = (v'*r)/(v'*u);
        end
        function x = update_x(obj, x, t, v)
            x = x + t*v;
        end
        function r = residual(obj, w)
            r = obj.b - w;
        end
        function s = orthogonalize(obj, v, u, r)
            s = -1 * (r'*u)/(v'*u);
        end

        function [x_min , flag, result, conv_iter] = main(obj)
            flag = 0;
            result = zeros(1, max(1, ceil(obj.max_iter / 1)));
            resultIndex = 1;
            % result_len = max(1, ceil(obj.max_iter ./ 50));
            % result = cell(1, result_len);
            x = obj.x0;
            w = obj.A * x;
            r = obj.b - w;
            
            for k = 1:obj.max_iter
                if k == 1
                    v = r;
                else
                    v = r + s*v;
                end 

                u = obj.A*v;
                t = obj.stepsize(v, u, r);
                x = obj.update_x(x, t, v);
                w = obj.A*x; 
                r = obj.residual(w);

                criteria = norm(r, obj.norm_type)/norm(obj.b, obj.norm_type);

                if mod(k, 1) == 0
                    result(resultIndex) = criteria;
                end

                if criteria < obj.tol
                    flag = 1;
                    conv_iter = k;
                    break;
                end

                s = orthogonalize(obj, v, u, r);
                % dips('iteration: ', k, 'residual: ', norm(r, '2')/norm(obj.b, '2'));
                
            end

            x_min = x;
            if flag == 0
                conv_iter = obj.max_iter+1;
            end
        end
    end
end
