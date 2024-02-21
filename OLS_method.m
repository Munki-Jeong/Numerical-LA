function driver
    % Data points
    data = [
        0,    1.0;
        0.15, 1.004;
        0.31, 1.031;
        0.5,  1.117;
        0.6,  1.223;
        0.75, 1.422
    ];

    coeff_store = cell(3, 1); % store coefficients for each degree
    errors = zeros(3, 1); % Store errors for each degree

    for n = 1:3
        coeff = main(n, data);
        disp(['degree: ', num2str(n)]);
        disp(coeff);
        coeff_store{n} = coeff; 
        errors(n) = calculateError(coeff, data);
    end

    plotPolynomials(data, coeff_store); % Call plotting function
    disp('Errors for polynomial degrees 1 to 3:');
    for n = 1:length(errors)
        disp(['Error for degree ', num2str(n), ': ', num2str(errors(n))]);
    end
end

function coeff = main(degree, data)
    % Calculate LHS and RHS matrices
    LHS_matrix = Cal_LHS(degree, data(:,1));
    RHS_matrix = Cal_RHS(degree, data);

    % Solve the normal equations
    coeff = LHS_matrix \ RHS_matrix;
end

function LHS_matrix = Cal_LHS(degree, x_data)
    n = degree;
    LHS_matrix = zeros(n+1, n+1);

    for row = 1:n+1
        for col = 1:n+1
            LHS_matrix(row, col) = sum(x_data .^ (row+col-2));
        end
    end
end

function RHS_matrix = Cal_RHS(degree, data)
    n = degree;
    RHS_matrix = zeros(n+1, 1);

    for iter = 0:n
        RHS_matrix(iter+1) = sum(data(:,2) .* (data(:,1) .^ iter));
    end
end


function error = calculateError(coeff, data)
    x = data(:,1);
    y = data(:,2);
    n = length(coeff) - 1; 
    poly_vals = polyval(flip(coeff), x);
    error = sum((y - poly_vals).^2); 
end

function plotPolynomials(data, coeff_store)
    x = data(:,1);
    y = data(:,2);
    
    figure; 
    scatter(x, y, 'filled'); 
    hold on;
    
    colors = ['r', 'g', 'b']; 
    
    for i = 1:length(coeff_store)
        coeffs = coeff_store{i};
        poly = @(x) polyval(flip(coeffs), x); 
        fplot(poly, [min(x), max(x)], colors(i)); 
    end
    
    hold off;
    xlabel('x');
    ylabel('y');
    title('Polynomial Fits of Different Degrees');
    legend('Data', 'Degree 1', 'Degree 2', 'Degree 3', 'Location', 'best');
end