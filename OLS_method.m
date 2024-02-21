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

    coeff_store = cell(3, 1); % Corrected to store coefficients for each degree

    for n = 1:3
        coeff = main(n, data);
        disp(['degree: ', num2str(n)]);
        disp(coeff);
        coeff_store{n} = coeff; % Corrected to properly store coefficients
    end

    plotPolynomials(data, coeff_store); % Call plotting function with correct arguments
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

function plotPolynomials(data, coeff_store)
    x = data(:,1);
    y = data(:,2);
    
    figure; % Start a new figure
    scatter(x, y, 'filled'); % Plot original data points
    hold on;
    
    colors = ['r', 'g', 'b']; % Define colors for different degrees
    
    for i = 1:length(coeff_store)
        coeffs = coeff_store{i};
        poly = @(x) polyval(flip(coeffs), x); % Define polynomial using coeffs
        fplot(poly, [min(x), max(x)], colors(i)); % Plot polynomial
    end
    
    hold off;
    xlabel('x');
    ylabel('y');
    title('Polynomial Fits of Different Degrees');
    legend('Data', 'Degree 1', 'Degree 2', 'Degree 3', 'Location', 'best');
end