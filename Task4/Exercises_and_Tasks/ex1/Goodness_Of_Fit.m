function [correlation_coefficient,a,b,N,...
    regression_string,correlation_string,N_string] = Goodness_Of_Fit(x,y,type)

% Calculation of the goodness of fit

nanindex = find(~isnan(x)&~isnan(y));

N = length(nanindex);
x = x(nanindex);
y = y(nanindex);

N_string = sprintf('%s%d','N = ',N);

% Calculation of fit: y2 = ay1 + b
fit_parameters = polyfit(x,y,1);
a = fit_parameters(1);
b = fit_parameters(2);

regression_string=sprintf('%s%.2f%s%.2f','y = ',a,'x + ',b);

    switch type
        case 'correlation'
            % Coefficient of determination of a correlation
            sigma_squared_xy = 1/(N-2) * sum((y-x).^2);  
            sigma_squared_yy = 1/(N-1) * sum((y-mean(y)).^2);
            correlation_coefficient = 1 - sigma_squared_xy/sigma_squared_yy;
            correlation_string = sprintf('%s%.3f','R^2 = ',correlation_coefficient);
        case 'linear'
            % Coefficient of determination of a linear fit
            sigma_squared_xy = 1/(N-2) * sum((y-(a*x+b)).^2);  
            sigma_squared_yy = 1/(N-1) * sum((y-mean(y)).^2);
            correlation_coefficient = 1 - sigma_squared_xy/sigma_squared_yy;
            correlation_string = sprintf('%s%.3f','R^2 = ',correlation_coefficient);
        case 'Pearson'
            % Pearson product-moment correlation coefficient
            sigma_matrix = sqrt(nancov(x,y));
            correlation_coefficient = (sigma_matrix(1,2)*sigma_matrix(2,1))/...
                                      (sigma_matrix(1,1)*sigma_matrix(2,2));
            correlation_string = sprintf('%s%.3f','\rho_{xy} = ',correlation_coefficient);                         
    end 

end

