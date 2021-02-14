function [x_range_optimum, max_R_squared] = linear_fit_finder_min_grad(x,y,n_points_fit_range,step_distance)

%  [x_range_optimum] = linear_fit_finder(x,y,n_points_fit_range,step_distance)
%
%finds optimum range in x values or steepest slope over which to take a linear fit
%   (when R^2 is greatest). requires the number of points to fit over and
%   how many points to move on by for each fit (defaults is 1, only
%   required to speed up function when number of points is high)



if ~exist('step_distance','var')
     % step_distance is not defined, defaulting to 1
      step_distance = 1;
else
 end

fit_range = n_points_fit_range-1;


min_P1 = inf; %minimum required value for the gradient
% min_R = -inf;
for int1 = 1:step_distance:(length(x)-fit_range)
    
    x_fit = x(int1:int1+fit_range);
    y_fit = y(int1:int1+fit_range);
    P1 = polyfit(x_fit,y_fit,1);
%     R_values = corrcoef(polyval(P1,x_fit),y_fit);
     R_values = corrcoef(x_fit,y_fit);
    
    R_squared = R_values(2).^2;
  
 
        if P1(1) <= min_P1 
            min_P1 = P1(1);
        max_R_squared = R_squared;
        x_range_optimum = int1:int1+fit_range;
    else
        end
        
    
    
    
    
end





end

