
clc
clear
format long

syms X Y;                                                                  %declaring the symbolic variables.
f = X.^2+6.*Y-2.*X.*Y;                                                     %declaring the function. 

x(1) = 1;                                                                  %declaring the X value.
y(1) = 1;                                                                  %declaring the Y value.
UL=1000;                                                                   %upper limit of the counter
tol = 0.01;                                                                %Convergence Criteria
i = 1;                                                                     %Iteration Counter
% Gradient and Hessian Computation:
df_dx = diff(f, X);                                                        %differentiating f with respect to X
df_dy = diff(f, Y);                                                        %differentiating f with respect to Y
J = [subs(df_dx,[X,Y], [x(1),y(1)]) subs(df_dy, [X,Y], [x(1),y(1)])];      %creating J matrix by substituting the values 
ddf_ddx = diff(df_dx,X);                                                   %differentiating df_dx with X                                                  
ddf_ddy = diff(df_dy,Y);                                                   %differentiating df_dy with Y
ddf_dxdy = diff(df_dx,Y);                                                  %differentiating df_dx with Y 
ddf_ddx_1 = subs(ddf_ddx, [X,Y], [x(1),y(1)]);                             %substituting the X,Y values on ddf_ddx  
ddf_ddy_1 = subs(ddf_ddy, [X,Y], [x(1),y(1)]);                             %substituting the X,Y values on ddf_ddy
ddf_dxdy_1 = subs(ddf_dxdy, [X,Y], [x(1),y(1)]);                           %substituting the X,Y values on ddf_dxdy
H = [ddf_ddx_1, ddf_dxdy_1; ddf_dxdy_1, ddf_ddy_1];                        % Creating the a matrix using above variables     
S = inv(H);                                                                %calculating the inverse of S
alpha=0.1;                                                                 %declaring the alpha value
% Optimization Condition:
for i = 1:UL                                                               %loop counter
    I = [x(i),y(i)]';                                                      %creating transpose matrix of initial coordinates
    x(i+1) = I(1)-alpha*S(1,:)*J';                                         %calculating x coordinates
    y(i+1) = I(2)-alpha*S(2,:)*J';                                         %calculating y coordinates
    i = i+1;                                                               %increment counter
    J = [subs(df_dx,[X,Y], [x(i),y(i)]) subs(df_dy, [X,Y], [x(i),y(i)])];  %Updated J matrix
    ddf_ddx_1 = subs(ddf_ddx, [X,Y], [x(i),y(i)]);                         %updating ddf_ddx_1
    ddf_ddy_1 = subs(ddf_ddy, [X,Y], [x(i),y(i)]);                         %updating ddf_ddy_1
    ddf_dxdy_1 = subs(ddf_dxdy, [X,Y], [x(i),y(i)]);                       %updating ddf_dxdy_1
    rel_error_x = norm(x(i) - x(i-1)); %calculating the error of present and previous values 
    rel_error_y = norm(y(i) - y(i-1)); %calculating the error of present and previous values
    if ((rel_error_x && rel_error_y) <tol)                                                %checking if the condition is statisified
                break                                                      %if true breaking out of the loop
     end
    H = [ddf_ddx_1, ddf_dxdy_1; ddf_dxdy_1, ddf_ddy_1];                    %Updating H Matrix
    S = inv(H);                                                            %inverse of H Matrix
end
% Result Table:`
Iter = 1:i;                                                                %iteration values                                                             
X_coordinate = x';                                                         %mapping x coordinate 
Y_coordinate = y';                                                         %mapping y coordinates
Iterations = Iter';                                                        %mapping iteration values
T = table(Iterations,X_coordinate,Y_coordinate);                           %creating tables
%Output:
fprintf('Initial Objective Function Value: %d\n\n',subs(f,[X,Y], [x(1),y(1)]));  %printing the results
if (norm(J) < tol)                                                         %logic testing for tolerence
    fprintf('Minimum succesfully obtained...\n\n');                        %printing
end
fprintf('Number of Iterations for Convergence: %d\n\n', i);                %printing maximum number of iteractions
fprintf('Point of Minima: [%d,%d]\n\n', x(i), y(i));                       %point of minima
fprintf('Objective Function Minimum Value after Optimization: %f\n\n', subs(f,[X,Y], [x(i),y(i)])); % function value at the minima point
disp(T)                                                                    %display table
