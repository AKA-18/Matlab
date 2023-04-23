clc, clearvars, close all % to clear all of workspace,variables,screen etc

W(1)=6;                                                                    %declaring the X value.
Q(1)=2;                                                                    %declaring the Y value.
lam(1)=1;                                                                  %declaring the lambda values.

syms x y lambda                                                            %declaring the symbolic variables.
f = x.^2+6.*y-x*y;                                                         %declaring the function.                                                                             
g = 3.*x+4.*y-25==0;                                                       %declaring the equality constraint.                                                                                                                             
L = f - lambda * lhs(g);                                                   %declaring the langrange equation.                                                                                  
dL_dx = (diff(L,x) );                                                      %differentiating the langrange equation with respect to x.                                                                                                           
dL_dy = (diff(L,y) );                                                      %differentiating the langrange equation with respect to y.
dL_dlambda = (diff(L,lambda) );                                            %differentiating the langrange equation with respect to lambda.
grad_fdx = matlabFunction(dL_dx);                                          %converting the dL_dx to matlab function.    
grad_fdy = matlabFunction(dL_dy);                                          %converting the dL_dy to matlab function.
grad_fdlamda = matlabFunction(dL_dlambda);                                 %converting the dL_dlambda to matlab function.

L11= (diff(dL_dx,x));                                                      %calculating the L11 element of L matrix.                                                                                                       
L12= (diff(dL_dx,y));                                                      %calculating the L12 element of L matrix.
L13= (diff(dL_dx,lambda));                                                 %calculating the L13 element of L matrix.
L21= (diff(dL_dy,x));                                                      %calculating the L21 element of L matrix.
L22= (diff(dL_dy,y));                                                      %calculating the L22 element of L matrix.
L23= (diff(dL_dy,lambda));                                                 %calculating the L23 element of L matrix.
L31= (diff(dL_dlambda,x));                                                 %calculating the L31 element of L matrix.
L32= (diff(dL_dlambda,y));                                                 %calculating the L32 element of L matrix.
L33= (diff(dL_dlambda,lambda));                                            %calculating the L33 element of L matrix.
L=[L11 L12 L13;L21 L22 L23;L31 L32 L33];                                   %declaring the L matrix.

H=inv(L);                                                                  %inverse of L matrix.
tol = 0.1;                                                                 %tolernce level of iteration.
alpha=0.1;
for i = 1:500
    a=grad_fdx(W(i),Q(i),lam(i));                                          %calculating the value of function dL_dx for the values W,Q,lamda.
    b=grad_fdy(W(i),lam(i));                                               %calculating the value of function dL_dy for the values W,Q,lamda.
    c=grad_fdlamda(W(i),Q(i));                                             %calculating the value of function dL_dlambda for the values W,Q,lamda.
    del_l=[a; b; c];                                                       %mapping the above values into a single matrix
    I=[W(i),Q(i),lam(i)];                                                  %mapping the x,y,z                                                                                                    
    W(i+1) = I(1)'-alpha*H(1,:)*del_l;                                     %calculating the x value                                                                          
    Q(i+1) = I(2)'-alpha*H(2,:)*del_l;                                     %calculating the y value
    lam(i+1) = I(3)'-0.1*H(3,:)*del_l;                                     %calculating the z value
    i = i+1;                                                               %increment counter
    rel_error_W = norm(W(i) - W(i-1));                                      %calculating the differnce between present and previous value 
    if ((rel_error_W) <tol)                                                %checking if the condition is statisified
                break                                                      %if true breaking out of the loop
    end
end

% Result Table:`
Iter = 1:i;                                                                %iteration values
X_coordinate = W';                                                         %mapping x coordinate
Y_coordinate = Q';                                                         %mapping y coordinates                                                        
Iterations = Iter';                                                        %mapping iteration values
T = table(Iterations,X_coordinate,Y_coordinate);                           %creating tables
%Output:
fprintf('Initial Objective Function Value: %d\n\n',subs(f,[x,y], [W(1),Q(1)])); %printing the results
if ((rel_error_W) <tol)                                                    %logic testing for tolerence
    fprintf('Minimum succesfully obtained...\n\n');                        %printing
end
fprintf('Number of Iterations for Convergence: %d\n\n', i )                %printing maximum number of iteractions
fprintf('Point of Minima: [%d,%d]\n\n', W(i), Q(i));                       %point of minima
fprintf('Objective Function Minimum Value after Optimization: %f\n\n', subs(f,[x,y], [W(i),Q(i)])); % function value at the minima point
disp(T)                                                                    %display table
