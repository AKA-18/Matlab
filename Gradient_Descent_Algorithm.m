% Gradient Descent to optimize the function f(x,y) = 2.*X.*Y+2*X-(X.^2)-2.*(Y.^2);
clc, clearvars, close all % to clear all of workspace,variables,screen etc

%declaring main function
f = @(x,y) x.^3+2.*(x-y).^2-3.*x;

%declaring symbolic main function
syms fg xg yg
fg = xg.^3+2.*(xg-yg).^2-3.*xg;

%creating the plot for performing gradient descent method
X = -2:0.01:3;
[X,Y] = meshgrid(X);
Z = f(X,Y);
surf(X,Y,Z,'FaceColor','c','FaceAlpha',0.3,'EdgeColor','none')
hold on
xlabel('x'),ylabel('y'),zlabel('z')

%finding partial derivatives (these are symbolic
%expression)
dfdx = diff(fg,xg);               
dfdy = diff(fg,yg);

%converting the symbolic expression to the format which can be used for
%logical purposes
grad_fdx = matlabFunction(dfdx);      
grad_fdy = matlabFunction(dfdy);

x(1) = 0; % initial value of x 
y(1) = 0; % initial value of y 
z(1) = f(x(1),y(1)); %initial value for z
alpha = 0.01;
 for i = 1:60
      x(i+1) = x(i) - alpha.*grad_fdx(x(i),y(i)); %x coordinate
      y(i+1)= y(i) - alpha.*grad_fdy(x(i),y(i)); %y coordinate
      z(i+1) = f(x(i+1),y(i+1));
     %z(i+1)=2.*x(i+1).*y(i+1)+2*x(i+1)-(x(i+1).^2)-2.*(y(i+1).^2);
 end
 plot3(x,y,z,'o','Markersize',3,'Color','red')
 hold off
 axis([min(x),max(x), min(y),max(y), min(z), max(z)]); 
 %view(150,10)
 title('Gradient Descent Optimization for f(x,y) = 2xy+2x-x^2-2*y^2 ')
 xlabel x-axis; ylabel y-axis; zlabel z-axis;
