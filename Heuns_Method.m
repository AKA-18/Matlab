%%
%Heunâ€™s Method
%clearing all variable and workspaces
clc;
clear all;
clear workspace;
clear all variables;

dxdt =@(t,x,y,z) 10*(y-x);                                                 %declaring first function
dydt =@(t,x,y,z) x*(28-z)-y;                                               %declaring second function
dzdt =@(t,x,y,z) x*y-8/3*z                                                 %declaring third function% 
a = 0;                                                                     %starting value
tf= 50;                                                                    %inputing final iteration time
n = 5000;                                                                  %subintervals bewteeen starting point and ending point
x0 = -8;                                                                   %first starting x-coordinate
y0 = 5;                                                                    %first starting y-coordinate
z0 = 20;                                                                   %inputing first starting z-coordinate
delt=0.01;                                                                  %step size
dx=[];                                                                     %declaring empty array for dx
dy=[];                                                                     %declaring empty array for dy
dz=[];                                                                     %declaring empty array for dz
t=[];                                                                      %declaring empty array for t
npoints=(tf)/delt;                                                         %calculating npoints for counter
k=1;                                                                       %declaring k variable
b=tf;                                                                      %mapping b variable to k
h = (b-a)/n;                                                               %calculating h variable
t(k)=0;dx(k)=x0;dy(k)=y0;dz(k)=z0;                                         %mapping the starting variable
while(k <= npoints)                                                        %while loop counter                                             
     k1x=dxdt(t(k),dx(k),dy(k),dz(k));                                     %calculte the first x value
     k1y=dydt(t(k),dx(k),dy(k),dz(k));                                     %calculte the first y value
     k1z=dzdt(t(k),dx(k),dy(k),dz(k));                                     %calculte the first z value
     t(k+1)=t(k)+delt;                                                     %calculte the first t value
     k2x=dxdt(t(k+1),dx(k)+delt*k1x,dy(k)+delt*k1x,dz(k)+delt*k1x);        %calculating the predictor x value
     k2y=dydt(t(k+1),dx(k)+delt*k1y,dy(k)+delt*k1y,dz(k)+delt*k1y);        %calculating the predictor y value
     k2z=dzdt(t(k+1),dx(k)+delt*k1z,dy(k)+delt*k1z,dz(k)+delt*k1z);        %calculating the predictor z value
     dx(k+1)=dx(k)+(delt/2)*(k1x+k2x);                                     %calculating the corrected x value
     dy(k+1)=dy(k)+(delt/2)*(k1y+k2y);                                     %calculating the corrected y value 
     dz(k+1)=dz(k)+(delt/2)*(k1z+k2z);                                     %calculating the corrected z value
     k=k+1;                                                                %incrementing k value
 end
 
%curves  
figure(1)
plot3(dx,dy,dz)                                                            %plotting
title('Figure(1) for Lorentz Plot using Heuns Method')                     %giving titles
% 
%%
%ode45 method
w0 = [x0;y0;z0];                                                           %column vector of intial values                                                         %range of t
[t,y] = ode45(@ODEsystem,a:h:tf,w0);                                       %solving the system of equation with ode45
figure(2)
plot3(y(:,1),y(:,2),y(:,3))                                                %plotting
title('Figure(2) for Lorentz Plot using ode45 Method')                     %giving titles             
function dydt = ODEsystem(t,y)                                             %declaring the function for ODE System
dydt = zeros(3,1);                                                         %declaring the intial function
dydt(1) = 10*(y(2)-y(1));                                                  %declaring first function
dydt(2) = y(1)*(28-y(3))-y(2);                                             %declaring second function 
dydt(3) = y(1)*y(2)-8/3*y(3);                                              %decalring third function
end
