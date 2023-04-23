close all
clear all
clc
%% Code for finding the 2D steady state error using manual method
tic                                                                        %function used to calculate the time
n=16;                                                                      %declaring the grid size
timesim = 100;                                                             %declaring the simulation time   
T=zeros(n);                                                                %declaring the T variable
TT=900;                                                                    %declaring boundary temperature
TL=200;                                                                    %declaring boundary temperature
TR=200;                                                                    %declaring boundary temperature
TB=200;                                                                    %declaring boundary temperature
T(1,1:n)=TT ;                                                              %top side temperature
T(n,1:n)=TB ;                                                              %bottom side temperature
T(1:n,1)=TL ;                                                              %left side temperature
T(1:n,n)=TR ;                                                              %right side temperature
T(1,1)=(TL+TT)/2;                                                          %top left cell value                                                          
T(n,1)=(TR+TB)/2;                                                          %bottom left cell value
T(1,n)=(TR+TT)/2;                                                          %top right cell value
T(n,n)=(TR+TB)/2;                                                          %bottom right cell value                                                                 %
iteration = 0;                                                             %counter variable one
figure(1)                                                                  %figure one
while iteration < timesim                                                  %while loop counter
    if (mod(iteration,10)==0)                                              %if condition 
        imagesc(T);                                                        %showing the image of T(the 16*16 grid)
        colorbar;                                                          %showing color bar
        title('Temperature (Steady State) by using manual function'),xlabel('width'),ylabel('height'); %labeling, titles, etc
        drawnow;                                                           %plotting
    end             
    iteration = iteration+1;                                               %increasing the counter value
    Told=T;                                                                %mapping the T value to Told
        for i=2:n-1;                                                       %loop counter 1                                                     
            for j=2:n-1;                                                   %loop counter 2
                T(i,j) = .25 *(T(i,j-1) + T(i-1,j) + T(i+1,j) + T(i,j+1)); %laplace equation for the temperature at various point within the 16*16 grid by iteractive method
            end
        end
end
fprintf('The time taken to complete the operation using manual logic ')    %printing the time taken for the execution
toc
%%
tic                                                                        %function used to calculate the time
n=16;                                                                      %declaring the grid size
TD=zeros(n);                                                               %declaring the TD variable
TTD=900;                                                                   %declaring boundary temperature
TLD=200;                                                                   %declaring boundary temperature
TRD=200;                                                                   %declaring boundary temperature
TBD=200;                                                                   %declaring boundary temperature
TD(1,1:n)=TTD ;                                                            %top side temperature
TD(n,1:n)=TBD ;                                                            %bottom side temperature
TD(1:n,1)=TLD ;                                                            %left side temperature
TD(1:n,n)=TRD ;                                                            %right side temperature
TD(1,1)=(TLD+TTD)/2;                                                       %top left cell value 
TD(n,1)=(TRD+TBD)/2;                                                       %bottom left side temperature
TD(1,n)=(TRD+TTD)/2;                                                       %top right cell value
TD(n,n)=(TRD+TBD)/2;                                                       %bottom right cell value
iteration=0;                                                               %counter variable
figure(2)                                                                  %figure two
while iteration < timesim                                                  %while loop counter
    if (mod(iteration,10)==0)                                              %if condition  
        imagesc(TD);                                                       %showing the image of TD(the 16*16 grid)
        colorbar;                                                          %showing color bar
        title('Temperature (Steady State) plot with using del2 function'),xlabel('width'),ylabel('height');   %labeling, titles, etc
        drawnow;                                                           %plotting
    end
    iteration = iteration+1;                                               %increasing the counter value
    TDold=TD;                                                              %mapping the TD value to Told
    Lu = del2(TD);                                                         %Laplacian of TD
    TD(2:n-1,2:n-1) = TD(2:n-1,2:n-1)+Lu(2:n-1,2:n-1);                     %laplace equation for the temperature at various point within the 16*16 grid by using del2 operation
end
fprintf('The time taken to complete the operation using del2 function ')   %printing the time taken for the execution
toc
