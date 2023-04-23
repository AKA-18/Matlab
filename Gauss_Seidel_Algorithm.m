% clearing all variable and workspaces
clc;
clear all;
clear workspace;
clear all variables;

%declaring the base quations for performing the Jacobi Iteration
syms x y z
eqns = [4*x-y+z == 7,
        4*x-8*y+z == -21,
        -2*x+y+5*z == 15];
[A,b] = equationsToMatrix(eqns);

dia_domi=all((2*abs(diag(A)))- sum(abs(A),2)>=0); %checking for diagional dominance
TF = ismatrix(A) %checking if the input variable A is a Matrix
while TF == 1
    if dia_domi==1 
        disp('The Obtained Matrix A is Diagional Dominant Matrix')
        f=zeros(3,1); %guess values
        %f=[1;2;2]
        g1(1)=f(1);g2(1)=f(2);g3(1)=f(3)%guess values remaping
        tol=0.00000001; % declaring the tolernace

        m=inf;%declaring maximum number of iterations
        k=1; %declaring the starting number of iterations
        str_vle=k %locking starting values
        while k<=m
            g1(k+1)=(b(1)-(A(1,2)*g2(k))-(A(1,3)*g3(k)))/A(1,1);
            g2(k+1)=(b(2)-(A(2,1)*g1(k+1))-(A(2,3)*g3(k)))/A(2,2);
            g3(k+1)=(b(3)-(A(3,1)*g1(k+1))-(A(3,2)*g2(k+1)))/A(3,3);
            g=[g1;g2;g3];
            g
            rel_error = norm(g(:,k+1) - g(:,k));
            if rel_error<=tol %checking if the condition is statisified
                break %if true breaking out of the loop
            end
            k=k+1
            end_vle=k % locking ending values
    end
    % Logic for printing the final iteration and coordinate value where
    % converging is happening
    fprintf('Total iteration required to completed the convergence : %i ',end_vle-str_vle) 
    disp(' ')
    disp('The Converging Coordinates are :')
    fprintf(' %s ',g(:,k))
    else
       disp('The Obtained Matrix A is not Diagional Matrix, Kindly provide diagional dominant matrix')
       break
    end
    break
end

