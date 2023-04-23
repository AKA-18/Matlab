clc;
clear all;
clear variables;
%% Inferneces : we would like to compared to the inbuilt function efficiency is far greater than that of user defined function.
               %% There is however two downsides:
               %% The vectorization will consume a lot of memory (since a vector N*N have to be created)
               %% Hence it won't be faster than the built-in function.
%% Without using inbuilt function

% Performing DFT without using inbuilt function
x = [2 3 -1 4 5 8 9 2 8 9 6 3 4 1 2 5 23 6 5 8 7 1 2 6 4 8 9 3 2];         %univariate data (can be changed)
N = length(x);                                                             %finding the length of the univariate data
X = zeros(N,1);                                                            %creating a new zero matrix with N rows and 1 columns
for k = 0:N-1                                                              %for loop counter
    for n = 0:N-1                                                          %for loop counter
        X(k+1) = X(k+1) + x(n+1)*exp(-1i*2*pi*(k)*(n)/N );                 %DFT expression
    end                                                                    %end
end  

%Plot one - Input Signal
t = 0:N-1;                                                                 %declaring time variable
figure(1)                                                                  %declaring figure(1)
subplot(511)                                                               %first subplot
stem(t,x);                                                                 %plotting time vs input univariate data
xlabel('Time (s)');                                                        %xlabel
ylabel('Amplitude');                                                       %ylabel
title('Time domain - Input sequence - User Defined Function')              %title


%Plot Two - Magnitude response
subplot(512)                                                               %second subplot
stem(t,abs(X))                                                             %plotting time vs DFT transformed data
xlabel('Frequency');                                                       %xlabel
ylabel('|X(k)|');                                                          %ylabel
title('Frequency domain - Magnitude response - User Defined Function')     %title

%Plot Three - Phase response
subplot(513)                                                               %3rd subplot
stem(t,angle(X))                                                           %plotting time vs angle of DFT transformed data
xlabel('Frequency');                                                       %xlabel
ylabel('Phase');                                                           %ylabel
title('Frequency domain - Phase response - User Defined Function')         %title

% Reconstructing the Input Univariate Data from the Output Data using Inverse Discrete Fourier Transform
rs_X=zeros(1,N)
for k = 0:N-1                                                              %for loop counter
    for n = 0:N-1                                                          %for loop counter
        rs_X(k+1) = rs_X(k+1) + (1/N)*X(n+1)*exp(1i*2*pi*(k)*(n)/N)        %IDFT expression
    end                                                                    %end
end   

%Plot Four - Reconstructed Signal
subplot(514)                                                               %subplot
stem(t,rs_X)                                                               %plotting the reconstructed signal
xlabel('Time');                                                            %xlabel
ylabel('Amplitude');                                                       %ylabel
title('Time domain - Reconstructed Signal - User Defined Function')  

% Reconstructing the Input Univariate Data from the Output Data using
% Inverse Discrete Fourier Transform with lower number of transforms

Y = X;                                                                     %mapping output DFT signal to new variable                                                                                                                             
Y(ceil(N/2)+1:end) = 0;                                                    %ceiling the lower values to zero                                                 

% Perform the IDFT on the truncated DFT coefficients
Y_r = zeros(1,N);                                                          %creating empty matrix of size N
for k= 0:N-1                                                               %for loop counter
    for n = 0:N-1                                                          %for loop counter
        Y_r(k+1) = Y_r(k+1) + (2/N)*Y(n+1)*exp(1i*2*pi*k*n/N);             %IFFT with truncated coefficients
    end                                                                    %end
end                                                                        %end

%Plot Fifa - Reconstructed Signal using lesser transform count
subplot(515)                                                               %subplot
stem(t,Y_r)                                                                %plotting the reconstructed signal
xlabel('Time');                                                            %xlabel
ylabel('Amplitude');                                                       %ylabel
title('Time domain - Reconstructed Signal - User Defined Function using lesser transform count')  %title

%% DFT and IDFT for univariate data using inbuilt function

%Plot one - Input Signal
t=0:N-1;                                                                   %declaring time variable                                                  
figure(2)                                                                  %figure(2)
subplot(411)                                                               %first subplot
stem(t,x);                                                                 %plotting time vs input univariate data                 
xlabel('Time (s)');                                                        %xlabel
ylabel('Amplitude');                                                       %ylabel
title('Input sequence -Inbuilt Function')                                  %title

%Plot Two - Magnitude response
subplot(412);                                                              %second subplot
stem(0:N-1,abs(fft(x)));                                                   %plotting time vs abs value fft transformed data using inbuilt function
xlabel('Frequency');                                                       %xlabel
ylabel('|X(k)|');                                                          %ylabel
title('Magnitude Response-Inbuilt Function');                               %title

%Plot Three - Phase response
subplot(413);                                                              %3rd subplot
stem(0:N-1,angle(fft(x)));                                                 %plotting time vs angle value fft transformed data using inbuilt function
xlabel('Frequency');                                                       %xlabel
ylabel('Phase');                                                           %ylabel
title('Phase Response-Inbuilt Function');                                   %title

%Plot Four - Reconstructed Signal
rs = real(ifft(fft(x)));                                                   %reconstructing original signal using inbuilt function 
subplot(414)                                                               %subplot
stem(t,rs)                                                                 %plotting the reconstructed signal
xlabel('Time');                                                            %xlabel
ylabel('Amplitude');                                                       %ylabel
title('Time domain - Reconstructed Signal In-built Function')               %title
%%
