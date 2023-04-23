clc
close all
clear all variables
 
f =rgb2gray(imread('car2.jpg'));                                            %loading the image
 
[M, N] = size(f);                                                          %calculating the rows and columns of image
wM        = zeros(M, M);                                                   %creating zero matrix with M size
wN        = zeros(N, N);                                                   %creating zero matrix with N size

%fourier transform of row data
for u = 0 : (M - 1)                                                        %for loop counter
    for x = 0 : (M - 1)                                                    %for loop counter
        wM(u+1, x+1) = exp(-2 * pi * 1i / M * x * u);                      %DFT expression                
    end                                                                    %end
end                                                                        %end
 
%fourier transform of column data
for v = 0 : (N - 1)                                                        %for loop counter
    for y = 0 : (N - 1)                                                    %for loop counter
        wN(y+1, v+1) = exp(-2 * pi * 1i / N * y * v);                      %DFT expression
    end                                                                    %end
end                                                                        %end

 
F = wM * im2double(f) * wN;                                                %mutliplying the original image values with the Fourier transformed value of wM and wN                                             


[rM, rN] = size(F);                                                        %calculating the rows and columns of image
wrM        = zeros(rM, rM);                                                %creating zero matrix with M size
wrN        = zeros(rN, rN);                                                %creating zero matrix with N size
 
for u = 0 : (rM - 1)                                                       %for loop counter
    for x = 0 : (rM - 1)                                                   %for loop counter
        wrM(u+1, x+1) = (1/rM)*exp(2 * pi * 1i* x * u / rM );              %IDFT expression 
    end                                                                    %end
end                                                                        %end
 
for v = 0 : (rN - 1)                                                       %for loop counter                                                      
    for y = 0 : (rN - 1)                                                   %for loop counter
        wrN(y+1, v+1) = (1/rN)*exp(2 * pi * 1i* y * v / rN );              %IDFT expression
    end                                                                    %end
end                                                                        %end

Fr = wrM * im2double(F) * wrN;                                             %recreating the reconstructed signal using DFT variable F with IDFT variable wrM and wrN
i = imread('car2.jpg');                                                     %reading the image                                        
x = double(rgb2gray(i));                                                   %converting the image to gray scale and double
N = length(x);                                                             %finding the length of the univariate data
X = zeros(N,1) ;                                                           %creating a new zero matrix with N rows and 1 columns

%%
%Plot one - Input Signal
t=0:N-1;                                                                   %declaring time variable                                                  
figure(1)                                                                  %figure(2)
subplot(511)                                                               %first subplot
stem(t,x);                                                                 %plotting time vs input univariate data                 
xlabel('Time (s)');                                                        %xlabel
ylabel('Amplitude');                                                       %ylabel
title('Input sequence')                                                    %title


%Plot Two - Magnitude response
subplot(512);                                                              %second subplot
stem(0:N-1,abs(F));                                                        %plotting time vs abs value fft transformed data using inbuilt function
xlabel('Frequency');                                                       %xlabel
ylabel('|X(k)|');                                                          %ylabel
title('Magnitude Response');                                               %title

%Plot Three - Phase response
subplot(513);                                                              %3rd subplot
stem(0:N-1,angle(F));                                                      %plotting time vs angle value fft transformed data using inbuilt function
xlabel('Frequency');                                                       %xlabel
ylabel('Phase');                                                           %ylabel
title('Phase Response');                                                   %title

%Plot Four - Reconstructed Signal
subplot(514)                                                               %subplot
stem(t,Fr)                                                                 %plotting the reconstructed signal
xlabel('Time');                                                            %xlabel
ylabel('Amplitude');                                                       %ylabel
title('Time domain - Reconstructed Signal')                                %title

Y = F;                                                                     %mapping output DFT signal to new variable                                                                                                                             
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

%% Construction and Reconstruction of Image using user build function

hFigure =imtool(i)                                                         %showing original image
set(hFigure,'NumberTitle','off','Name','Original Image')                   %setting figure parameters


F = wM * im2double(f) * wN;                                                %mutliplying the original image values with the Fourier transformed value of wM and wN                                             
hFigure =imtool(F)                                                         %showing the origfianl figure
set(hFigure,'NumberTitle','off','Name','Original Image Converted with user build DFT Function') %setting figure parameters

g= wrM * im2double(F) * wrN;                                               %recreating the reconstructed signal using DFT variable F with IDFT variable wrM and wrN
hFigure_1 =imtool(g)                                                       %showing the origfianl figure
set(hFigure_1,'NumberTitle','off','Name','Original Image with user build Inverse DFT Function')  %setting figure parameters

%% Plot using inbuild fft2 and ifft2 function

%Plot one - Input Signal
t=0:N-1;                                                                   %declaring time variable                                                  
figure(2)                                                                  %figure(2)
subplot(411)                                                               %first subplot
stem(t,x);                                                                 %plotting time vs input univariate data                 
xlabel('Time (s)');                                                        %xlabel
ylabel('Amplitude');                                                       %ylabel
title('Input sequence -In-built Function')                                  %title

%Plot Two - Magnitude response
subplot(412);                                                              %second subplot
stem(0:N-1,abs(fft2(x)));                                                  %plotting time vs abs value fft transformed data using inbuilt function
xlabel('Frequency');                                                       %xlabel
ylabel('|X(k)|');                                                          %ylabel
title('Magnitude Response In-built Function');                               %title

%Plot Three - Phase response
subplot(413);                                                              %3rd subplot
stem(0:N-1,angle(fft2(x)));                                                %plotting time vs angle value fft transformed data using inbuilt function
xlabel('Frequency');                                                       %xlabel
ylabel('Phase');                                                           %ylabel
title('Phase Response In-built Function');                                   %title

%Plot Four - Reconstructed Signal
rs = real(ifft2(fft2(x)));                                                 %reconstructing original signal using inbuilt function 
subplot(414)                                                               %subplot
stem(t,rs)                                                                 %plotting the reconstructed signal
xlabel('Time');                                                            %xlabel
ylabel('Amplitude');                                                       %ylabel
title('Time domain - Reconstructed Signal In-built Function')               %title
