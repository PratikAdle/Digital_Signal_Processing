clc ;
clear ;
close all ;

Fs = 100 ;
B = 4 ;

f = -10 : 1/Fs : 10 ; 

% Rectangular Function is given by R
R = rectpuls(f,B) ;

figure('Name','Rectangular Function')
plot(f,R) , grid on ;
xlabel('Frequency') ;
ylabel('Magnitude') ;
axis([-10 10 0 1.5]) ;


%% GAUSSIAN FUNCTION

sigma = [0.3 0.2 0.1 0.05 0.0025] ;  

% Plotting Gaussian Function for Different Values of Sigma
% Gaussian Function is given by G
figure('Name','Gaussian Function for different Values of Sigma') ;

for i = 1:length(sigma)               
    
    numr = exp(-(f.^2)./(2*((sigma(i))^2)));
    denr = sigma(i)*((2*pi).^0.5);

    G = numr/denr ;
    
    plot(f,G) ;
    xlabel('Frequency') ;
    ylabel('Magnitude') ;
    grid on ; hold on ;
    axis([-2 2 0 10]) ;
    
end

legend(sprintf('sigma = %g',sigma(1)),sprintf('sigma = %g',sigma(2)),...
  sprintf('sigma = %g',sigma(3)),sprintf('sigma = %g',sigma(4)),...
  sprintf('sigma = %g',sigma(5))) 


%% DESIRED FREQUENCY RESPONSE for different Values of Sigma


% Plotting Desired Frequency Response for Different Values of Sigma
% D = convolution of R and G
figure('Name','Desired Frequency Response for Different Values of Sigma')

for i = 1:length(sigma)
    
    numr = exp(-(f.^2)./(2*((sigma(i))^2)));
    denr = sigma(i)*(sqrt(2*pi));

    G = numr./denr ;
     
    D = conv(R,G) ;

    f1 = -20:1/Fs:20 ;
    
    plot(f1/10,20*log10(D)) , grid on ;
    xlabel('Frequency') ;
    ylabel('Magnitude in dB') ;
    hold on ;
end

legend(sprintf('sigma = %g',sigma(1)),sprintf('sigma = %g',sigma(2)),...
  sprintf('sigma = %g',sigma(3)),sprintf('sigma = %g',sigma(4)),...
  sprintf('sigma = %g',sigma(5))) 


%% R,G and D in Time Domain

figure('Name','R , G , D in Time Domain')
t = -10 : 1/500 : 10 ;

% r is ifft of R
% r = sinc function in time domain

r = B.*sinc(B.*t) ;

subplot(3,1,1)
plot(t,r) ;
xlabel('Time') ;
ylabel('Amplitude') ;
title('Rectangular function r(t)') ;


subplot(3,1,2)
% g is ifft of R
% g = gaussian function in time domain
g = exp(-2*(pi*0.5.*t).^2) ;
plot(t,g) ;
xlabel('Time') ;
ylabel('Amplitude') ;
title('Gaussian function g(t)');


subplot(3,1,3)
d = r.*g ;
dc_t = B*sinc(B.*t).*exp(-(2*pi*0.5.*t).^2)  ;
plot(t,dc_t) ;
xlabel('Time') ;
ylabel('Amplitude') ;
title('Desired Impulse Response d(t)');


figure('Name','R , G , D in Time Domain')
plot(t,r,t,g,t,dc_t) , grid on ;
xlabel('Time') ;
ylabel('Amplitude') ;
axis([-5 5 -1 4]) ;
legend(sprintf('sinc'),sprintf('gaussian function'),...
  sprintf('impulse'))



%% Varying N (Order of FIR Filter) 

figure('Name','Frequency Response for sigma = 0.2 , B = 4 and diff values of N')
N = [20 40 60 80] ;

for j = 1:length(N)
    
    n = -N(j)/2 : N(j)/2 ;
    
    dc_n = 0.1.*B*sinc(B.*n.*0.1).*exp(-(2*pi*0.2.*n.*0.1).^2) ;
    [H, w] = freqz(dc_n) ;
    plot(20*log10(abs(H))) , grid on ;
    xlabel('Frequency') ;
    ylabel('Magnitude in dB') ;
    hold on ;
    
end

legend(sprintf('N = %g',N(1)),sprintf('N = %g',N(2)),...
  sprintf('N = %g',N(3)),sprintf('N = %g',N(4))) 

%%  Varying sigma

% Order = 80 and B = 4

M = 80 ;

if mod(M,2) == 0
    n = -M/2 : 1 : M/2 ;
else 
    n = -(M-1)/2 : (M-1)/2 ; 
end
   
dc_n = 0.1.*B*sinc(B.*n.*0.1).*exp(-(2*pi*0.5.*n.*0.1).^2) ;

figure('Name','Frequency Response for Order = 80 , B = 4 and diff values of sigma')

for i = 1:length(sigma)
    
    dc_n = 0.1.*B*sinc(B.*n.*0.1).*exp(-(2*pi*sigma(i).*n.*0.1).^2) ;
    [H, w] = freqz(dc_n) ;
    plot(20*log10(abs(H))) , grid on ;
    xlabel('Frequency') ;
    ylabel('Magnitude in dB') ;
    hold on ;
end

legend(sprintf('sigma = %g',sigma(1)),sprintf('sigma = %g',sigma(2)),...
  sprintf('sigma = %g',sigma(3)),sprintf('sigma = %g',sigma(4)),...
  sprintf('sigma = %g',sigma(5))) 


%% Varying N (Order of FIR Filter) 

figure('Name','Frequency Response for sigma = 0.0025 , B = 4 and diff values of N')
N = [1000 2000 3000 4000 5000 6500] ;

for j = 1:length(N)
    
    n = -N(j)/2 : N(j)/2 ;
    
    dc_n = 0.1.*B*sinc(B.*n.*0.1).*exp(-(2*pi*0.0025.*n.*0.1).^2) ;
    [H, w] = freqz(dc_n) ;
    plot(20*log10(abs(H))) ,grid on ;
    xlabel('Frequency') ;
    ylabel('Magnitude in dB') ;
    hold on ;
    
end

legend(sprintf('N = %g',N(1)),sprintf('N = %g',N(2)),...
  sprintf('N = %g',N(3)),sprintf('N = %g',N(4)),sprintf('N = %g',N(5)),...
  sprintf('N = %g',N(6))) 

%% sigma = 0.0225 ; order = 20 ;

sig = 0.0225 ;
ord = 20 ;

n = -ord/2:ord/2 ;
dc_n = 0.1.*B*sinc(B.*n.*0.1).*exp(-(2*pi*sig.*n.*0.1).^2) ;

figure('Name','Frequency Response for sigma = 0.0225 , N = 20')
[H, w] = freqz(dc_n) ;
plot(20*log10(abs(H))) ,grid on ;
xlabel('Frequency') ;
ylabel('Magnitude in dB') ;
