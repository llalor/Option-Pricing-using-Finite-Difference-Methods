%% Project
clc
clear

%% Input
S0=100;
K=110;
r=0.05;
T=3/12;
sigma=0.2;
M= 60;
N= 300;
dt=T/N; %T/N
%dt=1/1000;
xmin = log(S0)-2*sigma*sqrt(T);
xmax = log(S0)+2*sigma*sqrt(T);
dx=(xmax-xmin)/(M+1);
%dx=0.006557377049180; this is (xmax-xmin)/(M+1)
%dx=0.008;
%dx=0.0001;
matval = zeros(M+1,N+1); 
vetS = linspace(xmin,xmax,M+1)';
veti = 0:M;
vetj = 0:N;

% Setting up coefficients
a = ((0.5*(sigma^2/(dx^2)))-((r-.5*(sigma^2))/(2*dx)))*dt;
b = 1-(((sigma^2/(dx^2))+r)*dt);
c = ((0.5*(sigma^2/(dx^2)))+((r-.5*(sigma^2))/(2*dx)))*dt;

%% Run functions
[Call,Put] = blsprice(S0,K,r,T,sigma);
BlackScholesPut = Put;
EuPutVanilla = EuPutExpl1(S0,K,r,T,sigma,xmax,xmin,M,N,dx,dt);
AmPutVanilla = AmPutExpl1(S0,K,r,T,sigma,xmax,xmin,M,N,dx,dt);
EuPutExotic = EuPutExoExpl1(S0,K,r,T,sigma,xmax,xmin,M,N,dx,dt);
AmPutExotic = AmPutExoExpl1(S0,K,r,T,sigma,xmax,xmin,M,N,dx,dt);
table(BlackScholesPut, EuPutVanilla, AmPutVanilla, EuPutExotic, AmPutExotic)

%% 

%% Black Scholes
[Call,Put] = blsprice(S0,K,r,T,sigma);
BlackScholesPut = Put;
%% Explicit Method European Put

matval(:,N+1) = max(K-exp(vetS),0); %all the different maturity movements with respect to stock movements
matval(1,:) = (K-exp(xmin))*exp(-r*dt*(N-vetj)); %discounting it with respect to time
matval(M+1,:) = 0;%

for j=N:-1:1 %loop going backwards in time
   for i=2:M %moving forwards in terms of stock price change
      matval(i,j) = a*matval(i-1,j+1) + b*matval(i,j+1)+ ... %adding the difference
         c*matval(i+1,j+1);
   end
end
% return price, by linear interpolation outside the grid
EuPutVanillaPrice = interp1(exp(vetS), matval(:,1), S0);


figure1=figure();
plot(exp(vetS),matval(:,N+1),'x-');
hold on
title('European Put Explicit Method')
xlabel('Stock Price')
ylabel('Option Price')
hold off
grid on

%% Explicit Method American Put Option 

payoff = max(K-exp(vetS),0);
matval(:,N+1) = max(K-exp(vetS),0);
matval(1,:) = (K-exp(xmin))*exp(-r*dt*(N-vetj));
matval(M+1,:) = 0;

for j=N:-1:1 %loop going backwards in time
   for i=2:M %moving forwards in terms of stock price change
      matval(i,j) = max(a*matval(i-1,j+1) + b*matval(i,j+1)+ ... %adding the difference
         c*matval(i+1,j+1), payoff(i)); %K-exp(vetS(i))payoff(i)
   end
end
% return price, by linear interpolation outside the grid
AmPutVanillaPrice = interp1(exp(vetS),matval(:,1),S0);

figure2=figure();
plot(exp(vetS),matval(:,N+1),'+-');
hold on
title('American Put Explicit Method')
xlabel('Stock Price')
ylabel('Option Price')
grid on
hold off

%% Explicit method European Exotic Put Option

matval(:,N+1) = max((K-(100*(((exp(vetS))./S0).^(1./T)))),0); %1/0.25
matval(1,:) = ((K-(100*(((exp(xmin)./S0).^(1/T)))))).*exp(-r*dt*(N-vetj)); 
matval(M+1,:) = 0;%

for j=N:-1:1 %loop going backwards in time
   for i=2:M %moving forwards in terms of stock price change
      matval(i,j) = a*matval(i-1,j+1) + b*matval(i,j+1)+ ... %adding the difference
         c*matval(i+1,j+1);
   end
end
% return price, by linear interpolation outside the grid
EuPutExoticPrice = interp1(exp(vetS),matval(:,1),S0)

figure3=figure();
plot(exp(vetS), matval(:,end),'x-');
hold on
title('European Exotic Put Explicit Method')
xlabel('Stock Price')
ylabel('Option Price')
grid on
hold off


%% Explicit method American Exotic Put Option

payoff2 = max((K-(100*(((exp(vetS))./S0).^(1./T)))),0);
matval(:,N+1) = max((K-(100*(((exp(vetS))./S0).^(1./T)))),0); 
matval(1,:) = ((K-(100*(((exp(xmin)./S0).^(1/T)))))).*exp(-r*dt*(N-vetj)); 
matval(M+1,:) = 0;

for j=N:-1:1 %loop going backwards in time
   for i=2:M %moving forwards in terms of stock price change
      matval(i,j) = max(a*matval(i-1,j+1) + b*matval(i,j+1)+ ... %adding the difference
         c*matval(i+1,j+1),payoff2(i));
   end
end
% return price, by linear interpolation outside the grid
AmPutExoticPrice = interp1(exp(vetS),matval(:,1),S0);

figure4=figure();
plot(exp(vetS),matval(:,N+1),'x-');
hold on
title('American Exotic Put Explicit Method')
xlabel('Stock Price')
ylabel('Option Price')
grid on
hold off

%% Make Tables
table1 = table(BlackScholesPut, EuPutVanillaPrice, AmPutVanillaPrice, EuPutExoticPrice, AmPutExoticPrice)

%% Discretization steps: Option prices for changes in M 

for i=1:10
    M = i*10;
    dx = (xmax-xmin)/(M+1);
[EuPutVanillaPrice(i)] = EuPutExpl1(S0,K,r,T,sigma,xmax,xmin,M,N,dx,dt);
[AmPutVanillaPrice(i)] = AmPutExpl1(S0,K,r,T,sigma,xmax,xmin,M,N,dx,dt);
[EuPutExoticPrice(i)] = EuPutExoExpl1(S0,K,r,T,sigma,xmax,xmin,M,N,dx,dt);
[AmPutExoticPrice(i)] = AmPutExoExpl1(S0,K,r,T,sigma,xmax,xmin,M,N,dx,dt);   
        
    
end

M_Value=linspace(10,100,10)';
EuPutVanillaPrice = (EuPutVanillaPrice)';
AmPutVanillaPrice = (AmPutVanillaPrice)';
EuPutExoticPrice = (EuPutExoticPrice)';
AmPutExoticPrice = (AmPutExoticPrice)';
table(M_Value, EuPutVanillaPrice, AmPutVanillaPrice, EuPutExoticPrice, AmPutExoticPrice)

%% Discretization step: Option prices for changes in N

for i=1:10
    N = i*100;
    dt = T/N;
[EuPutVanillaPrice(i)] = EuPutExpl1(S0,K,r,T,sigma,xmax,xmin,M,N,dx,dt);
[AmPutVanillaPrice(i)] = AmPutExpl1(S0,K,r,T,sigma,xmax,xmin,M,N,dx,dt);
[EuPutExoticPrice(i)] = EuPutExoExpl1(S0,K,r,T,sigma,xmax,xmin,M,N,dx,dt);
[AmPutExoticPrice(i)] = AmPutExoExpl1(S0,K,r,T,sigma,xmax,xmin,M,N,dx,dt);   
        
    
end

N_Value=linspace(100,1000,10)';
EuPutVanillaPrice = (EuPutVanillaPrice)';
AmPutVanillaPrice = (AmPutVanillaPrice)';
EuPutExoticPrice = (EuPutExoticPrice)';
AmPutExoticPrice = (AmPutExoticPrice)';
table(N_Value, EuPutVanillaPrice, AmPutVanillaPrice, EuPutExoticPrice, AmPutExoticPrice)

