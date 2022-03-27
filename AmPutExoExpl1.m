function AmPutExoticPrice = AmPutExoExpl1(S0,K,r,T,sigma,xmax,xmin,M,N,dx,dt)

%dx=(xmax-xmin)/(M+1);
matval = zeros(M+1,N+1); 
vetS = linspace(xmin,xmax,M+1)';
veti = 0:M;
vetj = 0:N;

payoff2 = max((K-(100*(((exp(vetS))./S0).^(1./T)))),0);
matval(:,N+1) = max((K-(100*(((exp(vetS))./S0).^(1./T)))),0); 
matval(1,:) = ((K-(100*(((exp(xmin)./S0).^(1/T)))))).*exp(-r*dt*(N-vetj)); 
matval(M+1,:) = 0;

% Setting up coefficients
a = ((0.5*(sigma^2/(dx^2)))-((r-.5*(sigma^2))/(2*dx)))*dt;
b = 1-(((sigma^2/(dx^2))+r)*dt);
c = ((0.5*(sigma^2/(dx^2)))+((r-.5*(sigma^2))/(2*dx)))*dt;

for j=N:-1:1 %loop going backwards in time
   for i=2:M %moving forwards in terms of stock price change
      matval(i,j) = max(a*matval(i-1,j+1) + b*matval(i,j+1)+ ... %adding the difference
         c*matval(i+1,j+1),payoff2(i));
   end
end
% return price, by linear interpolation outside the grid
AmPutExoticPrice = interp1(exp(vetS),matval(:,1),S0)

%figure4=figure();
%plot(exp(vetS),matval(:,N+1),'x-');
%hold on
%title('American Exotic Put Explicit Method')
%xlabel('Stock Price')
%ylabel('Option Price')
%grid on
%hold off