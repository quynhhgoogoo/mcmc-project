%% Markov chain Monte Carlo and Gibbs sampling
 
% 1) P(beta|sigma^2,x,y)
% beta=exp{-1/2*[beta^T*(X^T*X/sigma^2+Ssigma0^-1 )*beta - 2*beta^T(X^T*Y)/sigma^2 + Ssigma0^-1*beta0)]}
% N(beta|mu,Ssigma)
% mu=A^-1*b   = (X^T*X/sigma^2+Ssigma0^-1 )^-1*(X^T*Y)/sigma^2 + Ssigma0^-1*beta0)
% Ssigma=A^-1 = (X^T*X/sigma^2+Ssigma0^-1 )^-1


% 2) P(sigma^2|x,y,beta)
% sigma^2= IG[(sigma^2)a+n/2 , b+1/2(y-x*beta)^T(y-x*beta)]
% a*= (sigma^2)a+n/2
% b*= b+1/2(y-x*beta)^T(y-x*beta)


close all; 
clear; 
home; 
format long g;


% Initial data
% Parameters
nn      = 100;             % Number of samples for examine the AC
N       = 10000;           % Number of samples (iterations)
burnin  = 10000;           % Number of runs until the chain approaches stationarity
lag     = 10;              % Thinning or lag period: storing only every lag-th point

% Storage
theta   = zeros(1,N);      % Samples drawn from the Markov chain (States)
acc     = 0;               % Accepted samples
sigma = 1;                 % the standard deviation
mu = 1;


proposal_PDF = @(x,mu) normpdf(x,mu,sigma);   
sample_from_proposal_PDF = @(mu) normrnd(mu,sigma);
%p = @(x,y) exp(-x.^2).*(2+sin(x*5)+sin(x*2)); 


y = @(x,y) normpdf(x,mu,sigma)


aa = -3;   bb = 3;       
t = 0.5;

proposal_PDF = @(mu,Ssigma)
      mu= (X.^T.*X./sigma.^2+Ssigma0^.-1 )^-1.*(X^T.*Y)./sigma.^2 + Ssigma0.^-1.*beta0);
      Ssigma= X.^T.*X./sigma.^2+Ssigma0.^-1).^-1;
        
% M-H routine
for i = 1:burnin    % First make the burn-in stage
    [t] = MH_routine(t, y, proposal_PDF, sample_from_proposal_PDF);
end
for i = 1:N         % Cycle to the number of samples
     for j = 1:lag   % Cycle to make the thinning
        [t a] = MH_routine(t, y, proposal_PDF, sample_from_proposal_PDF);
    end
    theta(i) = t;        % Samples accepted
    acc      = acc + a;  % Accepted ?
end
accrate = acc/N;     % Acceptance rate


% Histogram, target function and samples 
xx = aa:0.01:bb;   % x-axis (Graphs)
figure;
        
% Histogram and target dist
[n1 x1] = hist(theta, ceil(sqrt(N))); 
bar(x1, n1/(N*(x1(2)-x1(1))));   colormap summer;   hold on;  % Normalized histogram
plot(xx, y(xx)/trapz(xx,y(xx)), 'r-', 'LineWidth', 2);        % Normalized "PDF"
xlim([aa bb]); grid on; 
title('Distribution of samples', 'FontSize', 15);
ylabel('Probability density function', 'FontSize', 12);
text(aa+3,0.8,sprintf('Acceptace rate = %g', accrate),'FontSize',12);

% samples
figure;
plot(1:N, theta,  'b-');   xlim([0 N]); ylim([aa bb]);   grid on; 
xlabel('Iterations, N', 'FontSize', 12); 
ylabel('Location', 'FontSize', 12);


%% Metropolis-Adjusted Langevin Algorithm
close all; 
clear; 
home; 
format long g;

% beta(i)= N(.|beta(i-1)+delta_t/2*M(-Ssigma^-1()beta(i-1)-mu),M*delta_t)
% M = I

