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


y = @(x,y) normpdf(x,mu,sigma);


aa = -3;   bb = 3;       
t = 0.5;

% Note: In matlab to compute inverse matrix use the inv() command. (.^-1 is elementwise inverse) 
% To compute something like c = inv(A) * b effienctly, use c = linsolve(A, b)
% proposal_PDF = @(mu,Ssigma)
%      mu= (X.^T.*X./sigma.^2+Ssigma0^.-1 )^-1.*(X^T.*Y)./sigma.^2 + Ssigma0.^-1.*beta0);
%      Ssigma= X.^T.*X./sigma.^2+Ssigma0.^-1).^-1;
        
% M-H routine
for i = 1:burnin    % First make the burn-in stage
    %[t] = MH_routine(t, y, proposal_PDF, sample_from_proposal_PDF);
end
for i = 1:N         % Cycle to the number of samples
     for j = 1:lag   % Cycle to make the thinning
        %[t a] = MH_routine(t, y, proposal_PDF, sample_from_proposal_PDF);
    end
    %theta(i) = t;        % Samples accepted
    %acc      = acc + a;  % Accepted ?
end
%accrate = acc/N;     % Acceptance rate


% Histogram, target function and samples 
xx = aa:0.01:bb;   % x-axis (Graphs)
figure;
        
% Histogram and target dist
%[n1 x1] = hist(theta, ceil(sqrt(N))); 
%bar(x1, n1/(N*(x1(2)-x1(1))));   colormap summer;   hold on;  % Normalized histogram
%plot(xx, y(xx)/trapz(xx,y(xx)), 'r-', 'LineWidth', 2);        % Normalized "PDF"
%xlim([aa bb]); grid on; 
%title('Distribution of samples', 'FontSize', 15);
%ylabel('Probability density function', 'FontSize', 12);
%text(aa+3,0.8,sprintf('Acceptace rate = %g', accrate),'FontSize',12);

% samples
%figure;
%plot(1:N, theta,  'b-');   xlim([0 N]); ylim([aa bb]);   grid on; 
%xlabel('Iterations, N', 'FontSize', 12); 
%ylabel('Location', 'FontSize', 12);


%% Metropolis-Adjusted Langevin Algorithm
close all; 
clear; 
home; 
format long g;

% beta(i)= N(.|beta(i-1)+delta_t/2*M(-Ssigma^-1()beta(i-1)-mu),M*delta_t)
% M = I

% The first task is to simulate some data from the linear regression model

N_data = 100; % Number of datapoints in the dataset
beta_true = [0.5 1.5 -1.0 0.7]'; % The ground truth intercept and weights
sigma_true = 0.5; % Ground truth value for the standard deviation
d = length(beta_true); % Number of model parameters (excluding variance)
designs = rand(N_data, d-1); % Sample covariates from U(0, 1)
X = cat(2, ones(N_data, 1), designs); % Construct design matrix
y = normrnd(X * beta_true, sigma_true); % Sample outputs from our model

scatter(X(:, 2), y); % Visualize first convariate and outputs

% We also need to select (noninformative) prior distributions

% The prior for the variance is inverse gamma distribution with parameters 
% a0 and b0
a0 = 0.0001;
b0 = 0.0001;

% The prior for the weights and intercept is a multivariate Gaussian
% distribution with parameters mu_0 and Ssigma_0)
mu_0 = zeros(d, 1); % Prior mean for weights and intercept
sigma2_0 = 1000; % Prior variance for weights and intercept
Ssigma0 = eye(d) * sigma2_0; % Prior covariance matrix for weights and intercept

% MCMC parameters
nn      = 100;             % Number of samples for examine the AC
N       = 10000;           % Number of samples (iterations)
burnin  = 10000;           % Number of runs until the chain approaches stationarity
lag     = 10;              % Thinning or lag period: storing only every lag-th point

% Storage
samples   = zeros(N, d+1);      % Samples drawn from the Markov chain (States)
acc     = 0;                    % Accepted samples

% Initialize chain
sigma2 = 1;                 % variance parameter (sampled from inv. gamma)
beta = zeros(1, d);       % weight parameters (sampled with MALA)

% MALA parameters
delta_t = 0.1;
M = eye(d);

% Define helper functions

% proposal density is p(beta(i + 1) | beta(i), ...)
proposal_PDF_mala = @(beta_proposed, beta_current, mu, Ssigma) 0; % TODO

sample_from_proposal_PDF_mala = @(beta_current, mu, Ssigma) 0; % TODO

% target density is p(beta|sigma2, X, y)
target_density_mala = @(beta, mu, Ssigma) 0; % TODO

% Run the chain
for i = 1:burnin    % First make the burn-in stage
    
    % Compute a parameter of inverse gamma distribution
    a = a0 + N_data;
    
    % Compute b parameter of inverse gamma distribution
    b = b0 + 0.5 * (y - X * beta')' * (y - X * beta');
    
    % Sample sigma2 from p(sigma2 | beta, X, y) (inverse gamma)
    sigma2 = 1 / gamrnd(a, 1/b);
    
    % Compute mean of p(beta | sigma2, X, y)
    mu = 0; % TODO
    
    % Compute covariance matrix of p(beta | sigma2, X, y)
    Ssigma = 0; % TODO
    
    % Sample proposed beta from the MALA proposal distribution
    beta_proposed = 0; % TODO
    
    % Evaluate p(beta_proposed | sigma2, X, y)
    proposed_target_density = 0; % TODO
    
    % Evaluate p(beta_current | sigma2, X, y)
    current_target_density = 0; % TODO
    
    % Evaluate p(beta_proposed | beta_current, ...)
    proposal_density = 0; % TODO
    
    % Evaluate p(beta_current | beta_proposed, ...)
    reverse_proposal_density = 0; % TODO
    
    % Compute acceptance probability with MH acceptance rule
    alpha = 0; % TODO
    
    % Replace beta with proposal with probability alpha
    beta = beta; % TODO
end
for i = 1:N         % Cycle to the number of samples
     for j = 1:lag   % Cycle to make the thinning
         
        % Compute a parameter of inverse gamma distribution
        a = a0 + N_data;
        
        % Compute b parameter of inverse gamma distribution
        b = b0 + 0.5 * (y - X * beta')' * (y - X * beta');

        % Sample sigma2 from p(sigma2 | beta, X, y) (inverse gamma)
        sigma2 = 1 / gamrnd(a, 1/b);
        
        % Compute mean of p(beta | sigma2, X, y)
        mu = 0; % TODO

        % Compute covariance matrix of p(beta | sigma2, X, y)
        Ssigma = 0; % TODO

        % Sample proposed beta from the MALA proposal distribution
        beta_proposed = 0; % TODO

        % Evaluate p(beta_proposed | sigma2, X, y)
        proposed_target_density = 0; % TODO

        % Evaluate p(beta_current | sigma2, X, y)
        current_target_density = 0; % TODO

        % Evaluate p(beta_proposed | beta_current, ...)
        proposal_density = 0; % TODO

        % Evaluate p(beta_current | beta_proposed, ...)
        reverse_proposal_density = 0; % TODO

        % Compute acceptance probability with MH acceptance rule
        alpha = 0; % TODO

        % Replace beta with proposal with probability alpha
        beta = beta; % TODO
        
     end
    samples(i, 1) = sigma2;          % Store sigma2 sample
    samples(i, 2:d+1) = beta;        % Store beta sample
    acc      = acc + a;  % Accepted ?
end
accrate = acc/N;     % Acceptance rate

% TODO: Make parameter specific histograms and trace plots for sigma2,
% beta(1) (intercept) and beta(2) (first weight).

% TODO: Compute posterior mean (mean of beta samples)

% TODO: Simulate test data from the same model and compute MSE of the fit
% (i.e. measure average distance between predictions based on 
% beta mean to each y value in the test set).


