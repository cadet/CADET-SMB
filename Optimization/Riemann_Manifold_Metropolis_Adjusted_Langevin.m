function Riemann_Manifold_Metropolis_Adjusted_Langevin(params)

% =============================================================================
% Metropolis Adjusted Langevin algorithm under Riemann Manifold using Parallel tempering (PRML)

% PRML optimizes a problem by ...

%  Parameter:
%        - params. It is the specified parameters from the main function.
%        And let the main function informed that which parameters are need
%        to be optimized
% =============================================================================
    

    startTime = clock;

%   Get the sampler options
    opt = getOptions_PRML(params);

%   Initialization of the chains
    [states, MetricTensor] = initChains( (1 ./ opt.temperature), opt);       
    
%   Pre-allocation    
    chain       = zeros(opt.Nchain, opt.nCol+1, opt.nsamples);
    Population  = cell(opt.Nchain,1);
    sigmaChain  = zeros(opt.nsamples, opt.Nchain);
    Temperatures = zeros(opt.nsamples, opt.Nchain);  
    
    Temperatures(:, 1) = ones(opt.nsamples, 1);
    Temperatures(1, :) = opt.temperature;
    
%     sigmaSqu  = states(:,opt.nCol+1)/ (opt.nObserv - opt.nCol);
%     sigmaSqu0 = sigmaSqu; n0 = 1;

    sigmaSqu = 0.1;
%-----------------------------------------------------------------------------------------
%   The main loop     
    for i = 1:opt.nsamples

        beta = 1 ./ Temperatures(i, :);
        
        fprintf('Iter: %5d    done: %5d%%    accepted: %5d%% \n', i,...
            fix(i/opt.nsamples*100), fix(opt.accepted/i*100));
        
%       sampler
        [states, MetricTensor, opt] = PRML_sampler(states, MetricTensor, sigmaSqu, beta, opt);

%       SWAPS: at predetermined intervals
        if  mod(i, opt.Swapint) == 0
            
            a = ceil(rand*opt.Nchain);%  a = randi([1 opt.Nchain]);

            if a == opt.Nchain
                b = 1;
            else
                b = a+1;
            end
	
            SSa = states(a, opt.nCol+1);
            SSb = states(b, opt.nCol+1);	

%             rho = ( exp(-0.5 * (SSa-SSb)/sigmaSqu(a)) ) ^ (beta(b)-beta(a
            rho = ( exp(-0.5 * (SSa-SSb)/sigmaSqu) ) ^ (beta(b)-beta(a));

            if rand < min(rho, 1)
                temp         = states(a, :);
                states(a, :) = states(b, :);
                states(b, :) = temp;
                clear temp;
            end
            
        end
         
%       Append the newest states to the end of the 3-D matrix
        chain(:,:,i) = states;

        
%       Convergence diagnostics
        if mod(i, opt.convergint) == 0 && i > opt.convergint
            
            R = GelmanR_statistic(i, chain(:,:,1:i), opt);
            
            if all(R < 1.1) || i == opt.nsamples
                opt.nsamples = i;
                break
            end
            
        end
        
        
%       Variance of error distribution (sigma) was treated as a parameter to be estimated.
%         for k = 1:opt.Nchain
%             sigmaSqu(k)  = 1 ./ GammarDistrib(1,1,(n0+opt.nObserv)/2,2 ./(n0 .* sigmaSqu0(k) + chain(k,opt.nCol+1,i)));
%             sigmaChain(i,k) = sigmaSqu(k)';
%         end
        
%       Temperature dynamics       
        Temperatures = temperatureDynamics(i, states, Temperatures, sigmaSqu, opt); 
        
        
    end  % for i = 1:opt.nsamples
% ----------------------------------------------------------------------------------------    


%   Preparation for the plotting
    for kk = 1: opt.Nchain
        for ii = 1:opt.nCol+1
            for jj = 1:opt.nsamples
                Population{kk}.samples(jj,ii) = chain(kk, ii, jj);
            end
        end
        Population{kk}.beta = beta(kk);
    end
    
    id = floor(0.3*opt.nsamples);
    plotSamples = zeros(id, opt.nCol+1);    
    
    for k = 1:id
        for j = 1: opt.nCol+1
            plotSamples(k, j) = Population{1}.samples(k+(opt.nsamples-id), j);
        end
    end
    

%   Gather some useful information and store them
    [yValue, row]         = min(Population{1}.samples(:,opt.nCol+1));
    xValue                = Population{1}.samples(row,1:opt.nCol);
    
    result.NumChain       = opt.Nchain;
    result.optTime        = etime(clock,startTime)/3600;
    result.convergence    = R;
    result.accepted       = fix(opt.accepted/opt.nsamples*100);
    
    result.xValue         = xValue;
    result.yValue         = yValue;
    result.chain          = chain;
    result.index          = opt.nsamples;
    result.population     = Population;
    result.plotSamples    = plotSamples;
%     result.sigmaChain     = sigmaChain; 
    result.Temperatures   = Temperatures(1:opt.nsamples, :);
    
    save(sprintf('result_%3d.mat', fix(rand*1000)), 'result');
    fprintf('Markov chain simulation finished and the result saved as result.mat \n');
    
end   


function opt = getOptions_PRML(params)

% -----------------------------------------------------------------------------
%  The parameters for the Optimizer 

%  Parameter:
%       - params. It is the specified parameters from the main function.
%        And let the main function informed that which parameters are need
%        to be optimized.

%  Return:
%       - opt.
%           + Nchain. The number of the candidates (particles)
%           + dimension. The number of the optimized parameters
%           + bounds. The boundary limitation of the parameters
%           + nsamples. The maximum number of algorithm's iteration
%           + cr_Probability. The cross-over probability in DE's formula
%           + weight. The weight coefficients in DE's formula
%           + strategy. There are 1,2,3,4,5,6 strategies in DE algorithm to 
%               deal with the cross-over and mutation. As for detais, we
%               will refer your the original paper of Storn and Price
% -----------------------------------------------------------------------------


    opt = [];
    
    opt.temperature       = [1, 10, 100];
    opt.Nchain            = length(opt.temperature);
    opt.nCol              = length(fieldnames(params)); 
    opt.bounds            = log ( [0.20    0.30;
                                   130     230;
                                   9.0e-7  10e-7;
                                   0.9e-7  1.0e-7;
                                   1.0e-7  2.0e-7;
                                   1.0e-7  2.0e-7] )';
    opt.nsamples          = 500;
  
        
%   Check out the dimension of the set of parameters, and the boundary limitation
    [row, col] = size(opt.nCol);
    if row > 1 || col > 1true
        error('The initialized dimension of the set of parameters might be wrong');
    end
    
    [row, col] = size(opt.bounds);
    if col ~= opt.nCol && row ~= 2
        error('Please check your setup of the range of parameters');
    end
 
%   Those are the specific parameters for the Markov Chain Simulation
    opt.convergint        = 50;
    opt.burnin            = 500;
    opt.printint          = 100;
    opt.Swapint           = 100;
    opt.iterMax           = 1000000;
    opt.nObserv           = 1000;
    

    opt.epsilon           = 0.15;
    opt.iterationNum      = 0;
    opt.accepted          = 0;
    opt.continueIteration = true;
    opt.converged         = false;
    opt.saveMetricTensor  = true;

    if fix(opt.nsamples/opt.convergint) ~= opt.nsamples/opt.convergint
        error('please let the number of samples devided by the index of convergence diagnostics');
    end

end

function [initChain, MetricTensor] = initChains(beta, opt)

% -----------------------------------------------------------------------------
% The initilization of the chains

%  Parameter:
%       - opt.
%           + PopulSize. The number of the candidates (particles)
%           + IndivSize. The number of the optimized parameters
%           + IndivScope. The boundary limitation of the parameters
%           + loopCount. The maximum number of algorithm's iteration
%           + cr_Probability. The cross-over probability in DE's formula
%           + weight. The weight coefficients in DE's formula
%           + strategy. There are 1,2,3,4,5,6 strategies in DE algorithm to 
%               deal with the cross-over and mutation. As for detais, we
%               will refer your the original paper of Storn and Price

% Return:
%       - Population. The population of the particles, correspondingly the objective value
%       - OptPopul. The best fit found among the population.
% -----------------------------------------------------------------------------


    initChain = rand(opt.Nchain, opt.nCol+1);
    MetricTensor = cell(1, opt.Nchain);
    
    for k = 1: opt.nCol
        initChain(:,k) = opt.bounds(1,k) + initChain(:,k) .* (opt.bounds(2,k) - opt.bounds(1,k));
    end
 
    
    for j = 1: opt.Nchain
        
        [Res, Jac] = feval( @simulatedMovingBed, exp( initChain(j, 1:opt.nCol) ) );
        
        initChain(j,opt.nCol+1) = Res' * Res;
        Sigma2 = (initChain(j, opt.nCol+1) / (opt.nObserv - opt.nCol));
        
        MetricTensor{j}.G = beta(j) .* (Jac' * (1/Sigma2) * Jac);
        MetricTensor{j}.invG = inv(MetricTensor{j}.G + eye(opt.nCol)*1e-10);
        MetricTensor{j}.GradL = -Jac' * Res / Sigma2;
        
        [R, p] = chol( MetricTensor{j}.invG );
        if p == 0
            MetricTensor{j}.sqrtInvG = R;
        else
            [U,S,V] = svd( MetricTensor{j}.invG );
            S = diag(S); S(S<1e-10) = 1e-10;
            MetricTensor{j}.sqrtInvG = triu( U * diag(sqrt(S)) * V' );
        end  
        
    end
    
end

function [states, MetricTensor, opt] = PRML_sampler(states, MetricTensor, sigmaSqu, beta, opt)

%-----------------------------------------------------------------------------------------
%
%-----------------------------------------------------------------------------------------

   
    for j = 1: opt.Nchain        
        
%         proposal = states(j,1:opt.nCol) + 0.5 * opt.epsilon^2 * (MetricTensor{j}.G \...
%             MetricTensor{j}.GradL)'+ opt.epsilon * randn(1,opt.nCol) * chol(MetricTensor{j}.invG);
                
        proposal = states(j,1:opt.nCol) + 0.5 * opt.epsilon^2 * (MetricTensor{j}.G \...
            MetricTensor{j}.GradL)'+ opt.epsilon * randn(1,opt.nCol) * MetricTensor{j}.sqrtInvG;
            

        
        [newRes, newJac] = feval( @simulatedMovingBed, exp( proposal) );
        
        
        newSS = newRes' * newRes;               
        SS    = states(j,opt.nCol+1);
        
%         rho = (exp( -0.5*(newSS - SS) / sigmaSqu(j)))^beta(j);
        rho = (exp( -0.5*(newSS - SS) / sigmaSqu))^beta(j);
        
        if rand <= min(1, rho)
            
            states(j, 1:opt.nCol) = proposal;
            states(j, opt.nCol+1) = newSS;
            
%             MetricTensor{j}.G = beta(j) .* (newJac' * (1/sigmaSqu(j)) * newJac);
            MetricTensor{j}.G = beta(j) .* (newJac' * (1/sigmaSqu) * newJac);
            MetricTensor{j}.invG = inv(MetricTensor{j}.G + eye(opt.nCol)*1e-10);
%             MetricTensor{j}.GradL = -newJac' * newRes / sigmaSqu(j);
            MetricTensor{j}.GradL = -newJac' * newRes / sigmaSqu;
            
            [R, p] = chol( MetricTensor{j}.invG );
            if p == 0
                MetricTensor{j}.sqrtInvG = R;
            else
                [U,S,V] = svd( MetricTensor{j}.invG );
                S = diag(S); S(S<1e-10) = 1e-10;
                MetricTensor{j}.sqrtInvG = triu( U * diag(sqrt(S)) * V' );
            end
            
            if j == 1
                opt.accepted = opt.accepted + 1;
            end
            
        end
    end
    
end

function Temperatures = temperatureDynamics(t, states, Temperatures, sigmaSqu, opt)

%-----------------------------------------------------------------------------------------
%
%-----------------------------------------------------------------------------------------


    t0 = 1e3; nu = 100;
    beta = 1 ./ Temperatures(t, :);

    for i = 2: opt.Nchain
        
        b = i - 1;
        if i == opt.Nchain
            c = 1;
        else
            c = i + 1;
        end
	
        SSa = states(i, opt.nCol+1);
        SSb = states(b, opt.nCol+1);
        SSc = states(c, opt.nCol+1);

%         rho_ab = min(1, (exp(-0.5*(SSa-SSb)/sigmaSqu(i)))^(beta(b)-beta(i)) );
%         rho_ca = min(1, (exp(-0.5*(SSc-SSa)/sigmaSqu(c)))^(beta(i)-beta(c)) );
        rho_ab = min(1, (exp(-0.5*(SSa-SSb)/sigmaSqu))^(beta(b)-beta(i)) );
        rho_ca = min(1, (exp(-0.5*(SSc-SSa)/sigmaSqu))^(beta(i)-beta(c)) );

        differential = t0/(nu*(t+t0)) * (rho_ab - rho_ca);
    
        Temperatures(t+1, i) = Temperatures(t+1, i-1) + exp(log(Temperatures(t, i)-Temperatures(t, i-1)) + differential);

    end
    

end

function R = GelmanR_statistic(idx, chain, opt)

%-----------------------------------------------------------------------------------------
%
%-----------------------------------------------------------------------------------------


%   Split each chain into half and check all the resulting half-sequences
    index           = floor(0.5 * idx); 
    eachChain       = zeros(index, opt.dimension);
    betweenMean     = zeros(opt.Nchain, opt.dimension);
    withinVariance  = zeros(opt.Nchain, opt.dimension);   
    
%   Mean and variance of each half-sequence chain
    for i = 1: opt.Nchain
        for j = 1: opt.dimension
            for k = 1: index
                eachChain(k,j) = chain(i,j,k+index);
            end
        end
        
        betweenMean(i,:)    = mean(eachChain);
        withinVariance(i,:) = var(eachChain);
        
    end
 
%   Between-sequence variance
    Sum = 0;
    for i = 1: opt.Nchain
       Sum = Sum + (betweenMean(i,:) - mean(betweenMean)) .^ 2;
    end
    B = Sum ./ (opt.Nchain-1);

%   Within-sequence variance
    Sum = 0;
    for i = 1: opt.Nchain
        Sum = Sum + withinVariance(i,:);
    end
    W = Sum ./ opt.Nchain;

%   Convergence diagnostics
    R = sqrt(1 + B ./ W);
    
end
% =============================================================================
%  SMB - The Simulated Moving Bed Chromatography for separation of
%  target compounds, such as fructose and glucose.
%  
%  Author: QiaoLe He   E-mail: q.he@fz-juelich.de
%                                      
%  Institute: Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
%  
%  All rights reserved. Please see the license of CADET.
% =============================================================================