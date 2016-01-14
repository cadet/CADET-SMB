function Metropolis_Adjusted_Differential_Evolution(params)

% =============================================================================
% Metropolis Adjusted Differential Evolution algorithm (MADE)

% MADE optimizes a problem by combining the prominent features of Metropolis 
% Hastings algorithm and Differential Evolution algorithm. In the upper level, 
% each chain is accepted with the Metropolis probability, while in the lower 
% level, chains have an evolution with resort to heuristic method, Differential 
% Evolution algorithm.

% Unlike the algorithm, PSO and DE, the MADE obtain the parameter distributions
% rather than the single parameter set. It provides the confidential intervals 
% for each parameter, since it is based on the Markov Chain Monte Carlo (MCMC).

%  Parameter:
%        - params. It is the specified parameters from the main function.
%        And let the main function informed that which parameters are need
%        to be optimized
% =============================================================================
    
    
    startTime = clock;

%   Get the sampler options
    opt = getOptions_MADE(params);
    
%   Initialization of the chains
    states = initChains(opt);

%   Prellocation
    chain       = zeros(opt.Nchain, opt.dimension+1, opt.nsamples);
%     sigmaChain  = zeros(opt.nsamples, opt.Nchain);
    
%   Calculate the initial sigma values
%     sigmaSqu = states(:, opt.dimension+1) / (opt.nr - opt.dimension);
%     sigmaSqu0 = sigmaSqu; n0 = 1; 
    acpt = 0;
    sigmaSqu = 0.01;

%-----------------------------------------------------------------------------------------
%   The main loop
    for i = 1: opt.nsamples
        
        fprintf('Iter: %5d    average acceptance: %3d,    ', i, fix(acpt/(opt.Nchain*i)*100));
        
        
        [states, acpt] = MADE_sampler(states, sigmaSqu, acpt, opt);
        
%       Append the newest states to the end of the 3-D matrix
        chain(:,:,i) = states;
        
        if mod(i, opt.convergint) == 0 && i > opt.convergint
            
            R = GelmanR_statistic(i, chain(:,:,1:i), opt);
            
            if all(R < 1.01) || i == opt.nsamples
                indx = i;
                break
            end
            
        end
  
        
%       Variance of error distribution (sigma) was treated as a parameter to be estimated.
%         for k = 1:opt.Nchain
%             sigmaSqu(k)  = 1 ./ GammarDistribution(1, 1, (n0 + opt.nr)/2, 2 ./(n0 .* sigmaSqu0(k)...
%                 + chain(k, opt.dimension+1, i)));
%             sigmaChain(i,k) = sigmaSqu(k)';
%         end

%         if any(sigmaSqu < 0)
%             warning('There is something wrong in the sigma evolution, since sigma square must be nonegtive');
%         end
        
    end    
%----------------------------------------------------------------------------------------- 
  

%   Discard the former 50%, retain only the last 50 percentage of the chains.
    if indx < opt.nsamples
        id = floor(0.5*indx);
    else
        id = floor(0.5*opt.nsamples);
    end
    
    samplers = zeros(id * opt.Nchain, opt.dimension+1);
    for i = 1:opt.Nchain
        for k = 1:id
            for j = 1: opt.dimension+1
                samplers((i-1)*id+k, j) = chain(i,j,k+(indx-id));
            end
        end
    end
        
     
%   Thinning of the chains
    result.Population = zeros(id * opt.Nchain / opt.thinint, opt.dimension+1);
    for i = 1: id*opt.Nchain
        if mod(i, opt.thinint) == 0
            result.Population(i/opt.thinint, :) = samplers(i, :);
        end
    end
    
   
%   Gather some useful information and store them
    clear samplers
    [yValue, row]         = min(result.Population(:,opt.dimension+1));
    xValue                = result.Population(row,1:opt.dimension);  
    result.optTime        = etime(clock,startTime) / 3600;
    result.convergDiagno  = R;
    result.correlation    = corrcoef(result.Population(:,1:opt.dimension));
    result.NumChain       = opt.Nchain;
    result.Ieteration     = indx;
    result.acceptanceRate = fix( acpt / (opt.Nchain * indx) * 100 );

    result.chain          = chain;
%     result.sigmaChain     = sigmaChain; 
    result.xValue         = xValue; 
    result.yValue         = yValue;
    
    save(sprintf('result_%3d.mat', fix(rand*1000)), 'result');
    fprintf('Markov chain simulation finished and the result saved as result.mat \n');    
    
    
%   Figures plotting
    figure(1);clf
    for i = 1: opt.dimension
        
        subplot(opt.dimension,1,i,'Parent', figure(1));
        
        histfit(exp(result.Population(:,i)), 50, 'kernel');
        
        xStr = sprintf('Parameter x%d', i);
        yStr = sprintf('Frequency');
        
        xlabel(xStr, 'FontSize', 10);
        ylabel(yStr, 'FontSize', 10);
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
        
    end

    figure(2);clf
    for i = 1: opt.dimension
        for j = i: opt.dimension
            
            subplot(opt.dimension, opt.dimension, j+(i-1)*opt.dimension, 'Parent', figure(2));
            
            scatter(exp(result.Population(:,j)), exp(result.Population(:,i)));
            
            xStr = sprintf('Parameter x%d', j);
            yStr = sprintf('Parameter x%d', i);
            
            xlabel(xStr, 'FontName', 'Times New Roman', 'FontSize', 10);
            ylabel(yStr, 'FontName', 'Times New Roman', 'FontSize', 10);
            set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
        end
    end
    
end    
   
function opt = getOptions_MADE(params)

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
 
    opt.Nchain        = 20;
    opt.dimension     = length(fieldnames(params));
    opt.bounds        = log ( [0.20    0.30;
                               150     230;
                               8.0e-7  10e-7;
                               0.9e-7  2.0e-7;
                               0.7e-7  2.0e-7;
                               1.0e-7  2.0e-7] )';
    opt.nsamples      = 300; 
 
    
%   Check out the dimension of the set of parameters, and the boundary limitation
    [row, col] = size(opt.dimension);
    if row > 1 || col > 1
        error('The initialized dimension of the set of parameters might be wrong');
    end
    
    [row, col] = size(opt.bounds);
    if col ~= opt.dimension && row ~= 2
        error('Please check your setup of the range of parameters');
    end

%   Those are the specific parameters for the Markov Chain Simulation
    opt.thinint       = 2;
    opt.convergint    = 50;
    opt.burnin        = 300;
    opt.printint      = 10;
    opt.DR            = 1;
    opt.iterMax       = 10000;
    opt.nr            = 1000;
    
%   Those are the specific parameters for the DE algorithm
    opt.cr_Probability  = 0.5;
    opt.weight          = 0.3;
    opt.enablePlot      = 1;
    opt.strategy        = 5;
    opt.iter            = 0;
    opt.compValue       = 1e5;



end

function initChain = initChains(opt)

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


%   Initilization of the chains   
    initChain = rand(opt.Nchain, opt.dimension+1);

%   Use vectorization to speed up, it's actually equivalent to the for-loop 
    initChain(:,1:opt.dimension) = repmat(opt.bounds(1,:),opt.Nchain,1) + ...
        initChain(:,1:opt.dimension) .* repmat( (opt.bounds(2,:) - opt.bounds(1,:)),opt.Nchain,1 );    
  
    
%   Simulation of the sampled points, using subroutine simulatedMovingBed    
    initChain(:,opt.dimension+1) = arrayfun( @(idx) feval( @simulatedMovingBed, ...
        exp( initChain(idx,1:opt.dimension)) ), 1: opt.Nchain);
    
end

function [states, acpt] = MADE_sampler(states, sigmaSqu, acpt, opt)

%-----------------------------------------------------------------------------------------
%
%-----------------------------------------------------------------------------------------


%   The evolution of the chains in DE
    [tempPopulation, OptPopul] = DE_Evolution(states, opt);
    
%   Abstract best information so far from the population and display it
    XResult = exp(OptPopul(1, 1:opt.dimension));
    YResult = OptPopul(opt.Nchain+1, opt.dimension+1);
    
    fprintf('Minimum: %g,    Parameters:[ %g| %g| %g| %g| %g| %g ] \n', YResult, XResult);

    
%   For each chain, it is accepted in terms of the Metropolis probability
    for j = 1: opt.Nchain
        
        proposal = tempPopulation(j, 1:opt.dimension);
        
        if any(proposal < opt.bounds(1,:)) || any(proposal > opt.bounds(2,:))
            
            rho = 0;
        else            
            
            newSS = feval( @simulatedMovingBed, exp(proposal) );
            SS    = states(j, opt.dimension+1);
            
%             rho = exp( -0.5*(newSS - SS) / sigmaSqu(j));
            rho = exp( -0.5*(newSS - SS) / sigmaSqu);
        end
        
        if rand <= min(1, rho)
            states(j, 1:opt.dimension) = proposal;
            states(j, opt.dimension+1) = newSS;
            acpt = acpt + 1;
        end
        
        
    end
    
end

function [tempPop, OptPopul] = DE_Evolution(Population, opt)

% -----------------------------------------------------------------------------
% The evolution of population

% Parameters:
%       - Population. The population of the particles, correspondingly the objective value
%       - OptPopul. The best fit found among the population.
%       - opt. Please see the comments of the function, InitPopulation

% Return:
%       - Population. The population of the particles, correspondingly the objective value
%       - OptPopul. The best fit found among the population.
% -----------------------------------------------------------------------------
    
    
    R = opt.Nchain;
    C = opt.dimension;
    [minValue, row] = min(Population(:,opt.dimension+1));
    
    OptPopul = zeros(R+1, C+1);
    OptPopul(R+1, C+1) = minValue;
    
    OptPopul(1:R, 1:C) = repmat( Population(row, 1:C), R, 1 );
    clear row;
    
    indexArray = (0:1:R-1);
    index = randperm(4);
    
    cr_mutation = rand(R, C) < opt.cr_Probability;
    cr_old = cr_mutation < 0.5;

    ShuffRow1 = randperm(R);
    PopMutR1 = Population(ShuffRow1, 1:C);
    
    idxShuff = rem(indexArray + index(1), R);
    ShuffRow2 = ShuffRow1(idxShuff + 1);
    PopMutR2 = Population(ShuffRow2, 1:C);
    
    idxShuff = rem(indexArray + index(2), R);
    ShuffRow3 = ShuffRow2(idxShuff + 1);
    PopMutR3 = Population(ShuffRow3, 1:C);
    
    idxShuff = rem(indexArray + index(3), R);
    ShuffRow4 = ShuffRow3(idxShuff + 1);
    PopMutR4 = Population(ShuffRow4, 1:C);
    
    idxShuff = rem(indexArray + index(4), R);
    ShuffRow5 = ShuffRow4(idxShuff + 1);
    PopMutR5 = Population(ShuffRow5, 1:C);
    
    switch opt.strategy
        case 1
%       strategy 1
            tempPop = PopMutR3 + (PopMutR1 - PopMutR2) * opt.weight;
            tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;
        
        case 2
%       strategy 2
            tempPop = Population(1:R, 1:C) + opt.weight * (OptPopul(1:R, 1:C) - Population(1:R, 1:C)) + (PopMutR1 - PopMutR2) * opt.weight;
            tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;
        
        case 3
%       strategy 3
            tempPop = OptPopul(1:R, 1:C) + (PopMutR1 - PopMutR2) .* ((1 -0.9999) * rand(R, C) + opt.weight);
            tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;
        
        case 4
%       strategy 4
            f1 = (1-opt.weight) * rand(R, 1) + opt.weight;
            for i = 1: C
                PopMutR5(:, i) = f1;
            end
            
            tempPop = PopMutR3 + (PopMutR1 - PopMutR2) .* PopMutR5;
            tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;
        
        case 5
%       strategy 5
            f1 = (1-opt.weight) * rand + opt.weight;
            tempPop = PopMutR3 + (PopMutR1 - PopMutR2) * f1;
            tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;
        
        case 6
%       strategy 6
            if (rand < 0.5)
                tempPop = PopMutR3 + (PopMutR1 - PopMutR2) * opt.weight;
            else
                tempPop = PopMutR3 + 0.5 * (opt.weight + 1.0) * (PopMutR1 + PopMutR2 - 2 * PopMutR3);
            end
            
            tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;
    end
  

%   Check the boundary limitation
    loBound = repmat(opt.bounds(1,:), R, 1);
    [row, col] = find( (tempPop(1:R, 1:C) - loBound) < 0 );
    tempPop((col-1).*R + row) = loBound((col-1).*R + row);
    
    upBound = repmat(opt.bounds(2,:), R, 1);
    [row, col] = find( (tempPop(1:R, 1:C) - upBound) > 0 );
    tempPop((col-1).*R + row) = upBound((col-1).*R + row);
    clear row col;
    
    
end  

function y = GammarDistribution(m, n, a, b)

%-----------------------------------------------------------------------------------------
% GammarDistrib random deviates from gamma distribution
% 
%  GAMMAR_MT(M,N,A,B) returns a M*N matrix of random deviates from the Gamma
%  distribution with shape parameter A and scale parameter B:
%
%  p(x|A,B) = B^-A/gamma(A)*x^(A-1)*exp(-x/B)
%
%  Uses method of Marsaglia and Tsang (2000)

% G. Marsaglia and W. W. Tsang:
% A Simple Method for Generating Gamma Variables,
% ACM Transactions on Mathematical Software, Vol. 26, No. 3, September 2000, 363-372.
%-----------------------------------------------------------------------------------------

    if nargin < 4, b = 1; end
    
    y = zeros(m, n);
    for j=1: n
        for i=1: m
            y(i, j) = Gammar(a, b);
        end
    end
    
end

function y = Gammar(a, b)

    if a < 1
        
        y = Gammar(1+a, b) * rand(1) ^ (1/a);
        
    else
        
        d = a - 1/3;
        c = 1 / sqrt(9*d);
        
        while(1)
            
            while(1)
                x = randn(1);
                v = 1 + c*x;
                
                if v > 0, break, end
                
            end
            
            v = v^3; u = rand(1);

            if u < 1 - 0.0331*x^4, break, end

            if log(u) < 0.5 * x^2 + d * (1-v+log(v)), break, end
                
        end
        
        y = b * d * v;
        
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
