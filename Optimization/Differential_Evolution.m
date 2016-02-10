function Differential_Evolution(params)

% =============================================================================
% Differential Evolution algorithm (DE)
%
% DE optimizes a problem by maintaining a population of candidate solutions
% and creating new candidate solutions by combing existing ones according
% to its mathematical formula, and then keeping whichever candidate solution
% has the best score or fitness on the optimization problem at hand
%
%  Parameter:
%        - params. It is the specified parameters from the main function.
%        And let the main function informed that which parameters are need
%        to be optimized
% =============================================================================
    

    startTime = clock;

%   Get the optimizer options
    opt = getOptions_DE(params);
    
%   Initialization of the population
    [Population, OptPopul] = InitPopulation(opt);
        

% ----------------------------------------------------------------------------------------
%   The main loop
    for i = 1: opt.loopCount
        
%       The evolution of the particles in DE
        [Population, OptPopul] = DE_Evolution(Population, OptPopul, opt);
        
%       Abstract best information so far from the population and display it
        XResult = exp(OptPopul(1, 1:opt.IndivSize));
        YResult = OptPopul(opt.PopulSize+1, opt.IndivSize+1);
                
        fprintf('Iter = %3d,   Minimum: %g,   Parameters:[%g| %g| %g| %g| %g| %g] \n',...
            i, YResult, XResult);
        
        
%       The convergence criterion: when the standard deviation is smaller
%       than the specified tolerence
        if mod(i, 5) == 0
            
            delta = std(Population(:,1:opt.IndivSize)) ./ mean(Population(:,1:opt.IndivSize));

            if all(abs(delta) < 0.005) || i == opt.loopCount
                break
            end
            
        end
        
                
    end
%-----------------------------------------------------------------------------------------     
 
    
%   Gather some useful information and store them
    result.optTime        = etime(clock,startTime)/3600;
    result.convergence    = delta;
    result.correlation    = corrcoef(Population(:,1:opt.IndivSize));
    result.PopulationSize = opt.PopulSize;
    result.Iteration      = opt.loopCount;
    result.Population     = Population;
    result.xValue         = XResult;
    result.yValue         = YResult;
        
    save(sprintf('result_%2d.mat',fix(rand*100)),'result');
    fprintf('The results have been stored in the result.mat \n');
    
end

function opt = getOptions_DE(params)

% -----------------------------------------------------------------------------
%  The parameters for the Optimizer 
%
%  Parameter:
%       - params. It is the specified parameters from the main function.
%        And let the main function informed that which parameters are need
%        to be optimized.
%
%  Return:
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
% -----------------------------------------------------------------------------


    opt = [];
    
    opt.PopulSize      = 20;
    opt.IndivSize      = length(fieldnames(params));
    opt.IndivScope     = log ( [0.20    0.30;
                               	150     230;
                               	8.0e-7  10e-7;
                               	0.9e-7  2.0e-7;
                               	0.7e-7  2.0e-7;
                               	1.0e-7  2.0e-7] );
    opt.loopCount      = 300;
    
%   Check out the dimension of the set of parameters, and the boundary limiatation
    [row, col] = size(opt.IndivSize);
    if row > 1 || col > 1
        error('The initialized dimension of the set of parameters might be wrong');
    end
    
    [row, col] = size(opt.IndivScope);
    if row ~= opt.IndivSize && col ~= 2
        error('Please check your setup of the range of parameters');
    end
    
%   Those are the specific parameters for the DE algorithm   
    opt.cr_Probability = 0.5;
    opt.weight         = 0.3;
    opt.enablePlot     = 0;
    opt.strategy       = 5;
    opt.BoundConstr    = 1;
    opt.iter           = 0;
    opt.compValue      = 1e5;
    
end

function [Population, OptPopul] = InitPopulation(opt)

% -----------------------------------------------------------------------------
% The initilization of the population
%
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
%
% Return:
%       - Population. The population of the particles, correspondingly the objective value
%       - OptPopul. The best fit found among the population.
% -----------------------------------------------------------------------------

   
%   Initialization of the parameters
    Population = rand(opt.PopulSize, opt.IndivSize+1);
    
%   Use vectorization to speed up, it's actually equivalent to the for-loop  
    Population(:,1:opt.IndivSize) = repmat(opt.IndivScope(:,1)',opt.PopulSize,1) + ...
        Population(:,1:opt.IndivSize) .* repmat( (opt.IndivScope(:,2) - opt.IndivScope(:,1))',opt.PopulSize,1 );
    

%   Simulation of the sampled points, using subroutine simulatedMovingBed
    Population(:, opt.IndivSize+1) = arrayfun( @(idx) feval(@simulatedMovingBed, ...
        exp(Population(idx, 1:opt.IndivSize)) ), 1: opt.PopulSize );    
    
    
%   The statistics of the population
    OptPopul = zeros(opt.PopulSize+1, opt.IndivSize+1);
    [minValue, row] = min(Population(:, opt.IndivSize+1));
    
    OptPopul(opt.PopulSize+1, opt.IndivSize+1) = minValue;
    OptPopul(1:opt.PopulSize, 1:opt.IndivSize) = repmat( Population(row, 1:opt.IndivSize),opt.PopulSize,1);
    
    
end

function [Population, OptPopul] = DE_Evolution(Population, OptPopul, opt)

% -----------------------------------------------------------------------------
% The evolution of population
%
% Parameters:
%       - Population. The population of the particles, correspondingly the objective value
%       - OptPopul. The best fit found among the population.
%       - opt. Please see the comments of the function, InitPopulation
%
% Return:
%       - Population. The population of the particles, correspondingly the objective value
%       - OptPopul. The best fit found among the population.
% -----------------------------------------------------------------------------
    

    R = opt.PopulSize;
    C = opt.IndivSize;

    indexArray = (0:1:R-1);
    index = randperm(4);
    
    cr_mutation = rand(R, C) < opt.cr_Probability;
    cr_old = cr_mutation < 0.5;

    ShuffRow1 = randperm(R);
    idxShuff = rem(indexArray + index(1), R);
    ShuffRow2 = ShuffRow1(idxShuff + 1);
    idxShuff = rem(indexArray + index(2), R);
    ShuffRow3 = ShuffRow2(idxShuff + 1);
    idxShuff = rem(indexArray + index(3), R);
    ShuffRow4 = ShuffRow3(idxShuff + 1);
    idxShuff = rem(indexArray + index(4), R);
    ShuffRow5 = ShuffRow4(idxShuff + 1);
    
    PopMutR1 = Population(ShuffRow1, 1:C);
    PopMutR2 = Population(ShuffRow2, 1:C);
    PopMutR3 = Population(ShuffRow3, 1:C);
    PopMutR4 = Population(ShuffRow4, 1:C);
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
    loBound = repmat(opt.IndivScope(:,1)', R, 1);
    [row, col] = find( (tempPop(1:R, 1:C) - loBound) < 0 );
    tempPop((col-1).*R + row) = loBound((col-1).*R + row);
    
    upBound = repmat(opt.IndivScope(:,2)', R, 1);
    [row, col] = find( (tempPop(1:R, 1:C) - upBound) > 0 );
    tempPop((col-1).*R + row) = upBound((col-1).*R + row);
    clear row col;
 
    
%   Simulate the new population and compare their objective function values
    tempValue(:, 1) = arrayfun(@(idx) feval( @simulatedMovingBed, exp(tempPop(idx, 1:C)) ), 1:R ); 
    
    [row, ~] = find(tempValue < Population(:,C+1));
    
    Population(row, 1:C) = tempPop(row, 1:C);
    Population(row, C+1) = tempValue(row);
    clear row col;
        
    
%   Rank the objective function value to find the minimum
    [minValue, minRow] = min(arrayfun(@(idx) Population(idx, C+1), 1:R));
    
    OptPopul(R+1, C+1) = minValue;
    OptPopul(1:R, 1:C) = repmat( Population(minRow, 1:C), R, 1 );
    
    
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