
classdef OptAlgorithms < handle
 % =============================================================================
% This is the class of the functions of optimization algorithms.
%
% =============================================================================   


%   DE
    methods (Static = true, Access = 'public')

        function Differential_Evolution(params)
% -----------------------------------------------------------------------------
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
% -----------------------------------------------------------------------------


            startTime = clock;

%           Get the optimizer options
            opt = OptAlgorithms.getOptions_DE(params);

%           Initialization of the population
            [Population, OptPopul] = OptAlgorithms.InitPopulation(opt);


% ----------------------------------------------------------------------------------------
%           The main loop
            for i = 1: opt.loopCount

%               The evolution of the particles in DE
                [Population, OptPopul] = OptAlgorithms.DE_Evolution(Population, OptPopul, opt);

%               Abstract best information so far from the population and display it
                XResult = exp(OptPopul(1, 1:opt.IndivSize));
                YResult = OptPopul(opt.PopulSize+1, opt.IndivSize+1);

                fprintf('Iter = %3d,   Minimum: %g,   Parameters:[%g| %g| %g| %g| %g| %g] \n',...
                    i, YResult, XResult);


%               The convergence criterion: when the standard deviation is smaller
%               than the specified tolerence
                if mod(i, 5) == 0

                    delta = std(Population(:,1:opt.IndivSize)) ./ mean(Population(:,1:opt.IndivSize));

                    if all(abs(delta) < 0.002) || i == opt.loopCount
                        break
                    end

                end


            end
%-----------------------------------------------------------------------------------------     


%           Gather some useful information and store them
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

%           Check out the dimension of the set of parameters, and the boundary limiatation
            [row, col] = size(opt.IndivSize);
            if row > 1 || col > 1
                error('getOptions_DE: The initialized dimension of the set of parameters might be wrong');
            end

            [row, col] = size(opt.IndivScope);
            if row ~= opt.IndivSize || col ~= 2
                error('getOptions_DE: Please check your setup of the range of parameters');
            end

%           Those are the specific parameters for the DE algorithm   
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


%           Initialization of the parameters
            Population = rand(opt.PopulSize, opt.IndivSize+1);

%           Use vectorization to speed up, it's actually equivalent to the for-loop  
            Population(:,1:opt.IndivSize) = repmat(opt.IndivScope(:,1)',opt.PopulSize,1) + ...
                Population(:,1:opt.IndivSize) .* repmat( (opt.IndivScope(:,2) - opt.IndivScope(:,1))',opt.PopulSize,1 );


%           Simulation of the sampled points, using subroutine simulatedMovingBed
            Population(:, opt.IndivSize+1) = arrayfun( @(idx) feval(@simulatedMovingBed, ...
                exp(Population(idx, 1:opt.IndivSize)) ), 1: opt.PopulSize );    


%           The statistics of the population
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
%               strategy 1
                    tempPop = PopMutR3 + (PopMutR1 - PopMutR2) * opt.weight;
                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;

                case 2
%               strategy 2
                    tempPop = Population(1:R, 1:C) + opt.weight * (OptPopul(1:R, 1:C) - Population(1:R, 1:C)) + (PopMutR1 - PopMutR2) * opt.weight;
                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;

                case 3
%               strategy 3
                    tempPop = OptPopul(1:R, 1:C) + (PopMutR1 - PopMutR2) .* ((1 -0.9999) * rand(R, C) + opt.weight);
                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;

                case 4
%               strategy 4
                    f1 = (1-opt.weight) * rand(R, 1) + opt.weight;
                    for i = 1: C
                        PopMutR5(:, i) = f1;
                    end

                    tempPop = PopMutR3 + (PopMutR1 - PopMutR2) .* PopMutR5;
                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;

                case 5
%               strategy 5
                    f1 = (1-opt.weight) * rand + opt.weight;
                    tempPop = PopMutR3 + (PopMutR1 - PopMutR2) * f1;
                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;

                case 6
%               strategy 6
                    if (rand < 0.5)
                        tempPop = PopMutR3 + (PopMutR1 - PopMutR2) * opt.weight;
                    else
                        tempPop = PopMutR3 + 0.5 * (opt.weight + 1.0) * (PopMutR1 + PopMutR2 - 2 * PopMutR3);
                    end

                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;
            end

%           Check the boundary limitation
            loBound = repmat(opt.IndivScope(:,1)', R, 1);
            [row, col] = find( (tempPop(1:R, 1:C) - loBound) < 0 );
            tempPop((col-1).*R + row) = loBound((col-1).*R + row);

            upBound = repmat(opt.IndivScope(:,2)', R, 1);
            [row, col] = find( (tempPop(1:R, 1:C) - upBound) > 0 );
            tempPop((col-1).*R + row) = upBound((col-1).*R + row);
            clear row col;


%           Simulate the new population and compare their objective function values
            tempValue(:, 1) = arrayfun(@(idx) feval( @simulatedMovingBed, exp(tempPop(idx, 1:C)) ), 1:R ); 

            [row, ~] = find(tempValue < Population(:,C+1));

            Population(row, 1:C) = tempPop(row, 1:C);
            Population(row, C+1) = tempValue(row);
            clear row col;


%           Rank the objective function value to find the minimum
            [minValue, minRow] = min(arrayfun(@(idx) Population(idx, C+1), 1:R));

            OptPopul(R+1, C+1) = minValue;
            OptPopul(1:R, 1:C) = repmat( Population(minRow, 1:C), R, 1 );


        end 

    end


%   PSO    
    methods (Static = true, Access = 'public') 

        function Particle_Swarm_Optimization(params)
% -----------------------------------------------------------------------------
% Particle Swarm Optimization algorithm (PSO)
%
%  PSO optimizes a problem by having a population of candidates
%  (particles), and moving these particles around in the search space
%  according to mathematical formula ovet the particle's position and
%  velocity. Each particle's movement is influenced by its own local
%  best-known position but, is also guided toward the best-known positions
%  in the search space, which are updated as better positions are found by
%  other particles
% 
%  Parameter:
%        - params. It is the specified parameters from the main function.
%        And let the main function informed that which parameters are need
%        to be optimized
% -----------------------------------------------------------------------------


            startTime = clock;


%           Get the optimizer options
            options = OptAlgorithms.getOptions_PSO(params);  

%           Initialization of the population
            [ParSwarm, OptSwarm, ToplOptSwarm] = OptAlgorithms.InitSwarm(options);


%-----------------------------------------------------------------------------------------
%           The main loop
            for k = 1:options.loopCount

%               The evolution of the particles in PSO
                [ParSwarm, OptSwarm, ToplOptSwarm] = OptAlgorithms.ParticlesEvolution...
                    (ParSwarm, OptSwarm, ToplOptSwarm, k, options);

%               Abstract best information so far from the population and display it
                XResult = exp( OptSwarm(options.swarmSize+1, 1:options.particleSize) );
                YResult = OptSwarm(options.swarmSize+1, options.particleSize+1);        


                fprintf('Iter = %3d,  Minimum: %g,  Parameters:[ %g| %g| %g| %g| %g| %g ] \n'...
                    , k, YResult, XResult);


%               The convergence criterion: when the standard deviation is smaller
%               than the specified tolerance
                if mod(k, 5) == 0

                    delta = std(ParSwarm(:,1:options.particleSize)) ./ mean(ParSwarm(:,1:options.particleSize));

                    if all(abs(delta) < 0.002) || k == options.loopCount
                        break
                    end

                end


            end
%-----------------------------------------------------------------------------------------  

%           Gather some useful information and store them
            result.optTime        = etime(clock,startTime)/3600;
            result.convergence    = delta;
            result.correlation    = corrcoef(ParSwarm(:,1:options.particleSize));
            result.PopulationSize = options.swarmSize;
            result.Iteration      = options.loopCount;
            result.Population     = ParSwarm;
            result.xValue         = XResult;
            result.yValue         = YResult;

            save(sprintf('result_%2d.mat',fix(rand*100)),'result');
            fprintf('The results have been stored in the result.mat \n');

        end

        function opt = getOptions_PSO(params)
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
%           + swarmSize. The number of the candidates (particles)
%           + particleSize. The number of the optimized parameters
%           + paramsRange. The boundary limitation of the parameters
%           + loopCount. The maximum number of algorithm's iteration
%           + wMax; wMin. The boundary of the weight in PSO's formula
%           + accelCoeff. The accelarate coefficients in PSO's formula
%           + topology. The different schemes for the particles' communication
%               * Ring Topology. Under this scheme, only the adjacent
%                   particles exchange the information. So this is the slowest
%                   scheme to transfer the best position among particles.
%               * Random Topology. Under this scheme, the communication is
%                   isolated a little bit. This is the intermediate one.
%               * without Topology. Under this scheme, the the best position
%                   around the population is going to transfer immediately to
%                   the rest particles. This is the default one in the
%                   literatures. However in this case, it will results in the
%                   unmature convergence.
% -----------------------------------------------------------------------------


%           The number of population of the swarm intellectual (SI)
            opt.swarmSize      = 20;

%           The dimension of the optimized parameters
            opt.particleSize   = length(fieldnames(params));

%           the row represents the parameter, while the column denotes the upbound and lowerbound
            opt.paramsRange    = log( [0.10    0.20;
                                       220     300;
                                       2.5e-7  3.8e-7;
                                       0.9e-8  2.3e-8;
                                       1.5e-7  2.8e-7;
                                       1.0e-7  2.4e-7;
                                       4.0e-8  5.5e-8] );   
            opt.loopCount      = 300;

%           check out the dimension of the set of parameters, and the boundary limitation
            [row, col] = size(opt.particleSize);
            if row > 1 || col > 1
                error('getOptions_PSO: The initialized dimension of the set of parameters might be wrong');
            end

            [row, col] = size(opt.paramsRange);
            if row ~= opt.particleSize || col ~= 2
                error('getOptions_PSO: Please check your setup of the range of parameters');
            end

%           Those are the specific parameters for the PSO algorithm
            opt.wMax           = 0.9;
            opt.wMin           = 0.4;
            opt.accelCoeff     = [0.6, 0.4];
            opt.isPlot         = struct('enableCurPlot', '0', 'enableStatPlot', '0');
            opt.topology       = 'Random Topology';
            opt.randCount      = int32(opt.swarmSize * 0.6);
            opt.iter           = 0;
            opt.compValue      = 1e5;

        end

        function [ParSwarm,OptSwarm,ToplOptSwarm] = InitSwarm(opt)
% -----------------------------------------------------------------------------
% The initilization of the population
%
% Parameter:
%       - opt.
%           + swarmSize. The number of the candidates (particles)
%           + particleSize. The number of the optimized parameters
%           + paramsRange. The boundary limitation of the parameters
%           + loopCount. The maximum number of algorithm's iteration
%           + wMax; wMin. The boundary of the weight in PSO's formula
%           + accelCoeff. The accelarate coefficients in PSO's formula
%           + topology. The different schemes for the particles' communication
%               * Ring Topology. Under this scheme, only the adjacent
%                   particles exchange the information. So this is the slowest
%                   scheme to transfer the best position among particles.
%               * Random Topology. Under this scheme, the communication is
%                   isolated a little bit. This is the intermediate one.
%               * without Topology. Under this scheme, the the best position
%                   around the population is going to transfer immediately to
%                   the rest particles. This is the default one in the
%                   literatures. However in this case, it will results in the
%                   unmature convergence.
%
%  Return:
%       - ParSwarm. The swarm of the particles, correspondingly the objective value
%       - OptSwarm. The local optima that each particle ever encountered
%       - ToplOptSwarm. The global optima around the population that is
%           shared by the rest particles.
% -----------------------------------------------------------------------------


%           Check out the variables which are needed in this funciton
            if nargout < 3
                error('InitSwarm: There are not enough output in the subroutine InitSwarm');
            end


%           Initilization of the parameters, velocity of parameters   
            ParSwarm = rand(opt.swarmSize, 2 * opt.particleSize + 1);

%           Use vectorization to speed up, it's actually equivalent to the for-loop  
            ParSwarm(:,1:opt.particleSize) = repmat(opt.paramsRange(:,1)',opt.swarmSize,1) + ...
                ParSwarm(:,1:opt.particleSize) .* repmat( (opt.paramsRange(:,2) - opt.paramsRange(:,1))',opt.swarmSize,1 );    


%           Simulation of the sampled points, using subroutine simulatedMovingBed
            ParSwarm(: ,2 * opt.particleSize + 1) = arrayfun( @(idx) feval(@simulatedMovingBed,...
                exp(ParSwarm(idx, 1:opt.particleSize)) ), 1:opt.swarmSize );

%           The statistics of the population
            OptSwarm = zeros(opt.swarmSize + 1, opt.particleSize + 1);
            [minValue, row] = min(ParSwarm(:,2 * opt.particleSize + 1));

            OptSwarm(1:opt.swarmSize, 1:opt.particleSize) = ParSwarm(1: opt.swarmSize, 1: opt.particleSize);
            OptSwarm(1:opt.swarmSize, opt.particleSize+1) = ParSwarm(1: opt.swarmSize, 2* opt.particleSize+1);
            OptSwarm(opt.swarmSize+1, 1:opt.particleSize) = ParSwarm(row, 1:opt.particleSize);
            OptSwarm(opt.swarmSize+1, opt.particleSize+1) = minValue;

            ToplOptSwarm = OptSwarm(1:opt.swarmSize, 1:opt.particleSize);

        end

        function [ParSwarm,OptSwarm,ToplOptSwarm] = ParticlesEvolution(ParSwarm, OptSwarm, ToplOptSwarm, CurCount, opt)
% -----------------------------------------------------------------------------
% The evolution of particles, according to the local optima and the global optima.
%
% Parameters:
%       - ParSwarm. The swarm of the particles, correspondingly the objective value
%       - OptSwarm. The local optima that each particle ever encountered
%       - ToplOptSwarm. The global optima around the population that is
%           shared by the rest particles.
%       - CurCount. The iteration number in the main loop
%       - opt. Please see the comments of the function, InitSwarm
%
% Return:
%       - ParSwarm. The swarm of the particles, correspondingly the objective value
%       - OptSwarm. The local optima that each particle ever encountered
%       - ToplOptSwarm. The global optima around the population that is
%           shared by the rest particles.
% -----------------------------------------------------------------------------


            if nargin ~= 5
                error('ParticlesEvolution: There are error in the input of the function ParticlesEvolution')
            end
            if nargout ~= 3 
                error('ParticlesEvolution: There are error in the output of the function ParticlesEvolution')
            end


%           Different strategies of the construction of the weight
            weight = opt.wMax - CurCount * ((opt.wMax - opt.wMin) / opt.loopCount);
%             w=0.7;
%             w=(MaxW-MinW)*(CurCount/LoopCount)^2+(MinW-MaxW)*(2*CurCount/LoopCount)+MaxW;
%             w=MinW*(MaxW/MinW)^(1/(1+10*CurCount/LoopCount));


%           Different options of the accelaration coefficients
            c1 = opt.accelCoeff(1); 
            c2 = opt.accelCoeff(2);

%             c1=2.8;
%             c2=1.3;

%             con=1;
%             c1=4-exp(-con*abs(mean(ParSwarm(:,2*ParCol+1))-AdaptFunc(OptSwarm(ParRow+1,:))));
%             c2=4-c1;


%           LocalOptDiff: Each particle compares the difference of current position and the best
%           position it ever encountered
            [ParRow, ParCol] = size(ParSwarm);
            ParCol = (ParCol - 1) / 2;
            LocalOptDiff = OptSwarm(1:ParRow, 1:ParCol) - ParSwarm(:, 1:ParCol);    


            for row = 1:ParRow

%               GlobalOptDiff: Each particle compares the difference of current position and the best
%               position it ever encountered
                if strcmp(opt.topology, 'without Topology')
                    GlobalOptDiff = OptSwarm(ParRow+1, 1:ParCol) - ParSwarm(row, 1:ParCol);
                end

                if strcmp(opt.topology, 'Random Topology') || strcmp(opt.topology, 'Ring Topology')
                    GlobalOptDiff = ToplOptSwarm(row, 1:ParCol) - ParSwarm(row, 1:ParCol);
                end

%               The evolution of the velocity, according to the LocalOptDiff and GlobalOptDiff   
%               TempVelocity = weight .* ParSwarm(row,:) + c1 * unifrnd(0,1.0) .* LocalOptDiff(row,:)...
%                   + c2 * unifrnd(0,1.0) .* GlobalOptDiff;
                TempVelocity = weight .* ParSwarm(row, ParCol+1 : 2*ParCol) + c1 *LocalOptDiff(row, :) + c2 * GlobalOptDiff;


%               check the boundary limitation of the Velocity
                velocityBound = (opt.paramsRange(:,2)-opt.paramsRange(:,1) )';

                [~, veCol] = find( abs(TempVelocity(:, 1:ParCol)) > velocityBound  );
                TempVelocity(veCol) = velocityBound(veCol);

%               Value Assignment: store the updated Velocity into the matrix ParSwarm
                ParSwarm(row, ParCol+1 : 2*ParCol) = TempVelocity;



                step_size = 1;
%                 step_size = 0.729;

%               The evolution of the current positions, according to the Velocity and the step size
                ParSwarm(row, 1:ParCol) = ParSwarm(row, 1:ParCol) + step_size * TempVelocity;


%               Check the boundary limitation
                loBound = repmat(opt.paramsRange(:,1)', ParRow, 1);
                [loRow, loCol] = find( (ParSwarm(1:ParRow, 1:ParCol) - loBound) < 0 );
                ParSwarm((loCol-1).*ParRow + loRow) = loBound((loCol-1).*ParRow + loRow);

                upBound = repmat(opt.paramsRange(:,2)', ParRow, 1);
                [upRow, upCol] = find( (ParSwarm(1:ParRow, 1:ParCol) - upBound) > 0 );
                ParSwarm((upCol-1).*ParRow + upRow) = upBound((upCol-1).*ParRow + upRow);
                clear loRow loCol upRow upCol veCol;


%               Simulation of the sampled points, using subroutine simulatedMovingBed
                ParSwarm(row, 2*ParCol+1) = feval( @simulatedMovingBed, exp(ParSwarm(row, 1:ParCol)) );

%               if the updated position is better than the current position, the
%               particle flies to the updated positon; otherwise, keep still in current position
%               Update the LocalOpt for each particle
                if ParSwarm(row, 2*ParCol+1) < OptSwarm(row, ParCol+1)
                    OptSwarm(row, 1:ParCol) = ParSwarm(row, 1:ParCol);
                    OptSwarm(row, ParCol+1) = ParSwarm(row, 2*ParCol+1);
                end

%               Usint different type of the topology between particles, the global
%               optimal position is transferred to the swarm with different strategy
%               Update the GlobalOpt around the whole group
                if strcmp(opt.topology, 'Random Topology')

                    for i = 1: opt.randCount

                        rowtemp = randi(ParRow, 1);

                        if OptSwarm(row, ParCol+1) > OptSwarm(rowtemp, ParCol+1)
                            minrow = rowtemp;
                        else
                            minrow = row;
                        end

                    end

                    ToplOptSwarm(row,:) = OptSwarm(minrow,1:ParCol);

                elseif strcmp(opt.topology, 'Ring Topology')            

                    if row == 1
                        ValTemp2 = OptSwarm(ParRow, ParCol+1);
                    else
                        ValTemp2 = OptSwarm(row-1, ParCol+1);
                    end     

                    if row == ParRow
                        ValTemp3 = OptSwarm(1, ParCol+1);
                    else
                        ValTemp3 = OptSwarm(row+1, ParCol+1);
                    end

                    [~, mr] = sort([ValTemp2, OptSwarm(row,ParCol+1), ValTemp3]);
                    if  mr(1) == 3
                        if row == ParRow
                            minrow = 1;
                        else
                            minrow = row+1;
                        end
                    elseif mr(1) == 2
                        minrow = row;
                    else
                        if row == 1
                            minrow = ParRow;
                        else
                            minrow = row-1;
                        end
                    end

                    ToplOptSwarm(row,:) = OptSwarm(minrow, 1:ParCol);

                end

            end   


%           Statistics
            [minValue, row] = min(arrayfun(@(idx) OptSwarm(idx, ParCol+1), 1: ParRow));

            OptSwarm(ParRow+1, 1:ParCol) = OptSwarm(row,1:ParCol);
            OptSwarm(ParRow+1, ParCol+1) = minValue; 

        end

    end


%   MADE    
    methods (Static = true, Access = 'public')

        function Metropolis_Adjusted_Differential_Evolution(params)
% -----------------------------------------------------------------------------
% Metropolis Adjusted Differential Evolution algorithm (MADE)
%
% MADE optimizes a problem by combining the prominent features of Metropolis 
% Hastings algorithm and Differential Evolution algorithm. In the upper level, 
% each chain is accepted with the Metropolis probability, while in the lower 
% level, chains have an evolution with resort to heuristic method, Differential 
% Evolution algorithm.
%
% Unlike the algorithm, PSO and DE, the MADE obtain the parameter distributions
% rather than the single parameter set. It provides the confidential intervals 
% for each parameter, since it is based on the Markov Chain Monte Carlo (MCMC).
%
%  Parameter:
%        - params. It is the specified parameters from the main function.
%        And let the main function informed that which parameters are need
%        to be optimized
% -----------------------------------------------------------------------------


            startTime = clock;

%           Get the sampler options
            opt = OptAlgorithms.getOptions_MADE(params);

%           Initialization of the chains
            states = OptAlgorithms.initChains(opt);

%           Preallocation
            chain       = zeros(opt.Nchain, opt.dimension+1, opt.nsamples);
%             sigmaChain  = zeros(opt.nsamples, opt.Nchain);

%           Calculate the initial sigma values
%             sigmaSqu0 = sigmaSqu; n0 = 1; 
            acpt = 0;
            sigmaSqu = 0.01;

%-----------------------------------------------------------------------------------------
%           The main loop
            for i = 1: opt.nsamples

                fprintf('Iter: %5d    average acceptance: %3d,    ', i, fix(acpt/(opt.Nchain*i)*100));


                [states, acpt] = OptAlgorithms.MADE_sampler(states, sigmaSqu, acpt, opt);

%               Append the newest states to the end of the 3-D matrix
                chain(:,:,i) = states;

                if mod(i, opt.convergint) == 0 && i > opt.convergint

                    R = GelmanR_statistic(i, chain(:,:,1:i), opt);

                    if all(R < 1.01) || i == opt.nsamples
                        indx = i;
                        break
                    end

                end


%               Variance of error distribution (sigma) was treated as a parameter to be estimated.
%                 for k = 1:opt.Nchain
%                     sigmaSqu(k)  = 1 ./ OptAlgorithms.GammarDistribution(1, 1, (n0 + opt.nr)/2, 2 ./(n0 .* sigmaSqu0(k)...
%                         + chain(k, opt.dimension+1, i)));
%                     sigmaChain(i,k) = sigmaSqu(k)';
%                 end

%                 if any(sigmaSqu < 0)
%                     warning('There is something wrong in the sigma evolution, since sigma square must be nonegtive');
%                 end

            end    
%----------------------------------------------------------------------------------------- 


%           Discard the former 50%, retain only the last 50 percentage of the chains.
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


%           Thinning of the chains
            result.Population = zeros(id * opt.Nchain / opt.thinint, opt.dimension+1);
            for i = 1: id*opt.Nchain
                if mod(i, opt.thinint) == 0
                    result.Population(i/opt.thinint, :) = samplers(i, :);
                end
            end


%           Gather some useful information and store them
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
%             result.sigmaChain     = sigmaChain; 
            result.xValue         = xValue; 
            result.yValue         = yValue;

            save(sprintf('result_%3d.mat', fix(rand*1000)), 'result');
            fprintf('Markov chain simulation finished and the result saved as result.mat \n');    


%           Figures plotting
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
%
%  Parameter:
%       - params. It is the specified parameters from the main function.
%        And let the main function informed that which parameters are need
%        to be optimized.
%
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


%           Check out the dimension of the set of parameters, and the boundary limitation
            [row, col] = size(opt.dimension);
            if row > 1 || col > 1
                error('getOptions_MADE: The initialized dimension of the set of parameters might be wrong');
            end

            [row, col] = size(opt.bounds);
            if col ~= opt.dimension || row ~= 2
                error('getOptions_MADE: Please check your setup of the range of parameters');
            end

%           Those are the specific parameters for the Markov Chain Simulation
            opt.thinint       = 2;
            opt.convergint    = 50;
            opt.burnin        = 300;
            opt.printint      = 10;
            opt.DR            = 1;
            opt.iterMax       = 10000;
            opt.nr            = 1000;

%           Those are the specific parameters for the DE algorithm
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


%           The evolution of the chains in DE
            [tempPopulation, OptPopul] = OptAlgorithms.MADE_DE_Evolution(states, opt);

%           Abstract best information so far from the population and display it
            XResult = exp(OptPopul(1, 1:opt.dimension));
            YResult = OptPopul(opt.Nchain+1, opt.dimension+1);

            fprintf('Minimum: %g,    Parameters:[ %g| %g| %g| %g| %g| %g ] \n', YResult, XResult);


%           For each chain, it is accepted in terms of the Metropolis probability
            for j = 1: opt.Nchain

                proposal = tempPopulation(j, 1:opt.dimension);

                if any(proposal < opt.bounds(1,:)) || any(proposal > opt.bounds(2,:))

                    rho = 0;
                else            

                    newSS = feval( @simulatedMovingBed, exp(proposal) );
                    SS    = states(j, opt.dimension+1);

%                     rho = exp( -0.5*(newSS - SS) / sigmaSqu(j));
                    rho = exp( -0.5*(newSS - SS) / sigmaSqu);
                end

                if rand <= min(1, rho)
                    states(j, 1:opt.dimension) = proposal;
                    states(j, opt.dimension+1) = newSS;
                    acpt = acpt + 1;
                end


            end

        end

        function [tempPop, OptPopul] = MADE_DE_Evolution(Population, opt)
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
%               strategy 1
                    tempPop = PopMutR3 + (PopMutR1 - PopMutR2) * opt.weight;
                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;

                case 2
%               strategy 2
                    tempPop = Population(1:R, 1:C) + opt.weight * (OptPopul(1:R, 1:C) - Population(1:R, 1:C)) + (PopMutR1 - PopMutR2) * opt.weight;
                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;

                case 3
%               strategy 3
                    tempPop = OptPopul(1:R, 1:C) + (PopMutR1 - PopMutR2) .* ((1 -0.9999) * rand(R, C) + opt.weight);
                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;

                case 4
%               strategy 4
                    f1 = (1-opt.weight) * rand(R, 1) + opt.weight;
                    for i = 1: C
                        PopMutR5(:, i) = f1;
                    end

                    tempPop = PopMutR3 + (PopMutR1 - PopMutR2) .* PopMutR5;
                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;

                case 5
%               strategy 5
                    f1 = (1-opt.weight) * rand + opt.weight;
                    tempPop = PopMutR3 + (PopMutR1 - PopMutR2) * f1;
                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;

                case 6
%               strategy 6
                    if (rand < 0.5)
                        tempPop = PopMutR3 + (PopMutR1 - PopMutR2) * opt.weight;
                    else
                        tempPop = PopMutR3 + 0.5 * (opt.weight + 1.0) * (PopMutR1 + PopMutR2 - 2 * PopMutR3);
                    end

                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;
            end


%           Check the boundary limitation
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
                    y(i, j) = OptAlgorithms.Gammar(a, b);
                end
            end

        end

        function y = Gammar(a, b)

            if a < 1

                y = OptAlgorithms.Gammar(1+a, b) * rand(1) ^ (1/a);

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


%           Split each chain into half and check all the resulting half-sequences
            index           = floor(0.5 * idx); 
            eachChain       = zeros(index, opt.dimension);
            betweenMean     = zeros(opt.Nchain, opt.dimension);
            withinVariance  = zeros(opt.Nchain, opt.dimension);   

%           Mean and variance of each half-sequence chain
            for i = 1: opt.Nchain
                for j = 1: opt.dimension
                    for k = 1: index
                        eachChain(k,j) = chain(i,j,k+index);
                    end
                end

                betweenMean(i,:)    = mean(eachChain);
                withinVariance(i,:) = var(eachChain);

            end

%           Between-sequence variance
            Sum = 0;
            for i = 1: opt.Nchain
               Sum = Sum + (betweenMean(i,:) - mean(betweenMean)) .^ 2;
            end
            B = Sum ./ (opt.Nchain-1);

%           Within-sequence variance
            Sum = 0;
            for i = 1: opt.Nchain
                Sum = Sum + withinVariance(i,:);
            end
            W = Sum ./ opt.Nchain;

%           Convergence diagnostics
            R = sqrt(1 + B ./ W);

        end
        
    end


end
% =============================================================================
%  SMB - The Simulated Moving Bed Chromatography for separation of
%  target compounds, either binary or ternary.
% 
%      Copyright © 2008-2016: Eric von Lieres, Qiaole He
% 
%      Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
% 
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================