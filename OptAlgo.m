
classdef OptAlgo < handle
% =============================================================================
% This is the class of the functions of OPTimization ALGOrithms (OptAlgo).
%
% =============================================================================

    properties (Constant)
        swarm       = 20;              % population size
        sample      = 200;             % loop count for evolution (at least 100)
        dataPoint   = 1000;            % amount of data observation
        prior       = [];              % prior information (default empty)
        logScale    = true;            % log-scale for parameter domain (default true)
        FUNC        = @simulatedMovingBed;  % the objective function
    end


% Lower level algorithms, in charge of continuous decision variables optimization
% -----------------------------------------------------------------------------
%   DE
    methods (Static = true, Access = 'public')

        function [xValue, yValue] = Differential_Evolution(opt)
% -----------------------------------------------------------------------------
% Differential Evolution algorithm (DE)
%
% DE optimizes a problem by maintaining a population of candidate solutions
% and creating new candidate solutions by combing existing ones according
% to its mathematical formula, and then keeping whichever candidate solution
% has the best score or fitness on the optimization problem at hand
%
% Parameter:
%       - FUNC. The objective function from your field
%
% Returns:
%       - xValue. The optimal parameter set
%       - yValue. The value of objective function with optimal parameter set
% -----------------------------------------------------------------------------


            startTime = clock;

            if nargin < 1
                error('OptAlgo.DE: There are no enough inputs \n');
            end

            % Get the optimizer options
            opt = OptAlgo.getOptionsDE(opt);

            % Initialization of the population
            [Population, OptPopul] = OptAlgo.initChainDE(opt);

%----------------------------------------------------------------------------------------

% The main loop
% ---------------------------------------------------------------------------------------
            for i = 1:opt.nsamples

                % The evolution of each generation
                [Population, OptPopul] = OptAlgo.evolutionDE(Population, OptPopul, opt);

                % Abstract best information so far from the population and display it
                yValue = OptPopul(opt.Nchain+1, opt.Nparams+1);
                xValue = OptAlgo.pTransfer('exp', OptPopul(1, 1:opt.Nparams));
                Value = [xValue, yValue];
                save('chainDE.dat', 'Value', '-ascii', '-tabs', '-append');

                fprintf('Iter = %5d   ----------------   Minimum: %10.3g  ---------------- \n', i, yValue);
                fprintf('%10.3g | ', xValue); fprintf('\n');

                % The convergence criterion: when the standard deviation is smaller
                %   than the specified tolerence
                if mod(i, 5) == 0
                    delta = std(Population(:, 1:opt.Nparams)) ./ mean(Population(:, 1:opt.Nparams));

                    if all(abs(delta) < opt.criterion) || i == opt.nsamples
                        maxIter = i;
                        break
                    end
                end

            end % for i = 1:opt.nsamples

%-------------------------------------------------------------------------------

% Post-process
%-------------------------------------------------------------------------------
            % Gather some useful information and store them in the debug mode
            result.optTime        = etime(clock, startTime)/3600; % in hours
            result.convergence    = delta;
            result.correlation    = corrcoef(Population(:, 1:opt.Nparams));
            result.Iteration      = maxIter;
            result.xValue         = xValue;
            result.yValue         = yValue;

            fprintf('\n****************************************************** \n');
            save(sprintf('result_%2d.mat',fix(rand*100)),'result');
            fprintf('Time %5.3g hours elapsed after %5d iterations \n', result.optTime, result.Iteration);
            fprintf('Objective value: %10.3g \n', yValue);
            fprintf('Optimal set: ');
            fprintf(' %g |', xValue);
            fprintf('\nThe statistical information of DE is stored as result.mat \n');
            fprintf('****************************************************** \n');

        end % Differential_Evolution

        function opt = getOptionsDE(obj)
% -----------------------------------------------------------------------------
%  The parameters for the Optimizer
%
%  Return:
%       - opt.
%           + Nchain. The number of the candidates (particles)
%           + Nparams. The number of the optimized parameters
%           + bounds. The boundary limitation of the parameters
%           + nsamples. The maximum number of algorithm's iteration
%           + crossProb. The cross-over probability in DE's formula
%           + weight. The weight coefficients in DE's formula
%           + strategy. There are 1,2,3,4,5,6 strategies in DE algorithm to
%               deal with the cross-over and mutation. As for detais, we
%               will refer your the original paper of Storn and Price
% -----------------------------------------------------------------------------


            opt = [];

            opt.Nchain      = OptAlgo.swarm;
            opt.Nparams     = length(fieldnames(obj.params));
            opt.bounds      = OptAlgo.pTransfer('log', obj.paramBound);
            opt.nsamples    = OptAlgo.sample;
            opt.criterion   = 0.01;

            % Check out dimension of the parameter set, and the boundary limiatation
            [row, col] = size(opt.Nparams);
            if row > 1 || col > 1
                error('OptAlgo.getOptionsDE: The initialized dimension of the set of parameters might be wrong \n');
            end

            [row, col] = size(opt.bounds);
            if row ~= opt.Nparams || col ~= 2
                error('OptAlgo.getOptionsDE: Please check your setup of the parameter range \n');
            end

            if mod(opt.nsamples, 5) ~= 0
                error('OptAlgo.getOptionsDE: Please set the maximum iteration be divisible to 5 \n');
            end

            % Those are the specific parameters for the DE algorithm
            opt.crossProb   = 0.5;
            opt.weight      = 0.3;
            opt.strategy    = 5;

        end % getOptionsDE

        function [Population, OptPopul] = initChainDE(opt)
% -----------------------------------------------------------------------------
% The initilization of the population
%
%  Parameter:
%       - opt.
%           + Nchain. The number of the candidates (particles)
%           + Nparams. The number of the optimized parameters
%           + bounds. The boundary limitation of the parameters
%           + nsamples. The maximum number of algorithm's iteration
%           + crossProb. The cross-over probability in DE's formula
%           + weight. The weight coefficients in DE's formula
%           + strategy. There are 1,2,3,4,5,6 strategies in DE algorithm to
%               deal with the cross-over and mutation. As for detais, we
%               will refer your the original paper of Storn and Price
%
% Return:
%       - Population. The population of the particles, correspondingly the objective value
%       - OptPopul. The best fit found among the population.
% -----------------------------------------------------------------------------


            if nargin < 1
                error('OptAlgo.initChainDE: There are no enough input arguments \n');
            end

            % Initialization of the parameters
            Population = rand(opt.Nchain, opt.Nparams+1);

            % Use vectorization to speed up. Keep the parameters in the domain
            Population(:, 1:opt.Nparams) = repmat(opt.bounds(:,1)', opt.Nchain, 1) + ...
                Population(:, 1:opt.Nparams) .* repmat( (opt.bounds(:,2) - opt.bounds(:,1))', opt.Nchain, 1 );

            % Simulation of the sampled points
            Population(:, opt.Nparams+1) = arrayfun( @(idx) feval(OptAlgo.FUNC, ...
                OptAlgo.pTransfer('exp', Population(idx, 1:opt.Nparams))), 1:opt.Nchain );

            % if the parallel toolbox is available
%            value = zeros(1, opt.Nchain);
%            parameter = Population(1:opt.Nchain, 1:opt.Nparams);
%            parfor i = 1:opt.Nchain
%               value(i)= feval( OptAlgo.FUNC, OptAlgo.pTransfer('exp', parameter(i,:)) )
%            end
%            Population(:,opt.Nparams+1) = value';

            % The statistics of the population
            [minValue, row] = min(Population(:, opt.Nparams+1));

            % Preallocation and allocation
            OptPopul = zeros(opt.Nchain+1, opt.Nparams+1);
            OptPopul(opt.Nchain+1, opt.Nparams+1) = minValue;
            OptPopul(1:opt.Nchain, 1:opt.Nparams) = repmat( Population(row, 1:opt.Nparams), opt.Nchain, 1 );

        end % initChainDE

        function [Population, OptPopul] = evolutionDE(Population, OptPopul, opt)
% -----------------------------------------------------------------------------
% The evolution of population
%
% Parameters:
%       - Population. The population of the particles, correspondingly the objective value
%       - OptPopul. The best fit found among the population.
%       - opt. Please see the comments of the function, initChainDE
%
% Return:
%       - Population. The population of the particles, correspondingly the objective value
%       - OptPopul. The best fit found among the population.
% -----------------------------------------------------------------------------


            if nargin < 3
                error('OptAlgo.evolutionDE: There are no enough input arguments \n');
            end

            R = opt.Nchain;
            C = opt.Nparams;

            indexArray = (0:1:R-1);
            index = randperm(4);

            cr_mutation = rand(R, C) < opt.crossProb;
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
                % strategy 1
                    tempPop = PopMutR3 + (PopMutR1 - PopMutR2) * opt.weight;
                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;

                case 2
                % strategy 2
                    tempPop = Population(1:R, 1:C) + opt.weight * (OptPopul(1:R, 1:C) - ...
                        Population(1:R, 1:C)) + (PopMutR1 - PopMutR2) * opt.weight;
                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;

                case 3
                % strategy 3
                    tempPop = OptPopul(1:R, 1:C) + (PopMutR1 - PopMutR2) .* ((1 -0.9999) * rand(R, C) + opt.weight);
                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;

                case 4
                % strategy 4
                    f1 = (1-opt.weight) * rand(R, 1) + opt.weight;
                    PopMutR5 = repmat(f1, 1, C);

                    tempPop = PopMutR3 + (PopMutR1 - PopMutR2) .* PopMutR5;
                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;

                case 5
                % strategy 5
                    f1 = (1-opt.weight) * rand + opt.weight;
                    tempPop = PopMutR3 + (PopMutR1 - PopMutR2) * f1;
                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;

                case 6
                % strategy 6
                    if (rand < 0.5)
                        tempPop = PopMutR3 + (PopMutR1 - PopMutR2) * opt.weight;
                    else
                        tempPop = PopMutR3 + 0.5 * (opt.weight + 1.0) * (PopMutR1 + PopMutR2 - 2 * PopMutR3);
                    end

                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;
            end

            % Check the boundary limitation
            loBound = repmat(opt.bounds(:,1)', R, 1);
            upBound = repmat(opt.bounds(:,2)', R, 1);
            tempPop(tempPop < loBound) = loBound(tempPop < loBound);
            tempPop(tempPop > upBound) = upBound(tempPop > upBound);

            % Simulate the new population and compare their objective function values
            tempValue(:, 1) = arrayfun( @(idx) feval( OptAlgo.FUNC, ...
                OptAlgo.pTransfer('exp', tempPop(idx, 1:C)) ), 1:R );

            % if the parallel toolbox is available
%            parfor i = 1:R
%                tempValue(i,1)= feval( OptAlgo.FUNC, OptAlgo.pTransfer('exp', tempPop(i,:)) )
%            end

            Population(tempValue < Population(:, C+1), 1:C) = tempPop(tempValue < Population(:, C+1), 1:C);
            Population(tempValue < Population(:, C+1), C+1) = tempValue(tempValue < Population(:, C+1));

            % Rank the objective function value to find the minimum
            [minValue, minRow] = min(Population(:, C+1));

            OptPopul(R+1, C+1) = minValue;
            OptPopul(1:R, 1:C) = repmat( Population(minRow, 1:C), R, 1 );

        end % evolutionDE


    end % DE


%   PSO
    methods (Static = true, Access = 'public')

        function [xValue, yValue] = Particle_Swarm_Optimization(opt)
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
% Returns:
%       - xValue. The optimal parameter set
%       - yValue. The value of objective function with optimal parameter set
% -----------------------------------------------------------------------------


            startTime = clock;

            if nargin < 1
                error('OptAlgo.PSO: There are no enough inputs \n');
            end

            % Get the optimizer options
            opt = OptAlgo.getOptionsPSO(opt);

            % Initialization of the population
            [ParSwarm, OptSwarm, ToplOptSwarm] = OptAlgo.initChainPSO(opt);

%-------------------------------------------------------------------------------

% The main loop
%-------------------------------------------------------------------------------
            for i = 1:opt.nsamples

                % The evolution of the particles in PSO
                [ParSwarm, OptSwarm, ToplOptSwarm] = OptAlgo.evolutionPSO ...
                    (ParSwarm, OptSwarm, ToplOptSwarm, i, opt);

                % Abstract best information so far from the population and display it
                yValue = OptSwarm(opt.Nchain+1, opt.Nparams+1);
                xValue = OptAlgo.pTransfer('exp', OptSwarm(opt.Nchain+1, 1:opt.Nparams));
                Value = [xValue, yValue];
                save('chainPSO.dat', 'Value', '-ascii', '-tabs', '-append');

                fprintf('Iter = %5d   ----------------  Minimum: %10.3g  ---------------- \n', i, yValue);
                fprintf('%10.3g | ', xValue); fprintf('\n');

                % convergence criterion
                if mod(i, 5) == 0
                    delta = std(ParSwarm(:, 1:opt.Nparams)) ./ mean(ParSwarm(:, 1:opt.Nparams));

                    if all(abs(delta) < opt.criterion) || i == opt.nsamples
                        maxIter = i;
                        break
                    end
                end

            end % for i = 1:opt.nsamples

%-------------------------------------------------------------------------------

% Post-process
%-------------------------------------------------------------------------------
            % Gather some useful information and store them in the debug mode
            result.optTime        = etime(clock, startTime)/3600; % in hours
            result.convergence    = delta;
            result.correlation    = corrcoef(ParSwarm(:, 1:opt.Nparams));
            result.Iteration      = maxIter;
            result.xValue         = xValue;
            result.yValue         = yValue;

            fprintf('\n****************************************************** \n');
            save(sprintf('result_%2d.mat',fix(rand*100)),'result');
            fprintf('Time %5.3g hours elapsed after %5d iterations \n', result.optTime, result.Iteration);
            fprintf('Objective value: %10.3g \n', yValue);
            fprintf('Optimal set: ');
            fprintf(' %g |', xValue);
            fprintf('\nThe statistical information of PSO is stored as result.mat \n');
            fprintf('****************************************************** \n');

        end % Particle_Swarm_Optimization

        function opt = getOptionsPSO(obj)
% -----------------------------------------------------------------------------
%  The parameters for the Optimizer
%
%  Return:
%       - opt.
%           + Nchain. The number of the candidates (particles)
%           + Nparams. The number of the optimized parameters
%           + bounds. The boundary limitation of the parameters
%           + nsamples. The maximum number of algorithm's iteration
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


            opt = [];

            opt.Nchain      = OptAlgo.swarm;
            opt.Nparams     = length(fieldnames(obj.params));
            opt.bounds      = OptAlgo.pTransfer('log', obj.paramBound);
            opt.nsamples    = OptAlgo.sample;
            opt.criterion   = 0.01;

            % Check out the dimension of the set of parameters, and the boundary limitation
            [row, col] = size(opt.Nparams);
            if row > 1 || col > 1
                error('OptAlgo.getOptionsPSO: The initialized dimension of the parameter set might be wrong \n');
            end

            [row, col] = size(opt.bounds);
            if row ~= opt.Nparams || col ~= 2
                error('OptAlgo.getOptionsPSO: Please check your setup of the range of parameters \n');
            end

            if mod(opt.nsamples, 5) ~= 0
                error('OptAlgo.getOptionsPSO: Please let the maximum interation be divisible to 5 \n');
            end

            % Those are the specific parameters for the PSO algorithm
            opt.wMax           = 0.9;
            opt.wMin           = 0.4;
            opt.accelCoeff     = [0.6, 0.4];
            opt.topology       = 'Random'; % Null; Ring
            opt.randCount      = int32(opt.Nchain * 0.6);

        end % getOptionsPSO

        function [ParSwarm, OptSwarm, ToplOptSwarm] = initChainPSO(opt)
% -----------------------------------------------------------------------------
% The initilization of the population
%
% Parameter:
%       - opt.
%           + Nchain. The number of the candidates (particles)
%           + Nparams. The number of the optimized parameters
%           + bounds. The boundary limitation of the parameters
%           + nsamples. The maximum number of algorithm's iteration
%           + wMax; wMin. The boundary of the weight in PSO's formula
%           + accelCoeff. The accelarate coefficients in PSO's formula
%           + topology. The different schemes for the particles' communication
%               * Ring. Under this scheme, only the adjacent
%                   particles exchange the information. So this is the slowest
%                   scheme to transfer the best position among particles.
%               * Random. Under this scheme, the communication is
%                   isolated a little bit. This is the intermediate one.
%               * Null. Under this scheme, the the best position
%                   around the population is going to transfer immediately to
%                   the rest particles. This is the default one in the
%                   literatures. However in this case, it will results in the
%                   unmature convergence.
%
%  Return:
%       - ParSwarm. The swarm of the particles, correspondingly the objective value
%       - OptSwarm. The local optima that each particle ever encountered
%       - ToplOptSwarm. The global optima that is shared by the rest particles.
% -----------------------------------------------------------------------------


            if nargout < 3
                error('OptAlgo.initChainPSO: There are not enough output \n');
            end

            % Initilization of the parameters, and velocity of parameters
            ParSwarm = rand(opt.Nchain, 2*opt.Nparams+1);

            % Use vectorization to speed up. Keep the parameters in the domain
            ParSwarm(:, 1:opt.Nparams) = repmat(opt.bounds(:,1)', opt.Nchain, 1) + ...
                ParSwarm(:, 1:opt.Nparams) .* repmat( (opt.bounds(:,2) - opt.bounds(:,1))', opt.Nchain, 1 );

            % Simulation of the sampled points
            ParSwarm(:, 2*opt.Nparams+1) = arrayfun( @(idx) feval( OptAlgo.FUNC, ...
                OptAlgo.pTransfer('exp', ParSwarm(idx, 1:opt.Nparams)) ), 1:opt.Nchain );

            % The statistics of the population
            OptSwarm = zeros(opt.Nchain+1, opt.Nparams+1);
            [minValue, minRow] = min(ParSwarm(:, 2*opt.Nparams+1));

            OptSwarm(1:opt.Nchain, 1:opt.Nparams) = ParSwarm(1:opt.Nchain, 1:opt.Nparams);
            OptSwarm(1:opt.Nchain, opt.Nparams+1) = ParSwarm(1:opt.Nchain, 2*opt.Nparams+1);
            OptSwarm(opt.Nchain+1, 1:opt.Nparams) = ParSwarm(minRow, 1:opt.Nparams);
            OptSwarm(opt.Nchain+1, opt.Nparams+1) = minValue;

            ToplOptSwarm = OptSwarm(1:opt.Nchain, 1:opt.Nparams);

        end % initChainPSO

        function [ParSwarm, OptSwarm, ToplOptSwarm] = evolutionPSO(ParSwarm, OptSwarm, ToplOptSwarm, iter, opt)
% -----------------------------------------------------------------------------
% The evolution of particles, according to the local optima and the global optima.
%
% Parameters:
%       - ParSwarm. The swarm of the particles, correspondingly the objective value
%       - OptSwarm. The local optima that each particle ever encountered
%       - ToplOptSwarm. The global optima that is shared by the rest particles.
%       - iter. The iteration number in the main loop
%       - opt. Please see the comments of the function, initChainPSO
%
% Return:
%       - ParSwarm. The swarm of the particles, correspondingly the objective value
%       - OptSwarm. The local optima that each particle ever encountered
%       - ToplOptSwarm. The global optima that is shared by the rest particles.
% -----------------------------------------------------------------------------


            if nargin < 5
                error('OptAlgo.evolutionPSO: There are no enough input arguments \n')
            end

            if nargout ~= 3
                error('OptAlgo.evolutionPSO: There are no enough output arguments \n')
            end

            R = opt.Nchain;
            C = opt.Nparams;

            % Different strategies of the construction of the weight
            weight = opt.wMax - iter * ((opt.wMax - opt.wMin) / opt.nsamples);
%            weight = 0.7;
%            weight = (opt.wMax - opt.wMin) * (iter / opt.nsamples)^2 ...
%                + (opt.wMin - opt.wMax) * (2 * iter / opt.nsamples) + opt.wMax;
%            weight = opt.wMin * (opt.wMax / opt.wMin)^(1 / (1 + 10 * iter / opt.nsamples));

            % The difference of current position and the best position the particle has encountered
            LocalOptDiff = OptSwarm(1:R, 1:C) - ParSwarm(1:R, 1:C);

            % The difference of current position and the best position the swarm has encountered
            if strcmp(opt.topology, 'Null')
                GlobalOptDiff = OptSwarm(R+1, 1:C) - ParSwarm(:, 1:C);
            elseif strcmp(opt.topology, 'Random') || strcmp(opt.topology, 'Ring')
                GlobalOptDiff = ToplOptSwarm(:, 1:C) - ParSwarm(:, 1:C);
            end

            % The evolution of the velocity matrix, according to LocalOptDiff and GlobalOptDiff
%            TempVelocity = weight .* ParSwarm(:,C+1:2*C) + ...
%               opt.accelCoeff(1) * unifrnd(0,1.0) .* LocalOptDiff + ...
%               opt.accelCoeff(2) * unifrnd(0,1.0) .* GlobalOptDiff;
            ParSwarm(:, C+1:2*C) = weight .* ParSwarm(:, C+1:2*C) + ...
                opt.accelCoeff(1) .* LocalOptDiff + opt.accelCoeff(2) .* GlobalOptDiff;

            % The evolution of the current positions, according to the Velocity and the step size
            stepSize = 1; % stepSize = 0.729;
            ParSwarm(:, 1:C) = ParSwarm(:, 1:C) + stepSize .* ParSwarm(:, C+1:2*C);

            % Check the boundary limitation
            loBound = repmat(opt.bounds(:,1)', R, 1);
            upBound = repmat(opt.bounds(:,2)', R, 1);
            ParSwarm(ParSwarm(:, 1:C) < loBound) = loBound(ParSwarm(:, 1:C) < loBound);
            ParSwarm(ParSwarm(:, 1:C) > upBound) = upBound(ParSwarm(:, 1:C) > upBound);

            % Simulation of the sampled points
            ParSwarm(:, 2*C+1) = arrayfun( @(idx) feval( OptAlgo.FUNC, ...
                OptAlgo.pTransfer('exp', ParSwarm(idx, 1:C)) ), 1:R );

            % Update the LocalOpt for each particle
            OptSwarm(ParSwarm(:, 2*C+1) < OptSwarm(1:R, C+1), 1:C) = ParSwarm(ParSwarm(:, 2*C+1) < OptSwarm(1:R, C+1), 1:C);
            OptSwarm(ParSwarm(:, 2*C+1) < OptSwarm(1:R, C+1), C+1) = ParSwarm(ParSwarm(:, 2*C+1) < OptSwarm(1:R, C+1), 2*C+1);

            for row = 1:R

                % Update the GlobalOpt around the whole group
                if strcmp(opt.topology, 'Random')

                    for i = 1:opt.randCount
                        rowtemp = randi(R, 1);

                        if OptSwarm(row, C+1) > OptSwarm(rowtemp, C+1)
                            minrow = rowtemp;
                        else
                            minrow = row;
                        end
                    end

                    ToplOptSwarm(row,:) = OptSwarm(minrow, 1:C);

                elseif strcmp(opt.topology, 'Ring')

                    if row == 1
                        ValTemp2 = OptSwarm(R, C+1);
                    else
                        ValTemp2 = OptSwarm(row-1, C+1);
                    end

                    if row == R
                        ValTemp3 = OptSwarm(1, C+1);
                    else
                        ValTemp3 = OptSwarm(row+1, C+1);
                    end

                    [~, mr] = sort([ValTemp2, OptSwarm(row, C+1), ValTemp3]);
                    if  mr(1) == 3
                        if row == R
                            minrow = 1;
                        else
                            minrow = row+1;
                        end
                    elseif mr(1) == 2
                        minrow = row;
                    else
                        if row == 1
                            minrow = R;
                        else
                            minrow = row-1;
                        end
                    end

                    ToplOptSwarm(row,:) = OptSwarm(minrow, 1:C);

                end % if strcmp(opt.topology)

            end % for row = 1:R

            % Statistics
            [minValue, minRow] = min(OptSwarm(:, C+1));

            OptSwarm(R+1, 1:C) = OptSwarm(minRow, 1:C);
            OptSwarm(R+1, C+1) = minValue;

        end % evolutionPSO


    end % PSO


%	MCMC
    methods (Static = true, Access = 'public')

        function [xValue, yValue] = Markov_Chain_Monte_Carlo(opt)
%------------------------------------------------------------------------------
% Markov Chain Monte Carlo (MCMC) simulation
%
% The central problem is that of determining the posterior probability for
% parameters given the data, p(thelta|y). With uninformative prior distribution
% and normalization constant p(y), the task is reduced to that of maximizing the
% likelihood of data, p(y|thelta).
%
% "Delayed rejection" was implemented in order to increase the acceptance ratio.
%
% If the Jacobian matrix can be obtained, it will be used to generate R
%
% Parameter:
%       - FUNC. The objective function from your field
%
% Returns:
%       - xValue. The optimal parameter set
%       - yValue. The value of objective function with optimal parameter set
%------------------------------------------------------------------------------


            startTime = clock;

            if nargin < 1
                error('OptAlgo.MCMC: There are no enough inputs \n');
            end

            % Get the MCMC options
            opt = OptAlgo.getOptionsMCMC(opt);

            % Preallocation
            accepted = 0; n0 = 1;
            lasti = 0; chaincov = []; chainmean = []; wsum = [];
            chain = zeros(opt.nsamples, opt.Nparams+1);

            % Inilization of a starting point
            [oldpar, SS, Jac] = OptAlgo.initPointMCMC(opt);

            % Construct the error standard deviation: sigma square
            if SS < 0, sumSquare = exp(SS); else, sumSquare = SS; end
            sigmaSqu = sumSquare / (opt.nDataPoint - opt.Nparams);
            sigmaSqu_0 = sigmaSqu;

            if ~isempty(Jac)

                % Construct the R = chol(covMatrix) for candidate generatio
                [~, S, V] = svd(Jac);

                % Truncated SVD
                if opt.Nparams == 1, S = S(1); else, S = diag(S); end
                S(S < 1e-3) = 1e-3;

                % Fisher information matrix
                % covariance matrix = sigmaSqu * inv(Fisher information matrix) = v'*s^{-2}*v
                covMatrix = V * diag(1./S.^2)* V' * sigmaSqu;

                % Construct the R = chol(covMatrix) for candidate generation
                R = chol(covMatrix .* 2.4^2 ./ opt.Nparams);

            else

                % If the Jacobian matrix cannot be obtained,
                %   a set of samples is used to generate R
                [R, oldpar, SS] = OptAlgo.burnInSamples(opt);

            end

%------------------------------------------------------------------------------

% Main loop
%------------------------------------------------------------------------------
            for j = 1:(opt.nsamples+opt.burn_in)

                if j > opt.burn_in
                    fprintf('Iter: %4d -------- Accept_ratio: %3d%% ---------- Minimum: %g ---------- \n',...
                        j, fix(accepted / j * 100), SS);
                    fprintf('%10.3g | ', OptAlgo.pTransfer('exp', oldpar)); fprintf('\n');
                end

                accept = false;

                % Generate the new proposal point with R
                newpar = oldpar + randn(1, opt.Nparams) * R;

                % Check the boundary limiation
                newpar( newpar < opt.bounds(1, :) ) = opt.bounds(1, newpar < opt.bounds(1, :));
                newpar( newpar > opt.bounds(2, :) ) = opt.bounds(1, newpar > opt.bounds(2, :));

                % Calculate the objective value of the new proposal
                newSS = feval( OptAlgo.FUNC, OptAlgo.pTransfer('exp', newpar) );

                % The Metropolis probability
                if OptAlgo.logScale
                    rho12 = exp( -0.5 *((newSS - SS) / sigmaSqu) + sum(newpar) - sum(oldpar) ) * ...
                        OptAlgo.priorPDF(newpar) / OptAlgo.priorPDF(oldpar);
                else
                    rho12 = exp( -0.5 * (newSS - SS) / sigmaSqu ) * ...
                        OptAlgo.priorPDF(newpar) / OptAlgo.priorPDF(oldpar);
                end

                % The new proposal is accepted with Metropolis probability
                if rand <= min(1, rho12)
                    accept   = true;
                    oldpar   = newpar;
                    SS       = newSS;
                    accepted = accepted + 1;
                end

                % If the poposal is denied, a Delayed Rejection procedure is adopted
                %   in order to increase the acceptance ratio
                if ~accept && opt.delayReject

                    % Shrink the searching domain by the factor 1/10
                    newpar2 = oldpar + randn(1, opt.Nparams) * (R ./ 10);

                    % Check the boundary limitation of the new generated point
                    newpar2(newpar2 < opt.bounds(1, :)) = opt.bounds(1, newpar2 < opt.bounds(1, :));
                    newpar2(newpar2 > opt.bounds(2, :)) = opt.bounds(1, newpar2 > opt.bounds(2, :));

                    % Calculate the objective value of the new proposal
                    newSS2 = feval( OptAlgo.FUNC, OptAlgo.pTransfer('exp', newpar2) );

                    if OptAlgo.logScale

                        rho32 = exp( -0.5 *((newSS - newSS2) / sigmaSqu) + sum(newpar) - sum(newpar2) ) * ...
                            OptAlgo.priorPDF(newpar) / OptAlgo.priorPDF(newpar2);

                        % The conventional version of calculation
%                        q2 = exp( -0.5 *((newSS2 - SS) / sigmaSqu) + sum(newpar2) - sum(oldpar) ) * ...
%                            OptAlgo.priorPDF(newpar2) / OptAlgo.priorPDF(oldpar);
%                        q1 = exp( -0.5 * (norm((newpar2 - newpar) * inv(R))^2 - norm((oldpar - newpar) * inv(R))^2) );

                        % The speed-up version of above calculation
                        q1q2 = exp( -0.5 *( (newSS2 - SS) / sigmaSqu + ...
                            (newpar2 - newpar) * (R \ (R' \ (newpar2' - newpar'))) - ...
                            (oldpar - newpar) * (R \ (R' \ (oldpar' - newpar'))) ) + ...
                            sum(newpar2) - sum(oldpar) ) * ...
                            OptAlgo.priorPDF(newpar2) / OptAlgo.priorPDF(oldpar);

                    else

                        rho32 = exp( -0.5 * (newSS - newSS2) / sigmaSqu ) * ...
                            OptAlgo.priorPDF(newpar) / OptAlgo.priorPDF(newpar2);

                        % The conventional version of calculation
%                        q2 = exp( -0.5 * (newSS2 - SS) / sigmaSqu ) * ...
%                            OptAlgo.priorPDF(newpar2) / OptAlgo.priorPDF(oldpar);
%                        q1 = exp( -0.5 * (norm((newpar2 - newpar) * inv(R))^2 - norm((oldpar - newpar) * inv(R))^2) );

                        % The speed-up version of above calculation
                        q1q2 = exp( -0.5 *( (newSS2 - SS) / sigmaSqu + ...
                            (newpar2 - newpar) * (R \ (R' \ (newpar2' - newpar'))) - ...
                            (oldpar - newpar) * (R \ (R' \ (oldpar' - newpar'))) ) ) * ...
                            OptAlgo.priorPDF(newpar2) / OptAlgo.priorPDF(oldpar);

                    end

                    rho13 = q1q2 * (1 - rho32) / (1 - rho12);

                    if rand <= min(1, rho13)
                        oldpar   = newpar2;
                        SS       = newSS2;
                        accepted = accepted + 1;
                    end

                end % if ~accept && delayReject

                % During the burn-in period, if the acceptance rate is extremly high or low,
                %   the R matrix is manually adjusted
                if j <= opt.burn_in
                    if mod(j, 50) == 0
                        if accepted/j < 0.05
                            fprintf('Acceptance ratio %3.2f smaller than 5 %%, scaled \n', accepted/j*100);
                            R = R ./ 5;
                        elseif accepted/j > 0.95
                            fprintf('Acceptance ratio %3.2f largeer than 95 %%, scaled \n', accepted/j*100);
                            R = R .* 5;
                        end
                    end
                end

                % After the burn-in period, the chain is stored
                if j > opt.burn_in
                    chain(j-opt.burn_in, 1:opt.Nparams) = oldpar;
                    chain(j-opt.burn_in, opt.Nparams+1) = SS;

                    temp = chain(j-opt.burn_in, :);
                    save('chainData.dat', 'temp', '-ascii', '-append');
                end

                if mod(j-opt.burn_in, opt.convergInt) == 0 && (j-opt.burn_in) > 0
                    % Updata the R according to previous chain
                    [chaincov, chainmean, wsum] = OptAlgo.covUpdate( chain(lasti+1:j-opt.burn_in, 1:opt.Nparams), ...
                        1, chaincov, chainmean, wsum );

                    lasti = j;
                    R = chol(chaincov + eye(opt.Nparams)*1e-7);

                    % Check the convergence condition
                    criterion = OptAlgo.Geweke( chain(1:j-opt.burn_in, 1:opt.Nparams) );
                    if all( abs(criterion) < opt.criterion ) || j == opt.nsamples+opt.burn_in
                        maxIter = j;
                        break
                    end
                end

                % Updata the sigma^2 according to the current objective value
                if SS < 0, sumSquare = exp(SS); else, sumSquare = SS; end
                sigmaSqu  = 1 / OptAlgo.GammarDistribution( 1, 1, (n0 + opt.nDataPoint)/2,...
                    2 / (n0 * sigmaSqu_0 + sumSquare) );

                save('sigmaSqu.dat', 'sigmaSqu', '-ascii', '-append');

            end % for j = 1:opt.nsamples

%------------------------------------------------------------------------------

% Post-process
%------------------------------------------------------------------------------
            clear chain;

            % Generate the population for figure plot
            Population = OptAlgo.conversionDataMCMC(maxIter, opt);

            OptAlgo.FigurePlot(Population, opt);

            [yValue, row] = min(Population(:, opt.Nparams+1));
            xValue = OptAlgo.pTransfer('exp', Population(row, 1:opt.Nparams));

            % Gather some useful information and store them
            result.optTime         = etime(clock, startTime) / 3600;
            result.covMatrix       = R' * R;
            result.Iteration       = maxIter;
            result.criterion       = criterion;
            result.accepted        = fix(accepted/maxIter * 100);
            result.xValue          = xValue;
            result.yValue          = yValue;
            result.sigma           = sqrt(sigmaSqu);

            fprintf('\n****************************************************** \n');
            save(sprintf('result_%2d.mat', fix(rand*100)), 'result');
            fprintf('Time %5.3g hours elapsed after %5d iterations \n', result.optTime, result.Iteration);
            fprintf('The minimal objective value found during sampling is: %10.3g \n', yValue);
            fprintf('The correspondingly optimal set is: ');
            fprintf(' %g |', xValue);
            fprintf('\nThe statistical information of MCMC is stored as result.mat \n');
            fprintf('The historical chain of MCMC is stored as Population.dat \n');
            fprintf('****************************************************** \n');

        end % Markov_Chain_Monte_Carlo

        function opt = getOptionsMCMC(obj)
%------------------------------------------------------------------------------
% The parameters for the Optimizer
%
% Return:
%       - opt.
%           + opt.Nparams. The number of optimized parameters
%           + opt.bounds. 2*nCol matrix of the parameter limitation
%           + opt.nsamples. The pre-defined maximal iteration
%           + opt.criterion. The error tolerance to stop the algorithm
%           + opt.burn_in. The burn-in period before adaptation begin
%           + opt.convergInt. The integer for checking of the convergence
%           + opt.rejectValue. This is specific used in the generation of R matrix
%               when Jacobian information is not available. In the generated sample,
%               the proposals whose objective value is larger than this will be rejected
%           + opt.nDataPoint. The number of observations in the measured data
%           + opt.deylayReject. By default DR= 1
%           + opt.Jacobian. If Jacobian matrix is available, set it to true
%------------------------------------------------------------------------------


            opt = [];

            opt.Nparams       = length(fieldnames(obj.params));
            opt.bounds        = OptAlgo.pTransfer('log', obj.paramBound)';
            opt.nsamples      = OptAlgo.sample;
            opt.criterion     = 0.0001;
            opt.burn_in       = 0;
            opt.convergInt    = 100;
            opt.rejectValue   = 1e5;
            opt.nDataPoint    = OptAlgo.dataPoint;
            opt.Jacobian      = false; % set it to true only when Jacobian matrix is available
            opt.delayReject   = true;

            [row, col] = size(opt.Nparams);
            if row > 1 || col > 1
                error('OptAlgo.getOptionsMCMC: The initialized dimension of the parameter set might be wrong \n');
            end

            [row, col] = size(opt.bounds);
            if col ~= opt.Nparams || row ~= 2
                error('OptAlgo.getOptionsMCMC: Please check your setup of the range of parameters \n');
            end

        end % getOptionsMCMC

        function [oldpar, SS, Jac]= initPointMCMC(opt)
%------------------------------------------------------------------------------
% Generate the initially guessing point for the MCMC algorithm
%------------------------------------------------------------------------------


            if nargin < 1
                error('OptAlgo.initPointMCMC: There are no enough input arguments \n');
            end

            oldpar = rand(1,opt.Nparams);

            oldpar(:, 1:opt.Nparams) = opt.bounds(1,:) + ...
                oldpar(:, 1:opt.Nparams) .* (opt.bounds(2,:) - opt.bounds(1,:));

            % Check the boundary limiation
%            oldpar(oldpar < opt.bounds(1, :)) = opt.bounds(1, oldpar < opt.bounds(1, :));
%            oldpar(oldpar > opt.bounds(2, :)) = opt.bounds(1, oldpar > opt.bounds(2, :));

            % Get the resudual value and Jacobian matrix of the guessing point
            Jac = [];
            if opt.Jacobian
                [SS, ~, Jac] = feval( OptAlgo.FUNC, OptAlgo.pTransfer('exp', oldpar) );
            else
                SS = feval( OptAlgo.FUNC, OptAlgo.pTransfer('exp', oldpar) );
            end

        end % initPointMCMC

        function [R, oldpar, SS] = burnInSamples(opt)
%------------------------------------------------------------------------------
% It is used for generating samples when Jocabian matrix is not available
%------------------------------------------------------------------------------


            if nargin < 1
                error('OptAlgo.burnInSamples: There are no enough input arguments \n');
            end

            Swarm = OptAlgo.swarm;

            % Generate random sample with swarm scale
            ParSwarm = rand(Swarm, opt.Nparams+1);

            % Keep all the swarm in the searching domain
            ParSwarm(:, 1:opt.Nparams) = repmat(opt.bounds(1,:), Swarm, 1) + ...
                ParSwarm(:, 1:opt.Nparams) .* repmat( (opt.bounds(2,:) - opt.bounds(1,:)), Swarm, 1 );

            % All the proposals are simulated
            ParSwarm(:, opt.Nparams+1) = arrayfun( @(idx) feval( OptAlgo.FUNC, ...
                OptAlgo.pTransfer('exp', ParSwarm(idx, 1:opt.Nparams)) ), 1:Swarm);

            % If the objective value is too big in the case, the proposals are deleted
            ParSwarm(ParSwarm(:, opt.Nparams+1) > opt.rejectValue, :) = [];
            [SS, minRow] = min(ParSwarm(:, opt.Nparams+1));
            oldpar = ParSwarm(minRow, 1:opt.Nparams);

            % Calculate the covariance matrix
            [chaincov, ~, ~] = OptAlgo.covUpdate(ParSwarm(:, 1:opt.Nparams), 1, [], [], []);

            if isempty(chaincov)
                error('OptAlgo.burnInSamples:\n %g randomly generated samples are all rejected with opt.rejectValue = %g \n Please change your rejection value in OptAlgo.getOptionsMCMC function', Swarm, opt.rejectValue);
            end

            % Cholesky decomposition
            R = chol( chaincov + eye(opt.Nparams) * 1e-7 );

        end % burnInSamples

        function [xcov, xmean, wsum] = covUpdate(x, w, oldcov, oldmean, oldwsum)
%------------------------------------------------------------------------------
% Recursive update the covariance matrix
%------------------------------------------------------------------------------


            [n, p] = size(x);

            if n == 0
                xcov = oldcov;
                xmean = oldmean;
                wsum = oldwsum;
                return
            end

            if nargin < 2 || isempty(w)
                w = 1;
            end

            if length(w) == 1
                w = ones(n,1) * w;
            end

            if nargin > 2 && ~isempty(oldcov)

                for i = 1:n
                    xi     = x(i,:);
                    wsum   = w(i);
                    xmeann = xi;
                    xmean  = oldmean + wsum / (wsum + oldwsum) * (xmeann - oldmean);

                    xcov =  oldcov + wsum ./ (wsum + oldwsum - 1) .* (oldwsum / (wsum + oldwsum) ...
                        .* ((xi - oldmean)' * (xi - oldmean)) - oldcov);
                    wsum    = wsum + oldwsum;
                    oldcov  = xcov;
                    oldmean = xmean;
                    oldwsum = wsum;
                end

            else

                wsum  = sum(w);
                xmean = zeros(1,p);
                xcov  = zeros(p,p);

                for i = 1:p
                    xmean(i) = sum(x(:,i) .* w) ./ wsum;
                end

                if wsum > 1
                    for i = 1:p
                        for j = 1:i
                            xcov(i,j) = (x(:,i) - xmean(i))' * ((x(:,j) - xmean(j)) .* w) ./ (wsum - 1);
                            if (i ~= j)
                                xcov(j,i) = xcov(i,j);
                            end
                        end
                    end
                end

            end

        end % covUpdate

        function Population = conversionDataMCMC(maxIter, opt)
%------------------------------------------------------------------------------
% Load and convert the data in ascii format to the mat format
% then generate the population for statistics
%------------------------------------------------------------------------------


            if nargin < 2
                error('OptAlgo.conversionDataMCMC: There are no enough input arguments \n');
            end

            load('chainData.dat');

            % Discard the former 50% chain, retain only the last 50%
            if maxIter < opt.nsamples + opt.burn_in
                idx = floor(0.5 * maxIter);
            else
                idx = floor(0.5 * (opt.nsamples + opt.burn_in));
            end

            if idx < length(chainData), idx = 0; end

            eval(sprintf('chainData(1:idx, :) = [];'));
            Population = chainData;

            save('population.dat', 'Population', '-ascii');

        end % conversionDataMCMC

        function z = Geweke(chain, a, b)
%------------------------------------------------------------------------------
% Geweke's MCMC convergence diagnostic
%------------------------------------------------------------------------------


            if nargin < 3
                a = 0.5;
                if nargin < 2
                    b = 0.8;
                end
            end

            n = length(chain);
            na = floor(a * n);
            nb = floor(b * n);

            if (n-na) / n >= 1
                error('OptAlgo.Geweke: Error with na and nb \n');
            end

            m1 = mean(chain(na:(nb-1),:));
            m2 = mean(chain(nb:end,:));

            % Spectral estimates for variance
            sa = OptAlgo.spectrum0(chain(na:(nb-1),:));
            sb = OptAlgo.spectrum0(chain(nb:end,:));

            z = (m1 - m2) ./ (sqrt( sa / (nb-na) + sb / (n-nb+1)));

        end % Geweke

        function s = spectrum0(x)
%------------------------------------------------------------------------------
% Spectral density at frequency zero
%------------------------------------------------------------------------------


            [m, n] = size(x);
            s = zeros(1,n);

            for i = 1:n
                spec = OptAlgo.spectrum(x(:,i),m);
                s(i) = spec(1);
            end

        end % spectrum0

        function [y, f] = spectrum(x, nfft, nw)
%------------------------------------------------------------------------------
% Power spectral density using Hanning window
%------------------------------------------------------------------------------


            if nargin < 3 || isempty(nw)
                nw = fix(nfft/4);
                if nargin < 2 || isempty(nfft)
                    nfft = min(length(x),256);
                end
            end

            noverlap = fix(nw/2);

            % Hanning window
            w = .5*(1 - cos(2*pi*(1:nw)' / (nw+1)));
            % Daniel
%            w = [0.5;ones(nw-2,1);0.5];
            n = length(x);

            if n < nw
                x(nw) = 0;
                n = nw;
            end

            x = x(:);

            k = fix((n - noverlap) / (nw - noverlap)); % no. of windows
            index = 1:nw;
            kmu = k * norm(w)^2; % Normalizing scale factor
            y = zeros(nfft,1);

            for i = 1:k
%                xw = w.*detrend(x(index),'linear');
                xw = w .* x(index);
                index = index + (nw - noverlap);
                Xx = abs(fft(xw,nfft)).^2;
                y = y + Xx;
            end

            y  = y * (1 / kmu); % normalize

            n2 = floor(nfft / 2);
            y  = y(1:n2);
            f  = 1 ./ n * (0:(n2-1));

        end % spectrum

	end % MCMC


%   MADE
    methods (Static = true, Access = 'public')

        function [xValue, yValue] = Metropolis_Adjusted_Differential_Evolution(opt)
% -----------------------------------------------------------------------------
% Metropolis Adjusted Differential Evolution algorithm (MADE)
%
% MADE optimizes a problem by combining the prominent features of Metropolis
% Hastings algorithm and Differential Evolution algorithm. In the upper level,
% each chain is accepted with the Metropolis probability, while in the lower
% level, chains have an evolution with resort to heuristic method, Differential
% Evolution algorithm.
%
% Unlike algorithms, PSO and DE, the MADE obtain the parameter distributions
% rather than the single parameter set. It provides the confidential intervals
% for each parameter, since it is based on the Markov Chain Monte Carlo (MCMC).
%
% Parameter:
%       - FUNC. The objective function from your field
%
% Returns:
%       - xValue. The optimal parameter set
%       - yValue. The value of objective function with optimal parameter set
% -----------------------------------------------------------------------------


            startTime = clock;

            if nargin < 2
                error('OptAlgo.MADE: There are no enough input arguments \n');
            end

            % Get the sampler options
            opt = OptAlgo.getOptionsMADE(opt);

            % Preallocation
            accepted = 0; n0 = 1;
            chain = zeros(opt.Nchain, opt.Nparams+1, opt.nsamples);

            % Initialization of the chains
            states = OptAlgo.initChainMADE(opt);

            % Calculate the initial sigma square values
            sumSquare = states(:, opt.Nparams+1);
            if any(sumSquare < 0)
                sumSquare(sumSquare < 0) = exp( sumSquare(sumSquare < 0) );
            end
            sigmaSqu = sumSquare ./ (opt.nDataPoint - opt.Nparams);
            sigmaSqu_0 = sigmaSqu;

%-----------------------------------------------------------------------------------------

% The main loop
%-----------------------------------------------------------------------------------------
            for i = 1:opt.nsamples

                fprintf('Iter: %4d ----- accept_ratio: %3d%%', i, fix(accepted/(opt.Nchain*i)*100));

                [states, accepted] = OptAlgo.samplerMADE(states, sigmaSqu, accepted, opt);

                % Append the newest states to the end of the 3-D matrix
                chain(:,:,i) = states;
                for j = 1:opt.Nchain
                    temp = states(j, :);
                    save(sprintf('chainData_%d.dat', j), 'temp', '-ascii', '-append');
                end

                % In each opt.convergInt interval, check the convergence condition
                if mod(i, opt.convergInt) == 0 && i >= opt.convergInt
                    criterion = OptAlgo.GelmanR(i, chain(:,:,1:i), opt);

                    if all(criterion < opt.criterion) || i == opt.nsamples
                        maxIter = i;
                        break
                    end
                end

                % Variance of error distribution (sigma) was treated as a parameter to be estimated.
                sumSquare = states(:, opt.Nparams+1);
                if any(sumSquare < 0)
                    sumSquare(sumSquare < 0) = exp( sumSquare(sumSquare < 0) );
                end

                % Gammar distribution for sigma evolution
                for k = 1:opt.Nchain
                    sigmaSqu(k)  = 1 ./ OptAlgo.GammarDistribution(1, 1, (n0 + opt.nDataPoint)/2, ...
                        2 / (n0 * sigmaSqu_0(k) + sumSquare(k)));
                end

            end % for i=1:opt.nsamples

%-----------------------------------------------------------------------------------------

% Post-process
%-----------------------------------------------------------------------------------------
            clear chain;

            % Generate the population for figure plot
            Population = OptAlgo.conversionDataMADE(maxIter, opt);

            % Plot the population
            OptAlgo.FigurePlot(Population, opt);

            [yValue, row]  = min(Population(:, opt.Nparams+1));
            xValue = OptAlgo.pTransfer('exp', Population(row, 1:opt.Nparams));

            % Gather some useful information and store them
            result.optTime        = etime(clock, startTime)/3600;
            result.Iteration      = maxIter;
            result.criterion      = criterion;
            result.accepted       = fix( accepted / (opt.Nchain * maxIter) * 100 );
            result.Nchain         = opt.Nchain;
            result.correlation    = corrcoef(Population(:, 1:opt.Nparams));
            result.population     = Population;
            result.xValue         = xValue;
            result.yValue         = yValue;

            fprintf('\n****************************************************** \n');
            save(sprintf('result_%2d.mat', fix(rand*100)), 'result');
            fprintf('Time %5.3g hours elapsed after %5d iterations \n', result.optTime, result.Iteration);
            fprintf('The minimal objective value found during sampling is: %10.3g \n', yValue);
            fprintf('The correspondingly optimal set is: ');
            fprintf(' %g |', xValue);
            fprintf('\nThe statistical information of MADE is stored as result.mat \n');
            fprintf('The historical chain of MADE is stored as Population.dat \n');
            fprintf('****************************************************** \n');

        end % Metropolis_Adjusted_Differential_Evolution

        function opt = getOptionsMADE(obj)
% -----------------------------------------------------------------------------
%  The parameters for the Optimizer
%
%  Return:
%       - opt.
%           + Nchain. The number of the candidates (particles)
%           + Nparams. The number of the optimized parameters
%           + bounds. The boundary limitation of the parameters
%           + nsamples. The maximum number of algorithm's iteration
%           + crossProb. The cross-over probability in DE's formula
%           + weight. The weight coefficients in DE's formula
%           + strategy. There are 1,2,3,4,5,6 strategies in DE algorithm to
%               deal with the cross-over and mutation. As for detais, we
%               will refer your the original paper of Storn and Price
% -----------------------------------------------------------------------------


            opt = [];

            opt.Nchain        = OptAlgo.swarm;
            opt.Nparams       = length(fieldnames(obj.params));
            opt.bounds        = OptAlgo.pTransfer('log', obj.paramBound)';
            opt.nsamples      = OptAlgo.sample;
            opt.criterion     = 1.01;
            opt.convergInt    = 100;
            opt.nDataPoint    = OptAlgo.dataPoint;

            % Check out the dimension of the set of parameters, and the boundary limitation
            [row, col] = size(opt.Nparams);
            if row > 1 || col > 1
                error('OptAlgo.getOptionsMADE: The initialized dimension of the set of parameters might be wrong \n');
            end

            [row, col] = size(opt.bounds);
            if col ~= opt.Nparams || row ~= 2
                error('OptAlgo.getOptionsMADE: Please check your setup of the range of parameters \n');
            end

            % Those are the specific parameters for the DE kernel
            opt.crossProb     = 0.5;
            opt.weight        = 0.3;
            opt.strategy      = 5;

        end % getOptionsMADE

        function initChain = initChainMADE(opt)
% -----------------------------------------------------------------------------
% The initilization of the chains
%
%  Parameter:
%       - opt.
%           + Nchain. The number of the candidates (particles)
%           + Nparams. The number of the optimized parameters
%           + bounds. The boundary limitation of the parameters
%           + nsamples. The maximum number of algorithm's iteration
%           + crossProb. The cross-over probability in DE's formula
%           + weight. The weight coefficients in DE's formula
%           + strategy. There are 1,2,3,4,5,6 strategies in DE algorithm to
%               deal with the cross-over and mutation. As for detais, we
%               will refer your the original paper of Storn and Price
%
% Return:
%       - Population. The population of the particles, correspondingly the objective value
%       - OptPopul. The best fit found among the population.
% -----------------------------------------------------------------------------


            if nargin < 1
                error('OptAlgo.initChainMADE: There are no enough input arguments \n');
            end

            % Initilization of the chains
            initChain = rand(opt.Nchain, opt.Nparams+1);

            % Keep the parameters in the domain
            initChain(:, 1:opt.Nparams) = repmat( opt.bounds(1,:), opt.Nchain, 1 ) + ...
                initChain(:, 1:opt.Nparams) .* repmat( (opt.bounds(2,:) - opt.bounds(1,:)), opt.Nchain, 1 );

            % Simulation of the sampled points
            initChain(:, opt.Nparams+1) = arrayfun( @(idx) feval( OptAlgo.FUNC, ...
                OptAlgo.pTransfer('exp', initChain(idx, 1:opt.Nparams)) ), 1:opt.Nchain);

        end % initChainMADE

        function [states, accepted] = samplerMADE(states, sigmaSqu, accepted, opt)
%-----------------------------------------------------------------------------------------
% The DE sampler for the Metropolis adjusted differential evolution method
%-----------------------------------------------------------------------------------------


            if nargin < 4
                error('OptAlgo.samplerMADE: There are no enough input arguments \n');
            end

            % The evolution of the chains with DE kernel
            [temp, OptPopul] = OptAlgo.evolutionMADE(states, opt);

            % Plot the best objective value and parameter set
            fprintf('----------------  Minimum: %g  ---------------- \n', OptPopul(opt.Nchain+1, opt.Nparams+1));
            fprintf('%10.3g | ', OptAlgo.pTransfer('exp', OptPopul(1, 1:opt.Nparams))); fprintf('\n');

            % In each chain, the proposal point is accepted in terms of the Metropolis probability
            for j = 1:opt.Nchain

                SS = states(j, opt.Nparams+1);

                proposal = temp(j, 1:opt.Nparams);

                if any(proposal < opt.bounds(1,:)) || any(proposal > opt.bounds(2,:))
                    newSS = Inf;
                else
                    newSS = feval( OptAlgo.FUNC, OptAlgo.pTransfer('exp', proposal) );
                end

                % The Metropolis probability
                if OptAlgo.logScale
                    rho = exp( -0.5 *((newSS - SS) / sigmaSqu(j)) + sum(newpar) - sum(oldpar) ) * ...
                        OptAlgo.priorPDF(proposal) / OptAlgo.priorPDF(states(j, 1:opt.Nparams));
                else
                    rho = exp( -0.5 * (newSS - SS) / sigmaSqu(j) ) * ...
                        OptAlgo.priorPDF(proposal) / OptAlgo.priorPDF(states(j, 1:opt.Nparams));
                end

                if rand <= min(1, rho)
                    states(j, 1:opt.Nparams) = proposal;
                    states(j, opt.Nparams+1) = newSS;
                    accepted = accepted + 1;
                end

            end % for j = 1:opt.Nchain

        end % samplerMADE

        function [tempPop, OptPopul] = evolutionMADE(Population, opt)
% -----------------------------------------------------------------------------
% The evolution of population
%
% Parameters:
%       - Population. The population of the particles, correspondingly the objective value
%       - OptPopul. The best fit found among the population.
%       - opt. Please see the comments of the function, initChainDE
%
% Return:
%       - Population. The population of the particles, correspondingly the objective value
%       - OptPopul. The best fit found among the population.
% -----------------------------------------------------------------------------


            if nargin < 2
                error('OptAlgo.evolutionMADE: There are no enough input arguments \n');
            end

            R = opt.Nchain;
            C = opt.Nparams;
            OptPopul = zeros(R+1, C+1);

            [minValue, row] = min(Population(:, opt.Nparams+1));

            OptPopul(R+1, C+1) = minValue;
            OptPopul(1:R, 1:C) = repmat( Population(row, 1:C), R, 1 );
            clear row;

            indexArray = (0:1:R-1);
            index = randperm(4);

            cr_mutation = rand(R, C) < opt.crossProb;
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
                % strategy 1
                    tempPop = PopMutR3 + (PopMutR1 - PopMutR2) * opt.weight;
                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;

                case 2
                % strategy 2
                    tempPop = Population(1:R, 1:C) + opt.weight * (OptPopul(1:R, 1:C) - Population(1:R, 1:C)) + (PopMutR1 - PopMutR2) * opt.weight;
                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;

                case 3
                % strategy 3
                    tempPop = OptPopul(1:R, 1:C) + (PopMutR1 - PopMutR2) .* ((1 -0.9999) * rand(R, C) + opt.weight);
                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;

                case 4
                % strategy 4
                    f1 = (1-opt.weight) * rand(R, 1) + opt.weight;
                    PopMutR5 = repmat(f1, 1, C);

                    tempPop = PopMutR3 + (PopMutR1 - PopMutR2) .* PopMutR5;
                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;

                case 5
                % strategy 5
                    f1 = (1-opt.weight) * rand + opt.weight;
                    tempPop = PopMutR3 + (PopMutR1 - PopMutR2) * f1;
                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;

                case 6
                % strategy 6
                    if (rand < 0.5)
                        tempPop = PopMutR3 + (PopMutR1 - PopMutR2) * opt.weight;
                    else
                        tempPop = PopMutR3 + 0.5 * (opt.weight + 1.0) * (PopMutR1 + PopMutR2 - 2 * PopMutR3);
                    end

                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;
            end

            % Check the boundary limitation
            loBound = repmat(opt.bounds(1,:), R, 1);
            upBound = repmat(opt.bounds(2,:), R, 1);
            tempPop(tempPop < loBound) = loBound(tempPop < loBound);
            tempPop(tempPop > upBound) = upBound(tempPop > upBound);

        end % evolutionMADE

        function Population = conversionDataMADE(maxIter, opt)
%------------------------------------------------------------------------------
% Load and convert the data in ascii format to the mat format
%   then generate the population for statistics
%------------------------------------------------------------------------------


            if nargin < 2
                error('OptAlgo.conversionDataMADE: There are no enough input arguments \n');
            end

            chain = []; Population = [];

            for i = 1:opt.Nchain

                load(sprintf('chainData_%d.dat', i));

                chain = [chain eval(sprintf('chainData_%d', i))];

            end

            save('chain.dat', 'chain', '-ascii');
            for  j = 1:opt.Nchain
                eval(sprintf('delete chainData_%d.dat',j));
            end

            % Discard the former 50% chain, retain only the last 50%
            if maxIter < opt.nsamples
                idx = floor(0.5 * maxIter);
            else
                idx = floor(0.5 * opt.nsamples);
            end

            if idx < length(chainData_1), idx = 0; end

            for k = 1:opt.Nchain
                eval(sprintf('chainData_%d(1:idx, :) = [];', k));
                Population = [Population; eval(sprintf('chainData_%d', k))];
            end

            save('population.dat', 'Population', '-ascii');

        end % conversionDataMADE


    end % MADE


%   PRIMAL
    methods (Static = true, Access = 'public')

        function [xValue, yValue] = Parallel_Riemann_Metropolis_Adjusted_Langevin(opt)
%------------------------------------------------------------------------------
% Parallel tempering RIemannian manifold Metropolis Adjusted Langevin (PRIMAL)
%
% The Langevin algorithm is combined with Metropolis probability, which leads to
% the MALA method in the literatures.
% Riemann manifold is introduced into the MALA method to enhance the searching capacity
% However, it still has the possibility to trap into local optima, so the parallel
% tempering technique is employed to tackle the multi-variable optimization
%
% Returns:
%       - xValue. The optimal parameter set
%       - yValue. The value of objective function with optimal parameter set
%------------------------------------------------------------------------------


            global accepted;

            startTime = clock;

            if nargin < 2
                error('OptAlgo.PRIMAL: There are no enough input arguments \n');
            end

            % Get the PRIMAL options
            opt = OptAlgo.getOptionsPRIMAL(opt);

            % Preallocation
            accepted = 0; n0 = 1;
            % The temperature matrix of the parallel tempering
            Temperatures       = zeros(opt.nsamples, opt.Nchain);
            Temperatures(:, 1) = ones(opt.nsamples, 1);
            Temperatures(1, :) = opt.temperature;

            chain = zeros(opt.Nchain, opt.Nparams+1, opt.nsamples);
%            sigmaChain = zeros(opt.nsamples, opt.Nchain);

            % Initialization
            [states, MetricTensor] = OptAlgo.initChainPRIMAL((1./opt.temperature), opt);

            % Calculate the initial sigma square values
            sumSquare = states(:, opt.Nparams+1);
            if any(sumSquare < 0)
                sumSquare(sumSquare < 0) = exp( sumSquare(sumSquare < 0) );
            end
            sigmaSqu  = sumSquare/ (opt.nDataPoint - opt.Nparams);
            sigmaSqu0 = sigmaSqu;

%------------------------------------------------------------------------------

% Main loop
%------------------------------------------------------------------------------
            for i = 1:(opt.nsamples+opt.burn_in)

                if i > opt.burn_in
                    fprintf('Iter: %3d ----- Accept_ratio: %2d%% ', i, fix(accepted/i*100));

                    % Abstract best information so far from the population and display it
                    [minValue, row] = min(states(1:opt.Nchain, opt.Nparams+1));

                    fprintf('----------------  Minimum: %3g  ---------------- \n', minValue);
                    fprintf('%10.3g | ', OptAlgo.pTransfer('exp', states(row, 1:opt.Nparams)) ); fprintf('\n');
                end

                % PRIMAL: Evolution of the chains
                [states, MetricTensor] = OptAlgo.samplerPRIMAL(states, MetricTensor, sigmaSqu, 1./Temperatures(i,:), opt);

                % Store the chain after burn-in period
                if i > opt.burn_in
                    chain(:,:,i) = states;

                    for j = 1:opt.Nchain
                        temp = states(j, :);
                        save(sprintf('chainData_%d.dat', j), 'temp', '-ascii', '-append');
                    end
                end

                % Implement the chain swap
                if mod(i, opt.swapInt) == 0
                    states = OptAlgo.chainSwap(states, sigmaSqu, 1./Temperatures(i,:), opt);
                end

                % Each opt.convergInt interval, check the convergence diagnostics
                if mod(i, opt.convergInt) == 0 && i >= opt.convergInt
                    criterion = OptAlgo.GelmanR(i, chain(:,:,1:i), opt);

                    if all(criterion < opt.criterion) || i == opt.nsamples
                        maxIter = i;
                        break
                    end
                end

                % Variance of error distribution (sigma) was treated as a parameter to be estimated.
                sumSquare = states(:, opt.Nparams+1);
                if any(sumSquare < 0)
                    sumSquare(sumSquare < 0) = exp( sumSquare(sumSquare < 0) );
                end

                % Gammar distribution for sigma evolution
                for k = 1:opt.Nchain
                    sigmaSqu(k)  = 1 ./ OptAlgo.GammarDistribution(1, 1, (n0 + opt.nDataPoint)/2, ...
                        2 / (n0 * sigmaSqu0(k) + sumSquare(k)));
%                    sigmaChain(i,k) = sigmaSqu(k)';
                end

                % Temperature dynamics
                Temperatures = OptAlgo.tpDynamics(i, states, Temperatures, sigmaSqu, opt);

            end  % for i = 1:opt.nsamples

%------------------------------------------------------------------------------

% Post-process
%------------------------------------------------------------------------------
            clear chain;

            % Generate the population for figure plot
            Population = OptAlgo.conversionDataPRIMAL(maxIter, opt);

            % Plot the population
            OptAlgo.FigurePlot(Population, opt)

            [yValue, row] = min(Population(:,opt.Nparams+1));
            xValue = OptAlgo.pTransfer('exp', Population(row,1:opt.Nparams));

            % Gather some useful information and store them
            result.optTime        = etime(clock, startTime)/3600;
            result.Iteration      = maxIter;
            result.criterion      = criterion;
            result.accepted       = fix(accepted/maxIter*100);
            result.NChain         = opt.Nchain;
            result.Temperatures   = Temperatures(1:opt.nsamples, :);
            result.correlation    = corrcoef(Population(1:opt.Nparams));
            result.population     = Population;
            result.xValue         = xValue;
            result.yValue         = yValue;

            fprintf('\n****************************************************** \n');
            save(sprintf('result_%2d.mat', fix(rand*100)), 'result');
            fprintf('Time %5.3g hours elapsed after %5d iterations \n', result.optTime, result.Iteration);
            fprintf('The minimal objective value found during sampling is: %10.3g \n', yValue);
            fprintf('The correspondingly optimal set is: ');
            fprintf(' %g |', xValue);
            fprintf('\nThe statistical information of MADE is stored as result.mat \n');
            fprintf('The historical chain of MADE is stored as Population.dat \n');
            fprintf('****************************************************** \n');

        end % Parallel_Riemann_Metropolis_Adjusted_Langevin

        function opt = getOptionsPRIMAL(obj)
%------------------------------------------------------------------------------
% The parameters for the optimizer
%
% Return:
%       - opt.
%           + temperature. The vector of the temperature of the parallel tempering
%           + Nchain. The number of the candidates
%           + Nparams. The number of the optimized parameters
%           + nsamples. The length of the sampling chain
%           + bounds. The boundary limitation of the sampling
%           + burn_in. The burn-in period that is discarded
%           + convergInt. The interval to check the convergence criterion
%           + swapInt. The interval to swap the N chains when sampling
%           + nDataPoint. The data point of the y
%           + criterion. The stopping tolerance
%           + epsilon. The epsilon value in the PRIMAL proposal formula
%------------------------------------------------------------------------------


            opt = [];

            opt.temperature       = [1];
            opt.Nchain            = length(opt.temperature);
            opt.Nparams           = length(fieldnames(obj.params));
            opt.nsamples          = OptAlgo.sample;
            opt.bounds            = OptAlgo.pTransfer('log', obj.paramBound)';

            opt.burn_in           = 0;
            opt.convergInt        = 100;
            opt.swapInt           = 100;
            opt.nDataPoint        = OptAlgo.dataPoint;
            opt.criterion         = 1.1;
            opt.epsilon           = 1e-3;

            if mod(opt.nsamples, opt.convergInt) ~= 0
                error('OptAlgo.getOptionsPRIMAL: Please set the samples be devisible to the convergInt \n');
            end

            [row, col] = size(opt.Nparams);
            if row > 1 || col > 1
                error('OptAlgo.getOptionsPRIMAL: The initialized dimension of the parameter set might be wrong \n');
            end

            [row, col] = size(opt.bounds);
            if col ~= opt.Nparams || row ~= 2
                error('OptAlgo.getOptionsPRIMAL: Please check your setup of the range of parameters \n');
            end

        end % getOptionsPRIMAL

        function [initChain, MetricTensor] = initChainPRIMAL(Beta, opt)
%------------------------------------------------------------------------------
% Generate initial chains for the PRIMAL algorithm
%
% parameter:
%       - Beta. The inverse of the temperature vector
%       - opt. The options of the PRIMAL
%
% Return:
%       - initChain. The Nchain x (Nparams+1) matrix
%       - MetricTensor.
%           + G. The Fisher information matrix (Metric tensor under Riemann manifold)
%           + invG. The inverse of the metric tensor
%           + GradL. The gradient of the log-density distribution L.
%           + sqrtInvG. The inverse of the Fisher information matrix
%------------------------------------------------------------------------------


            if nargin < 3
                error('OptAlgo.initChainPRIMAL: There are no enough input arguments \n');
            end

            % Initilization of the chains
            initChain = rand(opt.Nchain, opt.Nparams+1);
            MetricTensor = cell(1, opt.Nchain);

            % Use vectorization to speed up. Keep the parameters in the domain
            initChain(:, 1:opt.Nparams) = repmat(opt.bounds(1, :), opt.Nchain, 1) + ...
                initChain(:, 1:opt.Nparams) .* repmat( (opt.bounds(2,:)-opt.bounds(1,:)), opt.Nchain, 1 );

            for j = 1:opt.Nchain

                % Simulation
                [~, Res, Jac] = feval( OptAlgo.FUNC, initChain(j, 1:opt.Nparams) );

                initChain(j, opt.Nparams+1) = Res' * Res;
                SigmaSqu = (initChain(j, opt.Nparams+1) / (opt.nDataPoint - opt.Nparams));

                % Prepare the information for the PRIMAL proposal kernel
                %   metric tensor G, gradient vector, square root of inverse G
                G = Beta(j) .* (Jac' * (1/SigmaSqu) * Jac);
                MetricTensor{j}.invG = pinv( G + eye(opt.Nparams)*1e-10 );

                % square root of metrix invG = V * D * V'
                [V, D] = eig(MetricTensor{j}.invG);
                D = diag(D);
                for k = 1:opt.Nparams
                    if D(k) < 0
                        error('OptAlgo.PRIMAL: negative eigenvalue of metrix invG \nSet logScale to true in PRIMAL\n')
                    end
                    D(k) = sqrt( D(k) );
                end

                % sqrt(invG) = V * sqrt(D) * V'
                MetricTensor{j}.sqrtInvG = V * diag(D) * V';
                MetricTensor{j}.GradL = - Jac' * Res / SigmaSqu;

            end % for j = 1:opt.Nchain

        end % initChainPRIMAL

        function [states, MetricTensor] = samplerPRIMAL(states, MetricTensor, sigmaSqu, Beta, opt)
%------------------------------------------------------------------------------
% The Metropolis adjusted Langevin algorithms is used to generate new proposal
%
% Parameter:
%       - states. The Nchain x (Nparams+1) matrix
%       - MetricTensor. A struct that contains the Fisher information, gradient
%           of the log-density distribution, inverse of the metric tensor
%       - sigmaSqu. The sigma square matrix (the covariance matrix)
%       - Beta. Inverse of the temperature vector in parallel tempering
%       - opt. The PRIMAL options
%
% Return:
%       - states. Same
%       - MetricTensor. Same
%------------------------------------------------------------------------------


            if nargin < 6
                error('OptAlgo.samplerPRIMAL: There are no enough input arguments \n');
            end

            global accepted;

            for j = 1:opt.Nchain

                % New proposal formula based on the Metropolis adjusted Langevin algorithm
                proposal = states(j,1:opt.Nparams) + 0.5 * opt.epsilon^2 * (MetricTensor{j}.invG *...
                    MetricTensor{j}.GradL)'+ opt.epsilon * randn(1,opt.Nparams) * MetricTensor{j}.sqrtInvG;

                % Check the boundary limiation
                proposal(proposal < opt.bounds(1, :)) = opt.bounds(1, proposal < opt.bounds(1, :));
                proposal(proposal > opt.bounds(2, :)) = opt.bounds(2, proposal > opt.bounds(2, :));

                % Simulation of the new proposal
                [~, newRes, newJac] = feval( OptAlgo.FUNC, proposal );

                newSS = newRes' * newRes;
                SS    = states(j,opt.Nparams+1);

                % The Metropolis probability
                if OptAlgo.logScale
                    rho = ( exp( -0.5 *((newSS - SS) / sigmaSqu(j)) + sum(newpar) - sum(oldpar) ) * ...
                        OptAlgo.priorPDF(proposal) / OptAlgo.priorPDF(states(j, 1:opt.Nparams)) )^Beta(j);
                else
                    rho = ( exp( -0.5 * (newSS - SS) / sigmaSqu(j) ) * ...
                        OptAlgo.priorPDF(proposal) / OptAlgo.priorPDF(states(j, 1:opt.Nparams)) )^Beta(j);
                end

                % If the proposal is accepted
                if rand <= min(1, rho)

                    states(j, 1:opt.Nparams) = proposal;
                    states(j, opt.Nparams+1) = newSS;

                    % Prepare the information for the PRIMAL proposal kernel
                    %   metric tensor G, gradient vector, square root of inverse G
                    G = Beta(j) .* (newJac' * (1/sigmaSqu(j)) * newJac);
                    MetricTensor{j}.invG = pinv( G + eye(opt.Nparams)*1e-10 );

                    [V, D] = eig(MetricTensor{j}.invG);
                    D = diag(D);
                    for k = 1:opt.Nparams
                        if D(k) < 0
                            error('OptAlgo.PRIMAL:  negative eigenvalue of metrix invG \nSet logScale to true in PRIMAL\n')
                        end
                        D(k) = sqrt( D(k) );
                    end

                    MetricTensor{j}.sqrtInvG = V * diag(D) * V';
                    MetricTensor{j}.GradL = -newJac' * newRes / sigmaSqu(j);

                    if j == 1, accepted = accepted + 1; end

                end

            end % for j = 1:opt.Nchain

        end % samplerPRIMAL

        function states = chainSwap(states, sigmaSqu, Beta, opt)
%------------------------------------------------------------------------------
% Swap of the chains at pre-determined intervals
%
% Parameter:
%       - states. The Nchain x (Nparams+1) matrix
%       - sigmaSqu. The sigma square matrix (the covariance matrix)
%       - Beta. Inverse of the temperature vector in parallel tempering
%       - opt. The PRIMAL options
%
% Return:
%       - states. Same
%------------------------------------------------------------------------------


            if nargin < 4
                error('OptAlgo.chainSwap: There are no enough input arguments \n');
            end

            a = ceil(rand*opt.Nchain);

            if a == opt.Nchain
                b = 1;
            else
                b = a+1;
            end

            SSa = states(a, opt.Nparams+1);
            SSb = states(b, opt.Nparams+1);

            % Chains are swaped with certain Metropolis probability
            rho = ( exp(-0.5 * (SSa-SSb) / sigmaSqu(a)) )^(Beta(b) - Beta(a));

            if rand < min(1, rho)
                temp         = states(a, :);
                states(a, :) = states(b, :);
                states(b, :) = temp;
                clear temp;
            end

        end % chainSwap

        function Temperatures = tpDynamics(it, states, Temperatures, sigmaSqu, opt)
%------------------------------------------------------------------------------
% Temperature evolution in the parallel tempering
%
% parameter:
%       - it. The index of the current chain
%       - states. The Nchain x (Nparams+1) matrix
%       - Temperatures. The vector of the temperature in the parallel tempering
%       - simgaSqu. The sigma square vector (the covariance matrix)
%       - opt. The PRIMAL options
%
% Return:
%       - Temperatures. Same
%------------------------------------------------------------------------------


            if nargin < 5
                error('OptAlgo.tpDynamics: There are no enough input arguments \n');
            end

            t0 = 1e3; nu = 100;

            Beta = 1 ./ Temperatures(it, :);

            for i = 2: opt.Nchain

                b = i - 1;
                if i == opt.Nchain
                    c = 1;
                else
                    c = i + 1;
                end

                SSa = states(i, opt.Nparams+1);
                SSb = states(b, opt.Nparams+1);
                SSc = states(c, opt.Nparams+1);

                rho_ab = min( 1, (exp(-0.5*(SSa-SSb)/sigmaSqu(i)))^(Beta(b)-Beta(i)) );
                rho_ca = min( 1, (exp(-0.5*(SSc-SSa)/sigmaSqu(c)))^(Beta(i)-Beta(c)) );

                differential = t0 / (nu * (it + t0)) * (rho_ab - rho_ca);

                Temperatures(it+1, i) = Temperatures(it+1, i-1) + ...
                    exp( log(Temperatures(it, i) - Temperatures(it, i-1) ) + differential);

            end

        end % tpDynamics

        function Population = conversionDataPRIMAL(maxIter, opt)
%------------------------------------------------------------------------------
% Load and convert the data in ascii format to the mat format
%   then generate the population for statistics
%------------------------------------------------------------------------------


            if nargin < 2
                error('OptAlgo.conversionDataPRIMAL: There are no enough input arguments \n');
            end

            chain = []; Population = [];

            for i = 1:opt.Nchain

                load(sprintf('chainData_%d.dat', i));

                chain = [chain eval(sprintf('chainData_%d', i))];

            end

            save('chain.dat', 'chain', '-ascii');
            for  j = 1:opt.Nchain
                eval(sprintf('delete chainData_%d.dat',j));
            end

            % Discard the former 50% chain, retain only the last 50%
            if maxIter < opt.nsamples
                idx = floor(0.5 * maxIter);
            else
                idx = floor(0.5 * opt.nsamples);
            end

            for k = 1:opt.Nchain
                eval(sprintf('chainData_%d(1:idx, :) = [];', k));
                Population = [Population; eval(sprintf('chainData_%d', k))];
            end

            save('population.dat', 'Population', '-ascii');

        end % conversionDataPRIMAL


    end % PRIMAL


%   Miscellaneous
    methods (Static = true, Access = 'public')

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
%
% G. Marsaglia and W. W. Tsang:
% A Simple Method for Generating Gamma Variables,
% ACM Transactions on Mathematical Software, Vol. 26, No. 3, September 2000, 363-372.
%-----------------------------------------------------------------------------------------


            if nargin < 4, b = 1; end

            y = zeros(m, n);
            for j = 1:n
                for i=1: m
                    y(i, j) = Gammar(a, b);
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

            end % Gammar

        end % GammerDistribution

        function criterion = GelmanR(idx, chain, opt)
%-----------------------------------------------------------------------------------------
% Stopping criterion
%-----------------------------------------------------------------------------------------


            if nargin < 3
                error('OptAlgo.GelmanR: There are no enough input arguments \n');
            end

            % Split each chain into half and check all the resulting half-sequences
            index           = floor(0.5 * idx);
            eachChain       = zeros(index, opt.Nparams);
            betweenMean     = zeros(opt.Nchain, opt.Nparams);
            withinVariance  = zeros(opt.Nchain, opt.Nparams);

            % Mean and variance of each half-sequence chain
            for i = 1:opt.Nchain

                for j = 1:opt.Nparams
                    for k = 1:index
                        eachChain(k,j) = chain(i,j,k+index);
                    end
                end

                betweenMean(i,:)    = mean(eachChain);
                withinVariance(i,:) = var(eachChain);

            end

            % Between-sequence variance
            Sum = 0;
            for i = 1:opt.Nchain
               Sum = Sum + (betweenMean(i,:) - mean(betweenMean)) .^ 2;
            end
            B = Sum ./ (opt.Nchain-1);

            % Within-sequence variance
            Sum = 0;
            for i = 1:opt.Nchain
                Sum = Sum + withinVariance(i,:);
            end
            W = Sum ./ opt.Nchain;

            % Convergence diagnostics
            criterion = sqrt(1 + B ./ W);

        end % GelmanR

        function xLog = pTransfer(str, x)
%------------------------------------------------------------------------------
% If log-scale is enabled
%
% Parameters:
%       - str:
%           + 'exp'. X = exp(x)
%           + 'log'. X = log(x)
%       - x. Parameter set
%------------------------------------------------------------------------------


            if OptAlgo.logScale
                if strcmp('exp', str)
                    xLog = exp(x);
                elseif strcmp('log', str)
                    xLog = log(x);
                else
                    error('OptAlgo.LogTransfer \n');
                end
            else
                xLog = x;
            end

        end % logTransfer

        function prior = priorPDF(points)
%------------------------------------------------------------------------------
% This is a routine used for constructe the prior distribution for MCMC
%
% Parameters:
%       points. The estimated parameters
%
% Return:
%       prior. The prior possibility, prior = f(point).
%------------------------------------------------------------------------------


            if isempty(OptAlgo.prior)
                prior = 1;
                return;
            end

            % Load prior data file
            if OptAlgo.logScale
                data = OptAlgo.pTransfer('exp', OptAlgo.prior(:, 1:end-1));
            else
                data = OptAlgo.prior(:, 1:end-1);
            end
            [~, d] = size(data);

            % Get mesh grid and pdf of the prior
            if exist('tmp.mat', 'file') ~= 2
                OptAlgo.multivariatePrior(data, d);
            end

            proposedPoint = cell(1, d);
            for i = 1:d
                if OptAlgo.logScale
                    proposedPoint{i} = OptAlgo.pTransfer('exp', points(i));
                else
                    proposedPoint{i} = points(i);
                end
            end

            load('tmp.mat');
            % Evaluate the density value of the new proposal
            prior = interpn(fullAxisMesh{:}, pdfVal, proposedPoint{:});
            if isnan(prior), prior = 1; end

        end % priorPDF

        function multivariatePrior(data, d)
%------------------------------------------------------------------------------
% Build mesh grid and invoke mvkde (multivariate kernel density estimator) routine
%------------------------------------------------------------------------------


            % Preallocation
            temp = [];
            meshSize = 20;
            axisMesh = cell(1, d);
            fullAxisMesh = cell(1, d);

            maxVector = max(data, [], 1);
            minVector = min(data, [], 1);
            rangeVec  = maxVector - minVector;

            % Axial mesh
            for i = 1:d
                axisMesh{i} = minVector(i) : rangeVec(i)/(meshSize-1) : maxVector(i);
            end

            % Generate n-demensional grid
            [fullAxisMesh{:}] = ndgrid(axisMesh{:});

            % Reshape grid to call mvkde routine
            for i = 1:d
                temp = [temp, fullAxisMesh{i}(:)];
            end
            grids = reshape(temp, meshSize^d, d);

            % Invoke multivariate kernel density estimator
            pdfVal = OptAlgo.mvkde(data, grids);
            % Inverse reshape
            pdfVal = reshape(pdfVal, size(fullAxisMesh{1}));

            % Store the mesh and pdf information for further use
            save('tmp.mat', 'axisMesh', 'fullAxisMesh', 'pdfVal');

        end % multivariatePrior

        function FigurePlot(Population, opt)
%------------------------------------------------------------------------------
% Plot the histgram and scatter figures of the last population
%------------------------------------------------------------------------------


            if nargin < 2
                error('OptAlgo.FigurePlot: There are no enough input arguments \n');
            end

            figure(1);clf
            for i = 1: opt.Nparams

                subplot(round(opt.Nparams/2),2,i,'Parent', figure(1));

%                histfit(OptAlgo.pTransfer('exp', Population(:,i)), 50, 'kernel');
                [counts, centers] = hist( Population(:,i), 25 );
                bar(centers, counts, 'r');

                xlabel(sprintf('$\\ln(x_%d)$', i), 'FontSize', 16, 'Interpreter', 'latex');
                ylabel(sprintf('Frequency'), 'FontSize', 14, 'Interpreter', 'latex');
                set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
%                OptAlgo.tickLabelFormat(gca, 'x', '%0.2e');
%                OptAlgo.tickLabelFormat(gca, 'x', []);
%                set(gca, 'XTickLabel', num2str(get(gca, 'xTick')', '%g'));
%                OptAlgo.xtickLabelRotate([], 15, [], 'FontSize', 20, 'FontName', 'Times New Roman');
                set(gca, 'ygrid', 'on');

            end

            figure(2);clf
            for i = 1: opt.Nparams-1
                for j = i: opt.Nparams-1

                    subplot(opt.Nparams-1, opt.Nparams-1, j+(i-1)*(opt.Nparams-1), 'Parent', figure(2));

                    scatter(Population(:,j+1), Population(:,i), 8, 'o', 'MarkerEdgeColor', [0, 0.7, 0.7]);

                    xlabel(sprintf('$\\ln(x_%d)$', j+1), 'FontName', 'Times New Roman', 'FontSize', 16, 'Interpreter', 'latex');
                    ylabel(sprintf('$\\ln(x_%d)$', i), 'FontName', 'Times New Roman', 'FontSize', 16, 'Interpreter', 'latex');
                    set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
%                    OptAlgo.tickLabelFormat(gca, 'x', '%0.2e');
%                    set(gca, 'XTickLabel', num2str(get(gca, 'xTick')', '%g'));
%                    OptAlgo.xtickLabelRotate([], 15, [], 'FontSize', 20, 'FontName', 'Times New Roman');
                    grid on;

                end
            end

            figure(3);clf
            for i = 1: opt.Nparams

                subplot(round(opt.Nparams/2),2,i,'Parent', figure(3));

%                [y, x] = ksdensity(OptAlgo.pTransfer('exp', Population(:,i))); area(x, y, 'FaceColor', 'g');
                [y, x] = ksdensity( Population(:,i));
                area(x, y, 'FaceColor', 'g');

                xlabel(sprintf('$\\ln(x_%d)$', i), 'FontSize', 16, 'Interpreter', 'latex');
                ylabel(sprintf('$p(\\theta)$'), 'FontSize', 16, 'Interpreter', 'latex');
                set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
                set(gca, 'xgrid', 'on');
                set(gca, 'ygrid', 'on');

            end

        end % FigurePlot

        function [pdf, X1, X2] = mvkde(X, grid, gam)
%------------------------------------------------------------------------------
% adaptive kernel density estimation in high dimensions;
%
% INPUTS:   X  - data as a 'n' by 'd' vector;
%
%         grid - 'm' points of dimension 'd' over which pdf is computed;
%                default provided only for 2-dimensional data;
%                see example on how to construct it in higher dimensions
%
%          gam - cost/accuracy tradeoff parameter, where gam<n;
%                default value is gam=ceil(n^(1/2)); larger values
%                may result in better accuracy, but always reduce speed;
%                to speedup the code, reduce the value of "gam";
%
% OUTPUT: pdf   - the value of the estimated density at 'grid'
%         X1,X2 - grid only for 2 dimensional data
%
%
%  Reference:
%  Kernel density estimation via diffusion
%  Z. I. Botev, J. F. Grotowski, and D. P. Kroese (2010)
%  Annals of Statistics, Volume 38, Number 5, pages 2916-2957.
%------------------------------------------------------------------------------


            [n, d] = size(X);

            % begin scaling preprocessing
            MAX = max(X, [], 1);
            MIN = min(X, [], 1);
            scaling = MAX - MIN;

            MAX = MAX + scaling/10;
            MIN = MIN - scaling/10;
            scaling = MAX - MIN;

            X = bsxfun(@minus, X, MIN);
            X = bsxfun(@rdivide, X, scaling);

            % failing to provide grid
            if (nargin < 2) || isempty(grid)

                warning('Assuming data is 2 dimensional. For higher dimensions, provide a grid as in example.')

                % create meshgrid in 2-dimensions
                [X1, X2] = meshgrid( MIN(1):scaling(1)/(2^7-1):MAX(1), MIN(2):scaling(2)/(2^7-1):MAX(2) );

                % create grid for plotting
                grid=reshape([X1(:),X2(:)], 2^14, d);

            end

            mesh = bsxfun(@minus, grid, MIN);
            mesh = bsxfun(@rdivide, mesh, scaling);

            % failing to provide speed/accuracy tradeoff
            if nargin < 3
                gam = ceil(n^(1/2));
            end

            % algorithm initialization
            del = 0.1 / n^(d/(d+4));
            perm = randperm(n);
            mu = X(perm(1:gam), :);
            w = rand(1, gam);
            w = w / sum(w);
            Sig = bsxfun(@times, rand(d,d,gam), eye(d)*del);
            ent = -Inf;

            % begin algorithm
            for iter = 1:1500
                Eold = ent;

                % update parameters
                [w, mu, Sig, del, ent] = OptAlgo.regEM(w, mu, Sig, del, X);

                % stopping condition
                err = abs( (ent-Eold) / ent );
                if (err < 10^(-4)) || iter > 200, break, end
            end

            % now output density values at grid
            pdf = OptAlgo.probfun(mesh, w, mu, Sig) / prod(scaling); % evaluate density

            % adjust bandwidth for scaling
            del = del * scaling;

        end % mvkde

        function pdf = probfun(x, w, mu, Sig)
%------------------------------------------------------------------------------
%
%------------------------------------------------------------------------------


            [gam, d] = size(mu);

            pdf = 0;

            for k = 1:gam

                L = chol(Sig(:, :, k));
                s = diag(L);

                logpdf = -0.5 * sum( (bsxfun(@minus, x, mu(k, :)) / L).^2, 2 ) + log(w(k)) - ...
                    sum(log(s)) - d * log(2*pi) / 2;

                pdf = pdf + exp(logpdf);

            end

        end % probfun

        function [w, mu, Sig, del, ent] = regEM(w, mu, Sig, del, X)
%------------------------------------------------------------------------------
%
%------------------------------------------------------------------------------


            [gam, d] = size(mu);
            [n, d] = size(X);

            log_lh = zeros(n, gam);
            log_sig = log_lh;

            for i = 1:gam

                L = chol(Sig(:, :, i));

                Xcentered = bsxfun(@minus, X, mu(i,:));

                xRinv = Xcentered / L; xSig = sum((xRinv / L').^2,2) + eps;

                log_lh(:, i) =-0.5 * sum(xRinv.^2, 2) - sum(log(diag(L))) + ...
                    log(w(i)) - d * log(2*pi) / 2 - 0.5 * del^2 * trace((eye(d)/L)/L');

                log_sig(:, i) = log_lh(:, i) + log(xSig);

            end

            maxll = max (log_lh, [], 2);
            maxlsig = max (log_sig, [], 2);

            p = exp(bsxfun(@minus, log_lh, maxll));
            psig = exp(bsxfun(@minus, log_sig, maxlsig));

            density = sum(p, 2);
            psigd = sum(psig, 2);

            logpdf = log(density) + maxll;
            logpsigd = log(psigd) + maxlsig;

            p = bsxfun(@rdivide, p, density);

            ent = sum(logpdf);

            w = sum(p, 1);

            for i = find(w > 0)
                %compute mu's
                mu(i, :) = p(:, i)' * X / w(i);

                Xcentered = bsxfun(@minus, X, mu(i,:));
                Xcentered = bsxfun(@times, sqrt(p(:, i)), Xcentered);

                % compute sigmas
                Sig(:, :, i) = Xcentered' * Xcentered / w(i) + del^2 * eye(d);
            end

            % estimate curvature
            w = w / sum(w);

            curv = mean( exp(logpsigd - logpdf) );

            del = 1 / (4 * n * (4*pi)^(d/2) *curv)^(1 / (d + 2));

        end % regEM

        function tickLabelFormat(hAxes, axName, format)
%------------------------------------------------------------------------------
% Sets the format of the tick labels
%
% Syntax:
%    ticklabelformat(hAxes, axName, format)
%
% Input Parameters:
%    hAxes  - handle to the modified axes, such as returned by the gca function
%    axName - name(s) of axles to modify: 'x','y','z' or combination (e.g. 'xy')
%    format - format of the tick labels in sprintf format (e.g. '%.1f V') or a
%             function handle that will be called whenever labels need to be updated
%
%    Note: Calling TICKLABELFORMAT again with an empty ([] or '') format will revert
%    ^^^^  to Matlab's normal tick labels display behavior
%
% Examples:
%    ticklabelformat(gca,'y','%.6g V') - sets y axis on current axes to display 6 significant digits
%    ticklabelformat(gca,'xy','%.2f')  - sets x & y axes on current axes to display 2 decimal digits
%    ticklabelformat(gca,'z',@myCbFcn) - sets a function to update the Z tick labels on current axes
%    ticklabelformat(gca,'z',{@myCbFcn,extraData}) - sets an update function as above, with extra data
%
% Warning:
%    This code heavily relies on undocumented and unsupported Matlab functionality.
%    It works on Matlab 7+, but use at your own risk!
%
% Technical description and more details:
%    http://UndocumentedMatlab.com/blog/setting-axes-tick-labels-format
%
% Bugs and suggestions:
%    Please send to Yair Altman (altmany@gmail.com)
%------------------------------------------------------------------------------


            % Check # of args (we now have narginchk but this is not available on older Matlab releases)
            if nargin < 1
                help(mfilename)
                return;
            elseif nargin < 3
                error('OptAlgo.ticklabelformat: Not enough input arguments \n');
            end

            % Check input args
            if ~ishandle(hAxes) || (~isa(handle(hAxes),'axes') && ~isa(handle(hAxes),'matlab.graphics.axis.Axes'))
                error('OptAlgo.ticklabelformat: hAxes input argument must be a valid axes handle \n');
            elseif ~ischar(axName)
                error('OptAlgo.ticklabelformat: axName input argument must be a string \n');
            elseif ~isempty(format) && ~ischar(format) && ~isa(format,'function_handle') && ~iscell(format)
                error('OptAlgo.ticklabelformat: format input argument must be a string or function handle \n');
            end

            % normalize axes name(s) to lowercase
            axName = lower(axName);

            if strfind(axName,'x')
                install_adjust_ticklbl(hAxes,'X',format)
            elseif strfind(axName,'y')
                install_adjust_ticklbl(hAxes,'Y',format)
            elseif strfind(axName,'z')
                install_adjust_ticklbl(hAxes,'Z',format)
            end

            function install_adjust_ticklbl(hAxes,axName,format)
            % Install the new tick labels for the specified axes

                % If empty format was specified
                if isempty(format)

                    % Remove the current format (revert to default Matlab format)
                    set(hAxes,[axName 'TickLabelMode'],'auto')
                    setappdata(hAxes,[axName 'TickListener'],[])
                    return

                end

                % Determine whether to use the specified format as a
                % sprintf format or a user-specified callback function
                if ischar(format)
                    cb = {@adjust_ticklbl axName format};
                else
                    cb = format;
                end

                % Now install axis tick listeners to adjust tick labels
                % (use undocumented feature for adjustments)
                ha = handle(hAxes);
                propName = [axName 'Tick'];
                hp = findprop(ha,propName);

                try
                    % R2014a or older
                    hl = handle.listener(ha,hp,'PropertyPostSet',cb);

                    % Adjust tick labels now
                    % eventData.AffectedObject = hAxes;
                    % adjust_ticklbl([],eventData,axName,format)
                    set(hAxes,[axName 'TickLabelMode'],'manual')
                    set(hAxes,[axName 'TickLabelMode'],'auto')
                catch
                    % R2014b or newer
                    if iscell(cb)
                        cb = @(h,e) feval(cb{1},h,e,cb{2:end});
                    end
                    % hl(1) = addlistener(ha,propName,'PostSet',cb);
                    hl(1) = event.proplistener(ha,hp,'PostSet',cb);

                    % *Tick properties don't trigger PostSet events when updated automatically in R2014b - need to use *Lim
                    % addlistener(ha,[axName 'Lim'],'PostSet',cb);
                    hRuler = get(ha,[axName 'Ruler']);
                    % hl(2) = addlistener(hRuler,'MarkedClean',cb);
                    hl(2) = event.listener(hRuler,'MarkedClean',cb);

                    % Adjust tick labels now
                    eventData.AffectedObject = hAxes;
                    if ischar(format)
                        adjust_ticklbl([],eventData,axName,format)
                    else
                        hgfeval(format);
                    end
                    set(hAxes,[axName 'TickLabelMode'],'manual')
                    % set(hAxes,[axName 'TickLabelMode'],'auto')  % causes labels not to be updated in R2014b!
                end

                setappdata(hAxes,[axName 'TickListener'],hl)
                % drawnow;

                function adjust_ticklbl(hProp,eventData,axName,format)
                % Default tick labels update callback function (used if user did not specify their own function)

                    try
                        hAxes = eventData.AffectedObject;
                    catch
                        hAxes = ancestor(eventData.Source,'Axes');
                    end
                    tickValues = get(hAxes,[axName 'Tick']);
                    tickLabels = arrayfun(@(x)(sprintf(format,x)),tickValues,'UniformOutput',false);
                    set(hAxes,[axName 'TickLabel'],tickLabels);

                end

            end % install_adjust_ticklbl

        end % tickLabelFormat

        function hText = xtickLabelRotate(XTick, rot, varargin)
%------------------------------------------------------------------------------
% xtick label rotate
%
% Parameter:
%       - XTick. vector array of XTick positions & values (numeric) uses current
%           XTick values or XTickLabel cell array by default (if empty)
%       - rot. angle of rotation in degrees, 90° by default
%       - XTickLabel: cell array of label strings
%       - [var]. Optional. "Property-value" pairs passed to text generator
%           ex: 'interpreter','none', 'Color','m','Fontweight','bold'
%
% Return:
%       - hText. handle vector to text labels
%
% Example 1:  Rotate existing XTickLabels at their current position by 90°
%    xticklabel_rotate
%
% Example 2:  Rotate existing XTickLabels at their current position by 45° and change
% font size
%    xticklabel_rotate([],45,[],'FontSize',14)
%
% Example 3:  Set the positions of the XTicks and rotate them 90°
%    figure;  plot([1960:2004],randn(45,1)); xlim([1960 2004]);
%    xticklabel_rotate([1960:2:2004]);
%
% Example 4:  Use text labels at XTick positions rotated 45° without tex interpreter
%    xticklabel_rotate(XTick,45,NameFields,'interpreter','none');
%
% Example 5:  Use text labels rotated 90° at current positions
%    xticklabel_rotate([],90,NameFields);
%
% Example 6:  Multiline labels
%    figure;plot([1:4],[1:4])
%    axis([0.5 4.5 1 4])
%    xticklabel_rotate([1:4],45,{{'aaa' 'AA'};{'bbb' 'AA'};{'ccc' 'BB'};{'ddd' 'BB'}})
%
% This is a modified version of xticklabel_rotate90 by Denis Gilbert
% Modifications include Text labels (in the form of cell array)
%       Arbitrary angle rotation
%       Output of text handles
%       Resizing of axes and title/xlabel/ylabel positions to maintain same overall size
%          and keep text on plot
%          (handles small window resizing after, but not well due to proportional placement with
%           fixed font size. To fix this would require a serious resize function)
%       Uses current XTick by default
%       Uses current XTickLabel is different from XTick values (meaning has been already defined)
%
% Author: Brian FG Katz, bfgkatz@hotmail.com
%------------------------------------------------------------------------------


            % check to see if xticklabel_rotate has already been here (no other reason for this to happen)
            if isempty(get(gca,'XTickLabel'))
                error('OptAlgo.xtickLabelRotate: can not process, either xticklabel_rotate has already been run or XTickLabel field has been erased \n');
            end

            % Modified with forum comment by "Nathan Pust" allow the current text labels to be used and property value pairs to be changed for those labels
            if (nargin < 3 || isempty(varargin{1})) && (~exist('XTick') || isempty(XTick))

                xTickLabels = get(gca,'XTickLabel');

                if ~iscell(xTickLabels)
                    % remove trailing spaces if exist (typical with auto generated XTickLabel)
                    temp1 = num2cell(xTickLabels,2);
                    for loop = 1:length(temp1),
                        temp1{loop} = deblank(temp1{loop});
                    end
                    xTickLabels = temp1;
                end

                varargin = varargin(2:length(varargin));

            end

            % if no XTick is defined use the current XTick
            if (~exist('XTick') || isempty(XTick))
                XTick = get(gca,'XTick');
            end

            % Make XTick a column vector
            XTick = XTick(:);

            if ~exist('xTickLabels')
                % Define the xtickLabels
                % If XtickLabel is passed as a cell array then use the text
                if ~isempty(varargin) && (iscell(varargin{1}))
                    xTickLabels = varargin{1};
                    varargin = varargin(2:length(varargin));
                else
                    xTickLabels = num2str(XTick);
                end
            end

            if length(XTick) ~= length(xTickLabels)
                error('OptAlgo.xtickLabelRotate: must have same number of elements in "XTick" and "XTickLabel" \n');
            end

            % Set the Xtick locations and set XTicklabel to an empty string
            set(gca,'XTick',XTick,'XTickLabel','');

            if nargin < 2
                rot = 90 ;
            end

            % Determine the location of the labels based on the position of the xlabel
            hxLabel = get(gca,'XLabel');
            xLabelString = get(hxLabel,'String');

            set(hxLabel,'Units','data');
            xLabelPosition = get(hxLabel,'Position');
            y = xLabelPosition(2);

            % CODE below was modified following suggestions from Urs Schwarz
            y = repmat(y,size(XTick,1),1);
            % retrieve current axis' fontsize
            fs = get(gca,'fontsize');

            if ~iscell(xTickLabels)
                % Place the new xTickLabels by creating TEXT objects
                hText = text(XTick, y, xTickLabels,'fontsize',fs);
            else
                % Place multi-line text approximately where tick labels belong
                for cnt=1:length(XTick)
                    hText(cnt) = text(XTick(cnt),y(cnt),xTickLabels{cnt}, ...
                        'VerticalAlignment','top', 'UserData','xtick');
                end
            end

            % Rotate the text objects by ROT degrees
            % Modified with modified forum comment by "Korey Y" to deal with labels at top
            % Further edits added for axis position
            xAxisLocation = get(gca, 'XAxisLocation');
            if strcmp(xAxisLocation,'bottom')
                set(hText,'Rotation',rot,'HorizontalAlignment','right',varargin{:});
            else
                set(hText,'Rotation',rot,'HorizontalAlignment','left',varargin{:});
            end

            % Adjust the size of the axis to accomodate for longest label (like if they are text ones)
            % This approach keeps the top of the graph at the same place and tries to keep xlabel at the same place
            % This approach keeps the right side of the graph at the same place

            set(get(gca,'xlabel'),'units','data');
                labxorigpos_data = get(get(gca,'xlabel'),'position');
            set(get(gca,'ylabel'),'units','data');
                labyorigpos_data = get(get(gca,'ylabel'),'position');
            set(get(gca,'title'),'units','data');
                labtorigpos_data = get(get(gca,'title'),'position');

            set(gca,'units','pixel');
            set(hText,'units','pixel');
            set(get(gca,'xlabel'),'units','pixel');
            set(get(gca,'ylabel'),'units','pixel');

            origpos = get(gca,'position');

            % Modified with forum comment from "Peter Pan" to deal with case when only one XTickLabelName is given.
            x = get( hText, 'extent' );
            if iscell( x ) == true
                textsizes = cell2mat( x );
            else
                textsizes = x;
            end

            largest =  max(textsizes(:,3));
            longest =  max(textsizes(:,4));

            laborigext = get(get(gca,'xlabel'),'extent');
            laborigpos = get(get(gca,'xlabel'),'position');

            labyorigext = get(get(gca,'ylabel'),'extent');
            labyorigpos = get(get(gca,'ylabel'),'position');
            leftlabdist = labyorigpos(1) + labyorigext(1);

            % assume first entry is the farthest left
            leftpos = get(hText(1),'position');
            leftext = get(hText(1),'extent');
            leftdist = leftpos(1) + leftext(1);
            if leftdist > 0, leftdist = 0; end

            % Modified to allow for top axis labels and to minimize axis resizing
            if strcmp(xAxisLocation,'bottom')
                newpos = [origpos(1)-(min(leftdist,labyorigpos(1)))+labyorigpos(1) ...
                    origpos(2)+((longest+laborigpos(2))-get(gca,'FontSize')) ...
                    origpos(3)-(min(leftdist,labyorigpos(1)))+labyorigpos(1)-largest ...
                    origpos(4)-((longest+laborigpos(2))-get(gca,'FontSize'))];
            else
                newpos = [origpos(1)-(min(leftdist,labyorigpos(1)))+labyorigpos(1) ...
                    origpos(2) ...
                    origpos(3)-(min(leftdist,labyorigpos(1)))+labyorigpos(1)-largest ...
                    origpos(4)-(longest)+get(gca,'FontSize')];
            end
            set(gca,'position',newpos);

            % readjust position of text labels after resize of plot
            set(hText,'units','data');
            for loop= 1:length(hText)
                set(hText(loop),'position',[XTick(loop), y(loop)]);
            end

            % adjust position of xlabel and ylabel
            laborigpos = get(get(gca,'xlabel'),'position');
            set(get(gca,'xlabel'),'position',[laborigpos(1) laborigpos(2)-longest 0]);

            % switch to data coord and fix it all
            set(get(gca,'ylabel'),'units','data');
            set(get(gca,'ylabel'),'position',labyorigpos_data);
            set(get(gca,'title'),'position',labtorigpos_data);

            set(get(gca,'xlabel'),'units','data');
                labxorigpos_data_new = get(get(gca,'xlabel'),'position');
            set(get(gca,'xlabel'),'position',[labxorigpos_data(1) labxorigpos_data_new(2)]);

            % Reset all units to normalized to allow future resizing
            set(get(gca,'xlabel'),'units','normalized');
            set(get(gca,'ylabel'),'units','normalized');
            set(get(gca,'title'),'units','normalized');
            set(hText,'units','normalized');
            set(gca,'units','normalized');

            if nargout < 1
                clear hText
            end

        end % xtickLabelRotate

        function res = crashSaver(x, e)

            save('crashParams.dat', 'x', '-ascii', '-append');
            fprintf('%s The parameter set results in this crash has been stored in crashParams.dat\n', e.message);
            res = 68106800;

        end % crashSaver

    end % method


% Upper level algorithm, in charge of discrete structural optimization
% -----------------------------------------------------------------------------
    methods (Static = true, Access = 'public')

        function structure = discreteInit(opt)
% -----------------------------------------------------------------------------
% Generating the initial population of structures
%
% Encoding: each node between two adjacent columns are denoted by a number sequently
%           0 1 2 3 4 5 6 7 8 9 10 ...
% Here 0 represents the desorbent port all the time. As it is a loop, we always need a starting point
% And the sequence of ports are D E (ext_1 ext_2) F R constantly. The selective ranges of each pointer
% (E,F,R) are shown as follows in the binary scenario:
%           0 1 2 3 4 5 6 7 8 9 10 ...
%           D
%             E < ------- > E      : extract_pool
%               F < ------- > F    : feed_pool
%                 R < ------- > R  : raffinate_pool
% -----------------------------------------------------------------------------


            nodeIndex = opt.nColumn -1 ;

            % Preallocate of the structure matrix
            structure = zeros(opt.structNumber,opt.nZone+1);

            for i = 1:opt.structNumber

                if opt.nZone == 4
                    extract_pool = [1, nodeIndex-2];
                    structure(i,2) = randi(extract_pool);

                    feed_pool = [structure(i,2)+1, nodeIndex-1];
                    structure(i,3) = randi(feed_pool);

                    raffinate_pool = [structure(i,3)+1, nodeIndex];
                    structure(i,4) = randi(raffinate_pool);

                elseif opt.nZone == 5
                    extract1_pool = [1, nodeIndex-3];
                    structure(i,2) = randi(extract1_pool);

                    extract2_pool = [structure(i,2)+1, nodeIndex-2];
                    structure(i,3) = randi(extract2_pool);

                    feed_pool = [structure(i,3)+1, nodeIndex-1];
                    structure(i,4) = randi(feed_pool);

                    raffinate_pool = [structure(i,4)+1, nodeIndex];
                    structure(i,5) = randi(raffinate_pool);

                end
            end

        end

        function structID = structure2structID(opt,structure)
% -----------------------------------------------------------------------------
% This is the rountine that decode the structure into structID for simulation
%
% For instance, the structure is [0, 3, 5, 9] in a binary situation with column amount 10
% the structID is [3,2,4,1], there are three in the zone I, two in zone II, four in
% zone III, one in zone IV
% -----------------------------------------------------------------------------


            structID = zeros(1, opt.nZone);

            structID(1:opt.nZone-1) = structure(2:end) - structure(1:end-1);

            structID(end) = opt.nColumn - structure(end);

        end

        function mutant_struct = discreteMutation(opt, structure)
% -----------------------------------------------------------------------------
% This is the mutation part of the upper-level structure optimization algorithm
%
% First of all, two random selected structures are prepared;
% Then the optimal structure until now is recorded;
% Lastly, the mutant_struct = rand_struct_1 &+ rand_struct_2 &+ optima_struct
% -----------------------------------------------------------------------------


            % Record the optimal structure so far and select two random structures
            rand_struct_1 = structure(randi(structNumber), 1:opt.nZone);
            rand_struct_2 = structure(randi(structNumber), 1:opt.nZone);

            [~, id] = min(structure(:,opt.nZone+1));
            optima_struct = structure(id, 1:opt.nZone);

            % Preallocate of the mutation structure
            mutant_struct = zeros(1,nZone);

            for i = 1:opt.nZone
                if rand < 0.33
                    mutant_struct(i) = rand_struct_1(i);
                elseif (0.33 <= rand) && (rand<= 0.67)
                    mutant_struct(i) = rand_struct_2(i);
                else
                    mutant_struct(i) = optima_struct(i);
                end
            end

        end

        function trial_struct = discreteCrossover(opt, structure, mutant_struct)
% -----------------------------------------------------------------------------
% The crossover part of the upper-level structure optimization algorithm
%
% The generation of a new trial structure is achieved by randomly combining the original
% structure and mutation structure.
% -----------------------------------------------------------------------------


            trial_struct = zeros(1, opt.nZone);

            for i = 1:opt.nZone
                if rand < 0.5
                    trial_struct(i) = structure(i);
                else
                    trial_struct(i) = mutant_struct(i);
                end
            end

        end

    end % upper level


end % classdef
% =============================================================================
%  SMB - The Simulated Moving Bed Chromatography for separation of
%  target compounds, either binary or ternary.
%
%      Copyright © 2008-2017: Eric von Lieres, Qiaole He
%
%      Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
%
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
