function Particle_Swarm_Optimization(params)

% =============================================================================
% Particle Swarm Optimization algorithm (PSO)

%  PSO optimizes a problem by having a population of candidates
%  (particles), and moving these particles around in the search space
%  according to mathematical formula ovet the particle's position and
%  velocity. Each particle's movement is influenced by its own local
%  best-known position but, is also guided toward the best-known positions
%  in the search space, which are updated as better positions are found by
%  other particles

 
%  Parameter:
%        - params. It is the specified parameters from the main function.
%        And let the main function informed that which parameters are need
%        to be optimized
% =============================================================================


    startTime = clock;


%   Get the optimizer options
    options = getOptions_PSO(params);  
    
%   Initialization of the population
    [ParSwarm, OptSwarm, ToplOptSwarm] = InitSwarm(options);

        
%-----------------------------------------------------------------------------------------
%   The main loop
    for k = 1:options.loopCount
        
%       The evolution of the particles in PSO
        [ParSwarm, OptSwarm, ToplOptSwarm] = ParticlesEvolution(ParSwarm, OptSwarm, ToplOptSwarm, k, options);
      
%       Abstract best information so far from the population and display it
        XResult = exp( OptSwarm(options.swarmSize+1, 1:options.particleSize) );
        YResult = OptSwarm(options.swarmSize+1, options.particleSize+1);        
        
        
        fprintf('Iter = %3d,  Minimum: %g,  Parameters:[ %g| %g| %g| %g| %g| %g ] \n'...
            , k, YResult, XResult);
        
        
%       The convergence criterion: when the standard deviation is smaller
%       than the specified tolerance
        if mod(k, 5) == 0

            delta = std(ParSwarm(:,1:options.particleSize)) ./ mean(ParSwarm(:,1:options.particleSize));
            
            if all(abs(delta) < 0.005) || k == options.loopCount
                break
            end
            
        end
                     
        
    end
%-----------------------------------------------------------------------------------------  

%   Gather some useful information and store them
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

%  Parameter:
%       - params. It is the specified parameters from the main function.
%        And let the main function informed that which parameters are need
%        to be optimized.

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


%   The number of population of the swarm intellectual (SI)
    opt.swarmSize      = 20;
    
%   The dimension of the optimized parameters
    opt.particleSize   = length(fieldnames(params));
    
%   the row represents the parameter, while the column denotes the upbound and lowerbound
    opt.paramsRange    = log( [0.20    0.30;
                               150     230;
                               8.0e-7  10e-7;
                               0.9e-7  2.0e-7;
                               0.7e-7  2.0e-7;
                               1.0e-7  2.0e-7] );   
    opt.loopCount      = 300;
    
%   check out the dimension of the set of parameters, and the boundary limitation
    [row, col] = size(opt.particleSize);
    if row > 1 || col > 1
        error('The initialized dimension of the set of parameters might be wrong');
    end
    
    [row, col] = size(opt.paramsRange);
    if row ~= opt.particleSize && col ~= 2
        error('Please check your setup of the range of parameters');
    end

%   Those are the specific parameters for the PSO algorithm
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

%  Return:
%       - ParSwarm. The swarm of the particles, correspondingly the objective value
%       - OptSwarm. The local optima that each particle ever encountered
%       - ToplOptSwarm. The global optima around the population that is
%           shared by the rest particles.
% -----------------------------------------------------------------------------


%   Check out the variables which are needed in this funciton
    if nargout < 3
        error('There are not enough output in the subroutine InitSwarm');
    end

   
%   Initilization of the parameters, velocity of parameters   
    ParSwarm = rand(opt.swarmSize, 2 * opt.particleSize + 1);
    
%   Use vectorization to speed up, it's actually equivalent to the for-loop  
    ParSwarm(:,1:opt.particleSize) = repmat(opt.paramsRange(:,1)',opt.swarmSize,1) + ...
        ParSwarm(:,1:opt.particleSize) .* repmat( (opt.paramsRange(:,2) - opt.paramsRange(:,1))',opt.swarmSize,1 );    
   

%   Simulation of the sampled points, using subroutine simulatedMovingBed
    ParSwarm(: ,2 * opt.particleSize + 1) = arrayfun( @(idx) feval(@simulatedMovingBed,...
        exp(ParSwarm(idx, 1:opt.particleSize)) ), 1:opt.swarmSize );

%   The statistics of the population
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

% Parameters:
%       - ParSwarm. The swarm of the particles, correspondingly the objective value
%       - OptSwarm. The local optima that each particle ever encountered
%       - ToplOptSwarm. The global optima around the population that is
%           shared by the rest particles.
%       - CurCount. The iteration number in the main loop
%       - opt. Please see the comments of the function, InitSwarm

% Return:
%       - ParSwarm. The swarm of the particles, correspondingly the objective value
%       - OptSwarm. The local optima that each particle ever encountered
%       - ToplOptSwarm. The global optima around the population that is
%           shared by the rest particles.
% -----------------------------------------------------------------------------


    if nargin ~= 5
        error('There are error in the input of the function ParticlesEvolution')
    end
    if nargout ~= 3 
        error('There are error in the output of the function ParticlesEvolution')
    end
    

%   Different strategies of the construction of the weight
    weight = opt.wMax - CurCount * ((opt.wMax - opt.wMin) / opt.loopCount);
%   w=0.7;
%   w=(MaxW-MinW)*(CurCount/LoopCount)^2+(MinW-MaxW)*(2*CurCount/LoopCount)+MaxW;
%   w=MinW*(MaxW/MinW)^(1/(1+10*CurCount/LoopCount));


%   Different options of the accelaration coefficients
    c1 = opt.accelCoeff(1); 
    c2 = opt.accelCoeff(2);
    
%   c1=2.8;
%   c2=1.3;

%   con=1;
%   c1=4-exp(-con*abs(mean(ParSwarm(:,2*ParCol+1))-AdaptFunc(OptSwarm(ParRow+1,:))));
%   c2=4-c1;


%   LocalOptDiff: Each particle compares the difference of current position and the best
%   position it ever encountered
    [ParRow, ParCol] = size(ParSwarm);
    ParCol = (ParCol - 1) / 2;
    LocalOptDiff = OptSwarm(1:ParRow, 1:ParCol) - ParSwarm(:, 1:ParCol);    


    for row = 1:ParRow
        
%       GlobalOptDiff: Each particle compares the difference of current position and the best
%       position it ever encountered
        if strcmp(opt.topology, 'without Topology')
            GlobalOptDiff = OptSwarm(ParRow+1, 1:ParCol) - ParSwarm(row, 1:ParCol);
        end
        
        if strcmp(opt.topology, 'Random Topology') || strcmp(opt.topology, 'Ring Topology')
            GlobalOptDiff = ToplOptSwarm(row, 1:ParCol) - ParSwarm(row, 1:ParCol);
        end

%       The evolution of the velocity, according to the LocalOptDiff and GlobalOptDiff   
%       TempVelocity = weight .* ParSwarm(row,:) + c1 * unifrnd(0,1.0) .* LocalOptDiff(row,:)...
%           + c2 * unifrnd(0,1.0) .* GlobalOptDiff;
        TempVelocity = weight .* ParSwarm(row, ParCol+1 : 2*ParCol) + c1 *LocalOptDiff(row, :) + c2 * GlobalOptDiff;
        
        
%       check the boundary limitation of the Velocity
        velocityBound = (opt.paramsRange(:,2)-opt.paramsRange(:,1) )';
        
        [~, veCol] = find( abs(TempVelocity(:, 1:ParCol)) > velocityBound  );
        TempVelocity(veCol) = velocityBound(veCol);
        
%       Value Assignment: store the updated Velocity into the matrix ParSwarm
        ParSwarm(row, ParCol+1 : 2*ParCol) = TempVelocity;

        
        
        step_size = 1;
%       step_size = 0.729;

%       The evolution of the current positions, according to the Velocity and the step size
        ParSwarm(row, 1:ParCol) = ParSwarm(row, 1:ParCol) + step_size * TempVelocity;


%       Check the boundary limitation
        loBound = repmat(opt.paramsRange(:,1)', ParRow, 1);
        [loRow, loCol] = find( (ParSwarm(1:ParRow, 1:ParCol) - loBound) < 0 );
        ParSwarm((loCol-1).*ParRow + loRow) = loBound((loCol-1).*ParRow + loRow);
        
        upBound = repmat(opt.paramsRange(:,2)', ParRow, 1);
        [upRow, upCol] = find( (ParSwarm(1:ParRow, 1:ParCol) - upBound) > 0 );
        ParSwarm((upCol-1).*ParRow + upRow) = upBound((upCol-1).*ParRow + upRow);
        clear loRow loCol upRow upCol veCol;
    
        
%       Simulation of the sampled points, using subroutine simulatedMovingBed
        ParSwarm(row, 2*ParCol+1) = feval( @simulatedMovingBed, exp(ParSwarm(row, 1:ParCol)) );
 
%       if the updated position is better than the current position, the
%       particle flies to the updated positon; otherwise, keep still in
%       current position
%       Update the LocalOpt for each particle
        if ParSwarm(row, 2*ParCol+1) < OptSwarm(row, ParCol+1)
            OptSwarm(row, 1:ParCol) = ParSwarm(row, 1:ParCol);
            OptSwarm(row, ParCol+1) = ParSwarm(row, 2*ParCol+1);
        end

%       Usint different type of the topology between particles, the global
%       optimal position is transferred to the swarm with different
%       strategy
%       Update the GlobalOpt around the whole group
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

    
%   Statistics
    [minValue, row] = min(arrayfun(@(idx) OptSwarm(idx, ParCol+1), 1: ParRow));
    
    OptSwarm(ParRow+1, 1:ParCol) = OptSwarm(row,1:ParCol);
    OptSwarm(ParRow+1, ParCol+1) = minValue; 

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
