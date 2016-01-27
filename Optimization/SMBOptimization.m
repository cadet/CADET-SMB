function SMBOptimization()

% =============================================================================
% This is the main function of the optimization of the Simulated Moving Bed
% The optimized parameters in this case are 
%       - columnLength 
%       - switchTime 
%       - flowRates_recycle 
%       - flowRate_feed
%       - flowRate_desorbent
%       - flowRate_extrct
% 
%       theta = {L_c, t_s, Q_{re}, Q_F, Q_D, Q_E}
% 
% There are four types of algorithms are integrated into this code, either
% based on Heuristical theory or Deterministic theory, either optimization or sampling.
%       - Particle Swarm Optimizatio (PSO)
%       - Differential Evolution (DE)
%       - Metropolis Adjusted Differential Evolution (MADE)
%       - Riemann Manifold Metropolis Adjusted Langevin Algorithm (MLA)
% 
% =============================================================================



%   There are four optimization algorithms are availabe in this programme
    optimization_method = struct('Particle_Swarm_Optimization',[], 'Differential_Evolution',[],...
       'Metropolis_Adjusted_Differential_Evolution',[], 'Riemann_Manifold_Metropolis_Adjusted_Langevin',[],...
       'Deterministic_algorithm_fmincon',[]);

%   The set of the parameters which are optimized
    params = struct('columnLength',[], 'switch',[], 'recycle',[], 'feed',[], 'desorbent',[], 'extract',[]);

    
%   Select one method and make it true (correspondingly the rest methods false)
    optimization_method.Particle_Swarm_Optimization = true;
    optimization_method.Differential_Evolution = false;
    optimization_method.Metropolis_Adjusted_Differential_Evolution = false;
    optimization_method.Deterministic_algorithm_fmincon = false;


    if isfield(optimization_method, 'Particle_Swarm_Optimization') ...
            && optimization_method.Particle_Swarm_Optimization
        
        OptAlgorithms.Particle_Swarm_Optimization(params);
        
    elseif isfield(optimization_method, 'Differential_Evolution') ...
            && optimization_method.Differential_Evolution
        
        OptAlgorithms.Differential_Evolution(params);
        
    elseif isfield(optimization_method, 'Metropolis_Adjusted_Differential_Evolution') ...
            && optimization_method.Metropolis_Adjusted_Differential_Evolution
        
        OptAlgorithms.Metropolis_Adjusted_Differential_Evolution(params);
        
%     elseif isfield(optimization_method, 'Riemann_Manifold_Metropolis_Adjusted_Langevin') ...
%             && optimization_method.Riemann_Manifold_Metropolis_Adjusted_Langevin
%         
%         Riemann_Manifold_Metropolis_Adjusted_Langevin(params);

    elseif isfield(optimization_method, 'Deterministic_algorithm_fmincon') ...
            && optimization_method.Deterministic_algorithm_fmincon
        
        initParams = [0.25, 180, 9.62e-7, 0.98e-7, 1.96e-7, 1.54e-7];

        loBound = [0.20, 150, 8.0e-7, 0.9e-7, 0.7e-7, 1.0e-7];
        upBound = [0.30, 230, 10e-7,  2.0e-7, 2.0e-7, 2.0e-7];
        
        options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'iter',...
            'TolX',1e-6,'TolCon',1e-6,'TolFun',1e-6,'MaxIter',500);
        
        try
            [SMBparams, fval, exitflag, output, ~, grad] = fmincon( @simulatedMovingBed, ...
                initParams, [],[],[],[], loBound, upBound, [], options);
        catch exception
            disp('Errors in the MATLAB build-in optimizer: fmincon. \n Please check your input parameters and run again. \n');
            disp('The message from fmincon: %s \n', exception.message);
        end
        
        fprintf('Minimum: %g,   Parameters:[%g| %g| %g| %g| %g| %g] \n', fval, SMBparams);
        
    else
        
        warning('The method you selected is not provided in this programme');
        
    end
    
    
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
