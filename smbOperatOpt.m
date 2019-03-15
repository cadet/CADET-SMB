function [xValue, yValue] = smbOperatOpt(varargin)
% =============================================================================
% This is the main function of the optimization of the Simulated Moving Bed unit
% The optimized parameters in this case are
%       - columnLength
%       - switchTime
%       - flowRates_recycle
%       - flowRate_feed
%       - flowRate_desorbent
%       - flowRate_extract
%
%       theta = {L_c, t_s, Q_{re}, Q_F, Q_D, Q_E}
%
% In the FIVE-ZONE, the optimized parameters are
%       - columnLength
%       - switchTime
%       - flowRates_recycle
%       - flowRate_feed
%       - flowRate_desorbent
%       - flowRate_extract_1
%       - flowRate_extract_2
%
%       theta = {L_c, t_s, Q_{re}, Q_F, Q_D, Q_{E1}, Q_{E2}}
%
% There are four types of algorithms are integrated into this code, either
% based on Heuristical theory or Deterministic theory, either optimization or sampling.
%       - Particle Swarm Optimizatio (PSO)
%       - Differential Evolution (DE)
%       - Markov Chain Monte Carlo (MCMC)
%       - Metropolis Adjusted Differential Evolution (MADE)
%       - Parallel Riemann Metropolis Adjusted Langevin Algorithm (PRIMAL)
% =============================================================================


    if isempty(varargin)

        % The set of the parameters which are optimized
        params = struct('columnLength',[], 'switch',[], 'recycle',[], 'feed',[], 'desorbent',[], 'extract',[]); % binary scenario

        % There are four optimization algorithms availabe in this programme
        Algos = struct('PSO',[false], 'DE',[true], 'MCMC',[false], 'MADE',[false], 'PRIMAL',[false], 'fmincon',[false]);

        [opt,~,~] = getParameters( zeros(1,length(fieldnames(params))) );

        % The initial boundary of parameters: In the format of [x^1_min x^1_max; ...]
        opt.paramBound = [0.05 0.15; 250 350; 2.5e-7 3.8e-7; 1.5e-8 3.0e-8; 2.0e-7 3.5e-7; 2.0e-7 3.5e-7; 4.0e-8 5.5e-8];
        opt.params = params;
        % This is the demonstration case for the binary separation under FOUR-ZONE,
        %    in which 6 decision variables are optimized.
        initParams = [0.25, 180, 9.62e-7, 0.98e-7, 1.96e-7, 1.54e-7];

    else

        opt = varargin{1};
        Algos = varargin{2};
        initParams = opt.initParams;

    end


    if isfield(Algos, 'PSO') && Algos.PSO

        [xValue, yValue] = OptAlgo.Particle_Swarm_Optimization(opt);

    elseif isfield(Algos, 'DE') && Algos.DE

        [xValue, yValue] = OptAlgo.Differential_Evolution(opt);

    elseif isfield(Algos, 'MCMC') && Algos.MCMC

        [xValue, yValue] = OptAlgo.Markov_Chain_Monte_Carlo(opt);

    elseif isfield(Algos, 'MADE') && Algos.MADE

        [xValue, yValue] = OptAlgo.Metropolis_Adjusted_Differential_Evolution(opt);

    elseif isfield(Algos, 'fmincon') && Algos.fmincon

        loBound = opt.paramBound(:,1);
        upBound = opt.paramBound(:,2);

        options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'iter',...
            'TolX',1e-6,'TolCon',1e-6,'TolFun',1e-6,'MaxIter',500);

        try
            [xValue, yValue, exitflag, output, ~, grad] = fmincon( @SMB.simulatedMovingBed, ...
                initParams, [],[],[],[], loBound, upBound, [], options);
        catch exception
            disp('Errors in the MATLAB build-in optimizer: fmincon. \n Please check your input parameters and run again. \n');
            disp('The message from fmincon: %s \n', exception.message);
        end

    else

        error('smbOperatOpt: The method you selected is not provided in this programme \n');

    end

    if opt.enableDebug
        fprintf('----------------  Minimum: %10.3g  ---------------- \n', yValue);
        fprintf('%10.3g | ', xValue);
        fprintf('\n------------------------------------------------------ \n');
    end

end
% =============================================================================
%  SMB - The Simulated Moving Bed Chromatography for separation of
%  target compounds, either binary or ternary.
%
%      Copyright © 2008-2019: Eric von Lieres, Qiaole He
%
%      Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
%
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
