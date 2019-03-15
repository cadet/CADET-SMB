function result = smbStructOpt
% =============================================================================
% This is the main function in charge of the discrete structual optimization
% There are two levels in this algorithm
%       - upper level: in charge of the discrete structure evolution
%       - lower level: in charge of the decision variables optimization under certain
%       column configuration transferred from the upper level
%
% In the upper level, the algorithm is fixed right now, but in the lower there are four
% different types of algorithms are accessible, DE, PSO, MADE, fmincon
% =============================================================================


    global structID;

    % Assign the optimized decision variables under a certain column configuration
%    params = struct('columnLength',[], 'switch',[], 'recycle',[], 'feed',[], 'desorbent',[], 'extract',[]); % binary scenario
    params = struct('columnLength',[], 'switch',[], 'recycle',[], 'feed',[], 'desorbent',[], 'extract1',[], 'extract2',[]); % ternary scenario

    % There are four optimization algorithms availabe in the lower level programme
    Algos = struct('PSO',[false], 'DE',[true], 'MCMC',[false], 'MADE',[false], 'PRIMAL',[false], 'fmincon',[]);

    [opt,~,~] = getParameters( zeros(1,length(fieldnames(params))) );

    % Give the manual parameter boundaries, internal points searching
    opt.paramBound = [0.05 0.15; 250 350; 2.5e-7 3.8e-7; 1.5e-8 3.0e-8; 2.0e-7 3.5e-7; 2.0e-7 3.5e-7; 4.0e-8 5.5e-8];
    opt.params = params;
    if Algos.fmincon, opt.initParams = [0.25, 180, 9.62e-7, 0.98e-7, 1.96e-7, 1.54e-7]; end

    % Initilize the structure population and calculate the fitness value of each structure
    structure = OptAlgo.discreteInit(opt);
    decision_variables = zeros(opt.structNumber, length(fieldnames(params)));

    for i = 1:opt.structNumber

        structID = OptAlgo.structure2structID(opt, structure(i, 1:opt.nZone));
        [decision_variables(i, :), structure(i, opt.nZone+1)] = smbOperatOpt(opt, Algos);

    end


    maxIter = 20;
%-----------------------------------------------------------------------------------------
    % Main loop
    for k = 1:maxIter

        for i = 1:opt.structNumber

            % Mutation
            mutant_struct = OptAlgo.discreteMutation(opt, structure);

            % Crossover
            trial_struct = OptAlgo.discreteCrossover(opt, structure(i,1:opt.nZone), mutant_struct);

            % Simulation
            structID = OptAlgo.structure2structID(opt, trial_struct);
            [paramValue, trial_value] = OptAlgo.continuousUnitOptimization(opt, params, Algos);

            % Selection
            if trial_value < structure(i, opt.nZone+1)
                structure(i,1:opt.nZone) = trial_struct;
                decision_variables(i, :) = paramValue;
            end

        end


        [minValue, idStruct] = min(structure(:,opt.nZone+1));

        fprintf('================  Iter(Upper): %5d     Minimum: %10.3g  ================ \n', k, minValue);
        fprintf('Structure:'); fprintf('%d |',structure(idStruct, 1:opt.nZone));
        fprintf('\n---------------------------------------------------------------------- \n');
        fprintf('%10.3g | ', decision_variables(idStruct, :)); fprintf('\n');


        % Stopping criterion
        if mod(k, 5) == 0

            delta = std(structure(:, 1:opt.nZone)) ./ mean(structure(:, 1:opt.nZone));

            if all(abs(delta) < 0.01) || k == maxIter
                break
            end

        end

    end
%-----------------------------------------------------------------------------------------


    result.objective = minValue;
    result.structure = structure(idStruct, 1:opt.nZone);
    save(sprintf('result_%2d.mat',fix(rand*100)),'result');
    fprintf('The results have been stored in the result.mat \n');


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
