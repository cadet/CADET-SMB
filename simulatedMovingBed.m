function objective = simulatedMovingBed(varargin)
% =============================================================================
% Main function, which is in charge of switching to achieve CSS
%
% Binary separations (conventional four-zone scheme) and ternary separations
% (cascade, five-zone or eight zone schemes) are available in this program.
% As for the specific column configurations, I refer your to our paper.
% =============================================================================


    global string stringSet structID;

    tTotal = tic;

    % Generate alphabet strings to name columns
    stringSet = SMB.stringGeneration();

    % Read operating parameters of SMB unit
    [opt, interstVelocity, Feed, Desorbent] = getParameters(varargin{:});

    % Check interstitial velocities in SMB optimization
    % flag is true when anyone is negative, then assign a big value to objective function
    flag = SMB.interstVelocityCheck(interstVelocity, opt);
    if flag == 1, objective = 68106800; return; end

    if ~isfield(opt, 'structID'), opt.structID = structID; end

    if opt.enable_CSTR && opt.enable_DPFR
        error('It is not allowed have both the CSTR and DPFR in the simulation \n');
    end

    if opt.nColumn > length(stringSet)
        error('The simulation of %3g-column case in %3g-zone is not finished so far \n', opt.nColumn, opt.nZone);
    end

%   Preallocation
    [currentData, sequence, convergIndx, convergPrevious, initialState, initialState_DPFR, dyncData] ...
        = SMB.preallocation(Feed, opt);

%   Simulations
%----------------------------------------------------------------------------------------
    % Interactive plotting when debug mode is enabled
    SMB.UIplot('head', opt);

    % Main loop
    for i = 1:opt.nMaxIter

        % Switch implementation by means of attaching the analog column to corresponding position
        if mod(i, opt.nColumn) ~= 0
            k = string(mod(i, opt.nColumn));
        else
            k = string(opt.nColumn);
        end

        try
            % Single column simulation
            [outletProfile, lastState, initialState_DPFR] = SMB.colSolver(initialState,...
                initialState_DPFR, currentData, interstVelocity, Feed, Desorbent, opt, sequence, k, varargin{:});
        catch e
            fprintf('%s\n', e.message); objective = 68106800; return
        end

        % Store chromatogram and column state
        currentData{sequence.(k)}.outlet   = outletProfile.outlet;
        currentData{sequence.(k)}.colState = outletProfile.column;
        initialState = lastState;
        % If eight-zone scheme, store the intermediate outlet as the Feed2 inlet
        SMB.Feed2Connect(outletProfile, stringSet, opt, k);

        if opt.enableDebug
            % The collection of the dyncData for the trajectory plotting
            dyncData = SMB.dyncDataUpdate(dyncData, currentData, sequence, stringSet, opt, i);
        end

        % stopping criterion
        % flag is true when a stopping criterion is satisfied
        [convergPrevious, flag] = SMB.stoppingCriterion(convergPrevious, currentData, dyncData, convergIndx, opt, i);
        if flag == 1, break; end

    end % main loop

%   Post-process
%----------------------------------------------------------------------------------------
    % Compute the performance index, such Purity and Productivity
    Results = SMB.Purity_Productivity(currentData, opt);

    % Construct your own Objective Function and calculate the value
    objective = SMB.objectiveFunction(Results, opt);

    tTotal = toc(tTotal);
    % Store the final data into DATA.mat file when debug is active
    if opt.enableDebug
        fprintf('The time elapsed for reaching the Cyclic Steady State: %g sec \n', tTotal);
        SMB.concDataConvertToASCII(currentData, opt);
        SMB.trajDataConvertToASCII(dyncData, opt);
        c = clock;
        save(sprintf('Performance_%d%d.mat', c(3), c(4)), 'Results');
        fprintf('The results about concentration profiles and the trajectories have been stored \n');
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
