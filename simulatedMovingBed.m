function objective = simulatedMovingBed(varargin)
% =============================================================================
% Main function, which is in charge of switching to achieve CSS
%
% Binary separations (conventional four-zone scheme) and ternary separations
% (cascade, five-zone or eight zone schemes) are available in this program.
% As for the specific column configurations, I refer your to our paper.
% =============================================================================


    global string stringSet dummyProfile;

    tTotal = tic;

    % Generate alphabet strings to name columns
    stringSet = SMB.stringGeneration();

    % Read operating parameters of SMB unit
    [opt, interstVelocity, Feed, Desorbent] = getParameters(varargin{:});

    % If SMA isotherm, salt is considered as a component
    if strcmp('StericMassActionBinding', opt.BindingModel), opt.nComponents = opt.nComponents + 1; end

    % Check interstitial velocities in SMB optimization
    % flag is true when anyone is negative, then assign a big value to objective function
    flag = SMB.interstVelocityCheck(interstVelocity, opt);
    if flag == 1, objective = 68106800; return; end

    if opt.enable_CSTR && opt.enable_DPFR
        error('It is not allowed have both the CSTR and DPFR in the simulation \n');
    end

    if opt.nColumn > length(stringSet)
        error('The simulation of %3g-column case in %3g-zone is not finished so far \n', opt.nColumn, opt.nZone);
    end

%   Preallocation
    [currentData, tempData, transitionData, sequence, convergIndx, convergPrevious, ...
        dummyProfile, dyncData, plotData] = SMB.preallocation(Feed, opt);

%   Simulations
%----------------------------------------------------------------------------------------
    % Interactive plotting when debug mode is enabled
    SMB.UIplot('head', opt);

    % Main loop
    for i = 1:opt.nMaxIter

        % Switch implementation by means of attaching different column to corresponding position
        sequence = cell2struct( circshift( struct2cell(sequence),-1 ), stringSet(1:opt.nColumn) );

        simMex = cell(1, opt.nColumn);

        for j = 1:opt.nInterval

            % Do nColumn simulations in terms of string sequence
            for k = string'

                try
                    % Single column simulation
                    [outletProfile, lastState, currentData] = SMB.colSolver(currentData, ...
                        interstVelocity, Feed, Desorbent, opt, sequence, k, j, simMex, varargin{:});
                catch e
                    % If the simulation of single column crashes
                    fprintf('%s\n', e.message); objective = 68106800; return
                end

                % For the very first column, only the outletProfile, rather than the column state, is stored
                if ~strcmp(k, string(1))
                    currentData{sequence.(k)}.lastState = lastState;
                    if j == opt.nInterval
                        currentData{sequence.(k)}.colState = outletProfile.column;
                    end
                end
                currentData{sequence.(k)}.outlet = outletProfile.outlet;
                % Store the interval-wise concentration profile of all columns
                tempData{sequence.(k)}.concentration{j} = outletProfile.outlet.concentration;
                % If eight-zone case, store the intermediate outlet as the Feed2 inlet
                SMB.Feed2Connect(outletProfile, Feed, interstVelocity, stringSet, opt, k);

            end % for k = string'

            % Simulate the very first column again to update the actual column state, correspondingly the outletProfile is omitted
            [currentData{sequence.(string(1))}.colState, currentData{sequence.(string(1))}.lastState] = ...
                SMB.observerSimulation(sequence, interstVelocity, Feed, Desorbent, currentData{sequence.(string(end))}.outlet, ...
                currentData{sequence.(string(1))}.lastState, string(1), opt);

            if opt.enableDebug
                % The collection of the dyncData for the trajectory plotting
                dyncData = SMB.dyncDataUpdate(dyncData, currentData, sequence, stringSet, opt, i, j);
            end

            % Interactive plottings when debug mode is active
            SMB.UIplot('interval', opt, i, j);

            % Plot the dynamic trajectory
            if opt.enableDebug, SMB.plotDynamic(opt, dyncData(:,1:j+(i-1)*opt.nInterval), j+(i-1)*opt.nInterval); end

        end % for j = 1:opt.nInterval

        % Transfer the dummyProfile, which is the profile of last interval of last column
        dummyProfile = tempData{sequence.(string(end))};

        % Plot the outlet profile of each column in nColumn switches
        if opt.enableDebug, SMB.plotFigures(opt, currentData); end

        % Cat nInterval pieces of tempData into transitionData,
        % which is the internal column state after nColumn switches
        for ii = 1:opt.nColumn
            transitionData{ii}.outlet.concentration = cat( 1, tempData{ii}.concentration{:} );
        end

        % Store the instant outlet profiles of all columns in one iteration
        index = mod(i, opt.nColumn);
        if index == 0
            plotData(:,opt.nColumn) = transitionData';
        else
            plotData(:,index) = transitionData';
        end

        % stopping criterion
        % flag is true when a stopping criterion is satisfied
        [convergPrevious, flag] = SMB.stoppingCriterion(convergPrevious, transitionData, convergIndx, opt, i, j);
        if flag == 1, break; end

    end % main loop

%   Post-process
%----------------------------------------------------------------------------------------
    % Compute the performance index, such Purity and Productivity
    Results = SMB.Purity_Productivity(plotData, varargin{:});

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
