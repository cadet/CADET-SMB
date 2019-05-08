function objective = simulatedMovingBed(varargin)
% =============================================================================
% Main function, which is in charge of switching to achieve CSS
%
% Binary separations (conventional four-zone scheme) and ternary separations
% (cascade, five-zone or eight zone schemes) are available in this program.
% As for the specific column configurations, I refer your to our paper.
% =============================================================================


    global string stringSet structID dummyProfile;

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

    if ~isfield(opt, 'structID'), opt.structID = structID; end

    if opt.enable_CSTR && opt.enable_DPFR
        error('It is not allowed have both the CSTR and DPFR in the simulation \n');
    end

    if opt.nColumn > length(stringSet)
        error('The simulation of %3g-column case in %3g-zone is not finished so far \n', opt.nColumn, opt.nZone);
    end

%   Preallocation
    [currentData, sequence, convergIndx, convergPrevious, dummyProfile, dyncData, plotData] ...
        = SMB.preallocation(Feed, opt);

%   Simulations
%----------------------------------------------------------------------------------------
    % Interactive plotting when debug mode is enabled
    SMB.UIplot('head', opt);

    % Main loop
    for i = 1:opt.nMaxIter

        % Switch implementation by means of attaching different column to corresponding position
        sequence = cell2struct( circshift( struct2cell(sequence),-1 ), stringSet(1:opt.nColumn) );

        % Do nColumn simulations in terms of string sequence
        for k = string'

            try
                % Single column simulation
                [outletProfile, lastState, currentData] = SMB.colSolver(currentData, ...
                    interstVelocity, Feed, Desorbent, opt, sequence, k, varargin{:});
            catch e
                % If the simulation of single column crashes
                fprintf('%s\n', e.message); objective = 68106800; return
            end

            % Store chromatogram and column state
            currentData{sequence.(k)}.outlet     = outletProfile.outlet;
            currentData{sequence.(k)}.colState   = outletProfile.column;
            currentData{sequence.(k)}.lastState  = lastState;
            % If eight-zone case, store the intermediate outlet as the Feed2 inlet
            SMB.Feed2Connect(outletProfile, Feed, interstVelocity, stringSet, opt, k);

        end % for k = string'

        % Transfer the dymmyProfile, which is the profile of last column in terms of string
        dummyProfile = outletProfile.outlet;

        if opt.enableDebug
            % The collection of the dyncData for the trajectory plotting
            dyncData = SMB.dyncDataUpdate(dyncData, currentData, sequence, stringSet, opt, i);

            % Plot internal profile of columns after nColumn simulations
            SMB.plotFigures(opt, currentData);
            % Plot dynamic trajectories at withdrawn ports
            SMB.plotDynamic(opt, dyncData(:,1:i), i);
        end

        % Store the instant outlet profiles of all columns at each switching period
        index = mod(i, opt.nColumn);
        if index == 0
            plotData(:,opt.nColumn) = currentData';
        else
            plotData(:,index) = currentData';
        end

        % stopping criterion
        % flag is true when a stopping criterion is satisfied
        [convergPrevious, flag] = SMB.stoppingCriterion(convergPrevious, currentData, convergIndx, opt, i);
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
