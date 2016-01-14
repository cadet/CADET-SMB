function objective = simulatedMovingBed(ParSwarm)

% =============================================================================
% This is the main function which is charge of switching to reach the
% cyclic steady state. The layout of the columns and the configuration is
% listed as follow:

%              8-column SMB                                       4-column SMB
% Extract                            Feed       |    Extract                         Feed
%       \                            /          |         \                          /
%        --------Zone II(c/d)--------           |          --------Zone II(b)--------
%        |                          |           |          |                        | 
% Zone I(b/a)                    Zone III(e/f)  |     Zone I(a)               Zone III(c)
%        |                          |           |          |                        | 
%        --------Zone IV(h/g)--------           |          --------Zone IV(d)--------
%       /                            \          |         /                          \
% Desorbent                         Raffinate   |   Desorbent                       Raffinate

% Fluid phase goes from Zone I to Zone II to Zone III, while the ports switch direction
% is from Zone I to Zone IV to Zone III;
% =============================================================================
  

    global string;

    tTotal = tic;
    
    if length(ParSwarm) ~= 6
       error('The dimension of the parameter transferred from the optimizer does not match'); 
    end
    
    [opt, interstVelocity, Feed] = getParameters(ParSwarm);

%   Initialize the starting points, currentData
    currentData = cell(1, opt.nColumn);  
    for k = 1:opt.nColumn
        currentData{k}.outlet.time = linspace(0, opt.switch, opt.timePoints);
        currentData{k}.outlet.concentration = zeros(length(Feed.time), 2); 
        currentData{k}.lastState = []; 
    end
    
%   Numbered the columns for the sake of plotting
    if opt.nColumn == 4
        
        sequence.a = 4; sequence.b = 1; sequence.c = 2; sequence.d = 3;
        
        string = char('a','b','c','d');
        convergIndx = 3;
        
    elseif opt.nColumn == 8
        
        sequence.a = 8; sequence.b = 1; sequence.c = 2; sequence.d = 3;
        sequence.e = 4; sequence.f = 5; sequence.g = 6; sequence.h = 7;
        
        string = char('a','b','c','d','e','f','g','h');
        convergIndx = 5;
        
    else
        warning('The simulation of %3g_column case is not finished so far', opt.nColumn);
    end
        
%   prellocation
    plotData = cell(opt.nColumn,opt.nColumn);     
%   convergPrevious is used for stopping criterion
    convergPrevious = currentData{convergIndx}.outlet.concentration;
   
%-----------------------------------------------------------------------------------------
%   Main loop 
    for i = 1:opt.nMaxIter 

        if opt.nColumn == 4
            mid = sequence.a ; sequence.a = sequence.b; sequence.b = sequence.c; ...
                sequence.c = sequence.d; sequence.d = mid;
        elseif opt.nColumn == 8
            mid = sequence.a ; sequence.a = sequence.b; sequence.b = sequence.c; ...
                sequence.c = sequence.d; sequence.d = sequence.e; sequence.e = sequence.f;
                sequence.f = sequence.g; sequence.g = sequence.h; sequence.h = mid;
        end
        
%       The simulation of four columns by the sequence, say, 'a', 'b', 'c', 'd'
        for k = 1:opt.nColumn
            
           column = massConservation(currentData, interstVelocity, Feed, opt, sequence, string(k));
           [outletProfile, lastState] = secColumn(column.inlet, column.params, column.initialState, ParSwarm);
            
           currentData{eval(['sequence' '.' string(k)])}.outlet     = outletProfile;
           currentData{eval(['sequence' '.' string(k)])}.lastState  = lastState;
           
        end
      
        
%       Store the data, each one round (opt.nColumn switches), into plotData
%       plotData = column x position; 
        index = mod(i, opt.nColumn);
        if index == 0
            plotData(:,opt.nColumn) = currentData';
        else
            plotData(:,index) = currentData';
        end
        
        
%       convergence criterion was adopted in each nColumn iteration
%           ||( C(z, t) - C(z, t + 4 * t_s) ) / C(z, t)|| < tol for the column x
        if fix(i/opt.nColumn) == i/(opt.nColumn)
            
            diffNorm = norm( convergPrevious(:,1) - currentData{convergIndx}.outlet.concentration(:,1) ) + ...
                norm( convergPrevious(:,2) - currentData{convergIndx}.outlet.concentration(:,2) );
            
            stateNorm = norm( currentData{convergIndx}.outlet.concentration(:,1) ) + ...
                norm( currentData{convergIndx}.outlet.concentration(:,2));
            
            relativeDelta = diffNorm / stateNorm;
            if opt.enableDebug
                fprintf('---- Round: %3d    Switch: %4d    CSS_relError: %g \n', i/opt.nColumn, i, relativeDelta);
            end
                    
%           plot the outlet profile of each column in one round
            plotFigures(opt, plotData);
            
            if relativeDelta <= opt.tolIter
                break
            end
            
            convergPrevious = currentData{convergIndx}.outlet.concentration;
        end
    end
%-----------------------------------------------------------------------------------------

%   Compute the performance index, such Purity and Productivity
    Results = Purity_Productivity(plotData, opt);
    
%   Construct your own Objective Function and calculate the value
    objective = objectiveFunction(Results, opt);

	tTotal = toc(tTotal);
    if opt.enableDebug
        fprintf('Time elapsed for reaching the Cyclic Steady State: %g sec \n', tTotal);
    end
        
%   store the final data into DATA.mat file
%     save(sprintf('DATA_%2d.mat',fix(rand*100)),'Results');
%     fprintf('The results have been stored in the DATA.mat \n');
    
end

function Results = Purity_Productivity(plotData, opt)

%-----------------------------------------------------------------------------------------
% Calculation of the performance index of SMB, such Purity and Productivity

% 
%-----------------------------------------------------------------------------------------


    MolMass = opt.FructoseMolMass;
    Nominator = pi * (opt.columnDiameter/2)^2 * opt.columnLength * (1-opt.porosityColumn);

    if opt.nColumn == 4
%       using column 1 to calculate the integral of purity
%       Extract ports
        Purity_extract = trapz(plotData{1,1}.outlet.time, plotData{1,1}.outlet.concentration(:,2)) /...
            ( trapz(plotData{1,1}.outlet.time, plotData{1,1}.outlet.concentration(:,2)) +...
            trapz(plotData{1,1}.outlet.time, plotData{1,1}.outlet.concentration(:,1)) );

%       Raffinate ports  	
        Purity_raffinate = trapz(plotData{1,3}.outlet.time, plotData{1,3}.outlet.concentration(:,1)) / ...
            ( trapz(plotData{1,3}.outlet.time, plotData{1,3}.outlet.concentration(:,2)) +...
            trapz(plotData{1,3}.outlet.time, plotData{1,3}.outlet.concentration(:,1))	);	

%       per switching time, in the tank of extract port, such (unit: g/m^3) amount
%       of target component was collected.
        
        Productivity_extract = trapz(plotData{1,1}.outlet.time, plotData{1,1}.outlet.concentration(:,2))...
            * MolMass * opt.flowRate_recycle / Nominator;

        Productivity_raffinate = trapz(plotData{1,3}.outlet.time, plotData{1,3}.outlet.concentration(:,1))...
            * MolMass * opt.flowRate_recycle / Nominator;
    
    elseif opt.nColumn == 8
        
%       Extract ports
        Purity_extract = trapz(plotData{1,8}.outlet.time, plotData{1,8}.outlet.concentration(:,2)) /...
            ( trapz(plotData{1,8}.outlet.time, plotData{1,8}.outlet.concentration(:,2)) +...
            trapz(plotData{1,8}.outlet.time, plotData{1,8}.outlet.concentration(:,1)) );   
        
%       Raffinate ports  	
        Purity_raffinate = trapz(plotData{1,4}.outlet.time, plotData{1,4}.outlet.concentration(:,1)) / ...
            ( trapz(plotData{1,4}.outlet.time, plotData{1,4}.outlet.concentration(:,2)) +...
            trapz(plotData{1,4}.outlet.time, plotData{1,4}.outlet.concentration(:,1))	);	

%       per switch time, in the tank of extract port, such (unit: g/m^3) amount
%       of target component was collected.
        MolMass = opt.FructoseMolMass;
        Productivity_extract = trapz(plotData{1,8}.outlet.time, plotData{1,8}.outlet.concentration(:,2))...
            * MolMass * opt.flowRate_recycle / Nominator;

        Productivity_raffinate = trapz(plotData{1,4}.outlet.time, plotData{1,4}.outlet.concentration(:,1))...
            * MolMass * opt.flowRate_recycle / Nominator;
         
    end
    
    
    if opt.enableDebug
        fprintf('Purity (Extract): %g %% \n', Purity_extract * 100);
        fprintf('Purity (Raffinate): %g %% \n', Purity_raffinate * 100)
        fprintf('Productivity (Extract) in each switching time: %g g/m^3 \n', Productivity_extract);
        fprintf('Productivity (Raffinate) in each switching time: %g g/m^3 \n', Productivity_raffinate);
    end
    
    Results = struct('Purity_extract', Purity_extract, 'Purity_raffinate', Purity_raffinate, ...
        'Productivity_extract', Productivity_extract, 'Productivity_raffinate', Productivity_raffinate);
    
    Results.Data = plotData;

end

function objective = objectiveFunction(Results, opt)

%-----------------------------------------------------------------------------------------
% The objective function for the optimizers
% You can also define your own objective function here. The default function is: 

% Max Productivity_extract + Productivity_raffinate
% s.t. Purity_extract   >= 99% for more retained component
%      Purity_raffinate >= 99% for less retained component
%      other implicit constraints, such as upbound on Desorbent consumption
%-----------------------------------------------------------------------------------------


%   Construct the Penalty Function for the objective function
    penalty = abs( min(Results.Purity_extract - opt.Purity_extract_limit, 0) ) * opt.Penalty_factor ...
        + abs( min(Results.Purity_raffinate - opt.Purity_raffinate_limit, 0) ) * opt.Penalty_factor;

%   (-)Since in the optimizer, the defined programme is of optimization of minimum.    
    objective = -(Results.Productivity_extract + Results.Productivity_raffinate) + penalty;
    
    if opt.enableDebug
        fprintf('**** The objective value:  %g \n', objective);
    end
    
end

function plotFigures(opt, plotData)

%-----------------------------------------------------------------------------------------
%  This is the plot function 
%  The numbers in the figure() represent the number of the columns
%-----------------------------------------------------------------------------------------


    if nargin < 2
        disp('Error: there are no enough input data for the function, plotFigures');
    else
        if isempty(opt)
            disp('Error in plotFigures: the options of the parameters are missing');
        elseif isempty(plotData)
            disp('Error in plotFigures: the data for figures plotting are missing');
        end
    end
    
    if opt.enableDebug                        

        if opt.nColumn == 4

            for j = 1:opt.nColumn
                
                figure(j);clf

                y = [plotData{j,1}.outlet.concentration; plotData{j,2}.outlet.concentration;...
                    plotData{j,3}.outlet.concentration; plotData{j,4}.outlet.concentration];

                FigSet = plot(y); axis([0,opt.nColumn*opt.timePoints, 0,2e-3])
                ylabel('Concentration [Mol]', 'FontSize', 10);
                legend('comp 1', 'comp 2');

                set(FigSet, 'LineWidth', 2);
                set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
                set(gca, 'XTick', opt.timePoints/2:opt.timePoints:(opt.nColumn-0.5)*opt.timePoints);
                
                switch j
                    case 1
                        set(gca, 'XTickLabel', {'Zone I','Zone IV','Zone III','Zone II'});
                    case 2
                        set(gca, 'XTickLabel', {'Zone II','Zone I','Zone IV','Zone III'});
                    case 3
                        set(gca, 'XTickLabel', {'Zone III','Zone II','Zone I','Zone IV'});
                    case 4
                        set(gca, 'XTickLabel', {'Zone IV','Zone III','Zone II','Zone I'});
                end
                
                for i = 1: (opt.nColumn-1)
                    line([i*opt.timePoints,i*opt.timePoints], [0, 2e-3], 'color', 'k', 'LineStyle', ':');
                end
                
            end
            
        elseif opt.nColumn == 8
   
            for j = 1:opt.nColumn
                
                figure(j);clf

                y = [plotData{j,1}.outlet.concentration; plotData{j,2}.outlet.concentration;...
                    plotData{j,3}.outlet.concentration; plotData{j,4}.outlet.concentration;...
                    plotData{j,5}.outlet.concentration; plotData{j,6}.outlet.concentration;...
                    plotData{j,7}.outlet.concentration; plotData{j,8}.outlet.concentration];

                FigSet = plot(y); axis([0,opt.nColumn*opt.timePoints, 0,2e-3])
                ylabel('Concentration [Mol]', 'FontSize', 10);
                legend('comp 1', 'comp 2');

                set(FigSet, 'LineWidth', 2);
                set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
                
                switch j
                    case 1
                        set(gca, 'XTick', 0:2*opt.timePoints:opt.nColumn*opt.timePoints);
                        set(gca, 'XTickLabel', {'Zone I','Zone IV','Zone III','Zone II','Zone I'});
                    case 2
                        set(gca, 'XTick', opt.timePoints:2*opt.timePoints:(opt.nColumn-1)*opt.timePoints);
                        set(gca, 'XTickLabel', {'Zone I','Zone IV','Zone III','Zone II'});
                    case 3
                        set(gca, 'XTick', 0:2*opt.timePoints:opt.nColumn*opt.timePoints);
                        set(gca, 'XTickLabel', {'Zone II','Zone I','Zone IV','Zone III','Zone II'});
                    case 4
                        set(gca, 'XTick', opt.timePoints:2*opt.timePoints:(opt.nColumn-1)*opt.timePoints);
                        set(gca, 'XTickLabel', {'Zone II','Zone I','Zone IV','Zone III'});
                    case 5
                        set(gca, 'XTick', 0:2*opt.timePoints:opt.nColumn*opt.timePoints);
                        set(gca, 'XTickLabel', {'Zone III','Zone II','Zone I','Zone IV','Zone III'});
                    case 6
                        set(gca, 'XTick', opt.timePoints:2*opt.timePoints:(opt.nColumn-1)*opt.timePoints);
                        set(gca, 'XTickLabel', {'Zone III','Zone II','Zone I','Zone IV'});
                    case 7
                        set(gca, 'XTick', 0:2*opt.timePoints:opt.nColumn*opt.timePoints);
                        set(gca, 'XTickLabel', {'Zone IV','Zone III','Zone II','Zone I','Zone IV'});
                    case 8
                        set(gca, 'XTick', opt.timePoints:2*opt.timePoints:(opt.nColumn-1)*opt.timePoints);
                        set(gca, 'XTickLabel', {'Zone IV','Zone III','Zone II','Zone I'});
                end
                
                for i = 1: (opt.nColumn-1)
                    line([i*opt.timePoints,i*opt.timePoints], [0, 2e-3], 'color', 'k', 'LineStyle', ':');
                end
                
            end
            
        end  
        
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
