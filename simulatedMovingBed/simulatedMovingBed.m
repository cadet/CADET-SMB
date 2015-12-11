function simulatedMovingBed()

% =============================================================================
%  This is the main function which is charge of switching to reach the
%  cyclic steady state. The layout of the columns and the configuration is
%  listed as follow:

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

% fluid phase goes from Zone I to Zone II to Zone III, while the ports switch direction
% is from Zone I to Zone IV to Zone III;
% =============================================================================


    tTotal = tic;
    
    [opt, interstVelocity, Feed] = getParameters();

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
    inlet    = cell(1,opt.nColumn);
    for j = 1:opt.nColumn
        inlet{j}.time = linspace(0, opt.switch, opt.timePoints);
    end
 
    
%   convergPrevious is used for stopping criterion
    convergPrevious = currentData{convergIndx}.outlet.concentration;
   
    
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
           [outletProfile, lastState] = secColumn(column.inlet, column.params, column.initialState);
            
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
        
        
%       plot the outlet profile of each column in one round
        if index == 0
            plotFigures(opt, plotData);
        end
        
        
%       convergence criterion was adopted in each nColumn iteration
%           ||( C(z, t) - C(z, t + 4 * t_s) ) / C(z, t)|| < tol for the column x
        if fix(i/opt.nColumn) == i/(opt.nColumn)
            
            diffNorm = norm( convergPrevious(:,1) - currentData{convergIndx}.outlet.concentration(:,1) ) + ...
                norm( convergPrevious(:,2) - currentData{convergIndx}.outlet.concentration(:,2) );
            
            stateNorm = norm( currentData{convergIndx}.outlet.concentration(:,1) ) + ...
                norm( currentData{convergIndx}.outlet.concentration(:,2));
            
            relativeDelta = diffNorm / stateNorm;
            fprintf('---- Round: %3d    Switch: %4d    CSS_relError: %g \n', i/opt.nColumn, i, relativeDelta);
            
            if relativeDelta <= opt.tolIter
                break
            end
            
            convergPrevious = currentData{convergIndx}.outlet.concentration;
        end
    end
	

    tTotal = toc(tTotal);
    if opt.enableDebug
        fprintf('The time elapsed for reaching the Cyclic Steady State: %g sec \n', tTotal);
    end
    
%   store the final data into DATA.mat file
    save(sprintf('DATA_%2d.mat',fix(rand*100)),'currentData');
    fprintf('The results have been stored in the DATA.mat\n');
    
end


function plotFigures(opt, plotData)

%  This is the plot function 
%  The numbers in the figure() represent the number of the columns

    if opt.enableDebug                        

        if opt.nColumn == 4
%           column_1
            figure(01);clf

            subplot(1,opt.nColumn,1);
            FigSet1_1 = plot(plotData{1,1}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone I'), 'FontSize', 10, 'color', 'm');
            xlabel('Time [s]', 'FontSize', 10);
            ylabel('Concentration [Mol]', 'FontSize', 10);
            legend('comp 1', 'comp 2');
            set(FigSet1_1, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,2);
            FigSet1_2 = plot(plotData{1,2}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone IV'), 'FontSize', 10, 'color', 'g');
            set(FigSet1_2, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,3);
            FigSet1_3 = plot(plotData{1,3}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone III'), 'FontSize', 10, 'color', 'b');
            set(FigSet1_3, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,4);
            FigSet1_4 = plot(plotData{1,4}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone II'), 'FontSize', 10, 'color', 'r');	   
            set(FigSet1_4, 'LineWidth', 2);
            set(gca, 'FontSize', 10);
%           ------------------------------------------------------------------------------
%           column_2
            figure(02);clf

            subplot(1,opt.nColumn,1);
            FigSet2_1 = plot(plotData{2,1}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone II'), 'FontSize', 10, 'color', 'r');
            xlabel('Time [s]', 'FontSize', 10);
            ylabel('Concentration [Mol]', 'FontSize', 10);
            legend('comp 1', 'comp 2');
            set(FigSet2_1, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,2);
            FigSet2_2 = plot(plotData{2,2}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone I'), 'FontSize', 10, 'color', 'm');
            set(FigSet2_2, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,3);
            FigSet2_3 = plot(plotData{2,3}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone IV'), 'FontSize', 10, 'color', 'g');
            set(FigSet2_3, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,4);
            FigSet2_4 = plot(plotData{2,4}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone III'), 'FontSize', 10, 'color', 'b');	   
            set(FigSet2_4, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

%           ------------------------------------------------------------------------------
%           column_3
            figure(03);clf

            subplot(1,opt.nColumn,1);
            FigSet3_1 = plot(plotData{3,1}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone III'), 'FontSize', 10, 'color', 'b');
            xlabel('Time [s]', 'FontSize', 10);
            ylabel('Concentration [Mol]', 'FontSize', 10);
            legend('comp 1', 'comp 2');
            set(FigSet3_1, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,2);
            FigSet3_2 = plot(plotData{3,2}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone II'), 'FontSize', 10, 'color', 'r');
            set(FigSet3_2, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,3);
            FigSet3_3 = plot(plotData{3,3}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone I'), 'FontSize', 10, 'color', 'm');
            set(FigSet3_3, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,4);
            FigSet3_4 = plot(plotData{3,4}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone IV'), 'FontSize', 10, 'color', 'g');	   
            set(FigSet3_4, 'LineWidth', 2);
            set(gca, 'FontSize', 10);
            
%           ------------------------------------------------------------------------------
%           column_4
            figure(04);clf

            subplot(1,opt.nColumn,1);
            FigSet4_1 = plot(plotData{4,1}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone IV'), 'FontSize', 10, 'color', 'g');
            xlabel('Time [s]', 'FontSize', 10);
            ylabel('Concentration [Mol]', 'FontSize', 10);
            legend('comp 1', 'comp 2');
            set(FigSet4_1, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,2);
            FigSet4_2 = plot(plotData{4,2}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone III'), 'FontSize', 10, 'color', 'b');
            set(FigSet4_2, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,3);
            FigSet4_3 = plot(plotData{4,3}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone II'), 'FontSize', 10, 'color', 'r');
            set(FigSet4_3, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,4);
            FigSet4_4 = plot(plotData{4,4}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone I'), 'FontSize', 10, 'color', 'm');	   
            set(FigSet4_4, 'LineWidth', 2);
            set(gca, 'FontSize', 10);
            
        elseif opt.nColumn == 8
            
%           column_1
            figure(01);clf

            subplot(1,opt.nColumn,1);
            FigSet1_1 = plot(plotData{1,1}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone I'), 'FontSize', 10, 'color', 'm');
            xlabel('Time [s]', 'FontSize', 10);
            ylabel('Concentration [Mol]', 'FontSize', 10);
            legend('comp 1', 'comp 2');
            set(FigSet1_1, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,2);
            FigSet1_2 = plot(plotData{1,2}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone IV'), 'FontSize', 10, 'color', 'g');
            set(FigSet1_2, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,3);
            FigSet1_3 = plot(plotData{1,3}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone IV'), 'FontSize', 10, 'color', 'g');
            set(FigSet1_3, 'LineWidth', 2);
            set(gca, 'FontSize', 10);
            
            subplot(1,opt.nColumn,4);
            FigSet1_4 = plot(plotData{1,4}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone III'), 'FontSize', 10, 'color', 'b');
            set(FigSet1_4, 'LineWidth', 2);
            set(gca, 'FontSize', 10);
            
            subplot(1,opt.nColumn,5);
            FigSet1_5 = plot(plotData{1,5}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone III'), 'FontSize', 10, 'color', 'b');
            set(FigSet1_5, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,6);
            FigSet1_6 = plot(plotData{1,6}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone II'), 'FontSize', 10, 'color', 'r');
            set(FigSet1_6, 'LineWidth', 2);
            set(gca, 'FontSize', 10);
            
            subplot(1,opt.nColumn,7);
            FigSet1_7 = plot(plotData{1,7}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone II'), 'FontSize', 10, 'color', 'r');
            set(FigSet1_7, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,8);
            FigSet1_8 = plot(plotData{1,8}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone I'), 'FontSize', 10, 'color', 'm');	   
            set(FigSet1_8, 'LineWidth', 2);
            set(gca, 'FontSize', 10);
%           ------------------------------------------------------------------------------
%           column_2
            figure(02);clf

            subplot(1,opt.nColumn,1);
            FigSet2_1 = plot(plotData{2,1}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone I'), 'FontSize', 10, 'color', 'm');
            xlabel('Time [s]', 'FontSize', 10);
            ylabel('Concentration [Mol]', 'FontSize', 10);
            legend('comp 1', 'comp 2');
            set(FigSet2_1, 'LineWidth', 2);
            set(gca, 'FontSize', 10);
            
            subplot(1,opt.nColumn,2);
            FigSet2_2 = plot(plotData{2,2}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone I'), 'FontSize', 10, 'color', 'm');
            set(FigSet2_2, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,3);
            FigSet2_3 = plot(plotData{2,3}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone IV'), 'FontSize', 10, 'color', 'g');
            set(FigSet2_3, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,4);
            FigSet2_4 = plot(plotData{2,4}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone IV'), 'FontSize', 10, 'color', 'g');
            set(FigSet2_4, 'LineWidth', 2);
            set(gca, 'FontSize', 10);
            
            subplot(1,opt.nColumn,5);
            FigSet2_5 = plot(plotData{2,5}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone III'), 'FontSize', 10, 'color', 'b');
            set(FigSet2_5, 'LineWidth', 2);
            set(gca, 'FontSize', 10);
            
            subplot(1,opt.nColumn,6);
            FigSet2_6 = plot(plotData{2,6}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone III'), 'FontSize', 10, 'color', 'b');
            set(FigSet2_6, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,7);
            FigSet2_7 = plot(plotData{2,7}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone II'), 'FontSize', 10, 'color', 'r');
            set(FigSet2_7, 'LineWidth', 2);
            set(gca, 'FontSize', 10);
            
            subplot(1,opt.nColumn,8);
            FigSet2_8 = plot(plotData{2,8}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone II'), 'FontSize', 10, 'color', 'r');
            set(FigSet2_8, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

%           ------------------------------------------------------------------------------
%           column_3
            figure(03);clf

            subplot(1,opt.nColumn,1);
            FigSet3_8 = plot(plotData{3,1}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone II'), 'FontSize', 10, 'color', 'r');
            xlabel('Time [s]', 'FontSize', 10);
            ylabel('Concentration [Mol]', 'FontSize', 10);
            legend('comp 1', 'comp 2');
            set(FigSet3_8, 'LineWidth', 2);
            set(gca, 'FontSize', 10);
            
            subplot(1,opt.nColumn,2);
            FigSet3_2 = plot(plotData{3,2}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone I'), 'FontSize', 10, 'color', 'm');
            set(FigSet3_2, 'LineWidth', 2);
            set(gca, 'FontSize', 10);
            
            subplot(1,opt.nColumn,3);
            FigSet3_3 = plot(plotData{3,3}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone I'), 'FontSize', 10, 'color', 'm');
            set(FigSet3_3, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,4);
            FigSet3_4 = plot(plotData{3,4}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone IV'), 'FontSize', 10, 'color', 'g');
            set(FigSet3_4, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,5);
            FigSet3_5 = plot(plotData{3,5}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone IV'), 'FontSize', 10, 'color', 'g');
            set(FigSet3_5, 'LineWidth', 2);
            set(gca, 'FontSize', 10);
            
            subplot(1,opt.nColumn,6);
            FigSet3_6 = plot(plotData{3,6}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone III'), 'FontSize', 10, 'color', 'b');
            set(FigSet3_6, 'LineWidth', 2);
            set(gca, 'FontSize', 10);
            
            subplot(1,opt.nColumn,7);
            FigSet3_7 = plot(plotData{3,7}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone III'), 'FontSize', 10, 'color', 'b');
            set(FigSet3_7, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,8);
            FigSet3_8 = plot(plotData{3,8}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone II'), 'FontSize', 10, 'color', 'r');
            set(FigSet3_8, 'LineWidth', 2);
            set(gca, 'FontSize', 10);
            
            
%           ------------------------------------------------------------------------------
%           column_4
            figure(04);clf

            subplot(1,opt.nColumn,1);
            FigSet4_1 = plot(plotData{4,1}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone II'), 'FontSize', 10, 'color', 'r');
            xlabel('Time [s]', 'FontSize', 10);
            ylabel('Concentration [Mol]', 'FontSize', 10);
            legend('comp 1', 'comp 2');
            set(FigSet4_1, 'LineWidth', 2);
            set(gca, 'FontSize', 10);
            
            subplot(1,opt.nColumn,2);
            FigSet4_2 = plot(plotData{4,2}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone II'), 'FontSize', 10, 'color', 'r');
            set(FigSet4_2, 'LineWidth', 2);
            set(gca, 'FontSize', 10);
            
            subplot(1,opt.nColumn,3);
            FigSet4_3 = plot(plotData{4,3}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone I'), 'FontSize', 10, 'color', 'm');
            set(FigSet4_3, 'LineWidth', 2);
            set(gca, 'FontSize', 10);
            
            subplot(1,opt.nColumn,4);
            FigSet4_4 = plot(plotData{4,4}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone I'), 'FontSize', 10, 'color', 'm');
            set(FigSet4_4, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,5);
            FigSet4_5 = plot(plotData{4,5}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone IV'), 'FontSize', 10, 'color', 'g');
            set(FigSet4_5, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,6);
            FigSet4_6 = plot(plotData{4,6}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone IV'), 'FontSize', 10, 'color', 'g');
            set(FigSet4_6, 'LineWidth', 2);
            set(gca, 'FontSize', 10);
            
            subplot(1,opt.nColumn,7);
            FigSet4_7 = plot(plotData{4,7}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone III'), 'FontSize', 10, 'color', 'b');
            set(FigSet4_7, 'LineWidth', 2);
            set(gca, 'FontSize', 10);
            
            subplot(1,opt.nColumn,8);
            FigSet4_8 = plot(plotData{4,8}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone III'), 'FontSize', 10, 'color', 'b');
            set(FigSet4_8, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

        end  
        
    end

end
% =============================================================================
%  SMB - The Simulated Moving Bed Chromatography for separation of
%  target compounds, such as fructose and glucose.
%  
%  Author: QiaoLe He
%                                      
%  Institute: Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
%  
%  All rights reserved. Please see the license of CADET.
% =============================================================================