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
    
%   Initilize the starting points, currentData   
    currentData = cell(1, opt.nColumn);
    for k = 1:opt.nColumn
        currentData{k}.outlet.time = linspace(0, opt.switch, opt.timePoints);
        currentData{k}.outlet.concentration = zeros(length(Feed.time), 2);
        currentData{k}.outlet.time = linspace(0, opt.switch, opt.timePoints);
    end
    
    
    if opt.nColumn == 4
        
        sequence.a = 1; sequence.b = 2; sequence.c = 3; sequence.d = 4;
        
        string = char('c','b','a','d');
        convergIndx = 3;
        
    elseif opt.nColumn == 8
        
        sequence.a = 1; sequence.b = 2; sequence.c = 3; sequence.d = 4;
        sequence.e = 5; sequence.f = 6; sequence.g = 7; sequence.h = 8;
        
        string = char('e','d','c','b','a','h','g','f');
        convergIndx = 5;
        
    else
        warning('The simulation of %3g_column case is not finished so far', opt.nColumn);
    end
    
%   prellocation
    inlet       = cell(1,opt.nColumn);    
    for j = 1:opt.nColumn
        inlet{j}.time = linspace(0, opt.switch, opt.timePoints);
    end
   
    
%   convergPrevious is used for stopping criterion
    convergPrevious = currentData{convergIndx}.outlet.concentration;
    
    
%   Main loop
    for i = 1:opt.nMaxIter 

        if mod(i, opt.nColumn) ~= 0
            k = mod(i, opt.nColumn);
        else
            k = opt.nColumn;
        end
        
        if i == 1
            initialState = [];
        end
        
%       The simulation of four columns by the sequence, say, 'c', 'b', 'a', 'd'
        column = massConservation(currentData, interstVelocity, Feed, opt, sequence, string(k));
        [outletProfile, lastState] = secColumn(column.inlet, column.params, initialState);

        currentData{eval(['sequence' '.' string(k)])}.outlet     = outletProfile;
        initialState  = lastState;
        
        
        
%       plot the outlet profile of each column in once switch
        if mod(i, opt.nColumn) == 0
            plotFigures(opt, currentData);
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
    
%     store the final data into DATA.mat file
    save(sprintf('DATA_%2d.mat',fix(rand*100)),'currentData');
    fprintf('The results have been stored in the DATA.mat\n');
    
end


function plotFigures(opt, currentData)
%  This is the plot function 
%  The numbers in the figure() represent the number of the columns

    if opt.enableDebug                        

        if opt.nColumn == 4
            
            figure(01);clf

            subplot(1,opt.nColumn,1);
            FigSet1_1 = plot(currentData{4}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone IV'), 'FontSize', 10, 'color', 'g');
            xlabel('Time [s]', 'FontSize', 10);
            ylabel('Concentration [Mol]', 'FontSize', 10);
            legend('comp 1', 'comp 2');
            set(FigSet1_1, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,2);
            FigSet1_2 = plot(currentData{3}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone III'), 'FontSize', 10, 'color', 'b');
            set(FigSet1_2, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,3);
            FigSet1_3 = plot(currentData{2}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone II'), 'FontSize', 10, 'color', 'r');
            set(FigSet1_3, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,4);
            FigSet1_4 = plot(currentData{1}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone I'), 'FontSize', 10, 'color', 'm');	   
            set(FigSet1_4, 'LineWidth', 2);
            set(gca, 'FontSize', 10);
            
            
        elseif opt.nColumn == 8
            
            figure(01);clf

            subplot(1,opt.nColumn,1);
            FigSet3_8 = plot(currentData{8}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone IV'), 'FontSize', 10, 'color', 'g');
            xlabel('Time [s]', 'FontSize', 10);
            ylabel('Concentration [Mol]', 'FontSize', 10);
            legend('comp 1', 'comp 2');
            set(FigSet3_8, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,2);
            FigSet3_7 = plot(currentData{7}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone IV'), 'FontSize', 10, 'color', 'g');
            set(FigSet3_7, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,3);
            FigSet3_6 = plot(currentData{6}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone III'), 'FontSize', 10, 'color', 'b');
            set(FigSet3_6, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,4);
            FigSet3_5 = plot(currentData{5}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone III'), 'FontSize', 10, 'color', 'b');
            set(FigSet3_5, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,5);
            FigSet3_4 = plot(currentData{4}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone II'), 'FontSize', 10, 'color', 'r');
            xlabel('Time [s]', 'FontSize', 10);
            ylabel('Concentration [Mol]', 'FontSize', 10);
            legend('comp 1', 'comp 2');
            set(FigSet3_4, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,6);
            FigSet3_3 = plot(currentData{3}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone II'), 'FontSize', 10, 'color', 'r');
            set(FigSet3_3, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,7);
            FigSet3_2 = plot(currentData{2}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone I'), 'FontSize', 10, 'color', 'm');
            set(FigSet3_2, 'LineWidth', 2);
            set(gca, 'FontSize', 10);

            subplot(1,opt.nColumn,8);
            FigSet3_1 = plot(currentData{1}.outlet.concentration); axis([0,opt.timePoints, 0,2e-3]);
            title(sprintf('Zone I'), 'FontSize', 10, 'color', 'm');
            set(FigSet3_1, 'LineWidth', 2);
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