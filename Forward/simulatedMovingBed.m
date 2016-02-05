function simulatedMovingBed()

% =============================================================================
% This is the main function which is charge of switching to reach the
% cyclic steady state. The layout of the columns and the configuration is
% listed as follow:
%
%              4-column SMB                                       8-column SMB
% Extract                          Feed       |    Extract                           Feed
%       \                          /          |         \                            /
%        --------Zone II(b)--------           |          --------Zone II(c/d)--------
%        |                        |           |          |                          | 
% Zone I(a)                  Zone III(c)      |     Zone I(a/b)               Zone III(e/f)
%        |                        |           |          |                          | 
%        --------Zone IV(d)--------           |          --------Zone IV(h/g)--------
%       /                          \          |         /                            \
% Desorbent                       Raffinate   |   Desorbent                         Raffinate
%
%             12-column SMB                                       16-column SMB
% Extract                            Feed       |    Extract                         Feed
%       \                            /          |         \                          /
%        -- ----Zone II(d/e/f)-------           |          -----Zone II(e/f/g/h)-----
%        |                          |           |          |                        | 
% Zone I(c/b/a)                Zone III(g/h/i)  |  Zone I(a/b/c/d)           Zone III(i/j/k/l)
%        |                          |           |          |                        | 
%        -------Zone IV(l/k/j)-------           |          -----Zone IV(p/o/n/m)-----
%       /                            \          |         /                          \
% Desorbent                         Raffinate   |   Desorbent                       Raffinate
%
% Fluid phase goes from Zone I to Zone II to Zone III, while the ports switch direction
% is from Zone I to Zone IV to Zone III;
% =============================================================================
 

    global string;

    tTotal = tic;
    
    [opt, interstVelocity, Feed] = getParameters();
%     [opt, interstVelocity, Feed] = SMB.getParameters3(); % the default setup

%   Initialize the starting points, currentData
    currentData = cell(1, opt.nColumn);  
    for k = 1:opt.nColumn
        currentData{k}.outlet.time = linspace(0, opt.switch, opt.timePoints);
        currentData{k}.outlet.concentration = zeros(length(Feed.time), opt.nComponents); 
        currentData{k}.lastState = []; 
    end
    
%   Number the columns for the sake of plotting
    if opt.nColumn == 4
        
        sequence = cell2struct( [{4} {1} {2} {3}],{'a' 'b' 'c' 'd'},2 );
        string = char('a','b','c','d');
        convergIndx = 3;
        
    elseif opt.nColumn == 8
        
        sequence = cell2struct( [{8} {1} {2} {3} {4} {5} {6} {7}],{'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h'},2 );
        string = char('a','b','c','d','e','f','g','h');
        convergIndx = 5;
        
    elseif opt.nColumn == 12
        
        sequence = cell2struct( [{12} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11}],...
            {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l'},2 );
        string = char('a','b','c','d','e','f','g','h','i','j','k','l');
        convergIndx = 7;
        
    elseif opt.nColumn == 16
        
        sequence = cell2struct( [{16} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15}],...
            {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l' 'm' 'n' 'o' 'p'},2 );
        string = char('a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p');
        convergIndx = 9;
        
    else
        warning('The simulation of %3g_column case is not finished so far', opt.nColumn);
        
    end
        
%   preallocation
    plotData = cell(opt.nColumn,opt.nColumn);     
%   convergPrevious is used for stopping criterion
    convergPrevious = currentData{convergIndx}.outlet.concentration;
   

%-----------------------------------------------------------------------------------------
%   Main loop 
    for i = 1:opt.nMaxIter 

        if opt.nColumn == 4
            sequence = cell2struct( circshift( struct2cell(sequence),-1 ), ...
                {'a' 'b' 'c' 'd'} );
        elseif opt.nColumn == 8
            sequence = cell2struct( circshift( struct2cell(sequence),-1 ), ...
                {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h'} );
        elseif opt.nColumn == 12
            sequence = cell2struct( circshift( struct2cell(sequence),-1 ), ...
                {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l'} );
        elseif opt.nColumn == 16
            sequence = cell2struct( circshift( struct2cell(sequence),-1 ), ...
                {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l' 'm' 'n' 'o' 'p'} );
        end
        
        
%       The simulation of four columns by the sequence, say, 'a', 'b', 'c', 'd'
        for k = 1:opt.nColumn
            
           column = SMB.massConservation(currentData, interstVelocity, Feed, opt, sequence, string(k));
           [outletProfile, lastState] = SMB.secColumn(column.inlet, column.params, column.initialState);
            
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
            
            if opt.nComponents == 2
                diffNorm = norm( convergPrevious(:,1) - currentData{convergIndx}.outlet.concentration(:,1) ) + ...
                    norm( convergPrevious(:,2) - currentData{convergIndx}.outlet.concentration(:,2) );

                stateNorm = norm( currentData{convergIndx}.outlet.concentration(:,1) ) + ...
                    norm( currentData{convergIndx}.outlet.concentration(:,2));

            elseif opt.nComponents == 3
                diffNorm = norm( convergPrevious(:,1) - currentData{convergIndx}.outlet.concentration(:,1) ) + ...
                    norm( convergPrevious(:,2) - currentData{convergIndx}.outlet.concentration(:,2) ) + ...
                    norm( convergPrevious(:,3) - currentData{convergIndx}.outlet.concentration(:,3) );
            
                stateNorm = norm( currentData{convergIndx}.outlet.concentration(:,1) ) + ...
                    norm( currentData{convergIndx}.outlet.concentration(:,2) ) + ...
                    norm( currentData{convergIndx}.outlet.concentration(:,3) );
            end
            
            relativeDelta = diffNorm / stateNorm;

            fprintf('---- Round: %3d    Switch: %4d    CSS_relError: %g \n', i/opt.nColumn, i, relativeDelta);
                    
%           plot the outlet profile of each column in one round
            SMB.plotFigures(opt, plotData);
            
            if relativeDelta <= opt.tolIter
                break
            end
            
            convergPrevious = currentData{convergIndx}.outlet.concentration;
        end
    end  
%-----------------------------------------------------------------------------------------


%   Compute the performance index, such Purity and Productivity
    Results = SMB.Purity_Productivity(plotData, opt);
    
%   Construct your own Objective Function and calculate the value
    objective = SMB.objectiveFunction(Results, opt);

	tTotal = toc(tTotal);
    if opt.enableDebug
        fprintf('Time elapsed for reaching the Cyclic Steady State: %g sec \n', tTotal);
    end
        
%   store the final data into DATA.mat file
    save(sprintf('DATA_%2d.mat',fix(rand*100)),'Results');
    fprintf('The results have been stored in the DATA.mat\n');
    
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
