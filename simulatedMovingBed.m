function objective = simulatedMovingBed(varargin)

% =============================================================================
% This is the main function which is charge of switching to reach the
% cyclic steady state. The layout of the columns and the configuration is
% listed as follow:
%
%                                       FOUR-ZONE
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
%                                       FIVE-ZONE
%              5-column SMB                                       10-column SMB
%    Ext2                          Feed       |      Ext2                            Feed
%       \                          /          |         \                            /
%        --------Zone II(c)--------           |          --------Zone III(e/f)--------
%        |                        |           |          |                           | 
% Zone II(b)                      |           |     Zone II(d/c)                     |
%        |                        |           |          |                           |
% Ext1 --                    Zone IV(d)       |   Ext1 --                        Zone IV(g/h)
%        |                        |           |          |                           |
% Zone I(a)                       |           |     Zone I(b/a)                      |
%        |                        |           |          |                           | 
%        --------Zone V(e)---------           |          ---------Zone V(j/i)---------
%       /                          \          |         /                            \
% Desorbent                       Raffinate   |   Desorbent                         Raffinate
%
%             15-column SMB                                       20-column SMB
%    Ext2                            Feed       |      Ext2                              Feed
%       \                            /          |         \                              /
%        -------Zone II(g/h/i)-------           |          -------Zone III(i/g/k/l)-------
%        |                          |           |          |                             | 
% Zone II(f/e/d)                    |           | Zone II(h/g/f/e)                       |
%        |                          |           |          |                             |
% Ext1 --                    Zone IV(j/k/l)     |   Ext1 --                        Zone IV(m/n/o/p)
%        |                          |           |          |                             |
% Zone I(c/b/a)                     |           | Zone I(d/c/b/a)                        |
%        |                          |           |          |                             | 
%        -------Zone V(o/n/m)--------           |          -------Zone V(t/s/r/q)---------
%       /                            \          |         /                              \
% Desorbent                         Raffinate   |   Desorbent                           Raffinate
%
% Fluid phase goes from Zone I to Zone II to Zone III, while the ports switch direction
% is from Zone I to Zone IV to Zone III;
% =============================================================================


    global string;

    tTotal = tic;

    [opt, interstVelocity, Feed] = getParameters(varargin{:});

%   Initialize the starting points, currentData   
    currentData = cell(1, opt.nColumn);
    for k = 1:opt.nColumn
        currentData{k}.outlet.time = linspace(0, opt.switch, opt.timePoints);
        currentData{k}.outlet.concentration = zeros(length(Feed.time), opt.nComponents);
        currentData{k}.outlet.time = linspace(0, opt.switch, opt.timePoints);
    end

%   Numbered the columns for the sake of plotting 
%   Four-zone for binary separation, 1-1-1-1, 2-2-2-2, 3-3-3-3, and 4-4-4-4 configurations are available
    if opt.nZone == 4

        if opt.nColumn == 4

            sequence = cell2struct( [{1} {2} {3} {4}],{'a' 'b' 'c' 'd'},2 );
            string = char('c','b','a','d');
            convergIndx = 3;

        elseif opt.nColumn == 8

            sequence = cell2struct( [{1} {2} {3} {4} {5} {6} {7} {8}],{'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h'},2 );
            string = char('e','d','c','b','a','h','g','f');
            convergIndx = 5;

        elseif opt.nColumn == 12

            sequence = cell2struct( [{1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12}],...
                {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l'},2 );
            string = char('g','f','e','d','c','b','a','l','k','j','i','h');
            convergIndx = 7;

        elseif opt.nColumn == 16

            sequence = cell2struct( [{1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16}],...
                {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l' 'm' 'n' 'o' 'p'},2 );
            string = char('i','h','g','f','e','d','c','b','a','p','o','n','m','l','k','j');
            convergIndx = 9;

        else
            warning('The simulation of %3g_column case in %3g-zone is not finished so far', opt.nColumn, opt.nZone);

        end

%   Five-zone for binary separation, 1-1-1-1-1, 2-2-2-2-2, 3-3-3-3-3, and 4-4-4-4-4 configurations are available        
    elseif opt.nZone == 5

        if opt.nColumn == 5

            sequence = cell2struct( [{1} {2} {3} {4} {5}],{'a' 'b' 'c' 'd' 'e'},2 );
            string = char('d','c','b','a','e');
            convergIndx = 4;

        elseif opt.nColumn == 10

            sequence = cell2struct( [{1} {2} {3} {4} {5} {6} {7} {8} {9} {10}],...
                {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j'},2 );
            string = char('g','f','e','d','c','b','a','j','i','h');
            convergIndx = 7;

        elseif opt.nColumn == 15

            sequence = cell2struct( [{1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15}],...
                {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l' 'm' 'n' 'o' 'p'},2 );
            string = char('j','i','h','g','f','e','d','c','b','a','o','n','m','l','k');
            convergIndx = 10;

        elseif opt.nColumn == 20

            sequence = cell2struct( [{1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16} {17} {18} {19} {20}],...
                {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l' 'm' 'n' 'o' 'p' 'q' 'r' 's' 't'},2 );
            string = char('m','l','k','j','i','h','g','f','e','d','c','b','a','t','s','r','q','p','o','n');
            convergIndx = 13;

        else
            warning('The simulation of %3g_column case in %3g-zone is not finished so far', opt.nColumn, opt.nZone);

        end

    end

%   convergPrevious is used for stopping criterion
    convergPrevious = currentData{convergIndx}.outlet.concentration;


%----------------------------------------------------------------------------------------- 
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
        column = SMB.massConservation(currentData, interstVelocity, Feed, opt, sequence, string(k));
        [outletProfile, lastState] = SMB.secColumn(column.inlet, column.params, initialState, varargin{:});

        currentData{eval(['sequence' '.' string(k)])}.outlet = outletProfile;
        initialState  = lastState;


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
                    norm( currentData{convergIndx}.outlet.concentration(:,2)) + ...
                    norm( currentData{convergIndx}.outlet.concentration(:,3));

            elseif opt.nComponents == 4
                diffNorm = norm( convergPrevious(:,1) - currentData{convergIndx}.outlet.concentration(:,1) ) + ...
                    norm( convergPrevious(:,2) - currentData{convergIndx}.outlet.concentration(:,2) ) + ...
                    norm( convergPrevious(:,3) - currentData{convergIndx}.outlet.concentration(:,3) ) + ...
                    norm( convergPrevious(:,4) - currentData{convergIndx}.outlet.concentration(:,4) );

                stateNorm = norm( currentData{convergIndx}.outlet.concentration(:,1) ) + ...
                    norm( currentData{convergIndx}.outlet.concentration(:,2)) + ...
                    norm( currentData{convergIndx}.outlet.concentration(:,3)) + ...
                    norm( currentData{convergIndx}.outlet.concentration(:,4));
                
            end

            relativeDelta = diffNorm / stateNorm;

            if opt.enableDebug
                fprintf('---- Round: %3d    Switch: %4d    CSS_relError: %g \n', i/opt.nColumn, i, relativeDelta);
            end

%           plot the outlet profile of each column in one round
            SMB.plotFigures(opt, currentData);

            if relativeDelta <= opt.tolIter
                break
            end

            convergPrevious = currentData{convergIndx}.outlet.concentration;
        end
    end
%-----------------------------------------------------------------------------------------


%   Compute the performance index, such Purity and Productivity
    Results = SMB.Purity_Productivity(currentData, opt);

%   Construct your own Objective Function and calculate the value    
    objective = SMB.objectiveFunction(Results, opt);

    tTotal = toc(tTotal);
    if opt.enableDebug
        fprintf('The time elapsed for reaching the Cyclic Steady State: %g sec \n', tTotal);
    end

%   store the final data into DATA.mat file
    if opt.enableDebug
        save(sprintf('DATA_%2d.mat',fix(rand*100)),'Results');
        fprintf('The results have been stored in the DATA.mat \n');
    end

end
% =============================================================================
%  SMB - The Simulated Moving Bed Chromatography for separation of
%  target compounds, either binary or ternary.
%  
%  Author: QiaoLe He   E-mail: q.he@fz-juelich.de
%                                      
%  Institute: Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
%  
%  All rights reserved. Please see the license of CADET.
% =============================================================================