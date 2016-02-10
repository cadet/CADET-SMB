function column = massConservation(currentData, interstVelocity, Feed, opt, sequence, index)

% =============================================================================
% This is the function to calculate the concentration changes on each node.
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
% Fluid goes from Zone I to Zone II to Zone III, while the switch direction
% is from Zone I to Zone IV to Zone III;
%
% Parameters:
%       - currentData. Which includes each column's outlet concentration
%       (time-dependent), and the last state (which records every component's concentration 
%        in bulk phase and stationary phase, and used as the initial state for the next simulation).
%       - interstVelocity. The interstitial velocity of each column
%       - Feed. The initialied injection 
%       - opt. Options
%       - sequence. During switching, the structure used for storing the
%       sequence of columns
%       - index. It is a character. It tell this subroutine to calculate the specified column 
% 
% Returns: column
%   Preparation for next column simulation
%       - column.inlet. The new inlet concentration of each column, which is
%       obtained from mass conservation on each node.
%       - column.lastState. 
%       - column.params. Set the initial Mobile and Solid concentration to the
%       Simulator (if there is no lastState given), and also store the
%       interstitial velocity.
% =============================================================================


%   Time points
    column.inlet.time = linspace(0, opt.switch, opt.timePoints);

%   Get the interstitial velocity of each columns and initial condition
    params = getParams(sequence, interstVelocity, opt);

    
    if opt.nColumn == 4

        switch index
            
%           The calculation of the column in the Zone I
%           node DESORBENT (index a)
            case 'a' 

                column.params = params{sequence.a};

                %   C_a^in = Q_d * C_d^out / Q_a
                concentration = zeros(length(Feed.time), opt.nComponents);

                column.inlet.concentration = concentration .* params{sequence.d}.interstitialVelocity...
                    ./ params{sequence.a}.interstitialVelocity; 

                
%           The calculation of the column in the Zone II
%           node EXTRACT (index b)
            case 'b'

                column.params = params{sequence.b};

                %   C_b^in = C_a^out
                column.inlet.concentration = currentData{sequence.a}.outlet.concentration;

                
%           The calculation of the column in the Zone III
%           node FEED (index c)
            case 'c' 
                column.params = params{sequence.c};

                %   C_c^in = (Q_b * C_b^out + Q_F * C_F) / Q_c
                column.inlet.concentration = (currentData{sequence.b}.outlet.concentration .* ...
                    params{sequence.b}.interstitialVelocity + Feed.concentration .* interstVelocity.feed) ...
                    ./ params{sequence.c}.interstitialVelocity; 

                
%           The calculation of the column in the Zone IV
%           node RAFFINATE (index d)
            case 'd' 

                column.params = params{sequence.d};
                
                %   C_d^in = C_c^out
                column.inlet.concentration = currentData{sequence.c}.outlet.concentration;
        end
          
        
%     ------------------------------------------------------------------------------------    
    elseif opt.nColumn == 8
                
        switch index
        
%           The calculation of the column in the Zone I
%           node DESORBENT (index a)
            case 'a'  
                
                column.params = params{sequence.a};
                
                %   C_a^in = Q_h * C_h^out / Q_a
                concentration = zeros(length(Feed.time), opt.nComponents);
                
                column.inlet.concentration = concentration .* params{sequence.h}.interstitialVelocity...
                    ./ params{sequence.a}.interstitialVelocity;


%           node DESORBENT (index b)
            case 'b'  
                
                column.params = params{sequence.b};
                
                %   C_b^in = C_a^out
                column.inlet.concentration = currentData{sequence.a}.outlet.concentration;
    
                
%           The calculation of the column in the Zone II  
%           node EXTRACT (index c)
            case 'c'  
                
                column.params = params{sequence.c};
               
                %   C_c^in = C_b^out
                column.inlet.concentration = currentData{sequence.b}.outlet.concentration;
                
                
%           node EXTRACT (index d)
            case 'd'  
                
                column.params = params{sequence.d};
                
                %   C_d^in = C_c^out
                column.inlet.concentration = currentData{sequence.c}.outlet.concentration;
                
                
%           The calculation of the column in the Zone III
%           node FEED (index e)
            case 'e' 
                
                column.params = params{sequence.e};
                
                %   C_e^in = (Q_d * C_d^out + Q_F * C_F) / Q_e
                column.inlet.concentration = (currentData{sequence.d}.outlet.concentration .* ...
                	params{sequence.d}.interstitialVelocity + Feed.concentration .* interstVelocity.feed) ...
                	./ params{sequence.e}.interstitialVelocity;
                
            
%           node FEED (index f)
            case 'f' 
                
                column.params = params{sequence.f};
                
                %   C_f^in = C_e^out
                column.inlet.concentration = currentData{sequence.e}.outlet.concentration;
       
                
%           The calculation of the column in the Zone IV 
%           node RAFFINATE (index g)
            case 'g'  
                
                column.params = params{sequence.g};
                
                %   C_g^in = C_f^out
                column.inlet.concentration = currentData{sequence.f}.outlet.concentration;
                
                
%           node RAFFINATE (index h)
            case 'h' 
                
                column.params = params{sequence.h};
                
                %   C_h^in = C_g^out
                column.inlet.concentration = currentData{sequence.g}.outlet.concentration;
                
        end
        
%     ------------------------------------------------------------------------------------    
    elseif opt.nColumn == 12

        switch index

%           The calculation of the column in the Zone I
%           node DESORBENT (index a)
            case 'a'  

                column.params = params{sequence.a};

                %   C_a^in = Q_l * C_l^out / Q_a
                concentration = zeros(length(Feed.time), opt.nComponents);

                column.inlet.concentration = concentration .* params{sequence.l}.interstitialVelocity...
                    ./ params{sequence.a}.interstitialVelocity;


%           node DESORBENT (index b)
            case 'b'  

                column.params = params{sequence.b};

                %   C_b^in = C_a^out
                column.inlet.concentration = currentData{sequence.a}.outlet.concentration;


%           node DESORBENT (index c)
            case 'c'  

                column.params = params{sequence.c};

                %   C_c^in = C_b^out
                column.inlet.concentration = currentData{sequence.b}.outlet.concentration;


%           The calculation of the column in the Zone II  
%           node EXTRACT (index d)
            case 'd'  

                column.params = params{sequence.d};

                %   C_d^in = C_c^out
                column.inlet.concentration = currentData{sequence.c}.outlet.concentration;


%                   node EXTRACT (index e)
            case 'e'  

                column.params = params{sequence.e};

                %   C_e^in = C_d^out
                column.inlet.concentration = currentData{sequence.d}.outlet.concentration;


%           node EXTRACT (index f)
            case 'f'  

                column.params = params{sequence.f};

                %   C_f^in = C_e^out
                column.inlet.concentration = currentData{sequence.e}.outlet.concentration;


%           The calculation of the column in the Zone III
%           node FEED (index g)
            case 'g' 

                column.params = params{sequence.g};

                %   C_g^in = (Q_f * C_f^out + Q_F * C_F) / Q_g
                column.inlet.concentration = (currentData{sequence.f}.outlet.concentration .* ...
                	params{sequence.f}.interstitialVelocity + Feed.concentration .* interstVelocity.feed) ...
                	./ params{sequence.g}.interstitialVelocity;


%           node FEED (index h)
            case 'h' 

                column.params = params{sequence.h};

                %   C_h^in = C_g^out
                column.inlet.concentration = currentData{sequence.g}.outlet.concentration;


%           node FEED (index i)
            case 'i' 

                column.params = params{sequence.i};

                %   C_i^in = C_h^out
                column.inlet.concentration = currentData{sequence.h}.outlet.concentration;


%           The calculation of the column in the Zone IV 
%           node RAFFINATE (index j)
            case 'j'  

                column.params = params{sequence.j};

                %   C_j^in = C_i^out
                column.inlet.concentration = currentData{sequence.i}.outlet.concentration;


%           node RAFFINATE (index k)
            case 'k' 

                column.params = params{sequence.k};

                %   C_k^in = C_j^out
                column.inlet.concentration = currentData{sequence.j}.outlet.concentration;


%           node RAFFINATE (index l)
            case 'l' 

                column.params = params{sequence.l};

                %   C_l^in = C_k^out
                column.inlet.concentration = currentData{sequence.k}.outlet.concentration;

        end

%     ------------------------------------------------------------------------------------    
    elseif opt.nColumn == 16

        switch index

%           The calculation of the column in the Zone I
%           node DESORBENT (index a)
            case 'a'  

                column.params = params{sequence.a};

                %   C_a^in = Q_p * C_p^out / Q_a
                concentration = zeros(length(Feed.time), opt.nComponents);

                column.inlet.concentration = concentration .* params{sequence.p}.interstitialVelocity...
                    ./ params{sequence.a}.interstitialVelocity;


%           node DESORBENT (index b)
            case 'b'  

                column.params = params{sequence.b};

                %   C_b^in = C_a^out
                column.inlet.concentration = currentData{sequence.a}.outlet.concentration;


%           node DESORBENT (index c)
            case 'c'  

                column.params = params{sequence.c};

                %   C_c^in = C_b^out
                column.inlet.concentration = currentData{sequence.b}.outlet.concentration;


%           node DESORBENT (index d)
            case 'd'  

                column.params = params{sequence.d};

                %   C_d^in = C_c^out
                column.inlet.concentration = currentData{sequence.c}.outlet.concentration;


%           The calculation of the column in the Zone II  
%           node EXTRACT (index e)
            case 'e'  

                column.params = params{sequence.e};

                %   C_e^in = C_d^out
                column.inlet.concentration = currentData{sequence.d}.outlet.concentration;


%           node EXTRACT (index f)
            case 'f'  

                column.params = params{sequence.f};

                %   C_f^in = C_e^out
                column.inlet.concentration = currentData{sequence.e}.outlet.concentration;


%           node EXTRACT (index g)
            case 'g'  

                column.params = params{sequence.g};

                %   C_g^in = C_f^out
                column.inlet.concentration = currentData{sequence.f}.outlet.concentration;


%           node EXTRACT (index h)
            case 'h'  

                column.params = params{sequence.h};

                %   C_h^in = C_g^out
                column.inlet.concentration = currentData{sequence.g}.outlet.concentration;


%           The calculation of the column in the Zone III
%           node FEED (index i)
            case 'i' 

                column.params = params{sequence.i};

                %   C_i^in = (Q_h * C_h^out + Q_F * C_F) / Q_i
                column.inlet.concentration = (currentData{sequence.h}.outlet.concentration .* ...
                	params{sequence.h}.interstitialVelocity + Feed.concentration .* interstVelocity.feed) ...
                	./ params{sequence.i}.interstitialVelocity;


%           node FEED (index j)
            case 'j' 

                column.params = params{sequence.j};

                %   C_j^in = C_i^out
                column.inlet.concentration = currentData{sequence.i}.outlet.concentration;


%           node FEED (index k)
            case 'k' 

                column.params = params{sequence.k};

                %   C_k^in = C_j^out
                column.inlet.concentration = currentData{sequence.j}.outlet.concentration;


%           node FEED (index l)
            case 'l' 

                column.params = params{sequence.l};

                %   C_l^in = C_k^out
                column.inlet.concentration = currentData{sequence.k}.outlet.concentration;


%           The calculation of the column in the Zone IV 
%           node RAFFINATE (index m)
            case 'm'  

                column.params = params{sequence.m};

                %   C_m^in = C_l^out
                column.inlet.concentration = currentData{sequence.l}.outlet.concentration;


%           node RAFFINATE (index n)
            case 'n' 

                column.params = params{sequence.n};

                %   C_n^in = C_m^out
                column.inlet.concentration = currentData{sequence.m}.outlet.concentration;


%           node RAFFINATE (index o)
            case 'o' 

                column.params = params{sequence.o};

                %   C_o^in = C_n^out
                column.inlet.concentration = currentData{sequence.n}.outlet.concentration;


%           node RAFFINATE (index p)
            case 'p' 

                column.params = params{sequence.p};

                %   C_p^in = C_o^out
                column.inlet.concentration = currentData{sequence.o}.outlet.concentration;

        end
        
    end

end

function params = getParams(sequence, interstVelocity, opt)

%-----------------------------------------------------------------------------------------
% After each swtiching, the value of velocities and initial conditions are
% changed 
%-----------------------------------------------------------------------------------------
	

	global string;
    
    params = cell(1, opt.nColumn);
    for k = 1:opt.nColumn
        params{k} = struct('initMobilCon', [], 'initSolidCon', [], 'interstitialVelocity', []);
    end
    
    for j = 1:opt.nColumn
%       set the initial conditions to the solver, but when lastState is used, this setup will be ignored        
        params{eval(['sequence' '.' string(j)])}.initMobilCon = zeros(1,opt.nComponents);
	    params{eval(['sequence' '.' string(j)])}.initSolidCon = zeros(1,opt.nComponents);
    end
    
    if opt.nColumn == 4
        
        for i = 1: opt.nColumn
%           Interstitial velocity of each ZONE
            if strcmp('a', string(i))
				params{eval(['sequence' '.' string(i)])}.interstitialVelocity = interstVelocity.recycle;
            elseif strcmp('b', string(i))
        		params{eval(['sequence' '.' string(i)])}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract;
            elseif strcmp('c', string(i))
				params{eval(['sequence' '.' string(i)])}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract + interstVelocity.feed;
            elseif strcmp('d', string(i))
        		params{eval(['sequence' '.' string(i)])}.interstitialVelocity = interstVelocity.recycle - interstVelocity.desorbent;
            end
        end

    elseif opt.nColumn == 8
        
        for i = 1: opt.nColumn
%           Interstitial velocity of each ZONE
            if strcmp('a', string(i)) || strcmp('b', string(i))
				params{eval(['sequence' '.' string(i)])}.interstitialVelocity = interstVelocity.recycle;
            elseif strcmp('c', string(i)) || strcmp('d', string(i))
        		params{eval(['sequence' '.' string(i)])}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract;
            elseif strcmp('e', string(i)) || strcmp('f', string(i))
				params{eval(['sequence' '.' string(i)])}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract + interstVelocity.feed;
            elseif strcmp('g', string(i)) || strcmp('h', string(i))
        		params{eval(['sequence' '.' string(i)])}.interstitialVelocity = interstVelocity.recycle - interstVelocity.desorbent;
            end
        end
        
    elseif opt.nColumn == 12

        for i = 1: opt.nColumn
%                   Interstitial velocity of each ZONE
            if strcmp('a', string(i)) || strcmp('b', string(i)) || strcmp('c', string(i))
                params{eval(['sequence' '.' string(i)])}.interstitialVelocity = interstVelocity.recycle;
            elseif strcmp('d', string(i)) || strcmp('e', string(i)) || strcmp('f', string(i))
                params{eval(['sequence' '.' string(i)])}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract;
            elseif strcmp('g', string(i)) || strcmp('h', string(i)) || strcmp('i', string(i))
                params{eval(['sequence' '.' string(i)])}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract + interstVelocity.feed;
            elseif strcmp('j', string(i)) || strcmp('k', string(i)) || strcmp('l', string(i))
                params{eval(['sequence' '.' string(i)])}.interstitialVelocity = interstVelocity.recycle - interstVelocity.desorbent;
            end
        end

    elseif opt.nColumn == 16

        for i = 1: opt.nColumn
%                   Interstitial velocity of each ZONE
            if strcmp('a', string(i)) || strcmp('b', string(i)) || strcmp('c', string(i)) || strcmp('d', string(i))
                params{eval(['sequence' '.' string(i)])}.interstitialVelocity = interstVelocity.recycle;
            elseif strcmp('e', string(i)) || strcmp('f', string(i)) || strcmp('g', string(i)) || strcmp('h', string(i))
                params{eval(['sequence' '.' string(i)])}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract;
            elseif strcmp('i', string(i)) || strcmp('j', string(i)) || strcmp('k', string(i)) || strcmp('l', string(i))
                params{eval(['sequence' '.' string(i)])}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract + interstVelocity.feed;
            elseif strcmp('m', string(i)) || strcmp('n', string(i)) || strcmp('o', string(i)) || strcmp('p', string(i))
                params{eval(['sequence' '.' string(i)])}.interstitialVelocity = interstVelocity.recycle - interstVelocity.desorbent;
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