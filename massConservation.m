function column = massConservation(currentData, interstVelocity, Feed, opt, sequence, index)

% =============================================================================
% This is the function to calculate the concentration changes on each node.

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

% fluid goes from Zone I to Zone II to Zone III, while the switch direction
% is from Zone I to Zone IV to Zone III;

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
                column.initialState = currentData{sequence.a}.lastState;

                %   C_a^in = Q_d * C_d^out / Q_a
                concentration = zeros(length(Feed.time), 2);

                column.inlet.concentration = concentration .* params{sequence.d}.interstitialVelocity...
                    ./ params{sequence.a}.interstitialVelocity; 

                
%           The calculation of the column in the Zone II
%           node EXTRACT (index b)
            case 'b'

                column.params = params{sequence.b};
                column.initialState = currentData{sequence.b}.lastState;

                %   C_b^in = C_a^out
                column.inlet.concentration = currentData{sequence.a}.outlet.concentration;

                
%           The calculation of the column in the Zone III
%           node FEED (index c)
            case 'c' 
                column.params = params{sequence.c};
                column.initialState = currentData{sequence.c}.lastState;

                %   C_c^in = (Q_b*C_b^out + Q_F * C_F) / Q_c
                column.inlet.concentration = (currentData{sequence.b}.outlet.concentration .* ...
                    params{sequence.b}.interstitialVelocity + Feed.concentration .* interstVelocity.feed) ...
                    ./ params{sequence.c}.interstitialVelocity; 

                
%           The calculation of the column in the Zone IV
%           node RAFFINATE (index d)
            case 'd' 

                column.params = params{sequence.d};
                column.initialState = currentData{sequence.d}.lastState;
                
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
                column.initialState = currentData{sequence.a}.lastState;
                
                %   C_a^in = Q_d * C_d^out / Q_a
				concentration = zeros(length(Feed.time), 2);

                column.inlet.concentration = concentration .* params{sequence.h}.interstitialVelocity...
                    ./ params{sequence.a}.interstitialVelocity;


%           node DESORBENT (index b)
            case 'b'  
                
                column.params = params{sequence.b};
                column.initialState = currentData{sequence.b}.lastState;
                
                %   C_b^in = C_a^out
                column.inlet.concentration = currentData{sequence.a}.outlet.concentration;
    
                
%           The calculation of the column in the Zone II  
%           node EXTRACT (index c)
            case 'c'  
                
                column.params = params{sequence.c};
                column.initialState = currentData{sequence.c}.lastState;
               
                %   C_c^in = C_b^out
                column.inlet.concentration = currentData{sequence.b}.outlet.concentration;
                
                
%           node EXTRACT (index d)
            case 'd'  
                
                column.params = params{sequence.d};
                column.initialState = currentData{sequence.d}.lastState;
                
                %   C_d^in = C_c^out
                column.inlet.concentration = currentData{sequence.c}.outlet.concentration;
                
                
%           The calculation of the column in the Zone III
%           node FEED (index e)
            case 'e' 
                
                column.params = params{sequence.e};
                column.initialState = currentData{sequence.e}.lastState;
                
                %   C_e^in = (Q_d*C_d^out + Q_F * C_F) / Q_e
                column.inlet.concentration = (currentData{sequence.d}.outlet.concentration .* ...
                params{sequence.d}.interstitialVelocity + Feed.concentration .* interstVelocity.feed) ...
                ./ params{sequence.e}.interstitialVelocity;
                
            
%           node FEED (index f)
            case 'f' 
                
                column.params = params{sequence.f};
                column.initialState = currentData{sequence.f}.lastState;
                
                %   C_f^in = C_e
                column.inlet.concentration = currentData{sequence.e}.outlet.concentration;
       
                
%           The calculation of the column in the Zone IV 
%           node RAFFINATE (index g)
            case 'g'  
                
                column.params = params{sequence.g};
                column.initialState = currentData{sequence.g}.lastState;
                
                %   C_g^in = C_f^out
                column.inlet.concentration = currentData{sequence.f}.outlet.concentration;
                
                
%           node RAFFINATE (index h)
            case 'h' 
                
                column.params = params{sequence.h};
                column.initialState = currentData{sequence.h}.lastState;
                
                %   C_h^in = C_g^out
                column.inlet.concentration = currentData{sequence.g}.outlet.concentration;
                
        end
        
    end

end

function params = getParams(sequence, interstVelocity, opt)

    if opt.nColumn == 4
%       set the initial conditions to the solver, but when lastState is used, this setup will be ignored
        params{sequence.a}.initMobilCon = [0, 0];
        params{sequence.a}.initSolidCon = [0, 0];
        params{sequence.d} = params{sequence.a};
        params{sequence.c} = params{sequence.a};
        params{sequence.b} = params{sequence.a};

%       Interstitial velocity of each ZONE
        params{sequence.a}.interstitialVelocity = interstVelocity.recycle;
        params{sequence.b}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract;
        params{sequence.c}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract + interstVelocity.feed;
        params{sequence.d}.interstitialVelocity = interstVelocity.recycle - interstVelocity.desorbent;

    elseif opt.nColumn == 8
        
        params{sequence.a}.initMobilCon = [0, 0];
        params{sequence.a}.initSolidCon = [0, 0];
        params{sequence.b} = params{sequence.a};
        params{sequence.c} = params{sequence.a};
        params{sequence.d} = params{sequence.a};
        params{sequence.e} = params{sequence.a};
        params{sequence.f} = params{sequence.a};
        params{sequence.g} = params{sequence.a};
        params{sequence.h} = params{sequence.a};

        params{sequence.a}.interstitialVelocity = interstVelocity.recycle;
        params{sequence.b}.interstitialVelocity = interstVelocity.recycle;
        params{sequence.c}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract;
        params{sequence.d}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract;
        params{sequence.e}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract + interstVelocity.feed;
        params{sequence.f}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract + interstVelocity.feed;
        params{sequence.g}.interstitialVelocity = interstVelocity.recycle - interstVelocity.desorbent;
        params{sequence.h}.interstitialVelocity = interstVelocity.recycle - interstVelocity.desorbent;
              
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
