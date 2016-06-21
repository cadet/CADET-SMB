
% =============================================================================
% This is a simple script to invoke the four-zone cascade simulation
% =============================================================================


% nominate the connection between sub-unit one and sub-unit two
intermediate_feed = 'raffinate';

% the sequential simulation of two sub-units
for i = 1:2
    
    objective = simulatedMovingBed(i, intermediate_feed);
%     pause;

end

% =============================================================================
%  SMB - The Simulated Moving Bed Chromatography for separation of
%  target compounds, either binary or ternary.
% 
%      Copyright © 2008-2016: Eric von Lieres, Qiaole He
% 
%      Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
% 
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
