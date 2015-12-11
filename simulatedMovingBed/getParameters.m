function [opt, interstVelocity, Feed] = getParameters()

% =============================================================================
%     This is the function to input all the necessary data for simulation

% Returns: 
%       1. opt stands for options, which involves the parameter settings
%       for the algorithm, the binding isotherm, and the model equations

%       2. interstVelocity is calculated from flowrate of each column and inlet. 
%       interstitial_velocity = flow_rate / (across_area * porosity_Column)

%       3. Feed initializes the injection concentration
% =============================================================================


%   The parameter setting for simulator
    opt.tolIter         = 1e-4;
    opt.nMaxIter        = 1000;
    opt.nThreads        = 4;
    opt.nCellsColumn    = 40;
    opt.nCellsParticle  = 1;
    opt.switch          = 180;
    opt.timePoints      = 1000;
    opt.ABSTOL          = 1e-10;
    opt.INIT_STEP_SIZE  = 1e-14;
    opt.MAX_STEPS       = 5e6;
    
    opt.enableDebug = true;
    opt.nColumn = 4; % opt.nColumn = 4;

%   Binding: Linear Binding isotherm
    opt.nComponents = 2;
    opt.KA = [5.72 7.7];
    opt.KD = [1, 1];
    
%   Transport
    opt.dispersionColumn          = 3.8148e-20;     %
    opt.filmDiffusion             = [100 100];      % unknown 
    opt.diffusionParticle         = [1.6e4 1.6e4];  % unknown
    opt.diffusionParticleSurface  = [0.0 0.0];

%   Geometry
    opt.columnLength        = 0.25;      % m
    opt.particleRadius      = 0.0005;    % m % user-defined one in this case
    opt.porosityColumn      = 0.83;
    opt.porosityParticle    = 0.000001;  % unknown

%   Parameter units transformation
%   The flow rate of Zone I was defined as the recycle flow rate
    crossArea = pi * (0.02/2)^2;
    flowRate.recycle    = 9.62e-7;      % m^3/s  
    flowRate.feed       = 0.98e-7;      % m^3/s
    flowRate.raffinate  = 1.40e-7;      % m^3/s
    flowRate.desorbent  = 1.96e-7;      % m^3/s
    flowRate.extract    = 1.54e-7;      % m^3/s
    
%   Interstitial velocity = flow_rate / (across_area * opt.porosityColumn)
    interstVelocity.recycle   = flowRate.recycle / (crossArea*opt.porosityColumn);      % m/s 
    interstVelocity.feed      = flowRate.feed / (crossArea*opt.porosityColumn);         % m/s
    interstVelocity.raffinate = flowRate.raffinate / (crossArea*opt.porosityColumn);    % m/s
    interstVelocity.desorbent = flowRate.desorbent / (crossArea*opt.porosityColumn);    % m/s
    interstVelocity.extract   = flowRate.extract / (crossArea*opt.porosityColumn);      % m/s
   
    concentrationFeed = 0.55;   % g/m^3
    FructoseMolMass = 262.1535; % g/mol
    GlucoseMolMass  = 262.1535; % g/mol

    
%   Feed concentration setup   
    Feed.time = linspace(0, opt.switch, opt.timePoints);
    Feed.concentration = zeros(length(Feed.time), 2);

    Feed.concentration(1:end,1) = (concentrationFeed / FructoseMolMass );
    Feed.concentration(1:end,2) = (concentrationFeed / GlucoseMolMass);

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
