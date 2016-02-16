function updateAvailable = isSMBupdateAvailable

%==============================================================================
% This is the function to check if there is a new version SMB available.
%==============================================================================


    localPath = fileparts(mfilename('fullpath'));
    fileID = fopen([localPath filesep 'version.txt']);
    
    currentVersion = [];    
    try
        currentVersion = fgets(fileID);
    end
    
    if isempty(currentVersion)
        fprintf('The local version.txt file has been deleted. \n');
        return;
    end
    
    splitted = textscan(currentVersion, '%s', 'delimiter', '.');
    currentVersion = splitted{1};

    stableVersion = [];
    try
        stableVersion = urlread('https://raw.githubusercontent.com/modsim/CADET-SMB/master/version.txt'); 
    end
    
    if isempty(stableVersion)
        fprintf('The internet connection is not available now. \n');
        return;
    end
    
    stableVersionString = stableVersion;
    splitted = textscan(stableVersion, '%s', 'delimiter', '.');
    stableVersion = splitted{1};

    for i = 1:min(length(stableVersion), length(currentVersion))
        if str2double(stableVersion(i)) > str2double(currentVersion(i))
            updateAvailable = true;
            fprintf('There is new version available in the GitHub: https://github.com/modsim/CADET-SMB/releases. \n');
            break;
        end
    end
    
    fprintf('The newest version is now installed, Version %s', stableVersionString);

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