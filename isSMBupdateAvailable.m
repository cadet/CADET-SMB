function updateAvailable = isSMBupdateAvailable

%==============================================================================
% This is the function to check if there is a new version SMB available.
%==============================================================================


    localPath = fileparts(mfilename('fullpath'));
    fprintf('Adding %s to MATLAB PATH\n', localPath);
    path(localPath, path);
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
%      Copyright © 2008-2016: Eric von Lieres, Qiaole He
% 
%      Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
% 
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================