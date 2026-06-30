function varargout = findSettingsFile
    % Find the MPQC settings file and return as a string
    %
    % function settingsFile = mpqc.settings.findSettingsFile()
    %
    % Purpose
    % There are three places the MPQC settings file could be located. This function
    % Looks through this in the following order:
    % 1. In the user's home folder in a folder called MPQC_Settings
    % 2. In C:\MPQC_Settings
    % 3. In the mpqc\SETTINGS folder
    %
    % If there are settings files in more than one location, it is the first one to be found
    % that is read. The rest are ignored.
    %
    % If no files are found, the user is prompted where they want an empty file to be created.
    %
    %
    % Inputs
    % none
    %
    % Outputs
    % settingsFile - [optional] Full path to settings file
    % backupDirectory - [optional] Full path to the directory that contains the backup settings
    % possibleSettingsFileLocations - Cell array containing a list of possible locations for
    %                                the settings file.
    % availableSettingsFiles - Cell array containing the list of all settings file in all valid paths
    %
    %
    % What if you wish to move your settings file?
    % 1. Find the location of the settings file right now using output argument 1.
    % 2. Find the list of paths defining where the settings could could be using output argument 3.
    % 3. Choose a directory path from the list. Make the directory if necessary.
    % 4. Move your settings file using Windows
    %
    %
    % Rob Campbell, SWC AMF, initial commit 2022


    settingsDirs = mpqc.settings.possibleSettingsLocations; % Potential locations for the settings file
    settingsFname = mpqc.settings.returnMPQC_SettingsFileName; % The name of the settings file


    % Look for a the settings file
    settingsFileLocations = {};
    backupDirLocations = {};

    for ii=1:length(settingsDirs)

        if ~settingsDirs(ii).settingsLocation_exists
            continue
        end

        tSettings = fullfile(settingsDirs(ii).settingsLocation,settingsFname);
        if exist(tSettings,'file')
            settingsFileLocations{end+1} = tSettings;
            backupDirLocations{end+1} = settingsDirs(ii).backupSettingsLocation;
        end

    end

    if ~isempty(settingsFileLocations)
        settingsFilePath = settingsFileLocations{1}; % Choose first on the list
        backupDir = backupDirLocations{1};
    end


    % No settings files found so we prompt the user where to set up
    if isempty(settingsFileLocations)
        clc
        fprintf('\n\n *** No MPQC settings file found. Where do you want to store these? *** \n\n')

        for ii=1:length(settingsDirs)

            switch settingsDirs(ii).locationType
            case 'homefolder'
                fprintf('%d. At "%s"\n', ii, settingsDirs(ii).settingsLocation)
                fprintf(['   Choose this option if you always use the same user account and ', ...
                        'are not updating MPQC via Git.\n\n'])
            case 'C'
                fprintf('%d. At "%s"\n', ii, settingsDirs(ii).settingsLocation)
                fprintf(['    RECOMMENDED: Choose this if you change user accounts and are ', ...
                        'not updating MPQC via Git.\n\n'])
            case 'mpqc'
                fprintf('%d. At "%s"\n', ii, settingsDirs(ii).settingsLocation)
                fprintf(['   ADVANCED: Choose this option if you like your settings files along with the code ', ...
                        'and are updating MPQC via Git.\n\n'])
            end
        end

        % Handle user input
        while true

            reply = input(sprintf('[1 ... %d]?  ', length(settingsDirs)),'s');
            reply = str2num(reply);
            if isempty(reply)
                continue
            end

            if reply > length(settingsDirs) || reply < 1
                continue
            end

            % Make folders
            S = settingsDirs(reply);
            fprintf('\nCreating settings directories and files -- \n')

            if ~S.settingsLocation_exists
                fprintf('Making directory %s\n', S.settingsLocation);
                mkdir(S.settingsLocation)
            end
            if ~S.backupSettingsLocation_exists
                fprintf('Making directory %s\n', S.backupSettingsLocation);
                mkdir(S.backupSettingsLocation)
            end
            break

        end

        settingsFilePath = fullfile(S.settingsLocation,settingsFname);
        backupDir = S.backupSettingsLocation;
    end % if isempty



    DEFAULT_SETTINGS = mpqc.settings.default_settings;
    if ~exist(settingsFilePath)
        fprintf('\nCan not find system settings file: making empty default file at %s\n', settingsFilePath)
        mpqc.yaml.WriteYaml(settingsFilePath,DEFAULT_SETTINGS);
    end


    % Just print to screen if no output args are specified.
    if nargout==0
        fprintf('%s\n', settingsFilePath)
        return
    end

    if nargout>0
        varargout{1} = settingsFilePath;
    end

    if nargout>1
        varargout{2} = backupDir;
    end

    if nargout>2
        varargout{3} = {settingsDirs.settingsLocation};
    end

    if nargout>3
        varargout{4} = settingsFileLocations;
    end
