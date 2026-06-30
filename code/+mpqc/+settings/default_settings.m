function [settings,setTests] = default_settings
    % Return a set of default system settings to write to a file in the settings directory
    %
    % function [settings, setTests] = mpqc.settings.default_settings
    %
    % Purpose
    % Defines the default settings of MPQC. These are what are written to YAML when there
    % is no settings file. If the settings file has invalid values, these are corrected with
    % values found here. The function also defines a structure with identical fields but with
    % function handles defining tests for each item. These functions are defined as static
    % methods in the class settingsValuesTests
    %
    % Inputs
    % none
    %
    % Outputs
    % settings - default_settings
    % setTests - structure with tests required for each field as a cell array of function
    %               handles. The functions are static methods of the class settingsValuesTests.
    %
    % Rob Campbell, SWC AMF, initial commit 2023


    % Import the functions we use for checking that stuff is valid
    import mpqc.settings.settingsValuesTests.*


    %% NOTE!
    % In each case, every default setting value is defined on the first line then a test for its validity
    % is on the line after it.


    % General microscope details
    settings.microscope.name = 'ENTER_MICROSCOPE_NAME';  % Default value
    setTests.microscope.name = {@check_ischar}; % Tests of the default value

    settings.microscope.roomNumber = 'ENTER_ROOM_NUMBER';  % Default value
    setTests.microscope.roomNumber = {}; % Tests of the default value


    % Objective
    settings.objective.name = 'ENTER_MANUFACTURER 10x NA 0.8';
    setTests.objective.name = {@check_ischar};

    settings.objective.serialNumber = 'ENTER_OR_NONE';
    setTests.objective.serialNumber = {}; % No test


    % PMTs
    settings.PMT_1.model = 'ENTER_MANUFACTURER_AND_MODEL';
    setTests.PMT_1.model = {@check_ischar};
    settings.PMT_1.serialNumber = 'ENTER_VALUE';
    setTests.PMT_1.serialNumber = {};
    settings.PMT_1.microscopeChannelName = 'ENTER_FRIENDLY_CHAN_NAME e.g. Chan 1 Green';
    setTests.PMT_1.microscopeChannelName = {@check_ischar};
    settings.PMT_1.bandPassFilter = 'ENTER_PART_NUMBER_MUST_BE_STRING';
    setTests.PMT_1.bandPassFilter = {@check_ischar};

    settings.PMT_2.model = '';
    setTests.PMT_2.model = {@check_ischar};
    settings.PMT_2.serialNumber = '';
    setTests.PMT_2.serialNumber = {};
    settings.PMT_2.microscopeChannelName = '';
    setTests.PMT_2.microscopeChannelName = {@check_ischar};
    settings.PMT_2.bandPassFilter = '';
    setTests.PMT_2.bandPassFilter = {@check_ischar};


    settings.PMT_3.model = '';
    setTests.PMT_3.model = {@check_ischar};
    settings.PMT_3.serialNumber = '';
    setTests.PMT_3.serialNumber = {};
    settings.PMT_3.microscopeChannelName = '';
    setTests.PMT_3.microscopeChannelName = {@check_ischar};
    settings.PMT_3.bandPassFilter = '';
    setTests.PMT_3.bandPassFilter = {@check_ischar};


    settings.PMT_4.model = '';
    setTests.PMT_4.model = {@check_ischar};
    settings.PMT_4.serialNumber = '';
    setTests.PMT_4.serialNumber = {};
    settings.PMT_4.microscopeChannelName = '';
    setTests.PMT_4.microscopeChannelName = {@check_ischar};
    settings.PMT_4.bandPassFilter = '';
    setTests.PMT_4.bandPassFilter = {@check_ischar};


    %% Lasers
    settings.imagingLaser_1.model = 'ENTER_MANUFACTURER_AND_MODEL';
    setTests.imagingLaser_1.model = {@check_ischar};
    settings.imagingLaser_1.serialNumber = 'ENTER_SERIAL_NUMBER';
    setTests.imagingLaser_1.serialNumber = {};

    settings.imagingLaser_2.model = '';
    setTests.imagingLaser_2.model = {@check_ischar};
    settings.imagingLaser_2.serialNumber = '';
    setTests.imagingLaser_2.serialNumber = {};

    settings.imagingLaser_3.model = '';
    setTests.imagingLaser_3.model = {@check_ischar};
    settings.imagingLaser_3.serialNumber = '';
    setTests.imagingLaser_3.serialNumber = {};

    %% QC tools
    settings.QC.sourceIDs={};
    setTests.QC.sourceIDs = {@check_isCellArrayOfStrings};

end % default_settings
