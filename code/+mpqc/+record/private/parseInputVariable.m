function out = parseInputVariable(varargin)
    % Parse power and wavelength, and if required, depth and step size
    %
    % out = parseInputVariable('param1', val1, ...)
    %
    % Purpose
    % Recording functions need the power and wavelength, and sometimes depth and step
    % size, supplied by the user. The user may either do this by supplying input args as
    % parameter/value pairs or, if they do not, as interactive inputs that have a default
    % value. This function handles this. This function is called by the mpqc.record.
    % functions such as mpqc.record.uniform_slide, .lens_paper, and .PSF.
    % The behavior of this function depends on the function that called it. e.g. it
    % will not return depthMicrons and stepSize if it was called by the lens paper
    % recording function.
    %
    %
    % Inputs (optional param/val pairs)
    %  'wavelength' - Excitation wavelength of the laser. Defined in nm.
    %  'power' - Power at the sample. Defined in mW.
    %  'depthMicrons' - The number of microns over which to take a Z stack. Defined in um.
    %  'stepSize' - Step size between optical planes in a z-stack. Defined in um.
    %  'beamIndex' - The index of beam to be calibrated.
    %
    % Notes:
    %  1. All inputs apart from stepSize are rounded to the nearest whole number.
    %  2. Any param/val pairs other than the above are returned in the output structure
    %     without any processing or checks.
    %
    %
    % Outputs
    % out - A structure containing the user's choices.
    % e.g.
    %  out.wavelength = 920;
    %  out.power = 20;
    %
    %
    % Isabell Whiteley, SWC AMF, initial commit 2024


    % Make the inputParser object
    params = inputParser;
    params.CaseSensitive = false; % So we do not have to be case sensitive
    params.KeepUnmatched = true;
    % add parameters
    params.addParameter('wavelength', [], @(x) isnumeric(x));
    params.addParameter('power', [], @(x) isnumeric(x));
    params.addParameter('depthMicrons', [], @(x) isnumeric(x));
    params.addParameter('stepSize', [], @(x) isnumeric(x));
    params.addParameter('beamIndex', [], @(x) isnumeric(x));

    % Parse the input arguments
    params.parse(varargin{:});

    % Extract the variables (both the ones defined above and any the user happens to add)
    tFields = fields(params.Results);
    for ii=1:length(tFields)
        out.(tFields{ii}) = params.Results.(tFields{ii});
    end

    tFields = fields(params.Unmatched);
    for ii=1:length(tFields)
        out.(tFields{ii}) = params.Unmatched.(tFields{ii});
    end




    % Used to determine the identity of the calling function
    dStack = dbstack;
    if length(dStack)>1
        callerFile = dStack(2).file;
    else
        callerFile = '';
    end

    % Interactively handle each input argument if it was not supplied as a param/val pair
    if isempty(params.Results.depthMicrons) && strcmp(callerFile,'PSF.m')
        default=20;
        txt = sprintf('Please enter depth (um) [%d]: ',default);
        out.depthMicrons =  round(parseResponse(txt,default));
    end

    if isempty(params.Results.stepSize) && strcmp(callerFile,'PSF.m')
        default=0.25;
        txt = sprintf('Please enter step size (um) [%0.3f]: ',default);
        out.stepSize = parseResponse(txt,default);
    end

    if isempty(params.Results.wavelength)
        default=780;
        txt = sprintf('Please enter wavelength (nm) [%d]: ',default);
        out.wavelength = round(parseResponse(txt,default));
    end

    if isempty(params.Results.beamIndex) && strcmp(callerFile,'power.m')
        if evalin('caller', 'numBeams') == 1
            out.beamIndex = 1; % Only one beam is configured in ScanImage 
        else
            default=1;
            txt = sprintf('Please enter index of beam for calibration [%d]: ',default);
            out.beamIndex = round(parseResponse(txt,default));
        end
    end

    if isempty(params.Results.power) && ~strcmp(callerFile,'power.m')
        default=75;
        txt = sprintf('Please enter power (mW) [%d]: ',default);
        out.power = round(parseResponse(txt,default));
    end



function response = parseResponse(promptString,default)
    % Conducts an interactive prompt to help the user choose a value

    response = [];
    while isempty(response)
        response = input(promptString,'s');
        if isempty(response)
            response = default;
        else
            response = str2num(response);
        end
   end

