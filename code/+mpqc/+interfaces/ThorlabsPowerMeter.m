classdef ThorlabsPowerMeter < matlab.mixin.Copyable
    % ThorlabsPowerMeter Matlab class to control Thorlabs power meters
    %
    %   Interface class for ThorLabs power meters. This is a 'wrapper' to control
    %   Thorlabs devices via the Thorlabs .NET DLLs.
    %
    %   User Instructions:
    %       1. Download the Optical Power Monitor from the Thorlabs website:
    %       https://www.thorlabs.com/software_pages/ViewSoftwarePage.cfm?Code=OPM
    %       [The latest version is 4.0.4100.700 - Accessed on 01 SEP 2022]
    %       
    %       !! This must be done prior to connecting the powermeter or
    %       incorrect drivers are likely installed. If the code gives this error
    %       "Device seems not to be connected: deviceNET property is empty"
    %       incorrect driver is the likely cause. Disconnect powermeter, 
    %       delete powermeter drivers, and reconnect powermeter to fix.
    %
    %       2. Read the manual in the installation folder or the sofware help page
    %       https://www.thorlabs.com/software/MUC/OPM/v3.0/TL_OPM_V3.0_web-secured.pdf
    %
    %       3. Following the instructions in section 9: Write Your Own Application
    %       The common path of the *.dll files on Windows is:
    %       C:\Program Files\IVI Foundation\VISA\VisaCom64\Primary Interop Assemblies\Thorlabs.TLPM_64.Interop.dll
    %
    %       4. This class needs only the .NET wrapper dll, so follow the instruction for C#/.Net
    %
    %       5. Edit MOTORPATHDEFAULT below to point to the location of the DLLs
    %
    %       6. Connect your Power Meter with sensor to the PC USB port and power it on.
    %
    %       7. Please refer to the examples provided
    %
    %   For developers:
    %   The definition for all the classes can be found in the C# example
    %   provided by ThorLabs. (Shipped together with the software.)
    %
    %
    %
    %   EXAMPLES
    %
    %   Connecting 01:
    %   meter_list=mpqc.interfaces.ThorlabsPowerMeter;              % Initiate the meter_list
    %   DeviceDescription=meter_list.listdevices;               	% List available device(s)
    %   test_meter=meter_list.connect(DeviceDescription);           % Connect single/the first devices
    %
    %   Connecting 02:
    %   test_meter=meter_list.connect(DeviceDescription,1);         % Connect single/the first devices
    %   test_meter.setWaveLength(635);                              % Set sensor wavelength
    %   test_meter.setDispBrightness(0.3);                          % Set display brightness
    %   test_meter.setAttenuation(0);                               % Set Attenuation
    %   test_meter.sensorInfo;                                      % Retrive the sensor info
    %
    %   Setting auto-range:
    %   test_meter.setPowerAutoRange(1);                            % Set Autorange
    %   test_meter.setPowerRange(0.01);                            % Set manual range
    %
    %   Setting other values:
    %   test_meter.setAverageTime(0.01);                            % Set average time for the measurement
    %   test_meter.setTimeout(1000);                                % Set timeout value
    %
    %   PMT400 only:
    %   test_meter.darkAdjust;
    %   test_meter.getDarkOffset;
    %
    %   Reading power:
    %   test_meter.readPower
    %   test_meter.readVoltage
    %
    %
    %   Disconnect and release:
    %   test_meter.disconnect
    %
    %   -----------------
    %
    %   Author: Zimo Zhao
    %   Dept. Engineering Science, University of Oxford, Oxford OX1 3PJ, UK
    %   Email: zimo.zhao@eng.ox.ac.uk (please email issues and bugs)
    %   Website: https://eng.ox.ac.uk/smp/
    %   GitHub: https://github.com/Tinyblack/Matlab-Driver-for-Thorlabs-power-meter
    %
    %   Initially Developed On:
    %       Optical Power Monitor
    %           Application 3.1.3778.562
    %           TLPM__32 5.1.3754.327
    %       Matlab
    %           2020b
    %
    %   Test pass:
    %       Optical Power Monitor
    %           Application 4.0.4100.700
    %           TLPMX__32 5.3.4101.525
    %       Matlab
    %           2022a
    %
    %   Version History:
    %   1.00 ----- 21 May 2021 ----- Initial Release
    %   1.01 ----- 17 Aug 2021 ----- Clarify the way of utilizing *.dll files
    %   2.00 ----- 27 Aug 2021 ----- Support connection of multiple power meters
    %   2.01 ----- 26 Sep 2021 ----- Add force connection function to bypass the device availability check.
    %   3.00 ----- 01 Feb 2022 ----- Add functions: setPowerRange, setPowerAutoRange, setTimeout, setAverageTime, updateReading_V
    %   3.10 ----- 01 SEP 2022 ----- Test the script on latest TLPM driver and MATLAB. Some bugs are corrected as well
    %   3.11 ----- 16 JUN 2025 ----- Fail gracefully if DLL not installed. Minor tidy to docs. [RAAC]


    properties (Hidden)
        % Path to .net *.dll files (edit as appropriate)
        % pwd --- Current working directory of this file
        % (depending on the location where you put this file)
        % This line points to folder 'Thorlabs_DotNet_dll' under the same directory
        % Comment out this line and uncomment next line to use customized dll file directory
        % METERPATHDEFAULT=[pwd '\Thorlabs_DotNet_dll\'];
        METERPATHDEFAULT='C:\Program Files (x86)\Microsoft.NET\Primary Interop Assemblies\';

        %   *.dll files to be loaded
        %
        % NOTE
        %   No significant difference was noticed between
        %   "Thorlabs.TLPMX_64.Interop.dll" and "Thorlabs.TLPM_64.Interop.dll"
        %   But if you are going to use "Thorlabs.TLPMX_64.Interop.dll", please change
        %
        %         TLPMDLL='Thorlabs.TLPM_64.Interop.dll';
        %         TLPMCLASSNAME='Thorlabs.TLPM_64.Interop.TLPM';
        %
        %     into
        %
        %         TLPMDLL='Thorlabs.TLPMX_64.Interop.dll';
        %         TLPMCLASSNAME='Thorlabs.TLPM_64.Interop.TLPMX';
        %
        TLPMDLL='Thorlabs.TLPM_64.Interop.dll';
        TLPMCLASSNAME='Thorlabs.TLPM_64.Interop.TLPM';
    end

    properties
        % These properties are within Matlab wrapper
        isConnected=false;          % Flag set to device index in list if device connected
        deviceList                  % Structure that is the output of listdevices
        resourceNameConnected;      % USB resource name
        sensorName;                 % Sensor name
        sensorSerialNumber;         % Sensor serial number
        sensorCalibrationMessage;   % Sensor calibration information
        sensorType;                 % Sensor type
        sensorSubType;              % Sensor subtype
        sensorFlags;                % Sensor flag
        DarkOffset_Voltage;         % (PM400 ONLY) Dark offset voltage
        DarkOffset_Voltage_Unit;    % (PM400 ONLY) Dark offset voltage unit
        meterPowerReading;          % Last power reading
        meterPowerUnit;             % Power reading unit
        meterVoltageReading;        % Last voltage reading
        meterVoltageUnit;           % Voltage reading unit
    end

    properties (Hidden)
        % These are properties within the .NET environment.
        deviceNET;                  % Device object within .NET
    end

    methods
        function obj = ThorlabsPowerMeter(deviceList)
            %ThorlabsPowerMeter Construct an instance of this class
            %   This function first loads the dlls from the path and then
            %   list all the device available. It will return a list of all
            %   the available device(s).


            success=obj.loaddlls;
            if ~success
                return
            end

            if nargin>0
                obj.deviceList = deviceList;
            else
                obj.listdevices;
            end

            if isempty(obj.deviceList)
                obj.isConnected=false;
                warning('No Resource is found, please check the connection.');
            else
                numberOfResources=length(obj.deviceList);
                fprintf('Found the following %d device(s):\r',numberOfResources);
                for ii=1:length(obj.deviceList)
                    fprintf('\t\t%d) %s\r',ii,obj.deviceList(ii).resourceName);
                end
            end

            % Read the sensor head information
            obj.sensorInfo
        end

        function delete(obj)
            %DELETE Deconstruct the instance of this class
            %   Usage: obj.delete;
            %   This function disconnects the device and exits.
            if obj.isConnected
                try
                    obj.disconnect;
                catch ME
                    warning('Failed to release the device:\n%s',ME.message);
                end
            else 
                % Cannot disconnect because device is not connected
            end

        end

        function connect(obj,resource_index,ID_Query,Reset_Device)
            %CONNECT Connect to the specified resource.
            %   Usage: obj.connect(resource);
            %   By default, it will connect to the first resource on the
            %   list [resource_index=1] with ID query [ID_Query=1] and
            %   reset [Reset_Device=1];
            %   Use
            %   obj.connect(resource,ID_Query,Reset_Device,resource_index)
            %   to modify the default values.
            arguments
                obj
                resource_index (1,1) {mustBeNumeric} = 1 % (default) First resource
                ID_Query (1,1) {mustBeNumeric} = 1 % (default) Query the ID
                Reset_Device (1,1) {mustBeNumeric} = 1 % (default) Reset
            end

            if ~obj.isConnected && ~isempty(obj.deviceList)
                try
                    % The core method to create the power meter instance
                    resource = obj.deviceList(resource_index).resourceName;
                    obj.deviceNET=Thorlabs.TLPM_64.Interop.TLPM(resource,logical(ID_Query),logical(Reset_Device));
                    fprintf('Successfully connected to:\r\t\t%s\r',resource);
                    obj.resourceNameConnected=resource(resource_index,:);
                    obj.isConnected=resource_index;
                    obj.deviceList(resource_index).DeviceAvailable=0;
                catch ME
                    error('Failed to connect the device:\n%s', ME.message);
                end
            else
                if obj.isConnected==1
                    warning('Device is connected.');
                end
                if obj.deviceList(resource_index).DeviceAvailable==0
                    warning('Device is not available.');
                end
            end
        end

        function disconnect(obj)
            %DISCONNECT Disconnect the specified resource.
            %   Usage: obj.disconnect;
            %   Disconnect the specified resource.
            if obj.isConnected && ~isempty(obj.deviceNET)
                fprintf('\tDisconnecting ... %s\r',obj.resourceNameConnected);
                try
                    obj.deviceNET.Dispose();  %Disconnect the device
                    obj.isConnected=false;
                    fprintf('\tDevice Released Properly.\r\r');
                catch
                    warning('Unable to disconnect device.');
                end
            else % Cannot disconnect because device not connected
                warning('Device not connected.')
            end
        end

        function setAverageTime(obj,AverageTime)
            %SETAVERAGETIME Set the sensor average time.
            %   Usage: obj.setAverageTime(AverageTime);
            %   Set the sensor average time. This method will check the input
            %   and force it to a vaild value if it is out of the range.
            %
            % NOTE: setting this to over about 0.5 seems to cause VISA errors and 
            %       the device has to be reset. 

            if ~obj.isDeviceNetConnected
                return
            end

            if AverageTime>0.5
                warning('Setting averaging time to a value over 0.5s; this might cause VISA errors')
            end

            [~,AverageTime_MIN]=obj.deviceNET.getAvgTime(1);
            [~,AverageTime_MAX]=obj.deviceNET.getAvgTime(2);
            if (AverageTime_MIN<=AverageTime && AverageTime<=AverageTime_MAX)
                obj.deviceNET.setAvgTime(AverageTime);
                fprintf('\tSet integration time to %.3f s\r',AverageTime);
            else
                if AverageTime_MIN>AverageTime
                    warning('Exceed minimum average time! Force to the minimum.');
                    obj.deviceNET.setAvgTime(AverageTime_MIN);
                    fprintf('\tSet integration time to %.3f s\r',AverageTime_MIN);
                end
                if AverageTime>AverageTime_MAX
                    warning('Exceed maximum average time! Force to the maximum.');
                    obj.deviceNET.setAvgTime(AverageTime_MAX);
                    fprintf('\tSet integration time to %.3f s\r',AverageTime_MAX);
                end
            end
        end

        function setTimeout(obj,Timeout)
            %SETTIMEOUT Set the power meter timeout value.
            %   Usage: obj.setTimeout(Timeout);
            %   Set the sensor timeout value.

            if ~obj.isDeviceNetConnected
                return
            end

            obj.deviceNET.setTimeoutValue(Timeout);
            fprintf('\tSet Timeout Value to %.4fms\r',Timeout);
        end

        function setWaveLength(obj,wavelength)
            %SETWAVELENGTH Set the sensor wavelength.
            %   Usage: obj.setWaveLength(wavelength);
            %   Set the sensor wavelength. This method will check the input
            %   and force it to a vaild value if it is out of the range.

            if ~obj.isDeviceNetConnected
                return
            end

            [~,wavelength_MIN]=obj.deviceNET.getWavelength(1);
            [~,wavelength_MAX]=obj.deviceNET.getWavelength(2);
            if (wavelength_MIN<=wavelength && wavelength<=wavelength_MAX)
                obj.deviceNET.setWavelength(wavelength);
                fprintf('\tSet wavelength to %.4f\r',wavelength);
            else
                if wavelength_MIN>wavelength
                    warning('Exceed minimum wavelength! Force to the minimum.');
                    obj.deviceNET.setWavelength(wavelength_MIN);
                    fprintf('\tSet wavelength to %.4f\r',wavelength_MIN);
                end
                if wavelength>wavelength_MAX
                    warning('Exceed maximum wavelength! Force to the maximum.');
                    obj.deviceNET.setWavelength(wavelength_MAX);
                    fprintf('\tSet wavelength to %.4f\r',wavelength_MAX);
                end
            end
        end

        function setPowerAutoRange(obj,enable)
            obj.deviceNET.getPowerRange(enable);
        end

        function setPowerRange(obj,range)
            %SETPOWERRANGE Set the sensor power range.
            %   Usage: obj.setPowerRange(range);
            %   Set the sensor power range. This method will check the input
            %   and force it to a vaild value if it is out of the range.

            if ~obj.isDeviceNetConnected
                return
            end

            [~,range_MIN]=obj.deviceNET.getPowerRange(1);
            [~,range_MAX]=obj.deviceNET.getPowerRange(2);
            if (range_MIN<=range && range<=range_MAX)
                obj.deviceNET.setPowerRange(range);
                fprintf('\tSet range to %.4f\r',range);
            else
                if range_MIN>range
                    warning('Exceed minimum range! Force to the minimum.');
                    obj.deviceNET.setPowerRange(range_MIN);
                    fprintf('\tSet range to %.4f\r',range_MIN);
                end
                if range>range_MAX
                    warning('Exceed maximum range! Force to the maximum.');
                    obj.deviceNET.setPowerRange(range_MAX);
                    fprintf('\tSet range to %.4f\r',range_MIN);
                end
            end
        end

        function setDispBrightness(obj,Brightness)
            %SETDISPBRIGHTNESS Set the display brightness.
            %   Usage: obj.setDispBrightness(Brightness);
            %   Set the display brightness. This method will check the input
            %   and force it to a vaild value if it is out of the range.

            if ~obj.isDeviceNetConnected
                return
            end

            if (0.0<=Brightness && Brightness<=1.0)
                obj.deviceNET.setDispBrightness(Brightness);
            else
                if 0.0>Brightness
                    warning('Exceed minimum brightness! Force to the minimum.');
                    Brightness=0.0;
                    obj.deviceNET.setDispBrightness(Brightness);
                end
                if Brightness>1.0
                    warning('Exceed maximum brightness! Force to the maximum.');
                    Brightness=1.0;
                    obj.deviceNET.setDispBrightness(Brightness);
                end
            end
            fprintf('Set Display Brightness to %d%%\r',Brightness*100);
        end

        function setAttenuation(obj,Attenuation)
            %SETATTENUATION Set the attenuation.
            %   Usage: obj.setAttenuation(Attenuation);
            %   Set the attenuation.

            if ~obj.isDeviceNetConnected
                return
            end

            if any(strcmp(obj.modelName,{'PM100A', 'PM100D', 'PM100USB', 'PM200', 'PM400'}))
                [~,Attenuation_MIN]=obj.deviceNET.getAttenuation(1);
                [~,Attenuation_MAX]=obj.deviceNET.getAttenuation(2);
                if (Attenuation_MIN<=Attenuation && Attenuation<=Attenuation_MAX)
                    obj.deviceNET.setAttenuation(Attenuation);
                else
                    if Attenuation_MIN>Attenuation
                        warning('Exceed minimum Attenuation! Force to the minimum.');
                        Attenuation=Attenuation_MIN;
                        obj.deviceNET.setAttenuation(Attenuation);
                    end
                    if Attenuation>Attenuation_MAX
                        warning('Exceed maximum Attenuation! Force to the maximum.');
                        Attenuation=Attenuation_MAX;
                        obj.deviceNET.setAttenuation(Attenuation);
                    end
                end
                fprintf('Set Attenuation to %.4f dB, %.4fx\r',Attenuation,10^(Attenuation/20));
            else
                warning('This command is not supported on %s.',obj.modelName);
            end
        end

        function sensorInfo=sensorInfo(obj)
            %SENSORINFO Retrive the sensor information.
            %   Usage: obj.sensorInfo;
            %   Read the information of sensor connected and store it in
            %   the properties of the object.

            if ~obj.isDeviceNetConnected
                return
            end

            for ii=1:3
                descr{ii}=System.Text.StringBuilder;
                descr{ii}.Capacity=1024;
            end
            [~,type,subtype,sensor_flag]=obj.deviceNET.getSensorInfo(descr{1}, descr{2}, descr{3});
            obj.sensorName=char(descr{1}.ToString);
            obj.sensorSerialNumber=char(descr{2}.ToString);
            obj.sensorCalibrationMessage=char(descr{3}.ToString);
            switch type
                case 0x00
                    obj.sensorType='No sensor';
                    switch subtype
                        case 0x00
                            obj.sensorSubType='No sensor';
                        otherwise
                            warning('Unknown sensor.');
                    end
                case 0x01
                    obj.sensorType='Photodiode sensor';
                    switch subtype
                        case 0x01
                            obj.sensorSubType='Photodiode adapter';
                        case 0x02
                            obj.sensorSubType='Photodiode sensor';
                        case 0x03
                            obj.sensorSubType='Photodiode sensor with integrated filter identified by position';
                        case 0x12
                            obj.sensorSubType='Photodiode sensor with temperature sensor';
                        otherwise
                            warning('Unknown sensor.');
                    end
                case 0x02
                    obj.sensorType='Thermopile sensor';
                    switch subtype
                        case 0x01
                            obj.sensorSubType='Thermopile adapter';
                        case 0x02
                            obj.sensorSubType='Thermopile sensor';
                        case 0x12
                            obj.sensorSubType='Thermopile sensor with temperature sensor';
                        otherwise
                            warning('Unknown sensor.');
                    end
                case 0x03
                    obj.sensorType='Pyroelectric sensor';
                    switch subtype
                        case 0x01
                            obj.sensorSubType='Pyroelectric adapter';
                        case 0x02
                            obj.sensorSubType='Pyroelectric sensor';
                        case 0x12
                            obj.sensorSubType='Pyroelectric sensor with temperature sensor';
                        otherwise
                            warning('Unknown sensor.');
                    end
                otherwise
                    warning('Unknown sensor.');
            end
            tag=rem(sensor_flag,16);
            switch tag
                case 0x0000

                case 0x0001
                    obj.sensorFlags=[obj.sensorFlags,'Power sensor '];
                case 0x0002
                    obj.sensorFlags=[obj.sensorFlags,'Energy sensor '];
                otherwise
                    warning('Unknown flag.');
            end
            sensor_flag=sensor_flag-tag;
            tag=rem(sensor_flag,256);
            switch tag
                case 0x0000

                case 0x0010
                    obj.sensorFlags=[obj.sensorFlags,'Responsivity settable '];
                case 0x0020
                    obj.sensorFlags=[obj.sensorFlags,'Wavelength settable '];
                case 0x0040
                    obj.sensorFlags=[obj.sensorFlags,'Time constant settable '];
                otherwise
                    warning('Unknown flag.');
            end
            sensor_flag=sensor_flag-tag;
            tag=rem(sensor_flag,256*16);
            switch tag
                case 0x0000

                case 0x0100
                    obj.sensorFlags=[obj.sensorFlags,'With Temperature sensor '];
                otherwise
                    warning('Unknown flag.');
            end
            sensorInfo.Type=obj.sensorType;
            sensorInfo.SubType=obj.sensorSubType;
            sensorInfo.Flags=obj.sensorFlags;
        end

        function powerReading = readPower(obj)
            % readPower - Return power incident on sensor in W or dB
            %
            %   Usage: P = obj.readPower;
            %
            %   Details:
            %    Return power incident on sensor and also store the last read power
            %    value in the "meterPowerReading" property. The unit the meter reads
            %    in can be found in the meterPowerUnit property.
            %

            if ~obj.isDeviceNetConnected
                return
            end

            [~,powerReading]=obj.deviceNET.measPower;
            obj.meterPowerReading = powerReading;
            [~,meterPowerUnit_]=obj.deviceNET.getPowerUnit;
            switch meterPowerUnit_
                case 0
                    obj.meterPowerUnit='W';
                case 1
                    obj.meterPowerUnit='dBm';
                otherwise
                    warning('Unknown');
            end
        end

        function voltageReading = readVoltage(obj)
            % readVoltage - return the voltage value generated by the meter
            %
            %   obj.readVoltage;
            %
            %  Details:
            %   Return the voltage value associated with the power reading. Also stores
            %   the last read value in the property "meterVoltageReading". Works only on:
            %   Only for PM100D, PM100A, PM100USB, PM160T, PM200, PM400 with certain
            %   sensors.

            if ~obj.isDeviceNetConnected
                return
            end

            if any(strcmp(obj.modelName,{'PM100D', 'PM100A', 'PM100USB', 'PM160T', 'PM200', 'PM400'}))
                try
                    [~,obj.meterVoltageReading]=obj.deviceNET.measVoltage;
                    voltageReading = obj.meterVoltageReading;
                    obj.meterVoltageUnit='V';
                catch
                    warning('Wrong sensor type for this operation');
                end
            else
                voltageReading = [];
            end
        end

        function darkAdjust(obj)
            %DARKADJUST (PM400 Only) Initiate the Zero value measurement.
            % UPDATE: confirmed working on PM100D
            %   Usage: obj.darkAdjust;
            %   Start the measurement of Zero value.

            if ~obj.isDeviceNetConnected
                return
            end
              obj.deviceNET.startDarkAdjust;
                [~,DarkState]=obj.deviceNET.getDarkAdjustState;
                while DarkState
                    [~,DarkState]=obj.deviceNET.getDarkAdjustState;
                end
        end

        function [DarkOffset_Voltage,DarkOffset_Voltage_Unit]=getDarkOffset(obj)
            %GETDARKOFFSET (PM400 Only) Read the Zero value from powermeter.
            %   Usage: [DarkOffset_Voltage,DarkOffset_Voltage_Unit]=obj.getDarkOffset;
            %   Retrive the Zero value from power meter and store it in the
            %   properties of the object

            if ~obj.isDeviceNetConnected
                return
            end

            if any(strcmp(obj.modelName,'PM400'))||any(strcmp(obj.modelName,'PM100D'))
                [~,DarkOffset_Voltage]=obj.deviceNET.getDarkOffset;
                DarkOffset_Voltage_Unit='V';
                obj.DarkOffset_Voltage=DarkOffset_Voltage;
                obj.DarkOffset_Voltage_Unit=DarkOffset_Voltage_Unit;
            else
                warning('This command is not supported on %s.',obj.modelName);
            end
        end

        function success = loaddlls(obj) % Load DLLs
            %LOADDLLS Load needed dll libraries.
            %   Usage: success = obj.loaddlls;
            %
            %   returns true if successfully loads DLL. False otherwise.

            % For the DLLs to actually work we need to also add to the Windows path the
            % directory that contains "TLPM_64.dll". We assume that everyone is on 64 bit
            % so we hard-code this here:
            dllPath = 'C:\Program Files\IVI Foundation\VISA\Win64\Bin';
            if exist(dllPath,'dir') && ~contains(getenv('PATH'),dllPath)
                setenv('PATH', [dllPath ';' getenv('PATH')]);
            end

            % If needed change the path of dll to suit you application.
            fname = fullfile(obj.METERPATHDEFAULT,obj.TLPMDLL);
            success=false;
            if exist(fname,'file')
                try
                    fprintf('Importing DLL: %s\n', fname)
                    NET.addAssembly(fname);
                    success=true;
                catch % DLLs did not load
                    error('Unable to load .NET assemblies')
                 end
            else
                fprintf(['\nFailed to load ThorLabs power meter DLL!\n', ...
                    'Can not find file %s\n\n', ...
                    'You probably need to install the power meter GUI from: \n',...
                    'https://www.thorlabs.com/software_pages/ViewSoftwarePage.cfm?Code=OPM\n\n'], ...
                    fname)
            end
        end %loaddlls

        function isConnected = isDeviceNetConnected(obj)
            if isempty(obj.deviceNET)
                fprintf('Device seems not to be connected: deviceNET property is empty\n')
                isConnected = false;
            else
                isConnected = true;
            end
        end

        function deviceList=listdevices(obj)  % Read a list of resource names
            %LISTDEVICES List available resources.
            %   Usage: obj.listdevices;
            %   Retrive all the available devices and return as a structure:
            %
            %   resourceName          -- USB resource name
            %   modelName             -- Power meter model name
            %   serialNumber          -- Power meter serial number
            %   Manufacturer          -- Power meter manufacturer
            %   DeviceAvailable       -- Power meter avaliablity
            %
            % The structure is also copied to the deviceList property

            fprintf('Looking for devices...\n')
            findResource=Thorlabs.TLPM_64.Interop.TLPM(System.IntPtr);  % Build device list
            [~,count]=findResource.findRsrc; % Get device list

            for ii=1:4
                descr{ii}=System.Text.StringBuilder;
                descr{ii}.Capacity=2048;
            end

            if count>0
                n=1;
                for ii=0:count-1
                    findResource.getRsrcName(ii,descr{1});
                    [~,Device_Available]=findResource.getRsrcInfo(ii, descr{2}, descr{3}, descr{4});
                    deviceList(n).resourceName=char(descr{1}.ToString);
                    deviceList(n).modelName=char(descr{2}.ToString);
                    deviceList(n).serialNumber=char(descr{3}.ToString);
                    deviceList(n).Manufacturer=char(descr{4}.ToString);
                    deviceList(n).DeviceAvailable=Device_Available;
                    n=n+1;
                end
            else
                deviceList = [];
            end
            obj.deviceList=deviceList;
            findResource.Dispose();
        end

    end
end

