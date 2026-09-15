function MDFbackup
% Saves a backup of the machine data file.
%
% mpqc.record.MDFbackup
%
% Purpose: To create a back up or version controlled copy of the machine
% data file. Scanimage saves a new version with no back up produced each
% time it is opened. This can be prolematic if there are any issues with
% the NI DAQ or VDAQ and wiring placements are lost.
%
% Note: scanimage must be open
%
% Isabell Whiteley, SWC AMF 2026


 % Create 'diagnostic' directory in the user's desktop
    saveDir = mpqc.tools.makeTodaysDataDirectory;
    if isempty(saveDir)
        return
    end
    
% Find the MDF
MDFlocationPath = most.MachineDataFile.getInstance.fileName;

[filepath,name,ext] = fileparts(MDFlocationPath);
targetFile = fullfile(saveDir, sprintf('%s_%s%s', name, datestr(now,'yyyymmdd_HHMMSS'), ext));
copyfile(MDFlocationPath, targetFile);

end