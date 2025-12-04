function [imStack,metadata] = scanImage_stackLoad(fileName,subtractOffset)
% load a z-stack from ScanImage for analysis return also in a structure useful metadata
%
%
% function imStack = mpqc.tools.scanImage_stackLoad(fileName,subtractOffset)
%
% Purpose
% Return z stack and metadata. By default subtracts any offset if needed.
%
% Inputs (required)
% fileName - string defining the file name to load.
%
% Inputs (optional)
% subtractOffset - if true, subtract the offset from the image data. Default is true.
%              NOTE: setting to false only has an effect if the data were saved without
%              offset subtraction
%
% Outputs
% imStack - 3D stack
% metadata - useful metadata from scanimage header. Adds micsPerPixelXY.
%
% Example
% [imS,metadata] = mpqc.tools.scanImageStackLoad('Bead_*.tif')
%
%
%
% Rob Campbell, SWC AMF, initial commit 2022


    imStack = [];
    metadata = [];

    if ~exist(fileName,'file')
        fprintf('%s does not exist. scanImage_stackLoad will not load it.\n',fileName)
        return
    end

    if nargin<2 || isempty(subtractOffset)
        subtractOffset = true;
    end

    imStack = mpqc.tools.load3Dtiff(fileName);



    % Pull out the voxel size and other useful information
    metadata=sibridge.readTifHeader(fileName);
    if isempty(metadata)
        fprintf('\n\n *** TIFF header is missing ScanImage meta-data. Is this a ScanImage TIFF? *** \n\n')
        return
    end

    % Calculate the FOV and add this to the meta-data
    fov=diff(metadata.imagingFovUm(1:2));
    metadata.micsPerPixelXY = fov/metadata.linesPerFrame;


    if subtractOffset == false
        if any(metadata.channelSubtractOffset)
            fprintf('Data saved with offset subtracted. Returning offset subtracted data.\n')
        end
        return
    end

    % subtract the offset if needed
    savedChans = metadata.channelSave; % The acquired channels
    offset = metadata.channelOffset(savedChans); % Use acquired channels to index offsets
    offsetSubtracted = metadata.channelSubtractOffset;


    %expand offsets so they are the same length as the image data
    offset = repmat(offset,1,size(imStack,3)/length(offset));

    offsetSubtracted = repmat(offsetSubtracted,1,size(imStack,3)/length(savedChans));

    for ii=1:size(imStack,3)
        if offsetSubtracted(ii) == 0
            imStack(:,:,ii) = imStack(:,:,ii)-offset(ii);
        end
    end


