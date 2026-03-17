function OUT = calculate_PMT_SNR(dataDir)
	% Get the PMT SNR ** FUNCTION UNDER DEVELOPMENT **
	%
	% function OUT = calculate_PMT_SNR(dataDir)
	%
	% Purpose
	%
	%
	% Inputs
	% dataDir - if missing look in the current directory
	%
	% Output
	% Structure with extensive data on the recording and also the quantal size and offset.
	%
	%
	% Rob Campbell, SWC AMF, initial commit


	if nargin<1 || isempty(dataDir)
		dataDir = pwd;
	end

	if ~exist(dataDir,'dir')
		fprintf('Can not find directory %s\n',dataDir)
		return
	end



	ssFiles = getStandardSourceFiles(dataDir);
	drFiles = getDarkResponsedSourceFiles(dataDir);


	% Issue warning if the numbers do not match
	if length(ssFiles) ~= length(drFiles)
		fprintf('Warning: there are %d standard light source files but %d dark noise files\n', ...
			length(ssFiles), length(drFiles))
	end


	for ii=1:length(ssFiles)

		% Get the corresponding dark noise file name
		PMT_gain = mpqc.report.PMT_gain_from_fname(ssFiles{ii});
	    testName = sprintf('*dark_response_*%dV_*.tif', PMT_gain);

	    tDarkNoise = dir(fullfile(dataDir,testName));

	    if isempty(tDarkNoise)
	    	fprintf('Failed to find dark noise for standard source gain=%d\n', PMT_gain)
	    	continue
	    end

	    if length(tDarkNoise)>1
	    	fprintf('Found %d dark noise files for standard source gain=%d\n',...
		    	length(tDarkNoise),PMT_gain)
	    	fprintf('There should be just one. Skipping... \n')
	    	continue
	    end


	    % Load the standard source (light) data and the dark noise
		[ss_im,metadata]=loadImDeinterleaved(fullfile(dataDir,ssFiles{ii}));
		[dr_im,metadata]=loadImDeinterleaved(fullfile(dataDir,drFiles{ii}));

		% Subtract the offset from the standard source data
		offset = mean(dr_im,[1,2,4]);

		data = ss_im-offset;

		% calculate mean and variance of each channel
		t_mu = squeeze(mean(data,[1,2,4]));
		t_var = squeeze(var(data,[],[1,2,4]));
		t_std = squeeze(std(data,[],[1,2,4]));

		SNR = m;

		% Run the analysis

		%nChans = length(metadata.channelSave);

		% Fill in extra metadata
		OUT(ii).filename = ssFiles{ii};
		OUT(ii).gain = PMT_gain;
		OUT(ii).


		%%%% Find and fit standard source if present
		%% t_ssFiles = ssFiles(contains(ssFiles,sprintf('_%dV_',OUT(ii).gain)));

	end

end


function ss_files = getStandardSourceFiles(tDir)
	% Return cell array of standard source file names
	ss_files = dir(fullfile(tDir,'*_standard_light_source_*.tif'));
	ss_files = {ss_files(:).name};
end

function dr_files = getDarkResponsedSourceFiles(tDir)
	% Return cell array of dark response file names
	dr_files = dir(fullfile(tDir,'*_dark_response_*.tif'));
	dr_files = {dr_files(:).name};
end

function [im,metadata] = loadImDeinterleaved(pathToFile)
	% load scanimage data and return de-interleaved:
	% [rows,cols,channel,frames]
	% Also convert to double

	[im,metadata]=mpqc.tools.scanImage_stackLoad(fullfile(dataDir,ssFiles{ii}),false);

	im = double(im);
	imS = size(im);

	nChans = length(metadata.channelSave);

	im = reshape(im, imS(1), imS(2), nChans, []);


