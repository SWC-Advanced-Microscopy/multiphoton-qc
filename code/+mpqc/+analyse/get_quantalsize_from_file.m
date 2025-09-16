function [OUT,data] = get_quantalsize_quantalsize_from_file(fname,count_weight_gamma)
	% Get the quantal size and associated statistics from a file
	%
	% function OUT = mpqc.analyse.get_quantalsize_quantalsize_from_file(fname,count_weight_gamma)
	%
	% Purpose
	% Processes a file (likely a lens paper file) to extract the quantal size and associated
	% statistics from each recorded channel. If the channel with the structured target has
	% standard source data at the same gain, these are loaded and converted to a mean photon
	% count. The result is printed at the CLI and added to the output structure.
	%
	%
	% Inputs
	% fname - relate or absolute path to a file
	% count_weight_gamma - see help of  mpqc.analyse.compute_quantalsize
	%
	% Output
	% Structure with extensive data on the recording and also the quantal size and offset.
	% TODO -- explain this in more detail.
	%
	% See also
	% mpqc.analyse.compute_quantalsize
	% mpqc.analyse.plotPhotonFit
	%
	% Rob Campbell, SWC AMF, initial commit


	if nargin<1
		[fname,pathToFile] = uigetfile('*.tif');
		fname = fullfile(pathToFile,fname);
	end

	if nargin<2
		count_weight_gamma = [];
	end

	if ~exist(fname,'file')
		fprintf('Can not find file %s\n',fname)
		return
	end

	[im,metadata]=mpqc.tools.scanImage_stackLoad(fname,false); %Do not subtract the offset

	nChans = length(metadata.channelSave);

	pathToFile = fileparts(fname);

	ssFiles = getStandardSourceFiles(pathToFile);

	for ii=1:nChans
		% Run the analysis
		tChan = im(:,:,ii:nChans:end);
		[OUT(ii),data(ii)] = mpqc.analyse.compute_quantalsize(tChan,count_weight_gamma);

		% Fill in extra metadata
		OUT(ii).channel = metadata.channelSave(ii);
		OUT(ii).filename = fname;
		OUT(ii).gain = mpqc.report.PMT_gain_from_fname(fname);

		% Find and fit standard source if present
		t_ssFiles = ssFiles(contains(ssFiles,sprintf('_%dV_',OUT(ii).gain)));
		if ~isempty(t_ssFiles)
			OUT(ii).standard_source_results = convert_standardSource(t_ssFiles,OUT(ii));
		end
	end

end


function ss_files = getStandardSourceFiles(tDir)
	% Return cell array of standard source file names
	ss_files = dir(fullfile(tDir,'*_standard_light_source_*.tif'));
	ss_files = {ss_files(:).name};
end


function ssResults = convert_standardSource(t_ssFiles,OUT)


	for tFileInd = 1:length(t_ssFiles)
		fname = t_ssFiles{tFileInd};
		[im,metadata]=mpqc.tools.scanImage_stackLoad(fname,false);

		nChans = length(metadata.channelSave);
		n=1; %index

		for ii=1:nChans
			tChan = metadata.channelSave(ii);

			if OUT.channel ~= tChan
				continue
			end

			tData = im(:,:,ii:nChans:end);
			ssResults(n).fname = fname;
			ssResults(n).channel = tChan;
			ssResults(n).meanPhotonCount = mpqc.analyse.convertImageToPhotons(mean(tData(:)), OUT);

			fprintf('%s, Chan %d, mean photon count: %0.2f\n', ...
				fname, tChan, ssResults(n).meanPhotonCount)

			n = n+1;

		end % ii
	end % tFileInd

end
