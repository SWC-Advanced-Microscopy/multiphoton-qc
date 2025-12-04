function [result,dataForFit] = compute_quantalsize(frames, count_weight_gamma, min_count_proportion)
    % Compute photon quantal_size
    %
    % function result = mpqc.analyse.compute_quantalsize(frames, count_weight_gamma, min_count_proportion)
    %
    % Purpose
    % Compute the number of photons per pixel (the quantal size or system sensitivity) from
    % scanning microscope data. The input data (frames) should be a time series or Z-stack
    % that contains at least about 100 frames. The result can be affected by fluctuations in
    % laser power or contamination be electrical noise. Ensure the laser has stablised.
    % Clip beginning frames if necessary. Using moderately high PMT gains may be needed
    % if your system has appreciable noise after the amplification stage. e.g. with
    % multi-alkali PMTs set the gain to at least around 600 to 650V. Photons/pixel should
    % be constant over gain. You can generate a photons/pixel as a function of gain plot to
    % verify that your chosen gain is suitable.
    %
    % ** Likely the easiest function to use is mpqc.analyse.get_quantalsize_from_file **
    %
    %
    % Inputs (required)
    % frames: A time series or z stack in the format (height, width, time). Your data
    %   *must* contain a range of mean values and so ought to originate from imaging a
    %   structured sample. A self-luminous test target is not suitable, as all pixels will
    %   have the same mean value.
    %
    % Input (optional)
    %   count_weight_gamma: This setting influences the weighting of each count value in
    %     the robust regression. The default is 0.2. Use a weight of 0 (or a very small value)
    %     to weigh each intensity level equally. Use 1.0 to weigh each intensity in proportion
    %     to pixel counts. If unsure leave as default. You can verify the effect of this
    %     parameter by checking the fit with the plotPhotonFit function. NOTE: it seems
    %     it is necessary to set the coef to 2 in order to get a satisfactory fit with some
    %     datasets.
    %  min_count_proportion: This setting influences the magnitude of the smallest values
    %     that contribute to the fit. By default it is 0.07, meaning that the smallest values
    %     used for the fit will never be less than 7% of the maximum value used. If you see
    %     a peculiar shape in the distribution for low values with plotPhotonFit and this is
    %     throwing off the fit, you can increase the magnitude of this parameter.
    %
    %
    % Outputs
    %  result: A struct with fields:
    %         - model: The fitted regression model.
    %         - min_intensity: Minimum intensity used.
    %         - max_intensity: Maximum intensity used.
    %         - variance: Variances at intensity levels.
    %         - quantal_size: Sensitivity.
    %         - zero_level: X-intercept.
    %         - photons_per_pixel: The mean number of photons pixel in "frames"
    % dataForFit: The data that go into the fit
    %
    %
    % Example
    %  Run on a subset of frames then convert whole stack to photons. Note that the number
    %  of photons per pixel is typically very small. If taken from a multi-photon system,
    %  "Frames" is likely going to be int16 or possibly uint16. Once converted to photons
    %  the numbers will all be >=0 and likely very small. Data can therefore be converted
    %  to uint8. Note that to convert to photons we first subtract the zero (the X axis
    %  intercept) then divide by the quantal size.
    %
    %  >> frames = mpqc.tools.scanImage_stackLoad(fname);
    %  % We acquired two channels and so to analyse channel 2 we do:
    %  >> result = mpqc.analyse.compute_quantalsize(frames(:,:,2:2:end));
    %  >> mpqc.analyse.plotPhotonFit(OUT)
    %
    %
    %
    % Acknowledgements
    % This function is pretty much a direct copy of the Python work found here:
    % https://github.com/datajoint/anscombe-numcodecs by Dimitri Yatsenko
    % This function is based on work by Dimitri Yatsenko.
    %
    %
    % See also:
    % mpqc.analyse.get_quantalsize_quantalsize_from_file
    % mpqc.analyse.plotPhotonFit
    %
    %
    % Rob Campbell, SWC AMF, February 2025



    if nargin<2 || isempty(count_weight_gamma)
        count_weight_gamma = 0.2;
    end

    if nargin<3 || isempty(min_count_proportion)
        min_count_proportion = 0.05;
    end

    assert(ndims(frames) == 3, 'Input variable "frames" must have three dimensions.')

    %assert(count_weight_gamma>=0 && count_weight_gamma<=1, ...
    %    'count_weight_gamma should be between 0 and 1')

    % Esnure all values are positive
    min_pix_val = min(frames(:));
    frames = frames - min_pix_val;

    % Check if a large proportion of numbers are <0 and warn if so
    propUnderZero = sum(frames(:)<0)/length(frames(:));
    if propUnderZero>0.3
        fprintf(['%d%% of the raw values are negative. ', ...
            'Consider correcting: negative values are discarded.\n'], ...
            round(propUnderZero*100) )
    end

    % Replace all numbers lower than 0 with 0.
    frames = double(max(0, frames));


    % "intensity" will be used to determine mean values and "difference" will be
    % used to calculate the variance.
    intensity = floor((frames(:,:,1:end-1) + frames(:,:,2:end) + 1) / 2);
    difference = frames(:,:,1:end-1) - frames(:,:,2:end);

     % Convert to vectors
    intensity = intensity(:);
    difference = difference(:);

    % Determine the number of counts at each raw value
    % (histcounts) generates artifacts: do not use
    counts = accumarray(intensity(:) + 1, 1);
    mean_counts = mean(counts);

    % Calculate the longest stretch of values that are >1 of the
    % Note: since the number of bins in counts is equal to the number of possible
    % values, an index value of, say, 100 is equal to a count value of 100.
    [ind_start,ind_stop] = longest_run(counts > 0.01 * mean_counts);

    % Tweak the start index value so that it is never less than a defined proportion of
    % maximum value
    ind_start = max(floor(min_count_proportion * ind_stop), ind_start);

    assert((ind_stop - ind_start) > 100, ...
        ['The image does not have a sufficient range of intensities to compute the ', ...
        'noise transfer function.']);


    %  Keep only count values that are within range.
    if length(counts)<ind_stop
        ind_stop = length(counts);
    end
    counts = counts(ind_start+1:ind_stop);

    % Calculate the variance of each count value
    idx = (intensity >= ind_start) & (intensity < ind_stop);
    variance = accumarray(intensity(idx) - ind_start + 1, (difference(idx) .^ 2) / 2, [length(counts), 1]) ./ counts;



    %% Fit
    % We fit the variances using the mean-centered count values. If we do not
    % mean-center, MATLAB returns intercept values are wrong: forced to be zero.
    % We eventually return the non-centered intercept, so this process is transparent.
    % to the user. We use a robust fit and weight each point by the number of counts
    % multiplied by a scaling factor
    X = (ind_start:ind_stop-1)';
    W = counts .^ count_weight_gamma;
    useCentered = false;

    if useCentered
        Xvals = X - mean(X); % Must run on centered data or the intercept is 0
    else
        Xvals = X;
    end

    try
        coefs_raw = robustfit(Xvals, variance, @(r) W , []);
        weightedFit = true;
    catch
        fprintf('Weighted fit failed. Running fit with plain "huber" regression\n')
        weightedFit = false;
        coefs_raw = robustfit(Xvals, variance, 'huber');
    end



    % Uncenter
    if useCentered
        % Adjust the intercept by adding the mean of X
        coefs = [coefs_raw(1) - mean(X) * coefs_raw(2); coefs_raw(2:end)];
    else
        coefs = coefs_raw;
    end


    quantal_size = coefs(2);
    zero_level = -coefs(1) / quantal_size;

    % Some of the fields can be filled in by other routines. The "user" field is
    % a catch-all to put whatever the user wants
    result = struct(...
        'model', coefs, ...
        'count_weight_gamma', count_weight_gamma, ...
        'min_count_proportion', min_count_proportion, ...
        'weighted_fit', weightedFit, ...
        'counts', counts, ...
        'min_intensity', ind_start, ...
        'max_intensity', ind_stop, ...
        'variance', variance, ...
        'quantal_size', quantal_size, ...
        'zero_level', zero_level, ...
        'photons_per_pixel', mean(single(frames(:)-zero_level)/quantal_size), ...
        'image_min_pixel_value', min_pix_val,...
        'frame_xy_size', size(frames,[1,2]), ...
        'num_frames', size(frames,3), ...
        'gain', [], ...
        'filename', '', ...
        'channel', [], ...
        'standard_source_results', struct, ...
        'user', []);

    dataForFit.X = Xvals;
    dataForFit.Y = variance;
    dataForFit.weights =  W;


end

function [ind_start,ind_stop] = longest_run(bool_array)
    %    Find the longest contiguous segment of True values inside bool_array.
    %
    % Inputs
    %    bool_array: 1d boolean array.
    %
    % Returns
    %    ind_start - index at the start of the longest contiguous block of True values.
    %    ind_stop - index at the end of the longest contiguous block of True values.
    %
    %

    bool_array = [0,bool_array(:)',0];

    step = diff(bool_array);
    on = find(step == 1);
    off = find(step == -1);
    [~,ind]=max(off-on);

    ind_start = on(ind);
    ind_stop = off(ind);
end
