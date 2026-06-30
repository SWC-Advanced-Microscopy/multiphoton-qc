function im_p = convertImageToPhotons(im,STATS)
    % Convert image im to photons using coefficients in STATS
    %
    % im_p = mpqc.analyse convertImageToPhotons(im,STATS)
    %
    % Inputs
    % im - image stack, single frame, or average frame of an image time series.
    %
    % Outputs
    % im_p - the photon-convert image
    %
    %
    % Rob Campbell, SWC AMF, initial commit February 2025

    im_p = (im - STATS.image_min_pixel_value - STATS.zero_level) / STATS.quantal_size;

