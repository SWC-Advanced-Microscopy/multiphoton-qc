function varargout=generateMPQCreport(data_dir)
    % Makes all plots defined by the gen plot structure
    %
    % function generateMPQCreport(data_dir)
    %
    % Purpose
    % Makes all available plots from a data directory and tiles them over the screen.
    % If no input is provided, it looks in the current directory for data.
    % If you do not wish to generate a PDF report but wish to see all the figures, run
    % mpqc.report.plotAllBasic
    %
    % Inputs
    % data_dir - [optional] location of data directory.
    %
    % Outputs
    % none
    %
    % See also
    % mpqc.report.plotAllBasic
    %
    % Rob Campbell, SWC AMF, initial commit 2022


    % Bail out if the user lacks the report generator
    v = ver;
    if ~contains([v.Name] ,'MATLAB Report Generator')
        fprintf('\nThis function requires the MATLAB Report Generator to be installed\n')
        fprintf('If you can not install it, see mpqc.report.plotAllBasic\n\n')
        return
    end


    % Generate the structure we will use to make the plots
    if nargin<1
        data_dir = pwd;
    end

    GEN=mpqc.report.plot_functions_generator(data_dir);
    if isempty(GEN)
        return
    end

    % Java imports
    import mlreportgen.dom.*;
    import mlreportgen.report.*;


    % Preferences and useful instances of classes
    [~,t_dir]=fileparts(data_dir);
    t_dir = strrep(t_dir,' ', ''); %strip spaces
    pdfFileName = [t_dir,'_Session_Report'];
    rpt = Report(pdfFileName,'pdf');

    rpt.Layout.Landscape = true;
    imgStyle = {ScaleToFit(true)};
    br = PageBreak();


    %% Intro to the report by summarizing whatever information we can automatically generate
    chapter = Chapter('Title', 'Introduction');
    p1 = Paragraph(mpqc.report.generate_summary_text);


    add(chapter,p1)
    add(rpt,chapter)



    %% Electrical noise at each channel
    f=find(strcmp({GEN.type},'electrical_noise'));
    if ~isempty(f)
        chapter = Chapter('Title', 'Electrical Noise At Each Channel');

        for ii=1:length(f)
            GEN(f(ii)).plotting_func(GEN(f(ii)).full_path_to_data)

            fig = Figure();

            figImg = Image(getSnapshotImage(fig, rpt));
            figImg.Style = imgStyle;
            delete(gcf);

            add(chapter, figImg);
        end

        add(rpt,chapter)
    end


    %% Standard light source
    f=find(strcmp({GEN.type},'standard_source'));
    if ~isempty(f)
        chapter = Chapter('Title', 'Response to standard light source');

        for ii=1:length(f)
            GEN(f(ii)).plotting_func(GEN(f(ii)).full_path_to_data)

            fig = Figure();

            figImg = Image(getSnapshotImage(fig, rpt));
            figImg.Style = imgStyle;
            delete(gcf);

            add(chapter, figImg);
        end

        add(rpt,chapter)
    end


    %% Uniform slide image
    f=find(strcmp({GEN.type},'uniform_slide'));
    if ~isempty(f)
        chapter = Chapter('Title', 'Image uniformity');

        for ii=1:length(f)
            GEN(f(ii)).plotting_func(GEN(f(ii)).full_path_to_data);
            fig = Figure();

            figImg = Image(getSnapshotImage(fig, rpt));
            figImg.Style = imgStyle;
            delete(gcf);

            add(chapter, figImg);
        end

        add(rpt, chapter);
    end


    %% Image stability with time
    f=find(strcmp({GEN.type},'laser_stability'));
    if ~isempty(f)
        chapter = Chapter('Title', 'Laser stability');

        for ii=1:length(f)
            txt = GEN(f(ii)).plotting_func(GEN(f(ii)).full_path_to_data);

            p1 = Paragraph(txt);
            p1.Style = {FontSize('10pt')};
            fig = Figure();
            figImg = Image(getSnapshotImage(fig, rpt));
            figImg.Style = {Width('6.5in')};%imgStyle;

            delete(gcf);

            % Organise into a table
            tbl = Table({p1,figImg});
            tbl.Style={Border('none')};
            add(chapter,tbl)
        end

        add(rpt, chapter);
    end


    %% PSFs
    f=find(strcmp({GEN.type},'bead_psf'));
    if ~isempty(f)
        chapter = Chapter('Title', 'Bead PSF');

        for ii=1:length(f)
            GEN(f(ii)).plotting_func(GEN(f(ii)).full_path_to_data)

            fig = Figure();

            figImg = Image(getSnapshotImage(fig, rpt));
            figImg.Style = imgStyle;
            delete(gcf);

            add(chapter, figImg);
        end

        add(rpt, chapter);
    end


    %% Lens paper
    f=find(strcmp({GEN.type},'lens_paper'));
    if ~isempty(f)
        chapter = Chapter('Title', 'Gain from lens paper');

        for ii=1:length(f)
            [~,txt] = GEN(f(ii)).plotting_func(GEN(f(ii)).full_path_to_data,4);

            fig = Figure();

            figImg = Image(getSnapshotImage(fig, rpt));
            figImg.Style = imgStyle;
            delete(gcf);

            [~,fname] = fileparts(GEN(f(ii)).full_path_to_data);

            p1 = Paragraph(txt);

            add(chapter,p1)
            add(chapter, figImg);
            add(chapter,br)
        end

        add(rpt, chapter);
    end



    close(rpt);
    rptview(rpt);

    if nargout>0
        varargout{1} = rpt;
    end

    if nargout>1
        varargout{2} = GEN;
    end

