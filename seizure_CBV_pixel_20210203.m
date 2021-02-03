% seizure OIS analysis
% written by dychung@mgh.harvard.edu with Isra Tamim (isra.tamim@charite.de)
 
clear all; close all;
 
fs = 1; % put in actual frequency (fps) here
srate = 2; % this is the factor by which to reduce the data (filtered then decimated)
 
[selected_files image_directory] = uigetfile('*.jpg',...
    'Pick files (select multiple)','MultiSelect','on'); % select files
cd(image_directory); % go to the image folder
[start_file_number,was_it_ok]=listdlg('PromptString','Pick start file',...
    'ListString',selected_files,'ListSize',[200 500]); % pick start file

% name experiment
prompt = {'Enter experiment name'};
dlg_title = 'Experiment name';
defaultans = {'expt'};
expt_name = inputdlg(prompt,dlg_title,1,defaultans);
expt_name = expt_name{1};
 
% extract time data from file names
% these file names have hr.min.sec embeded in them, time of day
frame_time = single(zeros(length(selected_files),4));
n = 1;
while n <= length(selected_files)
    % fix issue where there are 2 times
    % look for (2) in the name of the file
    if contains(selected_files{n},'(2)')==1 % (2) comes before the non-(2) file in the file prompt, but depending on the OS matlab may read in correct order in selected_files
        % if (2) is in the name, then reorder selected files
        current_file_name = selected_files{n}; % this has (2) in it
        next_file_name = selected_files{n+1}; % this is normal
        
        % for the (2) time, derive time differently
        frame_time(n,2) = str2num(current_file_name(end-14:end-13)); % this pulls out the hour
        frame_time(n,3) = str2num(current_file_name(end-11:end-10)); % this pulls out the minutes
        frame_time(n,4) = str2num(current_file_name(end-8:end-7)); % this pulls out the seconds
        frame_time(n,1) = frame_time(n,4) + (60 * frame_time(n,3)) + ...
            (60 * 60 * frame_time(n,2)); % this converts time of day into seconds
        
        frame_time(n+1,:) = frame_time(n,:); % this is the same time
        
        % rearrange
        selected_files{n} = next_file_name;
        selected_files{n+1} = current_file_name;
        n = n + 2; % because we're doing 2 files here
%         n = n + 1;
    else
        current_file_name = selected_files{n};
        frame_time(n,2) = str2num(current_file_name(end-11:end-10)); % this pulls out the hour
        frame_time(n,3) = str2num(current_file_name(end-8:end-7)); % this pulls out the minutes
        frame_time(n,4) = str2num(current_file_name(end-5:end-4)); % this pulls out the seconds
        frame_time(n,1) = frame_time(n,4) + (60 * frame_time(n,3)) + ...
            (60 * 60 * frame_time(n,2)); % this converts time of day into seconds
        n = n + 1;
    end
end

% read one image to find the size and make mask
img = imread(selected_files{1});
img_g = img(:,:,2); % just look at 2nd green channel
[Ly0,Lx0] = size(img_g);

% reduce data with a 2D convolution, average pixels
ones_matrix = single(ones(srate));
get_reduction = conv2(img_g,ones_matrix,'valid');
img_g_reduced = get_reduction(1:srate:end,1:srate:end)/(srate^2);
%img_g_reduced = img_g(1:srate:end,1:srate:end);
[Ly,Lx] = size(img_g_reduced); % size of srate reduced frame

% make mask
brain_mask_fig = figure(1);
lowlim = 0; % set initial contrast parameters empirically
uplim = 255;
while 1
    brain_mask_img = imshow(img_g_reduced,[lowlim uplim]);
    axis image;
    set(brain_mask_fig,'Position',[800 500 700 500]);
    choice = questdlg('Is contrast OK?', ...
        'Contrast', ...
        'Yes','No','Yes');
    switch choice
        case 'Yes'
            break;
        case 'No'
            prompt = {'Enter Lower Contrast Limit:','Enter Upper Contrast Limit:'};
            dlg_title = 'Contrast';
            num_lines = 1;
            defaultans = {num2str(lowlim), num2str(uplim)};
            answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
            lowlim = str2double(answer(1));
            uplim = str2double(answer(2));
    end
end
get_roi = impoly;
pos_roi = getPosition(get_roi);
roi_poly = impoly(gca,pos_roi);
ROImask_get = createMask(roi_poly,brain_mask_img);
xROI_get = pos_roi(:,1); xROI_get = [xROI_get; pos_roi(1,1)]; % add extra row to close loop when plotting
yROI_get = pos_roi(:,2); yROI_get = [yROI_get; pos_roi(1,2)]; 
brain_mask = ROImask_get;
brain_mask = single(brain_mask);  
brain_mask(brain_mask==0) = NaN; % make zeros NaN in mask
xbrainmask = single(xROI_get);
ybrainmask = single(yROI_get);

Lz = length(selected_files);
Int = single(zeros(Ly,Lx,Lz));

% read in all images and preprocess them
for n = 1:Lz
    pre_int = imread(selected_files{n}); % imports all channels
    pre_int_2 = single(pre_int(:,:,2)); % just the green channel
    get_reduction = conv2(pre_int_2,ones_matrix,'valid'); % convolution to reduce data
    pre_int_3 = get_reduction(1:srate:end,1:srate:end)/(srate^2); % complete data reduction
    %pre_int_3 = pre_int_2(1:srate:end,1:srate:end); % decimate data by srate
    pre_int_4 =  pre_int_3; % green channel intensity signal over time
    Int(:,:,n) = times(brain_mask,pre_int_4); % mask everything other than brain
end
clear pre_int pre_int_2 pre_int_3 pre_int_4

% plot mean of the intensity signal over time
get_mean_int = nanmean(Int,1); % gets the mean over the 1st dimension
mean_int = nanmean(get_mean_int,2); % gets the mean over the 2nd dimension
mean_int = squeeze(mean_int); % get rid of the redundant dimension
clear get_mean_int;

% plot the mean of all green channel frames over time 
fig = figure(2);
set(fig,'Position',[200 80 1200 200]);
tracing = plot(mean_int); % plot the mean of the intensity signal over time
hold on;
title('Mean raw intensity over time'); xlabel('Frame #');
fig = gca;

% get the dHbt
mean_int_initial = mean(Int,3); % inital mean intensity for all pixels over time
dHbt_initial = log10(mean_int_initial./Int);

k=1;
fig_movie = figure(3);
proportion_range = [-.08 .08]; % initial proportion
change_step = 0.01; % this is the step change for contrast adjustment
% dialog instructions for navigation
instructions= msgbox({'Scroll with arrow keys.', 'Use r to re-define reference',...
    'Use k and l to move fast','Use f and d to adjust upper contrast.',...
    'Use v and c to adjust lower contrast.','Hit Enter when finished.'});
set(instructions,'Position',[130 450 210 120]);
while k < size(dHbt_initial, 3)
    %current_div_tc = imgaussfilt(current_div_tc_presmooth,1,'FilterSize',3); % make the image smooth
    %current_div_tc = imgaussfilt(current_div_tc_presmooth,1); % make the image pretty
    %current_div_tc = current_div_tc_presmooth; % don't mess with image
    imagesc(dHbt_initial(:,:,k),proportion_range);
    colormap jet;
    colorbar;
    title([num2str(k) ' / ' num2str(size(dHbt_initial, 3)) ]);
    was_a_key = waitforbuttonpress;
    if was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'leftarrow')
      k = max(1, k - 1);
    elseif was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'uparrow')
      k = max(1, k - 1);
    elseif was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'l')
      if k < (size(dHbt_initial, 3) - 100)
          k = k + 100;
      else
      end
    elseif was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'k')
        k = max(1, k - 100);
    elseif was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'return') % stops the loop and moves to ROI selection
        ROI_select_frame_number = k; % this is the frame where ROIs are selected
        k = size(dHbt_initial, 3); % stop loop
    elseif was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'r') % makes new working reference
        working_reference_one_over_tc = rdivide(1,tc(:,:,k));
        reference_frame_number = k;
    elseif was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'f') % increases upper contrast value
        proportion_range = [proportion_range(1) proportion_range(2)+change_step];
    elseif was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'd') % decreases upper contrast value
        if proportion_range(2) > (proportion_range(1) + change_step)
            proportion_range = [proportion_range(1) proportion_range(2)-change_step];
        else
        end
    elseif was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'v') % increase lower contrast value
        if (proportion_range(1) + change_step) < proportion_range(2)
            proportion_range = [proportion_range(1)+change_step proportion_range(2)];
        else
        end
    elseif was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'c') % decrease lower contrast value
        proportion_range(1) = proportion_range(1)-change_step;
    else
      k = k + 1;
    end
end
close(instructions); clear instructions;

% go back to figure 2
% ask if there are any timespans to be included
okay = 0;
mean_spans_all = [];
while okay == 0
    choice = questdlg('Are there any timespans to pick?', ...
        'Timespans to include?', ...
        'Yes','No','No');
    switch choice
        case 'Yes'
            instructions = helpdlg('Select timespans.');
            set(instructions, 'Position', [90 550 200 80]);
            fig = figure(2); hold on;
            mean_int_fig = figure(2);
            fig = gca;
            h = imrect;
            pos_h = getPosition(h);
            delete(h);
            x_width = pos_h;
            x_width(:,2)=[]; x_width(:,3)=[];
            integer_x_width = round(x_width); % format is x-coord and width
            if (integer_x_width(1) + integer_x_width(2)) > size(mean_int,1)
                integer_x_width(2) = size(mean_int,1) - integer_x_width(1);
            else
            end
            integer_x_width = [integer_x_width(1), integer_x_width(1)+integer_x_width(2)-1]; % changing format to actual first frame and actual final frame in span
            mean_spans_all = [mean_spans_all; integer_x_width]; % append rows to the end of the matrix, this as all the spans
            rect=rectangle('Position',[integer_x_width(1), fig.YLim(1), ...
                integer_x_width(2)-integer_x_width(1), fig.YLim(2)-fig.YLim(1)], ...
                'LineStyle','none','FaceColor',[0.8 0.8 0.8]); uistack(rect,'bottom');
        case 'No'
            instructions_3 = helpdlg('Select baseline.');
            set(instructions_3, 'Position', [90 550 200 80]);
            fig = figure(2); hold on;
            mean_int_fig = figure(2);
            fig = gca;
            h = imrect;
            pos_h = getPosition(h);
            delete(h);
            x_width = pos_h;
            x_width(:,2)=[]; x_width(:,3)=[];
            integer_x_width = round(x_width); % format is x-coord and width
            if (integer_x_width(1) + integer_x_width(2)) > size(mean_int,1)
                integer_x_width(2) = size(mean_int,1) - integer_x_width(1);
            else
            end
            integer_x_width = [integer_x_width(1), integer_x_width(1)+integer_x_width(2)-1]; % changing format to actual first frame and actual final frame in span
            baseline_span = integer_x_width; % selected baseline span
            rect=rectangle('Position',[integer_x_width(1), fig.YLim(1), ...
                integer_x_width(2)-integer_x_width(1), fig.YLim(2)-fig.YLim(1)], ...
                'LineStyle',':','FaceColor','r');
            uistack(tracing,'top');
            choice_2 = questdlg('Does everything look okay?', ...
                'Look okay?', 'Yes','No','Yes');
            switch choice_2
                case 'Yes'
                    okay = 1; % this means this step is finished and we can move on
                case 'No' % then the figure has to be reset and start over
                    mean_spans_all = [];
                    delete(fig);
                    fig = figure(2);
                    set(fig,'Position',[200 80 1200 200]);
                    plot(mean_int); % plot the mean of the intensity signal over time
                    hold on;
                    title('Mean raw intensity over time'); xlabel('Frame #');
                    fig = gca;
            end
    end
end
close(instructions); clear instructions;
clear instructions_3;
mean_spans_all = single(mean_spans_all);

% get new dHbt segments which we'll call xHbt because not actually Hbt
mean_int_baseline = nanmean(Int(:,:,baseline_span(1):baseline_span(2)),3); % mean intensity over selected baseline
for n = 1:size(mean_spans_all,1)
    mean_int_all(:,:,n) = nanmean(Int(:,:,mean_spans_all(n,1):mean_spans_all(n,2)),3); %  mean intensity for all pixels over timespan n
    xHbt{n} = log10(mean_int_all(:,:,n)./Int(:,:,mean_spans_all(n,1):mean_spans_all(n,2))); % get xHbt and put segments into a cell, mean from each segment
    xHbt_1_baseline{n} = log10(mean_int_baseline./Int(:,:,mean_spans_all(n,1):mean_spans_all(n,2))); % use a selected baseline
end

% stitch the segments together
xHbt_comb = zeros(Ly,Lx,mean_spans_all(end,end)); % premake variable to fill in below
xHbt_comb(xHbt_comb==0) = NaN;
xHbt_comb = single(xHbt_comb);
xHbt_comb_1_baseline = zeros(Ly,Lx,mean_spans_all(end,end)); % premake variable to fill in below
xHbt_comb_1_baseline(xHbt_comb_1_baseline==0) = NaN;
xHbt_comb_1_baseline = single(xHbt_comb_1_baseline);

for n = 1:size(mean_spans_all,1)
    current_xHbt = xHbt{n};
    current_span = mean_spans_all(n,:);
    xHbt_comb(:,:,current_span(1,1):current_span(1,2)) = current_xHbt;
    
    current_xHbt_1_baseline = xHbt_1_baseline{n};
    xHbt_comb_1_baseline(:,:,current_span(1,1):current_span(1,2)) = current_xHbt_1_baseline;
end
clear current_xHbt current_xHbt_1_baseline current_span;

% re-display images and adjust contrast
k=1;
fig_movie = figure(3);
proportion_range = [-.08 .08]; % initial proportion
change_step = 0.01; % this is the step change for contrast adjustment
% dialog instructions for navigation
instructions= msgbox({'Scroll with arrow keys.', 'Use r to re-define reference',...
    'Use k and l to move fast','Use f and d to adjust upper contrast.',...
    'Use v and c to adjust lower contrast.','Hit Enter when finished.'});
set(instructions,'Position',[130 450 210 120]);
while k < size(xHbt_comb, 3)
    %current_div_tc = imgaussfilt(current_div_tc_presmooth,1,'FilterSize',3); % make the image pretty
    %current_div_tc = imgaussfilt(current_div_tc_presmooth,1); % make the image pretty
    %current_div_tc = current_div_tc_presmooth; % don't mess with image
    %imagesc(xHbt_comb(:,:,k),proportion_range); % new mean for each segment
    imagesc(xHbt_comb_1_baseline(:,:,k),proportion_range); % use single baseline
    colormap jet;
    colorbar;
    title([num2str(k) ' / ' num2str(size(xHbt_comb, 3)) ]);
    was_a_key = waitforbuttonpress;
    if was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'leftarrow')
      k = max(1, k - 1);
    elseif was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'uparrow')
      k = max(1, k - 1);
    elseif was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'l')
      if k < (size(xHbt_comb, 3) - 100)
          k = k + 100;
      else
      end
    elseif was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'k')
        k = max(1, k - 100);
    elseif was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'return') % stops the loop and moves to ROI selection
        ROI_select_frame_number = k; % this is the frame where ROIs are selected
        break % stop loop
    elseif was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'r') % makes new working reference
        working_reference_one_over_tc = rdivide(1,tc(:,:,k));
        reference_frame_number = k;
    elseif was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'f') % increases upper contrast value
        proportion_range = [proportion_range(1) proportion_range(2)+change_step];
    elseif was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'd') % decreases upper contrast value
        if proportion_range(2) > (proportion_range(1) + change_step)
            proportion_range = [proportion_range(1) proportion_range(2)-change_step];
        else
        end
    elseif was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'v') % increase lower contrast value
        if (proportion_range(1) + change_step) < proportion_range(2)
            proportion_range = [proportion_range(1)+change_step proportion_range(2)];
        else
        end
    elseif was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'c') % decrease lower contrast value
        proportion_range(1) = proportion_range(1)-change_step;
    else
      k = k + 1;
    end
end
close(instructions); clear instructions;


% pick ROIs
% ask if you want to use an ROI from another analysis session
choice = questdlg('Do you want to use a prior ROI?','Use prior ROI?',...
        'Yes','No','Yes');
switch choice
    case 'Yes'
        % go over this, no checked yet
        [prior_ROI_file file_path] = uigetfile('*.mat','Select file containing prior ROI');
        load(fullfile(file_path,prior_ROI_file),'ROImask','xROI','yROI',...
            'current_div_tc','proportion_range','num_ROI','working_reference_one_over_tc');
        fig_movie = figure(3);
        imagesc(current_div_tc,proportion_range); colormap gray; colorbar;
        hold on;
        for n = 1:num_ROI
            plot(xROI{n}, yROI{n}, 'r', 'LineWidth',2); hold on; % draw ROI on figure
        end
        hold off;
    case 'No'
        % view tc proportion images
        % k=1; % CHANGE NUMBER TO START AT DIFFERENT FRAME - here, use k
        % from above
        
        %proportion_range = [0 2.5]; % initial proportion
        %change_step = 0.05; % this is the step change for contrast adjustment
        % dialog instructions for navigation
        instructions1= msgbox({'Scroll with arrow keys.', 'Use r to re-define reference',...
            'Use k and l to move fast','Use f and d to adjust upper contrast.',...
            'Use v and c to adjust lower contrast.','Hit Enter when finished.'});
        set(instructions1,'Position',[130 450 210 120]);
        instructions2= msgbox({'Pick 6 ROIs'});
        set(instructions2,'Position',[130 450 40 50]);
        while 1
            fig_movie = figure(3); % define ROIs
            roi_image = imshow(xHbt_comb(:,:,k),proportion_range);
            set(fig_movie,'Position',[800 500 700 500]);
            axis off;
            
            % ask for ROIs
            num_ROI = 6; % # of ROI
            for n = 1:num_ROI
                while k < Lz
                    fig_movie = figure(3);
                    %imagesc(xHbt_comb(:,:,k),proportion_range); % new mean for each segment
                    imagesc(xHbt_comb_1_baseline(:,:,k),proportion_range); % use single baseline
                    %colormap jet;
                    axis image off;
                    title(['Pick ROI ' num2str(k) ' / ' num2str(size(xHbt_comb, 3)) ]);
                    was_a_key = waitforbuttonpress;
                    if was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'leftarrow')
                      k = max(1, k - 1);
                    elseif was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'uparrow')
                      k = max(1, k - 1);
                    elseif was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'l')
                      if k < (size(xHbt_comb, 3) - 100)
                          k = k + 100;
                      else
                      end
                    elseif was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'k')
                        k = max(1, k - 100);
                    elseif was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'return') % stops the loop and moves to ROI selection
                        ROI_select_frame_number = k; % this is the frame where ROIs are selected
                        break; % stop loop
                    elseif was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'r') % makes new working reference
                        working_reference_one_over_tc = rdivide(1,tc(:,:,k));
                        reference_frame_number = k;
                        reference_frame_bp = reference_frame_number;
                    elseif was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'f') % increases upper contrast value
                        proportion_range = [proportion_range(1) proportion_range(2)+change_step];
                    elseif was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'd') % decreases upper contrast value
                        if proportion_range(2) > (proportion_range(1) + change_step)
                            proportion_range = [proportion_range(1) proportion_range(2)-change_step];
                        else
                        end
                    elseif was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'v') % increase lower contrast value
                        if (proportion_range(1) + change_step) < proportion_range(2)
                            proportion_range = [proportion_range(1)+change_step proportion_range(2)];
                        else
                        end
                    elseif was_a_key && strcmp(get(fig_movie, 'CurrentKey'), 'c') % decrease lower contrast value
                        if (proportion_range(1) - change_step) > 0
                            proportion_range = [proportion_range(1)-change_step proportion_range(2)];
                        else
                        end
                    else
                      k = k + 1;
                    end
                end
                fig_movie = figure(3);
                %roi_image = imshow(xHbt_comb(:,:,k),proportion_range); % new mean for each segment
                roi_image = imshow(xHbt_comb_1_baseline(:,:,k),proportion_range); % single baseline
                
                axis image;
                set(fig_movie,'Position',[800 500 700 500]);
                get_roi = impoly;
                pos_roi = getPosition(get_roi);
                roi_poly = impoly(gca,pos_roi);
                ROImask_get = createMask(roi_poly,roi_image);
                xROI_get = pos_roi(:,1); xROI_get = [xROI_get; pos_roi(1,1)]; % add extra row to close loop when plotting
                yROI_get = pos_roi(:,2); yROI_get = [yROI_get; pos_roi(1,2)]; 
                %[ROImask_get xROI_get yROI_get] = roipoly(xHbt_comb(:,:,k)); hold on; % use ROImask to get stuff
                hold on;
                plot(xROI_get, yROI_get, 'w','LineWidth',4); % draw ROI on figure
                plot(xROI_get, yROI_get, 'k','LineWidth',1);
                hold off; 
                ROImask{n} = ROImask_get;
                ROImask{n} = single(ROImask{n});  
                ROImask{n}(ROImask{n}==0) = NaN; % make zeros NaN in mask
                xROI{n} = xROI_get;
                yROI{n} = yROI_get;
                % ask which ROIS this is
                prompt = {'Which ROI is this?'};
                dlg_title = 'Which ROI?';
                defaultans = {''};
                name_ROI{n} = inputdlg(prompt,dlg_title,1,defaultans);
                pause(0.5);
                if n==1
%                     close(instructions2); clear instructions2;
                else
                end
            end
            
            % ask for lines
            num_line = 2; % number of lines to draw
            for n = 1:num_line
                fig_movie = figure(3);
                %roi_image = imshow(xHbt_comb(:,:,k),proportion_range); % new mean for each segment
                roi_image = imshow(xHbt_comb_1_baseline(:,:,k),proportion_range); % single baseline
                
                axis image;
                set(fig_movie,'Position',[800 500 700 500]);
                title('Pick Lines'); 
                get_roi = imline;
                pos_roi = getPosition(get_roi);
                roi_line = imline(gca,pos_roi);
                ROImask_get = createMask(roi_line,roi_image);
                xROI_get = pos_roi(:,1);
                yROI_get = pos_roi(:,2);
                hold on;
                plot(xROI_get, yROI_get, 'w','LineWidth',4); % draw ROI on figure
                plot(xROI_get, yROI_get, 'k','LineWidth',1);
                hold off;
                
                xLine{n} = xROI_get;
                yLine{n} = yROI_get;
                lengthLine(n) = sqrt((xLine{n}(2)-xLine{n}(1))^2+(yLine{n}(2)-yLine{n}(1))^2);
                % ask which ROIS this is
                prompt = {'Which Line is this?'};
                dlg_title = 'Which Line?';
                defaultans = {''};
                name_Line{n} = inputdlg(prompt,dlg_title,1,defaultans);
            end
            
            % plot ROIs and lines
            fig_roi = figure(4);
            %imagesc(xHbt_comb(:,:,k),proportion_range); % new mean for each segment
            imagesc(xHbt_comb_1_baseline(:,:,k),proportion_range); % use single baseline
            colormap jet; colorbar;
            hold on;
            for n = 1:num_ROI
                plot(xROI{n}, yROI{n}, 'w', 'LineWidth',4); hold on; % draw ROI on figure
                plot(xROI{n}, yROI{n}, 'k', 'LineWidth',1); hold on; % draw ROI on figure
            end
            for n = 1:num_line
                plot(xLine{n}, yLine{n}, 'w', 'LineWidth',4); hold on; % draw ROI on figure
                plot(xLine{n}, yLine{n}, 'k', 'LineWidth',1); hold on; % draw ROI on figure
            end
            hold off;
            ROI_check = questdlg('Are ROIs correct?','ROIs correct?',...
                'Yes','No','Yes');
            if strcmp(ROI_check,'Yes') == 1
                break;
            elseif strcmp(ROI_check,'No') == 1
                % move on
            else
            end
        end
        
end
close(instructions1); clear instructions1;

% get mean value in ROI
masked_xHbt_comb(:,:,num_ROI) = single(zeros(Ly,Lx)); % set aside memory
for n = 1:num_ROI
    ROImask{n} = single(ROImask{n}); % convert to single
    ROImask{n}(ROImask{n}==0) = NaN; % make zeros NaN in mask
    for k = 1:size(xHbt_comb,3)
        masked_xHbt_comb(:,:,n) = times(ROImask{n},xHbt_comb(:,:,k)); % mask, indiv means
        masked_xHbt_comb(:,:,n) = times(ROImask{n},xHbt_comb_1_baseline(:,:,k)); % mask, single baseline
        masked_xHbt_comb_avg(k,n) = nanmean(nanmean(masked_xHbt_comb(:,:,n)));
        
    end
end

% plot ROIs timecourses
fig = figure(5); set(fig,'Position',[500 150 1000 800]);
roi_fig_timecourses = figure(5);
for n = 1:num_ROI
%     subplot(num_ROI, 2, 2*n - 1);
%    subplot(num_ROI+1, 1, n);
    subplot(num_ROI, 1, n);
    plot(masked_xHbt_comb_avg(:,n));
    %ylim([0 1.5]);
    if n==1
        title('CBV changes in ROIs');
    elseif n==num_ROI
        xlabel('frame #');
    end
    ylabel(name_ROI{n});
    hold on;
    y_marker = fig.CurrentAxes.YLim(1);
end
% subplot(num_ROI+1, 1, num_ROI+1);
% plot(linspace(1,length(current_run),length(current_run(:,2))),...
% %     current_run(:,2)); % plotting MAP
% % ylim([50 120]);
% % ylabel('c');

% get values across lines
for n = 1:num_line
    X_length(n) = xLine{n}(2)-xLine{n}(1);
    Y_length(n) = yLine{n}(2)-yLine{n}(1);
    m_line(n) = (yLine{n}(2)-yLine{n}(1)) ./ (xLine{n}(2)-xLine{n}(1));
    m_line_y(n) = (xLine{n}(2)-xLine{n}(1)) ./ (yLine{n}(2)-yLine{n}(1));
    y_intersect(n) = yLine{n}(1) - m_line(n)*xLine{n}(1);
    x_intersect(n) = xLine{n}(1) - m_line_y(n)*yLine{n}(1);
    for X_step = 1:round(abs(X_length(n))) % x coordinate
        if X_length(n) < 0 % if this is negative
            X = (xLine{n}(1) - X_step); % change sign
        elseif X_length(n) > 0 % if positive
            X = (xLine{n}(1) + X_step);
        else
            % there is a problem or line is vertical
        end
        Y = m_line(n) * X + y_intersect(n); % y coordinate (likely non-integer)
        xHbt_line_x(X_step,:) = interp3(xHbt_comb,X,Y,1:size(xHbt_comb,3)); % have to interpolate given non-integers
        xHbt_line_x_all{n} = xHbt_line_x;
        
        % these don't calculate actual pixel distances
    end
    for Y_step = 1:round(abs(Y_length(n))) % y coordinate
        if X_length(n) < 0 % if this is negative
            Y = (yLine{n}(1) - Y_step); % change sign
        elseif X_length(n) > 0 % if positive
            Y = (yLine{n}(1) + Y_step);
        else
            % there is a problem or line is vertical
        end
        X = m_line_y(n) * Y + x_intersect(n); % x coordinate (likely non-integer)
        xHbt_line_y(Y_step,:) = interp3(xHbt_comb,X,Y,1:size(xHbt_comb,3)); % have to interpolate given non-integers
        xHbt_line_y_all{n} = xHbt_line_y;
        
        % these don't calculate actual pixel distances
    end
   
    for m = 1:size(xHbt_comb,3) % each frame
        xHbt_line_profile(:,m) = improfile(xHbt_comb(:,:,m),xLine{n},...
            yLine{n},lengthLine(n),'bilinear'); % indiv means
        xHbt_line_profile(:,m) = improfile(xHbt_comb_1_baseline(:,:,m),xLine{n},...
            yLine{n},lengthLine(n),'bilinear'); % single baseline
        
    end
    xHbt_line_profile_all{n} = single(xHbt_line_profile);
    clear xHbt_line_profile
end

% plot line timecourses
fig = figure(6); set(fig,'Position',[400 150 1000 800]);
line_fig_timecourses = figure(6);
upper_lower_lim = 0.04;
for n = 1:num_line
    subplot(num_line, 1, n);

    surf(xHbt_line_x_all{n},'EdgeColor','none','FaceColor','interp','FaceAlpha',1);
    colormap jet; colorbar;
    set(gca,'Ydir','reverse'); axis tight;    
    view(2);
    caxis([-upper_lower_lim upper_lower_lim]);
    if n==1
        title('CBV changes across line (X-step)');
    elseif n==num_ROI
        xlabel('frame #');
    end
    ylabel(name_Line{n});
    hold on;
    y_marker = fig.CurrentAxes.YLim(1);
end

% plot line timecourses
fig = figure(7); set(fig,'Position',[400 150 1000 800]);
line_fig_timecourses_y = figure(7);
upper_lower_lim = 0.04;
for n = 1:num_line
    subplot(num_line, 1, n);

    surf(xHbt_line_y_all{n},'EdgeColor','none','FaceColor','interp','FaceAlpha',1);
    colormap jet; colorbar;
    set(gca,'Ydir','reverse'); axis tight;    
    view(2);
    caxis([-upper_lower_lim upper_lower_lim]);
    if n==1
        title('CBV changes across line (Y-step)');
    elseif n==num_ROI
        xlabel('frame #');
    end
    ylabel(name_Line{n});
    hold on;
    y_marker = fig.CurrentAxes.YLim(1);
end

fig = figure(6); set(fig,'Position',[400 150 1000 800]);
line_fig_timecourses_profile = figure(6);
upper_lower_lim = 0.08;
for n = 1:num_line
    subplot(num_line, 1, n);
    surf(xHbt_line_profile_all{n},'EdgeColor','none','FaceColor','interp','FaceAlpha',1);
    colormap jet; colorbar;
    set(gca,'Ydir','reverse'); axis tight;    
    view(2);
    caxis([-upper_lower_lim upper_lower_lim]);
    if n==1
        title('CBV changes across line');
    elseif n==num_ROI
        xlabel('frame #');
    end
    ylabel([name_Line{n} ' (pixels)']);
    hold on;
    y_marker = fig.CurrentAxes.YLim(1);
end

%
% spectral power of ROIs

xlim_range = [0.0001 1];
ylim_range = [0 0.005];

spectra_roi_fig = figure(7);
L = length(masked_xHbt_comb_avg);
masked_xHbt_comb_avg_zeros = masked_xHbt_comb_avg; % convert NaNs to zeros so that fft works
masked_xHbt_comb_avg_zeros(isnan(masked_xHbt_comb_avg_zeros)==1) = 0;
for n = 1:num_ROI
    get_spect = fft(masked_xHbt_comb_avg_zeros(:,n));
    P2 = abs(get_spect/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:(L/2))/L;
    subplot(num_ROI, 1, n);
    semilogx(f,P1);
    set(gca,'FontSize',6);
    xlabel('f (Hz)','FontSize',10);
    ylabel(['|P1(f)|' name_ROI{n}],'FontSize',8);
    xlim(xlim_range);
    ylim(ylim_range);
    if n==1
        title('ROI frequency power spectrum','FontSize',10);
    else
    end
end


%%
% map power
band_freq = [0.01 0.1];
xHbt_comb_1_baseline_reshaped = reshape(xHbt_comb_1_baseline,...
    [size(xHbt_comb_1_baseline,1)*size(xHbt_comb_1_baseline,2) size(xHbt_comb_1_baseline,3)]); % make 1 xHbt 1 column per frame
xHbt_comb_1_baseline_reshaped = transpose(xHbt_comb_1_baseline_reshaped); % organize pixels into a row
xHbt_comb_1_baseline_reshaped(xHbt_comb_1_baseline_reshaped==Inf) = NaN; % turn Inf into NaNs
xHbt_comb_1_baseline_reshaped(isnan(xHbt_comb_1_baseline_reshaped)) = 0; % convert NaNs to zero
bandpower_xHbt_comb_1_baseline = bandpower(xHbt_comb_1_baseline_reshaped,...
    fs,band_freq);
bandpower_xHbt_comb_1_baseline = transpose(bandpower_xHbt_comb_1_baseline); % transpose back
bandpower_xHbt_comb_1_baseline = reshape(bandpower_xHbt_comb_1_baseline,...
    [size(xHbt_comb_1_baseline,1) size(xHbt_comb_1_baseline,2)]); % reshape from 1D to 2D
bandpower_xHbt_comb_1_baseline(bandpower_xHbt_comb_1_baseline==0) = NaN; % make zeros NaN, mostly in masked area

bandpower_fig=figure(8); set(bandpower_fig,'Position',[560 300 600 600]);
upper_range = 10e-5;
while 1
    bandpower_range = [0 upper_range];
    imagesc(bandpower_xHbt_comb_1_baseline,bandpower_range);
    colormap jet; axis image;
    choice = questdlg('Is contrast OK?', ...
        'Contrast', ...
        'Yes','No','Yes');
    switch choice
        case 'Yes'
            break
        case 'No'
            prompt = {'Enter Contrast Limit:'};
            dlg_title = 'Contrast';
            num_lines = 1;
            defaultans = {num2str(upper_range)};
            answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
            upper_range = str2double(answer);
            imagesc(bandpower_xHbt_comb_1_baseline,bandpower_range); % display new limit
            colormap jet; axis image;
    end
end
title(['Power (' num2str(band_freq(1)) '-' num2str(band_freq(2)) ' Hz)']);


%%
% save
saveas(brain_mask_fig,[expt_name '_1_brain_mask.jpg']);
saveas(mean_int_fig,[expt_name '_2_mean_intensities.jpg']);
saveas(fig_roi,[expt_name '_3_roi_fig.jpg']);
saveas(roi_fig_timecourses,[expt_name '_4_roi_timecourse.jpg']);
saveas(line_fig_timecourses_profile,[expt_name '_5_line_timecourse.jpg']);
saveas(spectra_roi_fig,[expt_name '_6_roi_spectra.jpg']);

save([expt_name '_analysis.mat'],'selected_files','image_directory',...
    'fs','srate','expt_name','Ly','Lx','Lz','frame_time','lowlim',...
    'xHbt_comb_1_baseline','mean_int_all','mean_int_baseline',...
    'mean_int','mean_spans_all','baseline_span','num_ROI','name_ROI',...
    'ROImask','xROI','yROI','masked_xHbt_comb_avg','num_line',...
    'name_Line','xLine','yLine','lengthLine','xHbt_line_profile_all');


%%
% save
saveas(brain_mask_fig,[expt_name '_1_brain_mask.jpg']);
saveas(mean_int_fig,[expt_name '_2_mean_intensities.jpg']);
saveas(fig_roi,[expt_name '_3_roi_fig.jpg']);
saveas(roi_fig_timecourses,[expt_name '_4_roi_timecourse.jpg']);
saveas(line_fig_timecourses_profile,[expt_name '_5_line_timecourse.jpg']);
saveas(spectra_roi_fig,[expt_name '_6_roi_spectra.jpg']);
saveas(bandpower_fig,[expt_name '_7_bandpower.jpg']);
% 
save([expt_name '_analysis.mat'],'selected_files','image_directory',...
    'fs','srate','expt_name','Ly','Lx','Lz','frame_time','lowlim',...
    'uplim','upper_lower_lim','proportion_range','brain_mask',...
    'xbrainmask','ybrainmask','Int','xHbt_comb',...
    'xHbt_comb_1_baseline','mean_int_all','mean_int_baseline',...
    'mean_int','mean_spans_all','baseline_span','num_ROI','name_ROI',...
    'ROImask','xROI','yROI','masked_xHbt_comb_avg','num_line',...
    'name_Line','xLine','yLine','lengthLine','xHbt_line_profile_all',...
    'xlim_range','ylim_range','bandpower_xHbt_comb_1_baseline',...
    'band_freq','bandpower_range');

%% make a movie (if not already made)
close(figure(3));
figure(3);
v = VideoWriter([expt_name '_movie.avi']);
open(v);
for n = 1:size(xHbt_comb_1_baseline,3)
    imagesc(xHbt_comb_1_baseline(:,:,n))
    colormap(jet); colorbar;
    axis image off;
%     view(90,90);
    rectangle('Position',[8 10 68 15],'FaceColor','k');
    figure(3);
%     axis on;
    upper_lower_lim = 0.08;
    caxis([-upper_lower_lim upper_lower_lim]);
    text(10,15,[num2str(frame_time(n,2)) ':' num2str(frame_time(n,3)) ':' num2str(frame_time(n,4))],...
        'HorizontalAlignment','left','FontSize',14,'Color',[1 1 1]);
    frame = getframe;
    writeVideo(v,frame);
end
close(v);







