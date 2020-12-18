% ======================= GENERAL INFO ====================================
%
% This is the Matlab code used for analyzing time on task effects on
% task-related and task-unrelated activity in time-frequency space of
% the EEG.
%
% Code by Stefan Arnau, November 2020
% Email: arnau@ifado.de
% GitHub: https://github.com/fischmechanik/theta_time_on_task
% OSF repository: https://osf.io/chmkj/
% OSF identifier: DOI 10.17605/OSF.IO/CHMKJ
%

% Remove residuals
clear all;

% ========================= PATH VARIABLES ===============================================================================================

% Path vars (there are some relative paths used in the script):
PATH_EEGLAB           = 'a_path';
PATH_FIELDTRIP        = 'a_path';
PATH_CUSTOM           = 'a_path';
PATH_RAW_DATA         = 'a_path';
PATH_ICSET            = 'a_path';
PATH_AUTOCLEANED      = 'a_path';
PATH_WAVELETS         = 'a_path';
PATH_TFDECOMP         = 'a_path';
PATH_PLOT             = 'a_path'; 
PATH_VEUSZ            = 'a_path';
PATH_REGRESSION       = 'a_path';
PATH_CLUSTSTATS       = 'a_path';
PATH_COMPONENT_DATA   = 'a_path'; 

% The subject list
subject_list =    {'Vp01', 'Vp02', 'Vp03',...
                   'Vp04', 'Vp05', 'Vp06',...
                   'Vp08', 'Vp09', 'Vp10',...
                   'Vp11', 'Vp12', 'Vp13',...
                   'Vp14', 'Vp40', 'Vp41',...
                   'Vp43', 'Vp44', 'Vp45'};


% SWITCH: Switch parts of script on/off
to_execute = {'part1'};

% ========================= PART 1: Event coding and preprocessing  =========================================================================================================

if ismember('part1', to_execute)

    % Init EEGlab
    addpath(PATH_EEGLAB);
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

    % Find chanlocfile
    channel_location_file = [PATH_EEGLAB 'plugins/dipfit3.3/standard_BESA/standard-10-5-cap385.elp'];

    % Iterating subject list
    for s = 1 : length(subject_list)

        % Current iteration subject
        subject = subject_list{s};
        id = str2num(subject(3 : 4));

        % Load data
        EEG = pop_loadset('filename', [subject '.set'], 'filepath', PATH_RAW_DATA, 'loadmode', 'all');
        EEG = pop_chanedit(EEG, 'lookup', channel_location_file);
        EEG.chanlocs_original = EEG.chanlocs;

        % Init event struct
        new_events = struct('latency', {},...
                            'oldtype', {},...
                            'type', {},...
                            'code', {},...
                            'id', {},...
                            'stimloc', {},... 
                            'corresp', {},...   
                            'saliency', {},...
                            'bl', {},...  
                            'sbl', {},...      
                            'sbl_total', {},...    
                            'stimnum', {},...     
                            'stimnum_sbl', {},...
                            'resploc', {},...    
                            'acc', {},...     
                            'rt', {},...         
                            'urevent', {},...
                            'duration', {}...
                            );

        % Sort events by latency
        [~, idx] = sort(cell2mat({EEG.event.latency}));
        EEG.event = EEG.event(idx);

        % Code experimental conditions
        stimnum = 0;
        stimnum_sbl = 0;
        sblnum = 0;
        for e = 1 : length(EEG.event)

            % If stim event
            if strcmpi(EEG.event(e).type(1), 'S') & ismember(str2num(EEG.event(e).type(2 : end)), [1 : 108])

                % Get eventcode and increase eventcount 
                enum = str2num(EEG.event(e).type(2 : end));
                stimnum = stimnum + 1;

                % Track change of sbl and increase sbl stimcount
                if ceil(enum / 12) > sblnum
                    sblnum = sblnum + 1;
                    stimnum_sbl = 0;
                end
                stimnum_sbl = stimnum_sbl + 1;

                % Code all the things
                new_events(stimnum).latency = EEG.event(e).latency;
                new_events(stimnum).oldtype = EEG.event(e).type;
                new_events(stimnum).type = 'stim';
                new_events(stimnum).code = 'stim';
                new_events(stimnum).id = id;
                if ismember(mod(enum, 4), [1, 0])
                    new_events(stimnum).stimloc = 'left';
                else
                    new_events(stimnum).stimloc = 'right';
                end
                if ismember(mod(enum, 4), [1, 2])
                    new_events(stimnum).corresp = 1; % corresponding
                else
                    new_events(stimnum).corresp = 2; % non-corresponding
                end
                if ismember(mod(enum, 12), [1, 2, 3, 4])
                    new_events(stimnum).saliency = 'low';
                elseif ismember(mod(enum, 12), [5, 6, 7, 8])
                    new_events(stimnum).saliency = 'mid';
                else
                    new_events(stimnum).saliency = 'high';
                end
                new_events(stimnum).bl = ceil(enum / 36);
                new_events(stimnum).sbl = ceil(ceil(enum / 12) / 3);
                new_events(stimnum).sbl_total = ceil(enum / 12);
                new_events(stimnum).stimnum = stimnum;
                new_events(stimnum).stimnum_sbl = stimnum_sbl;

                % Loop for response
                f = e + 1;
                stimcodes = {};
                for n = 1 : 108
                    if n < 10
                        stimcodes{end + 1} = ['S  ' num2str(n)];
                    elseif n < 100
                        stimcodes{end + 1} = ['S ' num2str(n)];
                    else
                        stimcodes{end + 1} = ['S' num2str(n)];
                    end
                end
                while ~strcmpi(EEG.event(f).type(1), 'R') &...
                      ~ismember(EEG.event(f).type, stimcodes) &...
                      f < length(EEG.event)
                    f = f + 1;
                end
                rnum = 0;
                if strcmpi(EEG.event(f).type(1), 'R')
                    rnum = str2num(EEG.event(f).type(2 : end));
                end

                if ismember(rnum, [3, 13])
                    new_events(stimnum).resploc = 'left';
                    new_events(stimnum).acc = 1;
                    new_events(stimnum).rt = (EEG.event(f).latency - EEG.event(e).latency) * 1000 / EEG.srate;
                elseif ismember(rnum, [4, 14])
                    new_events(stimnum).resploc = 'right';
                    new_events(stimnum).acc = 1;
                    new_events(stimnum).rt = (EEG.event(f).latency - EEG.event(e).latency) * 1000 / EEG.srate;
                elseif ismember(rnum, [5, 6, 9, 10, 11, 12, 15, 16])
                    new_events(stimnum).resploc = 'none';
                    new_events(stimnum).acc = 0;
                    new_events(stimnum).rt = NaN;
                elseif ismember(rnum, [7, 8])
                    new_events(stimnum).resploc = 'none';
                    new_events(stimnum).acc = 2;
                    new_events(stimnum).rt = NaN;
                else
                    new_events(stimnum).resploc = 'none';
                    new_events(stimnum).acc = NaN;
                    new_events(stimnum).rt = NaN;
                end

                new_events(stimnum).urevent = e;
                new_events(stimnum).duration = 1;

            end % End event check
        end % End eventit

        % Collect boundaries
        for e = 1 : length(EEG.event)
            if strcmpi(EEG.event(e).type, 'boundary')
                new_events(end + 1).latency = EEG.event(e).latency;
                new_events(end).oldtype = 'boundary';
                new_events(end).type = 'boundary';
                new_events(end).code = 'boundary';
                new_events(end).id = 0;
                new_events(end).stimloc = 'none';
                new_events(end).corresp = 0;
                new_events(end).saliency = 'none';
                new_events(end).bl = NaN;
                new_events(end).sbl = NaN;
                new_events(end).sbl_total = NaN;
                new_events(end).stimnum = NaN;
                new_events(end).stimnum_sbl = NaN;
                new_events(end).resploc = 'none';
                new_events(end).acc = NaN;
                new_events(end).rt = NaN;
                new_events(end).urevent = length(new_events);
                new_events(end).duration = 1;
            end
        end

        % Sort new events by latency
        [~, idx] = sort(cell2mat({new_events.latency}));
        new_events = new_events(idx);

        % Replace events by new events
        EEG.event = new_events;
        EEG = eeg_checkset(EEG, 'eventconsistency');

        % Bandpass filter data (ERPlab toolbox function) 
        EEG = pop_basicfilter(EEG, [1 : EEG.nbchan], 'Cutoff', [1, 30], 'Design', 'butter', 'Filter', 'bandpass', 'Order', 4, 'RemoveDC', 'on', 'Boundary', 'boundary');

        % Bad channel detection
        [EEG, i1] = pop_rejchan(EEG, 'elec', [1 : EEG.nbchan], 'threshold', 10, 'norm', 'on', 'measure', 'kurt');
        [EEG, i2] = pop_rejchan(EEG, 'elec', [1 : EEG.nbchan], 'threshold', 5, 'norm', 'on', 'measure', 'prob');
        EEG.chans_rejected = horzcat(i1, i2);
        EEG.chans_rejected_n = length(horzcat(i1, i2));

        % Reref common average
        EEG = pop_reref(EEG, []);

        % Resample data
        EEG = pop_resample(EEG, 200);

        % Epoch data 150 ms stim + jittered 2800ms ISI
        EEG = pop_epoch(EEG, {'stim'}, [-2.5, 1.5], 'newname', [num2str(id) '_seg'], 'epochinfo', 'yes');
        
        % Autoreject data before ICA
        EEG.segs_original_n = size(EEG.data, 3);
        [EEG, rejsegs] = pop_autorej(EEG, 'nogui', 'on', 'threshold', 1000, 'startprob', 5, 'maxrej', 5, 'eegplot', 'off');
        EEG.segs_rejected_before_ica = length(rejsegs);

        % Run ICA on randsample of 666 trials
        idx = randsample([1 : size(EEG.data, 3)], 666);
        ICA = pop_selectevent(EEG, 'epoch', idx, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
        ICA = pop_runica(ICA, 'extended', 1, 'interupt', 'on');
        EEG = pop_editset(EEG, 'icachansind', 'ICA.icachansind', 'icaweights', 'ICA.icaweights', 'icasphere', 'ICA.icasphere');

        % Save IC set
        EEG = pop_saveset(EEG, 'filename', [subject '_icset.set'], 'filepath', PATH_ICSET, 'check', 'on', 'savemode', 'twofiles');

    end % End subject loop

end% End part1

% ======================= PART2: Further data cleaning ====================================================================================

if ismember('part2', to_execute)

    % Init EEGlab
    addpath(PATH_EEGLAB);
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

    preprostats = [];

    % Iterating subject list
    cnt = 0;
    for s = 1 : length(subject_list)

        % Current iteration subject
        subject = subject_list{s};
        id = str2num(subject(3 : 4));

        % Things
        cnt = cnt + 1;
        preprostats(cnt, 1) = id;

        % Load data
        EEG = pop_loadset('filename', [subject '_icset.set'], 'filepath', PATH_ICSET, 'loadmode', 'all');
        preprostats(cnt, 2) = EEG.chans_rejected_n;

        % Run IClabel
        EEG = iclabel(EEG);
        EEG.ICout_IClabel = find(EEG.etc.ic_classification.ICLabel.classifications(:, 1) < 0.5);
        EEG = pop_subcomp(EEG, EEG.ICout_IClabel, 0);
        preprostats(cnt, 3) = length(EEG.ICout_IClabel);

        % Autoreject data again
        [EEG, rejsegs] = pop_autorej(EEG, 'nogui', 'on', 'threshold', 1000, 'startprob', 5, 'maxrej', 5);
        EEG.segs_rejected_after_ica = length(rejsegs);
        EEG.segs_rejected_overall_percentage = ((EEG.segs_rejected_before_ica + EEG.segs_rejected_after_ica) / EEG.segs_original_n) * 100;
        preprostats(cnt, 4) = EEG.segs_rejected_overall_percentage;

        % Interpol
        EEG = pop_interp(EEG, EEG.chanlocs_original, 'spherical');

        % Save number of trials
        preprostats(cnt, 5) = size(EEG.data, 3);

        % Save data
        EEG = pop_saveset(EEG, 'filename', [subject '_autocleaned.set'], 'filepath', PATH_AUTOCLEANED, 'check', 'on', 'savemode', 'twofiles');

    end % End subject loop

    % Save another thing
    dlmwrite([PATH_AUTOCLEANED 'preprostats.csv'], preprostats);

end % End part2

% Visual inspection of remaining ICs.
% ICs revoved after visual inspection (index):
% Vp04: 6
% Vp08: 2
% Vp11: 10
% Vp13: 24, 31
% Vp40: 18

% ======================= PART3: time frequency decomposition ===================================================================================================

if ismember('part3', to_execute)

    % Init EEGlab
    addpath(PATH_EEGLAB);
    addpath(PATH_CUSTOM);

    % Set complex Morlet wavelet parameters
    n_frq = 50;
    frqrange = [2, 25];
    tfres_range = [400, 100];
    EEG = pop_loadset('filename', [subject_list{1} '_autocleaned.set'], 'filepath', PATH_AUTOCLEANED, 'loadmode', 'info');

	% Set wavelet time
    wtime = -2 : 1 / EEG.srate : 2;
    
    % Determine fft frqs
	hz = linspace(0, EEG.srate, length(wtime));

    % Create wavelet frequencies and tapering Gaussian widths in temporal domain
    tf_freqs = linspace(frqrange(1), frqrange(2), n_frq);
    fwhmTs = logspace(log10(tfres_range(1)), log10(tfres_range(2)), n_frq);

	% Init matrices for wavelets
	cmw = zeros(length(tf_freqs), length(wtime));
	cmwX = zeros(length(tf_freqs), length(wtime));
    tlim = zeros(1, length(tf_freqs));
    
    % These will contain the wavelet widths as full width at 
    % half maximum in the temporal and spectral domain
	obs_fwhmT = zeros(1, length(tf_freqs));
	obs_fwhmF = zeros(1, length(tf_freqs));

	% Create the wavelets
	for frq = 1 : length(tf_freqs)

		% Create wavelet with tapering gaussian corresponding to desired width in temporal domain
		cmw(frq, :) = exp(2 * 1i * pi * tf_freqs(frq) .* wtime) .* exp((-4 * log(2) * wtime.^2) ./ (fwhmTs(frq) / 1000)^2);

		% Normalize wavelet
		cmw(frq, :) = cmw(frq, :) ./ max(cmw(frq, :));

		% Create normalized freq domain wavelet
		cmwX(frq, :) = fft(cmw(frq, :)) ./ max(fft(cmw(frq, :)));

		% Determine observed fwhmT
		midt = dsearchn(wtime', 0);
		cmw_amp = abs(cmw(frq, :)) ./ max(abs(cmw(frq, :))); % Normalize cmw amplitude
		obs_fwhmT(frq) = wtime(midt - 1 + dsearchn(cmw_amp(midt : end)', 0.5)) - wtime(dsearchn(cmw_amp(1 : midt)', 0.5));

		% Determine observed fwhmF
		idx = dsearchn(hz', tf_freqs(frq));
		cmwx_amp = abs(cmwX(frq, :)); 
		obs_fwhmF(frq) = hz(idx - 1 + dsearchn(cmwx_amp(idx : end)', 0.5) - dsearchn(cmwx_amp(1 : idx)', 0.5));

	end

	% Define time window of analysis
    pruned_segs = [-2000, 1000];
    EEG = pop_loadset('filename', [subject_list{1} '_autocleaned.set'], 'filepath', PATH_AUTOCLEANED, 'loadmode', 'info');
    tf_times = EEG.times(dsearchn(EEG.times', pruned_segs(1)) : dsearchn(EEG.times', pruned_segs(2)));

    % Save stuff
    dlmwrite([PATH_TFDECOMP 'pruned_segs.csv'], pruned_segs); 
    dlmwrite([PATH_TFDECOMP 'tf_times.csv'], tf_times); 
    dlmwrite([PATH_TFDECOMP 'tf_freqs.csv'], tf_freqs);
    dlmwrite([PATH_TFDECOMP 'fwhmT.csv'], obs_fwhmT); 
    dlmwrite([PATH_TFDECOMP 'fwhmF.csv'], obs_fwhmF); 

    % Iterating subject list
    for s = 1 : length(subject_list)

        subject = subject_list{s};
        id = str2num(subject(3 : 4));

        % Load clean data
        EEG = pop_loadset('filename', [subject '_autocleaned.set'], 'filepath', PATH_AUTOCLEANED, 'loadmode', 'all');
        d = double(EEG.data);

        % tf decomp
        for ch = 1 : size(d, 1)

            % Talk
            fprintf('\ntf decomp subject %i/%i | chan %i/%i...\n', s, numel(subject_list), ch, size(d, 1));

            % Pick channel data
            dch = squeeze(d(ch, :, :));

            % Set convolution length
            convlen = size(dch, 1) * size(dch, 2) + size(cmw, 2) - 1;

            % cmw to freq domain and scale
            cmwX = zeros(n_frq, convlen);
            for f = 1 : n_frq
                cmwX(f, :) = fft(cmw(f, :), convlen);
                cmwX(f, :) = cmwX(f, :) ./ max(cmwX(f, :));
            end

            % Get TF-power
            powcube = NaN(n_frq, size(dch, 1), size(dch, 2));
            tmp = fft(reshape(double(dch), 1, []), convlen);
            for f = 1 : n_frq
                as = ifft(cmwX(f, :) .* tmp); 
                as = as(((size(cmw, 2) - 1) / 2) + 1 : end - ((size(cmw, 2) - 1) / 2));
                as = reshape(as, size(dch, 1), size(dch, 2));
                powcube(f, :, :) = abs(as) .^ 2;          
            end

            % Cut edge artifacts
            powcube = powcube(:, dsearchn(EEG.times', pruned_segs(1)) : dsearchn(EEG.times', pruned_segs(2)), :);

            % Save single trial pows
            save([PATH_TFDECOMP subject '_powcube_chan_' num2str(ch)], 'powcube');

        end % End chanit

        % Build and save meta
        tmp = EEG.event(find(strcmpi({EEG.event.type}, 'stim')));
        meta = [cell2mat({tmp.id})',...
                cell2mat({tmp.bl})',...
                cell2mat({tmp.sbl})',...
                cell2mat({tmp.sbl_total})',...
                cell2mat({tmp.stimnum})',...
                cell2mat({tmp.stimnum_sbl})',...
                cell2mat({tmp.corresp})',...
                cell2mat({tmp.rt})',...
                cell2mat({tmp.acc})'];
        dlmwrite([PATH_TFDECOMP num2str(id) '_powcube_meta.csv'], meta); 

    end % End subject loop

end % End part3

% ===== PART4: DATA DESCRIPTIVES / PLOTS OF TF_SPACE AND FREQBAND TRACES & BEHAVIOR AT 4 TIME POINTS =====================================

if ismember('part4', to_execute)

    % Load info (chanlocs...)
    addpath(PATH_EEGLAB);
    eeglab;
    EEG = pop_loadset('filename', [subject_list{1} '_autocleaned.set'], 'filepath', PATH_AUTOCLEANED, 'loadmode', 'info');
    chanlocs = EEG.chanlocs;

    % Define stuff
    tf_times = dlmread([PATH_TFDECOMP 'tf_times.csv']); 
    tf_freqs = dlmread([PATH_TFDECOMP 'tf_freqs.csv']);
    
    % Set time windows as trial ranges for ersp calculation
    trialwins = {[1 : 200], [1241 : 1440], [2681 : 2880], [4121 : 4320]};

    % Define segment to analyze
    pruned_segs = [-1500, 1000];
    prune_time = tf_times(dsearchn(tf_times', pruned_segs(1)) : dsearchn(tf_times', pruned_segs(2)));

    % tf-power (ersp) at 4 time windows throughout the experiment for FCz and Pz with no, trial-general and trial-specific baseline
    raw_c = zeros(4, length(tf_freqs), length(prune_time));
    raw_nc = zeros(4, length(tf_freqs), length(prune_time));
    gen_c = zeros(4, length(tf_freqs), length(prune_time));
    gen_nc = zeros(4, length(tf_freqs), length(prune_time));

    % Collect data
    for s = 1 : length(subject_list)

        % Load subject meta
        subject = subject_list{s};
        id = str2num(subject(3 : 4));
        meta = dlmread([PATH_TFDECOMP num2str(id) '_powcube_meta.csv']);
   
        % Iterate channels
        for ch = 1 : EEG.nbchan

            % Talk. A little bit at least...
            fprintf('\nReading data. Sub: %i - Chan: %i...\n', s, ch);

            % Load data
            load([PATH_TFDECOMP subject '_powcube_chan_' num2str(ch)]);

            % Cut edge at beginning of trial (preceding trial activity...)
            powcube = powcube(:, dsearchn(tf_times', pruned_segs(1)) : dsearchn(tf_times', pruned_segs(2)), :);
            
            % Z-Standardize
            zcube = zeros(size(powcube));
            trial_bl_gen = [-1500, -200];
            tmp = powcube(:, prune_time >= trial_bl_gen(1) & prune_time <= trial_bl_gen(2), :);
            tmp = reshape(tmp, [size(tmp, 1), size(tmp, 2) * size(tmp, 3)]);
            mpow = repmat(mean(tmp, 2), [1, size(powcube, 2)]);
            stdpow = repmat(std(tmp, [], 2), [1, size(powcube, 2)]);
            for t = 1 : size(powcube, 3)
                zcube(:, :, t) = (squeeze(powcube(:, :, t)) - mpow) ./ stdpow;
            end

            % Sum up everything!
            for w = 1 : numel(trialwins)
                raw_c(w, :, :) = squeeze(raw_c(w, :, :)) + squeeze(mean(powcube(:, :, ismember(meta(:, 5), trialwins{w}) & meta(:, 7) == 1), 3));
                raw_nc(w, :, :) = squeeze(raw_nc(w, :, :)) + squeeze(mean(powcube(:, :, ismember(meta(:, 5), trialwins{w}) & meta(:, 7) == 2), 3));
                gen_c(w, :, :) = squeeze(gen_c(w, :, :)) + squeeze(mean(zcube(:, :, ismember(meta(:, 5), trialwins{w}) & meta(:, 7) == 1), 3));
                gen_nc(w, :, :) = squeeze(gen_nc(w, :, :)) + squeeze(mean(zcube(:, :, ismember(meta(:, 5), trialwins{w}) & meta(:, 7) == 2), 3));
            end
        end
    end

    % Average and save time-frequency representations for trialwins
    for w = 1 : numel(trialwins)
        raw_c(w, :, :) = squeeze(raw_c(w, :, :)) / (length(subject_list) * EEG.nbchan);
        dlmwrite([PATH_VEUSZ, 'descriptives_eeg_tf/', 'raw_c_win_' num2str(w) '.csv'], squeeze(raw_c(w, :, :)));
        raw_nc(w, :, :) = squeeze(raw_nc(w, :, :)) / (length(subject_list) * EEG.nbchan);
        dlmwrite([PATH_VEUSZ, 'descriptives_eeg_tf/', 'raw_nc_win_' num2str(w) '.csv'], squeeze(raw_nc(w, :, :)));
        gen_c(w, :, :) = squeeze(gen_c(w, :, :)) / (length(subject_list) * EEG.nbchan);
        dlmwrite([PATH_VEUSZ, 'descriptives_eeg_tf/', 'gen_c_win_' num2str(w) '.csv'], squeeze(gen_c(w, :, :)));
        gen_nc(w, :, :) = squeeze(gen_nc(w, :, :)) / (length(subject_list) * EEG.nbchan);
        dlmwrite([PATH_VEUSZ, 'descriptives_eeg_tf/', 'gen_nc_win_' num2str(w) '.csv'], squeeze(gen_nc(w, :, :)));
    end

    % Save frequency band lineplots and baseline spectrum
    for w = 1 : numel(trialwins)
        dlmwrite([PATH_VEUSZ, 'descriptives_eeg_lineplots/', 'c_theta_win_' num2str(w) '.csv'], mean(squeeze(raw_c(w, tf_freqs >= 4 & tf_freqs <= 7, :)), 1));
        dlmwrite([PATH_VEUSZ, 'descriptives_eeg_lineplots/', 'c_alpha_win_' num2str(w) '.csv'], mean(squeeze(raw_c(w, tf_freqs >= 8 & tf_freqs <= 12, :)), 1));
        dlmwrite([PATH_VEUSZ, 'descriptives_eeg_lineplots/', 'c_beta_win_' num2str(w) '.csv'], mean(squeeze(raw_c(w, tf_freqs >= 16 & tf_freqs <= 25, :)), 1));
        dlmwrite([PATH_VEUSZ, 'descriptives_eeg_lineplots/', 'c_spectrum_win_' num2str(w) '.csv'], mean(squeeze(raw_c(w, :, prune_time >= trial_bl_gen(1) & prune_time <= trial_bl_gen(2))), 2)');
        dlmwrite([PATH_VEUSZ, 'descriptives_eeg_lineplots/', 'nc_theta_win_' num2str(w) '.csv'], mean(squeeze(raw_nc(w, tf_freqs >= 4 & tf_freqs <= 7, :)), 1));
        dlmwrite([PATH_VEUSZ, 'descriptives_eeg_lineplots/', 'nc_alpha_win_' num2str(w) '.csv'], mean(squeeze(raw_nc(w, tf_freqs >= 8 & tf_freqs <= 12, :)), 1));
        dlmwrite([PATH_VEUSZ, 'descriptives_eeg_lineplots/', 'nc_beta_win_' num2str(w) '.csv'], mean(squeeze(raw_nc(w, tf_freqs >= 16 & tf_freqs <= 25, :)), 1));
        dlmwrite([PATH_VEUSZ, 'descriptives_eeg_lineplots/', 'nc_spectrum_win_' num2str(w) '.csv'], mean(squeeze(raw_nc(w, :, prune_time >= trial_bl_gen(1) & prune_time <= trial_bl_gen(2))), 2)');
    end

    % Save an x axis vectors for lineplots
    dlmwrite([PATH_VEUSZ, 'descriptives_eeg_lineplots/', 'xdat.csv'], [1 : size(raw_c, 3)]);
    dlmwrite([PATH_VEUSZ, 'descriptives_eeg_lineplots/', 'xdatfrq.csv'], [1 : size(raw_c, 2)]);

    % Collect behavioral data
    b_rt_c = zeros(length(subject_list), 4);
    b_ac_c = zeros(length(subject_list), 4);
    b_rt_nc = zeros(length(subject_list), 4);
    b_ac_nc = zeros(length(subject_list), 4);
    B = [];
    for s = 1 : length(subject_list)

        % Load subject meta
        subject = subject_list{s};
        id = str2num(subject(3 : 4));
        meta = dlmread([PATH_TFDECOMP num2str(id) '_powcube_meta.csv']);

         % Iterate trialwins
        for w = 1 : numel(trialwins)
            b_rt_c(s, w) = mean(meta(ismember(meta(:, 5), trialwins{w}) & meta(:, 9) == 1 & meta(:, 7) == 1, 8));
            b_ac_c(s, w) = sum(ismember(meta(:, 5), trialwins{w}) & meta(:, 9) == 1 & meta(:, 7) == 1) / sum(ismember(meta(:, 5), trialwins{w}) & meta(:, 7) == 1);
            b_rt_nc(s, w) = mean(meta(ismember(meta(:, 5), trialwins{w}) & meta(:, 9) == 1 & meta(:, 7) == 2, 8));
            b_ac_nc(s, w) = sum(ismember(meta(:, 5), trialwins{w}) & meta(:, 9) == 1 & meta(:, 7) == 2) / sum(ismember(meta(:, 5), trialwins{w}) & meta(:, 7) == 2);
        end

        % Piling up behavioral data...
        if s == 1
            B = meta;
        else
            B = [B; meta];
        end

    end % End subject iteration

    for s = 1 : length(subject_list)

        % Load subject meta
        subject = subject_list{s};
        id = str2num(subject(3 : 4));
        rt_c(s) = mean(B(B(:, 1) == id & B(:, 9) == 1 & B(:, 7) == 1, 8));
        rt_nc(s) = mean(B(B(:, 1) == id & B(:, 9) == 1 & B(:, 7) == 2, 8));
        ac_c(s) = (length(B(B(:, 1) == id & B(:, 9) == 1 & B(:, 7) == 1, 9)) / length(B(B(:, 1) == id & B(:, 7) == 1, 9))) * 100;
        ac_nc(s) = (length(B(B(:, 1) == id & B(:, 9) == 1 & B(:, 7) == 2, 9)) / length(B(B(:, 1) == id & B(:, 7) == 2, 9))) * 100;
    end

    % Table for plotting
    b_out = zeros(4, 8);
    for w = 1 : numel(trialwins)
        b_out(w, :) = [mean(b_rt_c(:, w)), std(b_rt_c(:, w)), mean(b_ac_c(:, w)) * 100, std(b_ac_c(:, w)) * 100,...
                       mean(b_rt_nc(:, w)), std(b_rt_nc(:, w)), mean(b_ac_nc(:, w)) * 100, std(b_ac_nc(:, w)) * 100];
    end
    dlmwrite([PATH_VEUSZ, 'descriptives_behavior/', 'behavior.csv'], b_out, 'delimiter', '\t');
    dlmwrite([PATH_VEUSZ, 'descriptives_behavior/', 'xax_wins.csv'], [1 : numel(trialwins)]);

    % Create tables
    varnames = {'id', 'bl' , 'sbl', 'sbl_total', 'stimnum', 'stimnum_sbl', 'corresp', 'rt', 'acc'};
    M = B(B(:, 9) == 1, :);
    tbl_rt = table(M(:, 1), M(:, 2), M(:, 3), M(:, 4), M(:, 5), M(:, 6), M(:, 7), M(:, 8), M(:, 9), 'VariableNames', varnames);
    M = B;
    tbl_ac = table(M(:, 1), M(:, 2), M(:, 3), M(:, 4), M(:, 5), M(:, 6), M(:, 7), M(:, 8), M(:, 9), 'VariableNames', varnames);

    % Cast vars
    tbl_rt.id = nominal(tbl_rt.id);
    tbl_ac.id = nominal(tbl_ac.id);
    tbl_rt.corresp = nominal(tbl_rt.corresp);
    tbl_ac.corresp = nominal(tbl_ac.corresp);

    % Compute LMEs
    lme_rt = fitlme(tbl_rt, 'rt ~ stimnum*corresp + (1|id)');
    lme_ac = fitlme(tbl_ac, 'acc ~ stimnum*corresp + (1|id)');


end % End part4

% ======================= PART5: Calculate linear regression on sensor space data ================================================

if ismember('part5', to_execute)

    % Define a trial baseline interval
    trial_bl = [-1500, -200];

    % Load info
    addpath(PATH_EEGLAB);
    eeglab;
    EEG = pop_loadset('filename', [subject_list{1} '_autocleaned.set'], 'filepath', PATH_AUTOCLEANED, 'loadmode', 'info');

    % Load tf params
    tf_times = dlmread([PATH_TFDECOMP 'tf_times.csv']); 
    tf_freqs = dlmread([PATH_TFDECOMP 'tf_freqs.csv']);

    % Define segment to analyze
    pruned_segs = [-1500, 1000];
    prune_time = tf_times(dsearchn(tf_times', pruned_segs(1)) : dsearchn(tf_times', pruned_segs(2)));

    % Iterating subject list
    for s = 1 : length(subject_list)

        % Load subject meta
        subject = subject_list{s};
        id = str2num(subject(3 : 4));
        meta = dlmread([PATH_TFDECOMP num2str(id) '_powcube_meta.csv']);

        % Split meta c and nc
        meta_c = meta(meta(:, 7) == 1, :);
        meta_nc = meta(meta(:, 7) == 2, :);

        % Regression design matrices
        desmat = [ones(size(meta, 1), 1), meta(:, 5)];
        desmat_c = [ones(size(meta_c, 1), 1), meta_c(:, 5)];
        desmat_nc = [ones(size(meta_nc, 1), 1), meta_nc(:, 5)];

        % Scale predictors
        desmat(:, 2) = desmat(:, 2) / max(abs(desmat(:, 2)));
        desmat_c(:, 2) = desmat_c(:, 2) / max(abs(desmat_c(:, 2)));
        desmat_nc(:, 2) = desmat_nc(:, 2) / max(abs(desmat_nc(:, 2)));

        % Result structs
        tt_c.powspctrm = zeros(EEG.nbchan, numel(tf_freqs), numel(prune_time));
        tt_c.dimord    = 'chan_freq_time';
        tt_c.label     = {EEG.chanlocs.labels};
        tt_c.freq      = tf_freqs;
        tt_c.time      = prune_time;

        tt_nc.powspctrm = zeros(EEG.nbchan, numel(tf_freqs), numel(prune_time));
        tt_nc.dimord    = 'chan_freq_time';
        tt_nc.label     = {EEG.chanlocs.labels};
        tt_nc.freq      = tf_freqs;
        tt_nc.time      = prune_time;

        tot_tru.powspctrm = zeros(EEG.nbchan, numel(tf_freqs), numel(prune_time));
        tot_tru.dimord    = 'chan_freq_time';
        tot_tru.label     = {EEG.chanlocs.labels};
        tot_tru.freq      = tf_freqs;
        tot_tru.time      = prune_time;

        tot_fak.powspctrm = zeros(EEG.nbchan, numel(tf_freqs), numel(prune_time));
        tot_fak.dimord    = 'chan_freq_time';
        tot_fak.label     = {EEG.chanlocs.labels};
        tot_fak.freq      = tf_freqs;
        tot_fak.time      = prune_time;

        tot_c_tru.powspctrm = zeros(EEG.nbchan, numel(tf_freqs), numel(prune_time));
        tot_c_tru.dimord    = 'chan_freq_time';
        tot_c_tru.label     = {EEG.chanlocs.labels};
        tot_c_tru.freq      = tf_freqs;
        tot_c_tru.time      = prune_time;

        tot_nc_tru.powspctrm = zeros(EEG.nbchan, numel(tf_freqs), numel(prune_time));
        tot_nc_tru.dimord    = 'chan_freq_time';
        tot_nc_tru.label     = {EEG.chanlocs.labels};
        tot_nc_tru.freq      = tf_freqs;
        tot_nc_tru.time      = prune_time;

        % regression
        for ch = 1 : EEG.nbchan

            % Talk
            fprintf('\nSingle trial regression subject %i/%i | chan %i/%i...\n', s, numel(subject_list), ch, EEG.nbchan);

            % Load data
            load([PATH_TFDECOMP subject '_powcube_chan_' num2str(ch)]);

            % Prune data
            powcube = powcube(:, dsearchn(tf_times', pruned_segs(1)) : dsearchn(tf_times', pruned_segs(2)), :);

            % Z-Standardize trials
            zcube = zeros(size(powcube));
            tmp = powcube(:, prune_time >= trial_bl(1) & prune_time <= trial_bl(2), :);
            tmp = reshape(tmp, [size(tmp, 1), size(tmp, 2) * size(tmp, 3)]);
            mpow = repmat(mean(tmp, 2), [1, size(zcube, 2)]);
            stdpow = repmat(std(tmp, [], 2), [1, size(zcube, 2)]);
            for t = 1 : size(powcube, 3)
                zcube(:, :, t) = (squeeze(powcube(:, :, t)) - mpow) ./ stdpow;
            end

            % Split data
            zcube_c = zcube(:, :, meta(:, 7) == 1);
            zcube_nc = zcube(:, :, meta(:, 7) == 2);

            % Save corresponding and non-corresponding power data for main effect
            tt_c.powspctrm(ch, :, :) = squeeze(mean(zcube_c, 3));
            tt_nc.powspctrm(ch, :, :) = squeeze(mean(zcube_nc, 3));

            % Reshape data (trials in rows, each tf point a column)
            d = reshape(zcube, numel(prune_time) * numel(tf_freqs), size(zcube, 3))';
            d_c = reshape(zcube_c, numel(prune_time) * numel(tf_freqs), size(zcube_c, 3))';
            d_nc = reshape(zcube_nc, numel(prune_time) * numel(tf_freqs), size(zcube_nc, 3))';

            % OLS fit
            tmp = (desmat' * desmat) \ desmat' * d;
            tot_tru.powspctrm(ch, :, :) = reshape(squeeze(tmp(2, :)), [numel(tf_freqs), numel(prune_time)]);

            tmp = (desmat_c' * desmat_c) \ desmat_c' * d_c;
            tot_c_tru.powspctrm(ch, :, :) = reshape(squeeze(tmp(2, :)), [numel(tf_freqs), numel(prune_time)]);

            tmp = (desmat_nc' * desmat_nc) \ desmat_nc' * d_nc;
            tot_nc_tru.powspctrm(ch, :, :) = reshape(squeeze(tmp(2, :)), [numel(tf_freqs), numel(prune_time)]);

            % Generate null hypothesis distribution data
            fakedesmat = desmat;
            fakedesmat(:, 2) = desmat(randperm(size(desmat, 1)), 2);
            tmp = (fakedesmat' * fakedesmat) \ fakedesmat' * d;
            tot_fak.powspctrm(ch, :, :) = reshape(squeeze(tmp(2, :)), [numel(tf_freqs), numel(prune_time)]);

        end % End chanit

        % Save regression data
        save([PATH_REGRESSION 'tt_c_' subject], 'tt_c');
        save([PATH_REGRESSION 'tt_nc_' subject], 'tt_nc');
        save([PATH_REGRESSION 'tot_c_tru_' subject], 'tot_c_tru');
        save([PATH_REGRESSION 'tot_nc_tru_' subject], 'tot_nc_tru');
        save([PATH_REGRESSION 'tot_tru_' subject], 'tot_tru');
        save([PATH_REGRESSION 'tot_fak_' subject], 'tot_fak');

    end % End subit

end % End part5

% ======================= PART6: Cluster based permutation tests ==============================

if ismember('part6', to_execute)

    % Load info
    addpath(PATH_EEGLAB);
    eeglab;
    EEG = pop_loadset('filename', [subject_list{1} '_autocleaned.set'], 'filepath', PATH_AUTOCLEANED, 'loadmode', 'info');

    % Build elec struct
    for ch = 1 : EEG.nbchan
        elec.label{ch} = EEG.chanlocs(ch).labels;
        elec.elecpos(ch, :) = [EEG.chanlocs(ch).X, EEG.chanlocs(ch).Y, EEG.chanlocs(ch).Z];
        elec.chanpos(ch, :) = [EEG.chanlocs(ch).X, EEG.chanlocs(ch).Y, EEG.chanlocs(ch).Z];
    end

    % Prepare layout
    cfg      = [];
    cfg.elec = elec;
    cfg.rotate = 90;
    layout = ft_prepare_layout(cfg);

    % Init ft
    rmpath(PATH_EEGLAB);
    addpath(PATH_FIELDTRIP);
    addpath(PATH_CUSTOM);
    ft_defaults;

    % Load data
    for s = 1 : length(subject_list)
        subject = subject_list{s};
        load([PATH_REGRESSION 'tot_tru_' subject]);
        d_tot_tru{s} = tot_tru;
        load([PATH_REGRESSION 'tot_fak_' subject]);
        d_tot_fak{s} = tot_fak;
        load([PATH_REGRESSION 'tot_c_tru_' subject]);
        d_tot_c_tru{s} = tot_c_tru;
        load([PATH_REGRESSION 'tot_nc_tru_' subject]);
        d_tot_nc_tru{s} = tot_nc_tru;
        load([PATH_REGRESSION 'tt_c_' subject]);
        d_tt_c{s} = tt_c;
        load([PATH_REGRESSION 'tt_nc_' subject]);
        d_tt_nc{s} = tt_nc;
    end

    % Calc beta GAs and save
    cfg = [];
    cfg.keepindividual = 'yes';
    GA_tot_tru = ft_freqgrandaverage(cfg, d_tot_tru{1, :});
    GA_tot_fak = ft_freqgrandaverage(cfg, d_tot_fak{1, :});
    save([PATH_CLUSTSTATS 'GA_tot_tru.mat'], 'GA_tot_tru');
    save([PATH_CLUSTSTATS 'GA_tot_fak.mat'], 'GA_tot_fak');
    GA_tot_c_tru = ft_freqgrandaverage(cfg, d_tot_c_tru{1, :});
    GA_tot_nc_tru = ft_freqgrandaverage(cfg, d_tot_nc_tru{1, :});
    save([PATH_CLUSTSTATS 'GA_tot_c_tru.mat'], 'GA_tot_c_tru');
    save([PATH_CLUSTSTATS 'GA_tot_nc_tru.mat'], 'GA_tot_nc_tru');
    GA_tt_c = ft_freqgrandaverage(cfg, d_tt_c{1, :});
    GA_tt_nc = ft_freqgrandaverage(cfg, d_tt_nc{1, :});
    save([PATH_CLUSTSTATS 'GA_tt_c.mat'], 'GA_tt_c');
    save([PATH_CLUSTSTATS 'GA_tt_nc.mat'], 'GA_tt_nc');

    % Define neighbours
    cfg                 = [];
    cfg.layout          = layout;
    cfg.feedback        = 'no';
    cfg.method          = 'triangulation'; 
    cfg.neighbours      = ft_prepare_neighbours(cfg, GA_tot_c_tru);
    neighbours          = cfg.neighbours;

    % Testparams
    testalpha   = 0.025;
    voxelalpha  = 0.01;
    nperm       = 1000;

    % Set config. Same for all tests
    cfg = [];
    cfg.tail             = 0;
    cfg.statistic        = 'depsamplesT';
    cfg.alpha            = testalpha;
    cfg.neighbours       = neighbours;
    cfg.minnbchan        = 2;
    cfg.method           = 'montecarlo';
    cfg.correctm         = 'cluster';
    cfg.clustertail      = 0;
    cfg.clusteralpha     = voxelalpha;
    cfg.clusterstatistic = 'maxsum';
    cfg.numrandomization = nperm;
    cfg.computecritval   = 'yes'; 
    cfg.ivar             = 1;
    cfg.uvar             = 2;
    cfg.design           = [ones(1, numel(subject_list)), 2 * ones(1, numel(subject_list)); 1 : numel(subject_list), 1 : numel(subject_list)];

    % The tests
    [stat_tot] = ft_freqstatistics(cfg, GA_tot_fak, GA_tot_tru);
    [stat_tot_interaction] = ft_freqstatistics(cfg, GA_tot_nc_tru, GA_tot_c_tru);
    [stat_tt] = ft_freqstatistics(cfg, GA_tt_nc, GA_tt_c);

    % Calculate and save effect sizes
    adjpetasq_tot = [];
    adjpetasq_tot_interaction = [];
    adjpetasq_tt = [];
    for ch = 1 : EEG.nbchan
        petasq = (squeeze(stat_tot.stat(ch, :, :)) .^ 2) ./ ((squeeze(stat_tot.stat(ch, :, :)) .^ 2) + (numel(subject_list) - 1));
        adj_petasq = petasq - (1 - petasq) .* (1 / (numel(subject_list) - 1));
        adjpetasq_tot(ch, :, :) = adj_petasq;
        petasq = (squeeze(stat_tot_interaction.stat(ch, :, :)) .^ 2) ./ ((squeeze(stat_tot_interaction.stat(ch, :, :)) .^ 2) + (numel(subject_list) - 1));
        adj_petasq = petasq - (1 - petasq) .* (1 / (numel(subject_list) - 1));
        adjpetasq_tot_interaction(ch, :, :) = adj_petasq;
        petasq = (squeeze(stat_tt.stat(ch, :, :)) .^ 2) ./ ((squeeze(stat_tt.stat(ch, :, :)) .^ 2) + (numel(subject_list) - 1));
        adj_petasq = petasq - (1 - petasq) .* (1 / (numel(subject_list) - 1));
        adjpetasq_tt(ch, :, :) = adj_petasq;
    end
    save([PATH_CLUSTSTATS 'effect_size_cube_tot.mat'], 'adjpetasq_tot');
    save([PATH_CLUSTSTATS 'effect_size_cube_tot_interaction.mat'], 'adjpetasq_tot_interaction');
    save([PATH_CLUSTSTATS 'effect_size_cube_tt.mat'], 'adjpetasq_tt');

    % Identify significant clusters
    clusts = struct();
    cnt = 0;
    stat_names = {'stat_tot', 'stat_tot_interaction', 'stat_tt'};
    for s = 1 : numel(stat_names)
        stat = eval(stat_names{s});
        if ~isempty(stat.negclusters)
            neg_idx = find([stat.negclusters(1, :).prob] < testalpha);
            for c = 1 : numel(neg_idx)
                cnt = cnt + 1;
                clusts(cnt).testlabel = stat_names{s};
                clusts(cnt).clustnum = cnt;
                clusts(cnt).time = stat.time;
                clusts(cnt).freq = stat.freq;
                clusts(cnt).polarity = -1;
                clusts(cnt).prob = stat.negclusters(1, neg_idx(c)).prob;
                clusts(cnt).idx = stat.negclusterslabelmat == neg_idx(c);
                clusts(cnt).stats = clusts(cnt).idx .* stat.stat * -1;
                clusts(cnt).chans_sig = find(logical(mean(clusts(cnt).idx, [2,3])));
            end
        end
        if ~isempty(stat.posclusters)
            pos_idx = find([stat.posclusters(1, :).prob] < testalpha);
            for c = 1 : numel(pos_idx)
                cnt = cnt + 1;
                clusts(cnt).testlabel = stat_names{s};
                clusts(cnt).clustnum = cnt;
                clusts(cnt).time = stat.time;
                clusts(cnt).freq = stat.freq;
                clusts(cnt).polarity = 1;
                clusts(cnt).prob = stat.posclusters(1, pos_idx(c)).prob;
                clusts(cnt).idx = stat.posclusterslabelmat == pos_idx(c);
                clusts(cnt).stats = clusts(cnt).idx .* stat.stat;
                clusts(cnt).chans_sig = find(logical(mean(clusts(cnt).idx, [2, 3])));
            end
        end
    end

    % Save cluster struct
    save([PATH_CLUSTSTATS 'significant_clusters.mat'], 'clusts');

    % Plot identified cluster
    clinecol = 'k';
    cmap = 'jet';
    chanlocs = EEG.chanlocs;
    for cnt = 1 : numel(clusts)

        figure('Visible', 'off'); clf;

        subplot(2, 2, 1)
        pd = squeeze(sum(clusts(cnt).stats, 1));
        contourf(clusts(cnt).time, clusts(cnt).freq, pd, 40, 'linecolor','none')
        hold on
        contour(clusts(cnt).time, clusts(cnt).freq, logical(squeeze(mean(clusts(cnt).idx, 1))), 1, 'linecolor', clinecol, 'LineWidth', 2)
        colormap(cmap)
        set(gca, 'xlim', [clusts(cnt).time(1), clusts(cnt).time(end)], 'clim', [-max(abs(pd(:))), max(abs(pd(:)))], 'YScale', 'lin', 'YTick', [4, 8, 12, 20])
        colorbar;
        title(['sum t across chans, plrt: ' num2str(clusts(cnt).polarity)], 'FontSize', 10)

        subplot(2, 2, 2)
        pd = squeeze(mean(clusts(cnt).idx, 1));
        contourf(clusts(cnt).time, clusts(cnt).freq, pd, 40, 'linecolor','none')
        hold on
        contour(clusts(cnt).time, clusts(cnt).freq, logical(squeeze(mean(clusts(cnt).idx, 1))), 1, 'linecolor', clinecol, 'LineWidth', 2)
        colormap(cmap)
        set(gca, 'xlim', [clusts(cnt).time(1), clusts(cnt).time(end)], 'clim', [-1, 1], 'YScale', 'lin', 'YTick', [4, 8, 12, 20])
        colorbar;
        title(['proportion chans significant'], 'FontSize', 10)

        subplot(2, 2, 3)
        pd = squeeze(sum(clusts(cnt).stats, [2, 3]));
        topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'on');
        colormap(cmap)
        set(gca, 'clim', [-max(abs(pd(:))), max(abs(pd(:)))])
        colorbar;
        title(['sum t per electrode'], 'FontSize', 10)

        subplot(2, 2, 4)
        pd = squeeze(mean(clusts(cnt).idx, [2, 3]));
        topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'on');
        colormap(cmap)
        set(gca, 'clim', [-1, 1])
        colorbar;
        title(['proportion tf-points significant'], 'FontSize', 10)

        saveas(gcf, [PATH_PLOT 'clustnum_' num2str(clusts(cnt).clustnum) '_' clusts(cnt).testlabel '.png']); 
    end

end % End part6

% ======================= PART7: Data for plotting ============================

if ismember('part7', to_execute)

    % Load info (chanlocs...)
    addpath(PATH_EEGLAB);
    eeglab;
    EEG = pop_loadset('filename', [subject_list{1} '_autocleaned.set'], 'filepath', PATH_AUTOCLEANED, 'loadmode', 'info');
    chanlocs = EEG.chanlocs;

    % Load tf params
    tf_times = dlmread([PATH_TFDECOMP 'tf_times.csv']); 
    tf_freqs = dlmread([PATH_TFDECOMP 'tf_freqs.csv']);
    pruned_segs = [-1500, 1000];
    prune_time = tf_times(dsearchn(tf_times', pruned_segs(1)) : dsearchn(tf_times', pruned_segs(2)));
    
    % Load sigificant cluster struct as'clusts'
    load([PATH_CLUSTSTATS 'significant_clusters.mat']);

    % Save Cluster contours
    dlmwrite([PATH_VEUSZ, 'inference_sensor_space/', 'cluster_tot_contour.csv'], logical(squeeze(mean(clusts(1).idx, 1))));
    dlmwrite([PATH_VEUSZ, 'inference_sensor_space/', 'cluster_tt_contour.csv'], logical(squeeze(mean(clusts(2).idx, 1))));

    % S-R correspondence cluster range in time
    clu_start = prune_time(min(find(sum(logical(squeeze(mean(clusts(2).idx, 1))), 1))));

    % Load effect sizes
    load([PATH_CLUSTSTATS 'effect_size_cube_tot.mat']);
    load([PATH_CLUSTSTATS 'effect_size_cube_tot_interaction.mat']);
    load([PATH_CLUSTSTATS 'effect_size_cube_tt.mat']);

    % Save effect sizes averaged across channels
    dlmwrite([PATH_VEUSZ, 'inference_sensor_space/' , 'effect_size_tot.csv'], squeeze(mean(adjpetasq_tot, 1)));
    dlmwrite([PATH_VEUSZ, 'inference_sensor_space/' , 'effect_size_tot_interaction.csv'], squeeze(mean(adjpetasq_tot_interaction, 1)));
    dlmwrite([PATH_VEUSZ, 'inference_sensor_space/' , 'effect_size_tt.csv'], squeeze(mean(adjpetasq_tt, 1)));

    % Load GA data
    load([PATH_CLUSTSTATS 'GA_tot_tru.mat']);
    load([PATH_CLUSTSTATS 'GA_tot_c_tru.mat']);
    load([PATH_CLUSTSTATS 'GA_tot_nc_tru.mat']);
    load([PATH_CLUSTSTATS 'GA_tt_c.mat']);
    load([PATH_CLUSTSTATS 'GA_tt_nc.mat']);

    % Save beta coefficients tot
    tot_betas = squeeze(mean(squeeze(mean(GA_tot_tru.powspctrm, 1)), 1));
    dlmwrite([PATH_VEUSZ, 'inference_sensor_space/', 'tot_betas.csv'], tot_betas);
    tot_betas_c = squeeze(mean(squeeze(mean(GA_tot_c_tru.powspctrm, 1)), 1));
    dlmwrite([PATH_VEUSZ, 'inference_sensor_space/', 'tot_betas_c.csv'], tot_betas_c);
    tot_betas_nc = squeeze(mean(squeeze(mean(GA_tot_nc_tru.powspctrm, 1)), 1));
    dlmwrite([PATH_VEUSZ, 'inference_sensor_space/', 'tot_betas_nc.csv'], tot_betas_nc);
    tt_ersp_c = squeeze(mean(squeeze(mean(GA_tt_c.powspctrm, 1)), 1));
    dlmwrite([PATH_VEUSZ, 'inference_sensor_space/', 'tt_ersp_c.csv'], tt_ersp_c);
    tt_ersp_nc = squeeze(mean(squeeze(mean(GA_tt_nc.powspctrm, 1)), 1));
    dlmwrite([PATH_VEUSZ, 'inference_sensor_space/', 'tt_ersp_nc.csv'], tt_ersp_nc);
    tt_ersp_c = squeeze(mean(squeeze(mean(GA_tt_c.powspctrm, 1)), 1));
    dlmwrite([PATH_VEUSZ, 'inference_sensor_space/', 'tt_ersp_c.csv'], tt_ersp_c);
    tt_ersp_diff = tt_ersp_nc - tt_ersp_c;
    dlmwrite([PATH_VEUSZ, 'inference_sensor_space/', 'tt_ersp_diff.csv'], tt_ersp_diff);
    
    % Topography tot beta
    pd = [];
    d = squeeze(mean(GA_tot_tru.powspctrm, 1));
    for ch = 1 : EEG.nbchan
        tmp = squeeze(d(ch, :, :)); % Channel betas
        pd(ch) = mean(tmp(:));
    end
    markercolor = 'k';
    cmap = 'jet';
    clim = [-0.2, 0.2];
    figure('Visible', 'off'); clf;
    topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'off', 'emarker2', {clusts(1).chans_sig, 'p', markercolor, 14, 1});
    colormap(cmap);
    caxis(clim);
    saveas(gcf, [PATH_VEUSZ, 'inference_sensor_space/', 'topo_tot_beta.png']);

    % Topography tot beta c
    pd = [];
    d = squeeze(mean(GA_tot_c_tru.powspctrm, 1));
    for ch = 1 : EEG.nbchan
        tmp = squeeze(d(ch, :, :)); % Channel betas
        pd(ch) = mean(tmp(:));
    end
    markercolor = 'k';
    cmap = 'jet';
    clim = [-0.2, 0.2];
    figure('Visible', 'off'); clf;
    topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'off', 'emarker2', {[], 'p', markercolor, 14, 1});
    colormap(cmap);
    caxis(clim);
    saveas(gcf, [PATH_VEUSZ, 'inference_sensor_space/', 'topo_tot_beta_c.png']);

    % Topography tot beta nc
    pd = [];
    d = squeeze(mean(GA_tot_nc_tru.powspctrm, 1));
    for ch = 1 : EEG.nbchan
        tmp = squeeze(d(ch, :, :)); % Channel betas
        pd(ch) = mean(tmp(:));
    end
    markercolor = 'k';
    cmap = 'jet';
    clim = [-0.2, 0.2];
    figure('Visible', 'off'); clf;
    topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'off', 'emarker2', {[], 'p', markercolor, 14, 1});
    colormap(cmap);
    caxis(clim);
    saveas(gcf, [PATH_VEUSZ, 'inference_sensor_space/', 'topo_tot_beta_nc.png']);

    % Topography ersp c
    pd = [];
    d = squeeze(mean(GA_tt_c.powspctrm, 1));
    for ch = 1 : EEG.nbchan
        tmp = squeeze(d(ch, :, :)); % Channel betas
        pd(ch) = mean(tmp(:));
    end
    markercolor = 'k';
    cmap = 'jet';
    clim = [-1, 1];
    figure('Visible', 'off'); clf;
    topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'off', 'emarker2', {clusts(2).chans_sig, 'p', markercolor, 14, 1});
    colormap(cmap);
    caxis(clim);
    saveas(gcf, [PATH_VEUSZ, 'inference_sensor_space/', 'topo_tt_ersp_c.png']);

    % Topography ersp nc
    pd = [];
    d = squeeze(mean(GA_tt_nc.powspctrm, 1));
    for ch = 1 : EEG.nbchan
        tmp = squeeze(d(ch, :, :)); % Channel betas
        pd(ch) = mean(tmp(:));
    end
    markercolor = 'k';
    cmap = 'jet';
    clim = [-1, 1];
    figure('Visible', 'off'); clf;
    topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'off', 'emarker2', {clusts(2).chans_sig, 'p', markercolor, 14, 1});
    colormap(cmap);
    caxis(clim);
    saveas(gcf, [PATH_VEUSZ, 'inference_sensor_space/', 'topo_tt_ersp_nc.png']);

    % Topography effect sizes tot
    pd = [];
    d = adjpetasq_tot;
    for ch = 1 : EEG.nbchan
        tmp = squeeze(d(ch, :, :)); % channel effect sizes
        pd(ch) = mean(tmp(:));
    end
    markercolor = 'k';
    cmap = 'jet';
    clim = [-0.5, 0.5];
    figure('Visible', 'off'); clf;
    topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'off', 'emarker2', {clusts(1).chans_sig, 'p', markercolor, 14, 1});
    colormap(cmap);
    caxis(clim);
    saveas(gcf, [PATH_VEUSZ, 'inference_sensor_space/', 'topo_adjpetasq_tot.png']);

    % Topography effect sizes interaction
    pd = [];
    d = adjpetasq_tot_interaction;
    for ch = 1 : EEG.nbchan
        tmp = squeeze(d(ch, :, :)); % channel effect sizes
        pd(ch) = mean(tmp(:));
    end
    markercolor = 'k';
    cmap = 'jet';
    clim = [-0.5, 0.5];
    figure('Visible', 'off'); clf;
    topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'off', 'emarker2', {[], 'p', markercolor, 14, 1});
    colormap(cmap);
    caxis(clim);
    saveas(gcf, [PATH_VEUSZ, 'inference_sensor_space/', 'topo_adjpetasq_tot_interaction.png']);

    % Topography effect sizes tt
    pd = [];
    d = adjpetasq_tt;
    for ch = 1 : EEG.nbchan
        tmp = squeeze(d(ch, :, :)); % channel effect sizes
        pd(ch) = mean(tmp(:));
    end
    markercolor = 'k';
    cmap = 'jet';
    clim = [-0.1, 0.1];
    figure('Visible', 'off'); clf;
    topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'off', 'emarker2', {clusts(2).chans_sig, 'p', markercolor, 14, 1});
    colormap(cmap);
    caxis(clim);
    saveas(gcf, [PATH_VEUSZ, 'inference_sensor_space/', 'topo_adjpetasq_tt.png']);

end % End part7

% ======================= PART8: Generalized eigendecomposition, apply spatial filter, time-frequency analysis and linear regression ==========================================================================

if ismember('part8', to_execute)

    % Init EEGlab
    addpath(PATH_EEGLAB);
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

	% Load info
    EEG = pop_loadset('filename', [subject_list{1} '_autocleaned.set'], 'filepath', PATH_AUTOCLEANED, 'loadmode', 'info');

	% Determine electrode distances based on cartesian coordinates (loosely adopted on elec_distance.m)
	dists = zeros(EEG.nbchan);
	cart_coords = [cell2mat({EEG.chanlocs.X})', cell2mat({EEG.chanlocs.Y})', cell2mat({EEG.chanlocs.Z})'];
	for ch1 = 1 : EEG.nbchan
		crds1 = cart_coords(ch1, :);
		len1 = sqrt(sum(crds1 .^ 2));
		for ch2 = 1 : EEG.nbchan
			crds2 = cart_coords(ch2, :);
			len2 = sqrt(sum(crds2 .^ 2));
			if ch1 == ch2
				dists(ch1, ch2) = 0;
			else
				r = (len1 + len2) / 2; % Estimate sphere radius from len1 and len2
				theta = acos(dot(crds1, crds2) / (len1 * len2)); % Angle between A & B in radians
				dists(ch1, ch2) = r * theta; % Arc length = radius * theta
			end
		end
	end

	% Create spatial filter template
	focuschan = 14;
	template_topo = dists(focuschan, :) / max(dists(focuschan, :)); % Normalize distances
	template_topo = ones(size(template_topo)) - template_topo; % Invert

	% Plot spatial filter map template
	figure('Visible', 'off'); clf;
	topoplot(template_topo, EEG.chanlocs, 'electrodes', 'off', 'numcontour', 0)
	title(['template topo'])
	set(gcf, 'PaperUnits', 'centimeters')
	set(gcf, 'PaperPosition', [0, 0, 10, 10])
    saveas(gcf, [PATH_VEUSZ, 'descriptives_spatial_filter/', 'filter_template.png']);

    % Load TF analysis parameters
    pruned_segs = dlmread([PATH_TFDECOMP 'pruned_segs.csv']); % [-2000, 1000]
    tf_times = dlmread([PATH_TFDECOMP 'tf_times.csv']); 

    % Set complex Morlet wavelet parameters
    n_frq = 50;
    frqrange = [2, 25];
    tfres_range = [400, 100];

    % Create wavelet frequencies and tapering Gaussian widths in temporal domain
    tf_freqs = linspace(frqrange(1), frqrange(2), n_frq);
    fwhmTs = logspace(log10(tfres_range(1)), log10(tfres_range(2)), n_frq);

    % Create wavelets
    wtime = -2 : 1 / EEG.srate : 2;
	for frq = 1 : length(tf_freqs)
		cmw(frq, :) = exp(2 * 1i * pi * tf_freqs(frq) .* wtime) .* exp((-4 * log(2) * wtime.^2) ./ (fwhmTs(frq) / 1000)^2);
		cmw(frq, :) = cmw(frq, :) ./ max(cmw(frq, :));
	end
    
    % Init ersp resmat for cmcomponentp data
    ersp_cmp_c = zeros(4, length(tf_freqs), length(dsearchn(tf_times', -500) : dsearchn(tf_times', 1000)));
    ersp_cmp_nc = zeros(4, length(tf_freqs), length(dsearchn(tf_times', -500) : dsearchn(tf_times', 1000)));

    % Iterating subject list
    for s = 1 : length(subject_list)

        subject = subject_list{s};
        id = str2num(subject(3 : 4));

        % Load clean data
        EEG = pop_loadset('filename', [subject '_autocleaned.set'], 'filepath', PATH_AUTOCLEANED, 'loadmode', 'all');
        d = double(EEG.data);

        % Prune for event-related analysis
        tidx = dsearchn(EEG.times', [pruned_segs(1), pruned_segs(2)]');
        d = d(:, tidx(1) : tidx(2), :);

        % Filter data for spatial filter determination
		fwhmf = 3;
		cfrq = 5.5;
		tmp = reshape(d, EEG.nbchan, []); % Trial after trial
		hz = linspace(0, EEG.srate, size(tmp, 2)); % Define frequency vector
		s  = fwhmf * (2 * pi - 1) / (4 * pi); % normalized width
		x  = hz - cfrq; % shifted frequencies
		fx = exp(-.5 * ( x / s) .^ 2); % Create gaussian in frq-domain
		fx = fx ./ max(fx); % normalize gaussian
		tmp = 2 * real(ifft(bsxfun(@times, fft(tmp, [], 2), fx), [], 2)); % Actually filter the data 
		d_filt = reshape(tmp, [EEG.nbchan, size(d, 2), size(d, 3)]); % Back to 3d

		% Find indices of time points for S & R selection
		tidx_S = dsearchn(tf_times', [300, 700]');
		tidx_R = dsearchn(tf_times', [300, 700]');

		% Calculate covariance matrices
		d_S = reshape(d_filt(:, tidx_S(1) : tidx_S(2), :), EEG.nbchan, []);
		d_S = bsxfun(@minus, d_S, mean(d_S, 2)); % Mean center data
		S = d_S * d_S' / diff(tidx_S(1 : 2));

		d_R = reshape(d(:, tidx_R(1) : tidx_R(2), :), EEG.nbchan, []);
		d_R = bsxfun(@minus, d_R, mean(d_R, 2)); % Mean center data
		R = d_R * d_R' / diff(tidx_R(1 : 2));

		% GED and sort eigenvalues and eigenvectors
		[evecs, evals] = eig(S, R);
		[evals, srtidx] = sort(diag(evals), 'descend'); % Sort eigenvalues
		evecs = evecs(:, srtidx); % Sort eigenvectors

		% Normalize eigenvectors
        evecs = bsxfun(@rdivide, evecs, sqrt(sum(evecs .^ 2, 1)));
        
        % Threshold eigenvalues
		thresh_eval_sd = 1; % In sd
		thresh_eigenvalue = median(evals) + std(evals) * thresh_eval_sd;
		suprathresh_eval_idx = find(evals > thresh_eigenvalue);  

        % Iterate components
        ged_maps = zeros(EEG.nbchan, EEG.nbchan);
        ged_tsrs = zeros(EEG.nbchan, length(tf_times), EEG.trials);
		for cmp = 1 : EEG.nbchan

			% Compute maps and flip sign if necessary
			ged_maps(cmp, :) = evecs(:, cmp)' * S;
			[~, idx] = max(abs(ged_maps(cmp, :)));
			ged_maps(cmp, :) = ged_maps(cmp, :) * sign(ged_maps(cmp, idx));

			% Compute time series data for components, i.e. apply spatial filters to unfiltered data
			tpts = evecs(:, cmp)' * reshape(d, EEG.nbchan, []);
			ged_tsrs(cmp, :, :) = reshape(tpts, [length(tf_times), size(d, 3)]);   

        end

        % Identify blink components [Fp1, Fp2, AFz -> 1, 2, 3] > [FC1, FCz, FC2 -> 13, 14, 15]
		blink_cmp = zeros(1, EEG.nbchan);
		for cmp = 1 : EEG.nbchan
			if mean(ged_maps(cmp, [1, 2, 3])) > mean(ged_maps(cmp, [13, 14, 15]))
				blink_cmp(cmp) = 1;
			end
		end

		% Find highest similarity in supra-threshold non-blink cmps
		cmpsim = 0;
		chosen_cmp = 0;
		for e = 1 : EEG.nbchan
			if ismember(e, suprathresh_eval_idx) & ~blink_cmp(e)
				tmp_cmp = ged_maps(e, :) / max(ged_maps(e, :)); % Normalize
				tmp = corrcoef(tmp_cmp, template_topo);
				if tmp(1, 2) * evals(e) > cmpsim
					cmpsim = tmp(1, 2) * evals(e);
					chosen_cmp = e;
				end
			end	
		end

		% Save filter topography
		figure('Visible', 'off'); clf;
		topoplot(ged_maps(chosen_cmp, :), EEG.chanlocs, 'electrodes', 'off', 'numcontour', 0)
        saveas(gcf, [PATH_VEUSZ, 'descriptives_spatial_filter/', 'filter_topo_' subject '.png']);
        
        % Get component signal
        cmp_sig = squeeze(ged_tsrs(chosen_cmp, :, :));

        % tf decomp of component
        convlen = size(cmp_sig, 1) * size(cmp_sig, 2) + size(cmw, 2) - 1;

        % cmw to freq domain and scale
        cmwX = zeros(length(tf_freqs), convlen);
        for f = 1 : length(tf_freqs)
            cmwX(f, :) = fft(cmw(f, :), convlen);
            cmwX(f, :) = cmwX(f, :) ./ max(cmwX(f, :));
        end

        % Get TF-power
        powcube = NaN(length(tf_freqs), size(cmp_sig, 1), size(cmp_sig, 2));
        tmp = fft(reshape(double(cmp_sig), 1, []), convlen);
        for f = 1 : length(tf_freqs)
            as = ifft(cmwX(f, :) .* tmp); 
            as = as(((size(cmw, 2) - 1) / 2) + 1 : end - ((size(cmw, 2) - 1) / 2));
            as = reshape(as, size(cmp_sig, 1), size(cmp_sig, 2));
            powcube(f, :, :) = abs(as) .^ 2;          
        end

        % Cut edges
        powcube = powcube(:, dsearchn(tf_times', -500) : dsearchn(tf_times', 1000), :);
        prune_time = tf_times(dsearchn(tf_times', -500) : dsearchn(tf_times', 1000));

        % Save component signal and single trial pows and pruned time vector
        save([PATH_COMPONENT_DATA subject '_signal_cmp'], 'cmp_sig');
        save([PATH_COMPONENT_DATA subject '_powcube_cmp'], 'powcube');
        save([PATH_COMPONENT_DATA 'tf_time_cmp'], 'prune_time');

        % Load subject meta
        meta = dlmread([PATH_TFDECOMP num2str(id) '_powcube_meta.csv']);

        % Regression design matrices
        meta_c = meta(meta(:, 7) == 1, :);
        meta_nc = meta(meta(:, 7) == 2, :);
        desmat = [ones(size(meta, 1), 1), meta(:, 5)];
        desmat(:, 2) = desmat(:, 2) / max(abs(desmat(:, 2))); % Scale
        desmat_c = [ones(size(meta_c, 1), 1), meta_c(:, 5)];
        desmat_c(:, 2) = desmat_c(:, 2) / max(abs(desmat_c(:, 2))); % Scale
        desmat_nc = [ones(size(meta_nc, 1), 1), meta_nc(:, 5)];
        desmat_nc(:, 2) = desmat_nc(:, 2) / max(abs(desmat_nc(:, 2))); % Scale

        % Z-Standardize trials
        zcube = zeros(size(powcube));

        % Apply single trial baseline
        blidx = dsearchn(prune_time', [-200, 0]');
        for t = 1 : size(powcube, 3)
            d_trial = squeeze(powcube(:, :, t)); % Get trial tfmat
            blvals = squeeze(mean(d_trial(:, blidx(1) : blidx(2)), 2)); % Get baseline
            blstd = std(d_trial(:, blidx(1) : blidx(2)), 0, 2);
            d_trial = bsxfun(@minus, d_trial, blvals);
            zcube(:, :, t) = bsxfun(@rdivide, d_trial, blstd);
        end

        % Set time windows as trial ranges for ersp calculation
        trialwins = {[1 : 200], [1241 : 1440], [2681 : 2880], [4121 : 4320]};

        % Save cmp ersps for timewins
        for w = 1 : numel(trialwins)
            ersp_cmp_c(w, :, :) = squeeze(ersp_cmp_c(w, :, :)) + squeeze(mean(zcube(:, :, ismember(meta(:, 5), trialwins{w}) & meta(:, 7) == 1), 3));
            ersp_cmp_nc(w, :, :) = squeeze(ersp_cmp_nc(w, :, :)) + squeeze(mean(zcube(:, :, ismember(meta(:, 5), trialwins{w}) & meta(:, 7) == 2), 3));
        end

        % Split data
        zcube_c = zcube(:, :, meta(:, 7) == 1);
        zcube_nc = zcube(:, :, meta(:, 7) == 2);

        % OLS fit all trials
        d = reshape(zcube, numel(prune_time) * numel(tf_freqs), size(zcube, 3))';
        tmp = (desmat' * desmat) \ desmat' * d;
        cmp_betas_tru = reshape(squeeze(tmp(2, :)), [numel(tf_freqs), numel(prune_time)]);
        save([PATH_COMPONENT_DATA subject '_cmp_betas_tru'], 'cmp_betas_tru');

        % Generate tot pseudo condition
        fakedesmat = desmat;
        fakedesmat(:, 2) = desmat(randperm(size(desmat, 1)), 2); % Permute tot column
        tmp = (fakedesmat' * fakedesmat) \ fakedesmat' * d;
        cmp_betas_fak = reshape(squeeze(tmp(2, :)), [numel(tf_freqs), numel(prune_time)]);
        save([PATH_COMPONENT_DATA subject '_cmp_betas_fak'], 'cmp_betas_fak');

        % OLS fit c trials
        d = reshape(zcube_c, numel(prune_time) * numel(tf_freqs), size(zcube_c, 3))';
        tmp = (desmat_c' * desmat_c) \ desmat_c' * d;
        cmp_betas_tru_c = reshape(squeeze(tmp(2, :)), [numel(tf_freqs), numel(prune_time)]);
        save([PATH_COMPONENT_DATA subject '_cmp_betas_tru_c'], 'cmp_betas_tru_c');

        % OLS fit nc trials
        d = reshape(zcube_nc, numel(prune_time) * numel(tf_freqs), size(zcube_nc, 3))';
        tmp = (desmat_nc' * desmat_nc) \ desmat_nc' * d;
        cmp_betas_tru_nc = reshape(squeeze(tmp(2, :)), [numel(tf_freqs), numel(prune_time)]);
        save([PATH_COMPONENT_DATA subject '_cmp_betas_tru_nc'], 'cmp_betas_tru_nc');

        % Save ersp of c and nc trials
        cmp_ersp_c = squeeze(mean(zcube_c, 3));
        cmp_ersp_nc = squeeze(mean(zcube_nc, 3));
        save([PATH_COMPONENT_DATA subject '_cmp_ersp_c'], 'cmp_ersp_c');
        save([PATH_COMPONENT_DATA subject '_cmp_ersp_nc'], 'cmp_ersp_nc');

    end % End subit

    % Average and save timewin cmp ersps
    for w = 1 : numel(trialwins)
        ersp_cmp_c(w, :, :) = squeeze(ersp_cmp_c(w, :, :)) / numel(subject_list);
        dlmwrite([PATH_VEUSZ, 'descriptives_spatial_filter/', 'ersp_cmp_c_win_' num2str(w) '.csv'], squeeze(ersp_cmp_c(w, :, :)));

        ersp_cmp_nc(w, :, :) = squeeze(ersp_cmp_nc(w, :, :)) / numel(subject_list);
        dlmwrite([PATH_VEUSZ, 'descriptives_spatial_filter/', 'ersp_cmp_nc_win_' num2str(w) '.csv'], squeeze(ersp_cmp_nc(w, :, :)));
    end

    % Load  data
    d_tru = zeros(length(subject_list), length(tf_freqs), length(prune_time));
    d_fak = zeros(length(subject_list), length(tf_freqs), length(prune_time));
    d_tru_c = zeros(length(subject_list), length(tf_freqs), length(prune_time));
    d_tru_nc = zeros(length(subject_list), length(tf_freqs), length(prune_time));
    d_ersp_c = zeros(length(subject_list), length(tf_freqs), length(prune_time));
    d_ersp_nc = zeros(length(subject_list), length(tf_freqs), length(prune_time));
    for s = 1 : length(subject_list)
        subject = subject_list{s};
        load([PATH_COMPONENT_DATA subject '_cmp_betas_tru']);
        load([PATH_COMPONENT_DATA subject '_cmp_betas_fak']);
        d_tru(s, :, :) = cmp_betas_tru;
        d_fak(s, :, :) = cmp_betas_fak;
        load([PATH_COMPONENT_DATA subject '_cmp_betas_tru_c']);
        load([PATH_COMPONENT_DATA subject '_cmp_betas_tru_nc']);
        d_tru_nc(s, :, :) = cmp_betas_tru_c;
        d_fak_nc(s, :, :) = cmp_betas_tru_nc;
        load([PATH_COMPONENT_DATA subject '_cmp_ersp_c']);
        load([PATH_COMPONENT_DATA subject '_cmp_ersp_nc']);
        d_ersp_c(s, :, :) = cmp_ersp_c;
        d_ersp_nc(s, :, :) = cmp_ersp_nc;
    end

    % Save average data
    dlmwrite([PATH_VEUSZ, 'inference_spatial_filter/', 'tot_beta.csv'], squeeze(mean(d_tru, 1)));
    dlmwrite([PATH_VEUSZ, 'inference_spatial_filter/', 'tot_beta_c.csv'], squeeze(mean(d_tru_c, 1)));
    dlmwrite([PATH_VEUSZ, 'inference_spatial_filter/', 'tot_beta_nc.csv'], squeeze(mean(d_tru_nc, 1)));
    dlmwrite([PATH_VEUSZ, 'inference_spatial_filter/', 'tt_ersp_c.csv'], squeeze(mean(d_ersp_c, 1)));
    dlmwrite([PATH_VEUSZ, 'inference_spatial_filter/', 'tt_ersp_nc.csv'], squeeze(mean(d_ersp_nc, 1)));
    dlmwrite([PATH_VEUSZ, 'inference_spatial_filter/', 'tt_ersp_diff.csv'], squeeze(mean(d_ersp_nc, 1)) - squeeze(mean(d_ersp_c, 1)));

    % Permtest params
    pval_voxel   = 0.01;
    pval_cluster = 0.025;
    n_perms      = 1000;

    % Test against zero c data
    d1 = d_tru;
    d2 = d_fak;
    n_freq = length(tf_freqs);
    n_time = length(prune_time);
    permuted_t = zeros(n_perms, n_freq, n_time);
    max_clust = zeros(n_perms, 2);
    desmat = [zeros(length(subject_list), 1), ones(length(subject_list), 1)];
    for p = 1 : n_perms
        fprintf('%i\n', p);
        toflip = find(round(rand(length(subject_list), 1)));
        d1_perm = d1;
        d1_perm(toflip, :, :) = d2(toflip, :, :);
        d2_perm = d2;
        d2_perm(toflip, :, :) = d1(toflip, :, :);
        tnum = squeeze(mean(d1_perm - d2_perm, 1));
        tdenum = squeeze(std(d1_perm - d2_perm, 0, 1)) / sqrt(length(subject_list));
        fake_t = tnum ./ tdenum;
        permuted_t(p, :, :) = fake_t;
        fake_t(abs(fake_t) < tinv(1 - pval_voxel, length(subject_list) - 1)) = 0;
        clusts = bwconncomp(fake_t);
        sum_t = [];
        for clu = 1 : numel(clusts.PixelIdxList)
            cidx = clusts.PixelIdxList{clu};
            sum_t(end + 1) = sum(fake_t(cidx));
        end
        max_clust(p, 1) = min([0, sum_t]);
        max_clust(p, 2) = max([0, sum_t]);      
    end
    tnum = squeeze(mean(d1 - d2, 1));
    tdenum = squeeze(std(d1 - d2, 0, 1)) / sqrt(length(subject_list));
    tmat = tnum ./ tdenum;
    tvals = tmat;
    tmat(abs(tmat) < tinv(1 - pval_voxel, length(subject_list) - 1)) = 0;
    threshtvals = tmat;
    clusts = bwconncomp(tmat);
    sum_t = [];
    for clu = 1 : numel(clusts.PixelIdxList)
        cidx = clusts.PixelIdxList{clu};
        sum_t(end + 1) = sum(tmat(cidx));
    end
    clust_thresh_lower = prctile(max_clust(:, 1), pval_cluster * 100);
    clust_thresh_upper = prctile(max_clust(:, 2), 100 - pval_cluster * 100);
    clust2remove = find(sum_t > clust_thresh_lower & sum_t < clust_thresh_upper);
    for clu = 1 : length(clust2remove)
        tmat(clusts.PixelIdxList{clust2remove(clu)}) = 0;
    end
    contourres = logical(tmat);

    % Save contour of effect
    dlmwrite([PATH_VEUSZ, 'inference_spatial_filter/', 'tot_contour.csv'], contourres);

    % Calculate and save effect sizes
    petasq = (tvals .^ 2) ./ ((tvals .^ 2) + (numel(subject_list) - 1));
    adj_petasq = petasq - (1 - petasq) .* (1 / (numel(subject_list) - 1));
    dlmwrite([PATH_VEUSZ, 'inference_spatial_filter/', 'tot_effect_sizes.csv'], adj_petasq);

    % Test tt
    d1 = d_ersp_c;
    d2 = d_ersp_nc;
    n_freq = length(tf_freqs);
    n_time = length(prune_time);
    permuted_t = zeros(n_perms, n_freq, n_time);
    max_clust = zeros(n_perms, 2);
    desmat = [zeros(length(subject_list), 1), ones(length(subject_list), 1)];
    for p = 1 : n_perms
        fprintf('%i\n', p);
        toflip = find(round(rand(length(subject_list), 1)));
        d1_perm = d1;
        d1_perm(toflip, :, :) = d2(toflip, :, :);
        d2_perm = d2;
        d2_perm(toflip, :, :) = d1(toflip, :, :);
        tnum = squeeze(mean(d1_perm - d2_perm, 1));
        tdenum = squeeze(std(d1_perm - d2_perm, 0, 1)) / sqrt(length(subject_list));
        fake_t = tnum ./ tdenum;
        permuted_t(p, :, :) = fake_t;
        fake_t(abs(fake_t) < tinv(1 - pval_voxel, length(subject_list) - 1)) = 0;
        clusts = bwconncomp(fake_t);
        sum_t = [];
        for clu = 1 : numel(clusts.PixelIdxList)
            cidx = clusts.PixelIdxList{clu};
            sum_t(end + 1) = sum(fake_t(cidx));
        end
        max_clust(p, 1) = min([0, sum_t]);
        max_clust(p, 2) = max([0, sum_t]);      
    end
    tnum = squeeze(mean(d1 - d2, 1));
    tdenum = squeeze(std(d1 - d2, 0, 1)) / sqrt(length(subject_list));
    tmat = tnum ./ tdenum;
    tvals = tmat;
    tmat(abs(tmat) < tinv(1 - pval_voxel, length(subject_list) - 1)) = 0;
    threshtvals = tmat;
    clusts = bwconncomp(tmat);
    sum_t = [];
    for clu = 1 : numel(clusts.PixelIdxList)
        cidx = clusts.PixelIdxList{clu};
        sum_t(end + 1) = sum(tmat(cidx));
    end
    clust_thresh_lower = prctile(max_clust(:, 1), pval_cluster * 100);
    clust_thresh_upper = prctile(max_clust(:, 2), 100 - pval_cluster * 100);
    clust2remove = find(sum_t > clust_thresh_lower & sum_t < clust_thresh_upper);
    for clu = 1 : length(clust2remove)
        tmat(clusts.PixelIdxList{clust2remove(clu)}) = 0;
    end
    contourres = logical(tmat);

    % Save contour of effect
    dlmwrite([PATH_VEUSZ, 'inference_spatial_filter/', 'tt_contour.csv'], contourres);

    % Calculate and save effect sizes
    petasq = (tvals .^ 2) ./ ((tvals .^ 2) + (numel(subject_list) - 1));
    adj_petasq = petasq - (1 - petasq) .* (1 / (numel(subject_list) - 1));
    dlmwrite([PATH_VEUSZ, 'inference_spatial_filter/', 'tt_effect_sizes.csv'], adj_petasq);

end % End part8

% ======================= PART9: SUBJECTIVE RATINGS =====================================

if ismember('part9', to_execute)

    % Read data
    sr = xlsread([PATH_RAW_DATA 'subjective_ratings.xlsx']);

    % Option to leave out first assessments after the breaks and at start
    if false
        sr(ismember(sr(:, 2), [1, 5, 9]), :) = [];
        sr(sr(:, 2) == 2, 2) = 1;
        sr(sr(:, 2) == 3, 2) = 2;  
        sr(sr(:, 2) == 4, 2) = 3;  
        sr(sr(:, 2) == 6, 2) = 4;  
        sr(sr(:, 2) == 7, 2) = 5;  
        sr(sr(:, 2) == 8, 2) = 6;  
        sr(sr(:, 2) == 10, 2) = 7;  
        sr(sr(:, 2) == 11, 2) = 8;  
        sr(sr(:, 2) == 12, 2) = 9;  
    end

    % Flip polarity of motivation
    sr(:, 4) = (sr(:, 4) * -1) + 10;

    % Create tables
    varnames = {'id', 'tot' , 'fat', 'mot'};
    tbl = table(sr(:, 1), sr(:, 2), sr(:, 3), sr(:, 4), 'VariableNames', varnames);
    
    % Cast vars
    tbl.id = nominal(tbl.id);

    % Compute LMEs
    lme_fat = fitlme(tbl, 'fat ~ tot + (1|id)');
    lme_mot = fitlme(tbl, 'mot ~ tot + (1|id)');

    % Calculate averages
    out = [];
    ts = unique(sr(:, 2));
    for t = 1 : length(ts)
        tmp = sr(sr(:, 2) == ts(t), 3);
        fat_m = mean(tmp);
        fat_s = std(tmp);
        tmp = sr(sr(:, 2) == ts(t), 4);
        mot_m = mean(tmp);
        mot_s = std(tmp);
        out(t, :) = [fat_m, fat_s, mot_m, mot_s];
    end

    dlmwrite([PATH_VEUSZ, 'descriptives_subjective/', 'ratings.csv'], out, 'delimiter', '\t');
    dlmwrite([PATH_VEUSZ, 'descriptives_subjective/', 'xax.csv'], [1 : length(ts)]);

end % End part9