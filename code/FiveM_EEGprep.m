function dataout = FiveM_EEGprep(varargin)

% Single function to import mff folder from a subject and export it
% following BIDS. After the export, a derivative folder is created and
% data are processed and basic power is returned allowing quick QC.
% WARNING: 20 min EEG for 20 min is massive, it is therefore essential to
% increase memory allocation: Home --> Preferences --> General --> Java
% Heap Memory and adjust there as needed
%
% FORMAT dataout = FiveM_EEGprep(folderin,subjectID,BIDSfolder)
%
% INPUTS if empty, user is prompted
%        folderin is the mff folder e.g. /.../Caludia_20250224_190613.mff
%        subject ID is a char to encode subject pseudoID e.g. '01'
%        BIDSfolder is the root folder where 5M data are stored
%
% OUTPUTS creates BIDSfolder/sub-subjectID/eeg/ ... folder and data
%         creates BIDSfolder/derivatives/sub-subjectID/eeg/ ... folder and data
%         creates BIDSfolder/derivatives/sub-subjectID/eeg/QC.tvs file
%         return created folders and QC as a structure
%                dataout.raw the BIDS raw data folder
%                dataout.derivatives the BIDS derivatives data folder
%                dataout.QC the average power spectrum in bands averaged
% %                         over some channels
%
% The preprossecing workflow is as follow:
% 1  - upsampling at 2K Hz
% 2  - MRI artefact correction
% 3  - Ballistocardiogram correction
% 4  - channel selection (no face electrodes, no impedences > 50kohms)
% 4  - downsampling at 250Hz
% 5  - 50Hz filtering
% 6  - bad channel removal and 1Hz high pass
% 7  - ICA and ICA label on 5Hz data and back projection on 1Hz data
%     (trying to minimize slow drift over 20 min acquistion)
% 8  - ASR to remove remaining bad segments
% 10 - basic stats on where we have power, and get average power spectrum
%
% Requires EEGLAB and plugins (plugins are installed automatically if not
% there, but requires internet connection)
%
% Cyril Pernet, Neurobiololy Research Unit, Copenhagen, DK
% Ilaria Mazzonetto,
%
% March 2025

%% start eeglab and get missing dependencies if needed
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
if ~exist('pop_importbids','file')
    plugin_askinstall('EEG-BIDS',[],1);
end
if ~exist('mff_import','file')
    plugin_askinstall('MFFMatlabIO',[],1);
end
if ~exist('eegplugin_fmrib.m','file')
    plugin_askinstall('fMRIb',[],1);
end
if ~exist('pop_zapline_plus','file')
    plugin_askinstall('zapline-plus',[],1);
end
if ~exist('eegplugin_eegstats.m','file')
    plugin_askinstall('eegstats',[],1);
end

%% check inputs
if isempty(varargin)
    folderin  = uigetdir(pwd,'get source eeg folder');
    if folderin == 0
        warning('selection aborded'); return
    end
    BIDSfolder  = uigetdir(pwd,'get destination folder');
    if BIDSfolder == 0
        warning('selection aborded'); return
    end

    subjectID = inputdlg2('Enter subject ID','PseudoID',1);
    if isempty(subjectID) || isempty(subjectID{1})
        warning('no name entered'); return
    end
    newname  = ['sub-' cell2mat(subjectID)];

else
    if ~exist(varargin{1},"dir")
        error('input folder %s not found',varargin{1})
    else
        folderin = varargin{1};
    end

    subjectID = varargin{2};
    if isnumeric(subjectID)
        subjectID = char(subjectID);
    end
    newname  = ['sub-' subjectID];

    if ~exist(varargin{3},"dir")
        error('output folder %s not found',varargin{3})
    else
        BIDSfolder = varargin{3};
    end
end

%% EEG data conversion to BIDS
bids_exportmff(folderin, subjectID, ...
    fullfile(BIDSfolder,[newname filesep 'eeg']));
dataout.raw = fullfile(BIDSfolder,[newname filesep 'eeg']);

%% Preprocess

% MR filter
EEG = pop_loadset(fullfile(BIDSfolder,[newname filesep 'eeg' filesep ...
    filesep newname '_task-rest_eeg.set']));
EEG.data   = double(EEG.data);
lpf        =  60; 
L          =  20; % following their suggestion: It is recommended that you up-sample the data to bring the sampling frequency to about 20 kHz)
window     = 10;% we tested different values, and this seemed to provide the best cleaning in terms of reducing the gradient artifact peak and its harmonics  );
Trigs      = round([EEG.event(find(strcmp({EEG.event(:).code},'TREV'))).latency]);
anc_chk    = 1;  
opt1       = 0.03; % no significant impact when adjusting this parameter
ECGchannel = find(arrayfun(@(x) strcmp(x.labels,'ECG'),EEG.chanlocs));
opt3       = 'auto';
EEG        = fmrib_fastr(EEG,lpf,L,window,Trigs,0,anc_chk,0,[],[],opt1,ECGchannel,opt3);
EEG        = eeg_checkset(EEG);

% Cardiac artefacts
EEG        = pop_fmrib_qrsdetect(EEG,ECGchannel,'qrs','no');
EEG        = pop_fmrib_pas(EEG,'qrs','obs',4); % OBS method with npc = 4.
EEG        = eeg_checkset(EEG);

% export to derivatives
destination  = fullfile(BIDSfolder,['derivatives' filesep newname filesep 'eeg']);
EEG.setname  = EEG.filename;
EEG.filename = [newname '_desc-ArtefactsCorr_task-rest'];
EEG.filepath = destination;
mkdir(destination)
pop_saveset( EEG, 'filename',[newname '_desc-ArtefactsCorr_task-rest'],'filepath',destination);
clear EEG

%% 'standard' preproc
% similar to https://pubmed.ncbi.nlm.nih.gov/33519362/ 
EEG = resting_eeg_preproc(fullfile(destination,[newname '_desc-ArtefactsCorr_task-rest']));

%% extract power in pre-defined freq bands for QC 
% freq as per COBIDAS MEEG https://www.nature.com/articles/s41593-020-00709-0
disp('extracting power ...')
EEG = pop_eegstats(EEG, 'thetarange',[4 7.9] ,'alpharange',[8 12.9] , ...
            'otherranges',[13 29.5; 30 49.5],'averagepower','off','channels','','winsize',2, ...
            'overlap',1,'iaf','on','iafminchan',1,'alphaasymmetry','off', ...
            'csvfile',fullfile(destination,[newname '_desc-power_eeg.csv']));
pop_saveset(EEG, 'filename',[newname '_desc-ArtefactsCorr-preproc_task-rest'],'filepath',destination);
dataout.derivatives = destination;

% re-export in a main table for easier QC
% from EEG.etc.eegstats average some channels of interests and make QC
% (TBD)
% dataout.QC.theta
% dataout.QC.alpha
% dataout.QC.beta
% dataout.QC.gamma

disp('data preprocessing and QC done')
end

%% ------------------------------------------------------------------------
%%                          ROUTINES
%% ------------------------------------------------------------------------
function bids_exportmff(eegfolder,subID, destination)

% routine to read an eeg mmf folder and export to bids
% assumes resting state data - ie no event code from stimuli

if ~exist(destination,"dir")
    mkdir(destination)
end

% metadata
if iscell(subID)
    subID = cell2mat(subID);
end
task_name          = 'rest';
export_task_events = 'on';

% read and write the data
EEG = pop_mffimport({eegfolder},{'code','mffkeys','sourcedevice'});
ECGchannel = find(arrayfun(@(x) strcmpi(x.labels,'ECG'),EEG.chanlocs));
if ~isempty(ECGchannel)
    EEG = pop_chanedit(EEG, 'changefield',{ECGchannel,'type','ECG'},'rplurchanloc',1); %bugs if no ECG channel is present
end
pop_saveset( EEG, 'filename',['sub-' subID '_task-' task_name],'filepath',fileparts(eegfolder));

% export to BIDS locally
file2export.file                = {fullfile(fileparts(eegfolder),['sub-' subID '_task-' task_name '.set'])};
file2export.task                = {task_name};
subinfo{1,1}                    = 'participant_id';
subinfo{2,1}                    = subID;
PL.TaskName                     = {task_name};
PL.PowerLineFrequency           = 50;
eventsinfo.onset.LongName       = 'Event onset';
eventsinfo.onset.Description    = 'Onset (in seconds) of the event measured from the beginning of the acquisition of the first volume in the corresponding task imaging data file';
eventsinfo.onset.Units          = 'second';
eventsinfo.duration.LongName    = 'Event duration';
eventsinfo.duration.Description = 'Duration of the event (measured from onset) in seconds';
eventsinfo.duration.Units       = 'second';
eventsinfo.sample.LongName      = 'Sample';
eventsinfo.sample.Description   = 'Onset of the event according to the sampling scheme of the recorded modality';
eventsinfo.value.LongName       = 'Event marker';
eventsinfo.value.Description    = 'Marker value associated with the event';
bids_export(file2export,'interactive','off','targetdir', ...
    fullfile(fileparts(eegfolder),'tmp_export'), ...
    'eInfoDesc', eventsinfo,'tInfo', PL, 'pInfo', subinfo);

% edit the electrodes.tsv with impedences
tmp_folder = fullfile(fileparts(eegfolder), ...
    ['tmp_export' filesep 'sub-' subinfo{2,1} filesep 'eeg']);
insides = readstruct([eegfolder filesep 'info1.xml']);
electrodes = readtable(fullfile(tmp_folder, ['sub-' subID '_electrodes.tsv']),'FileType','text');
for i = 1:length(insides.calibrations.calibration)
    imp     = arrayfun(@(x) x.Text,insides.calibrations.calibration(i).channels.ch)';
    impName = sprintf('impedendence_measurement_%g',i);
    if length(imp)<size(electrodes,1)
        imp = [imp;NaN(size(electrodes,1)-size(imp,1),1)];
    end
    electrodes = addvars(electrodes,imp,'NewVariableNames',impName);
end
writetable(electrodes, fullfile(tmp_folder, ['sub-' subID '_electrodes.tsv']), ...
    'FileType', 'text', 'Delimiter', '\t',"WriteVariableNames",true);

% move and clean-up
copyfile(tmp_folder,destination);
if ~exist(fullfile(fileparts(fileparts(destination)), 'dataset_description.json'),'file')
        copyfile(fullfile(fileparts(eegfolder), ['tmp_export' filesep 'dataset_description.json']),...
        fullfile(fileparts(fileparts(destination)), 'dataset_description.json'));
end
if ~exist(fullfile(fileparts(fileparts(destination)), 'task-rest_events.json'),'file')
        copyfile(fullfile(fileparts(eegfolder), ['tmp_export' filesep 'task-rest_events.json']),...
        fullfile(fileparts(fileparts(destination)), 'task-rest_events.json'));
end
% bids root info, update participants.tsv
if ~exist(fullfile(fileparts(fileparts(destination)), 'participants.tsv'),'file')
        copyfile(fullfile(fileparts(eegfolder), ['tmp_export' filesep 'participants.tsv']),...
        fullfile(fileparts(fileparts(destination)), 'participants.tsv'));
else
    participants = readtable(fullfile(fileparts(fileparts(destination)), 'participants.tsv'),'FileType','text');
    participants.participant_id{end+1,1} = ['sub-' subID];
    writetable(participants, fullfile(fileparts(fileparts(destination)), 'participants.tsv'), 'FileType', 'text', ...
        'Delimiter', '\t',"WriteVariableNames",true);
end

rmdir(fullfile(fileparts(eegfolder),'tmp_export'),'s');
if exist(file2export.file{1}, 'file')
    [root, name, ~] = fileparts(file2export.file{1});
    delete(fullfile(root, [name '.fdt']))
    delete(file2export.file{1});
end
end

function EEG = resting_eeg_preproc(fullfilename)
% preprocessing - note this is a little more agressive than usual
% the criterion use in bad channel, bad segments, etc are all higher
% than what would do with data collected outside the scanner
%
% fullfilename = fullfile(destination,[newname '_desc-ArtefactsCorr_task-rest'])

% some hard coded values (to change as nedded)
% see clean_rawdata, pop_runica
IMPTH              = 1.05; % threshold for impedences > 50kohm
ChannelCriterion   = 0.9; % interchannel correlation to remove bad channels
LineNoiseCriterion = 3 ; % devition from avg
ICAalg             = 'runica'; % which ICA method to use

% compute
if ~strcmp(fullfilename(end-3:end),'.set')
    fullfilename = [fullfilename '.set'];
end
EEG               = pop_loadset(fullfilename);
expected_chanlocs = EEG.chanlocs;

% remove unwanted channels
channellist = [];
ECGchannel  = find(arrayfun(@(x) strcmp(x.labels,'ECG'),EEG.chanlocs));
if ~isempty(ECGchannel)
    expected_chanlocs(ECGchannel) = [];
    channellist = ECGchannel;
end

% face channels, the conductivity profile is weird because of sinuses
channellist = [channellist, 237:256];

% to high impedences
root_folder = extractBefore(EEG.filepath,'derivatives');
subjectID   = cell2mat(extractBetween(EEG.filename,'sub-','_desc'));
impedences  = readtable(fullfile(root_folder,...
    ['sub-' subjectID filesep 'eeg' filesep 'sub-' subjectID '_electrodes.tsv']),...
    "FileType","text");
ECGchannel = find(arrayfun(@(x) strcmp(x,'ECG'),impedences.name));
if ~isempty(ECGchannel)
    impedences(ECGchannel,:) = [];
end
absoluteimp = sum(impedences(:,5:end)>IMPTH,2);
channellist = [channellist,find(table2array(absoluteimp))'];
EEG         = pop_select(EEG,'rmchannel',channellist);

% downsample
if EEG.srate ~= 250
    EEG = pop_resample(EEG, 250);
end

% 50Hz removal
EEG = pop_zapline_plus(EEG,'noisefreqs','line',...
    'coarseFreqDetectPowerDiff',4,'chunkLength',0,...
    'adaptiveNremove',1,'fixedNremove',1,'plotResults',0);

% remove bad channels 
EEG = pop_clean_rawdata( EEG,'FlatlineCriterion',5,'ChannelCriterion',ChannelCriterion, ...
    'LineNoiseCriterion',LineNoiseCriterion,'Highpass',[0.25 0.75], 'BurstCriterion','off',...
    'WindowCriterion','off','BurstRejection','off', 'Distance','Euclidian', ...
    'WindowCriterionTolerances','off' );

% interpolate missing channels and reference
[~,idx] = setdiff({expected_chanlocs.labels},{EEG.chanlocs.labels});
if ~isempty(idx)
    EEG = pop_interp(EEG, expected_chanlocs(idx), 'spherical');
end
EEG = pop_reref(EEG,[],'interpchan','off');

% ICA cleaning the usual way:
% EEG = pop_runica(EEG, 'icatype',ICAalg,'concatcond','on','options',{'pca',EEG.nbchan-1});
% EEG = pop_iclabel(EEG, 'default');
% EEG = pop_icflag(EEG,[NaN NaN;0.8 1;0.8 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN]);
% EEG = pop_subcomp(EEG,[],0);

% ICA cleaning but on 2Hz filtered data
% assuming slow drift with saline drying over 20 min acquisition
% this can bias the ICA, so compute on 2Hz filtered 
tmpeeg = pop_eegfiltnew(EEG,2,0);
tmpeeg = pop_runica(tmpeeg, 'icatype',ICAalg,'concatcond','on','options',{'pca',EEG.nbchan-1});
% project solution to the 0.5Hz filtered data
EEG.icasphere   = [];
EEG.icasphere   = tmpeeg.icasphere;
EEG.icaweights  = [];
EEG.icaweights  = tmpeeg.icaweights;
EEG.icachansind = [];
EEG.icachansind = tmpeeg.icachansind;
EEG.icaact      = [];
EEG.icawinv     = [];
EEG             = eeg_checkset(EEG); % re-compute EEG.icawinv
% clean data with ICs using IClabels
EEG = pop_iclabel(EEG, 'default');
EEG = pop_icflag(EEG,[NaN NaN;0.8 1;0.8 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN]);
EEG = pop_subcomp(EEG,[],0);
clear tmpeeg

% last remove remaining bad segments with ASR defaults
EEG = pop_clean_rawdata(EEG,'FlatlineCriterion','off','ChannelCriterion','off',...
    'LineNoiseCriterion','off','Highpass','off','BurstCriterion',20,...
    'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian',...
    'WindowCriterionTolerances',[-Inf 7] );
end
