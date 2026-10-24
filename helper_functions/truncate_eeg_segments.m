function eeg_segment_data = truncate_eeg_segments( ...
    matched_segments, eeg_input, aligned_eeg_channel_names, channel_overrides)
%TRUNCATE_EEG_SEGMENTS Extract ordered HD-EEG channels for each segment.
%   EEG_SEGMENT_DATA = TRUNCATE_EEG_SEGMENTS(MATCHED_SEGMENTS, EEG_INPUT,
%   ALIGNED_EEG_CHANNEL_NAMES) maps the requested aligned channel names to
%   HD-EEG channel labels and returns one channels-by-samples array per
%   matched segment.
%   Optional CHANNEL_OVERRIDES replaces selected default HD-EEG sources
%   using a two-column cell array, for example {'C3', 'LL5'}.
%
%   Aligned channel: Fp1 Fp2 F3  F4  C3  C4  O1   O2   M1  M2  VEOGL
%   HD-EEG channel:  L1  R1  LL2 RR2 LA2 RA2 LL11 RR11 LD6 RD6 VEOGL

if nargin < 4
    channel_overrides = {};
end
channel_mapping = get_hd_eeg_channel_mapping(channel_overrides);
requested_channel_names = aligned_eeg_channel_names;

[is_supported, mapping_index] = ismember( ...
    requested_channel_names, channel_mapping(:, 1));
if ~all(is_supported)
    unsupported_names = strjoin( ...
        requested_channel_names(~is_supported), ', ');
    error('truncate_eeg_segments:UnsupportedChannel', ...
        'Unsupported aligned EEG channel(s): %s.', unsupported_names)
end

requested_hd_eeg_names = channel_mapping(mapping_index, 2).';
eeg_channel_names = {eeg_input.EEG.chanlocs.labels};
[has_channel, eeg_channel_index] = ismember( ...
    requested_hd_eeg_names, eeg_channel_names);
if ~all(has_channel)
    missing_names = strjoin(requested_hd_eeg_names(~has_channel), ', ');
    error('truncate_eeg_segments:MissingChannel', ...
        'Required HD-EEG channel(s) not found: %s.', missing_names)
end

selected_eeg_data = eeg_input.EEG.data(eeg_channel_index, :);
eeg_segment_data = cell(height(matched_segments), 1);

for segment_i = 1:height(matched_segments)
    start_sample = matched_segments.eeg_start_sample(segment_i);
    end_sample = matched_segments.eeg_end_sample(segment_i);
    eeg_segment_data{segment_i} = ...
        selected_eeg_data(:, start_sample:end_sample);
end

end
