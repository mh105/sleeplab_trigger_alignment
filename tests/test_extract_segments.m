function tests = test_extract_segments
tests = functiontests(localfunctions);
end

function setupOnce(~)
project_path = fileparts(fileparts(mfilename('fullpath')));
addpath(project_path)
addpath(fullfile(project_path, 'helper_functions'))
end

function testAllMatchedChunksFormOneSegment(testCase)
[eeg, edf] = build_pair(complete_keep, complete_keep);
[result, eeg, edf] = match_triggers(eeg, edf);

report = extract_segments(result, eeg, edf);

verifyEqual(testCase, height(report), 1)
verifyEqual(testCase, ...
    string(report.Properties.VariableNames(1:4)), ...
    ["segment_index", "eeg_duration_sec", "edf_duration_sec", ...
    "edf_minus_eeg_sec"])
verifyFalse(testCase, ismember( ...
    "estimated_edf_trigger_Fs", string(report.Properties.VariableNames)))
verifyEqual(testCase, ...
    report{1, {'start_eeg_chunk', 'start_edf_chunk', 'start_interval', ...
    'end_eeg_chunk', 'end_edf_chunk', 'end_interval'}}, ...
    [1 1 1 5 5 49])
verifyEqual(testCase, report.eeg_start_sample, ...
    double(eeg.event_table.latency(1)))
verifyEqual(testCase, report.eeg_end_sample, ...
    result.eeg_terminal_anchor_latency)
verifyEqual(testCase, report.edf_start_sample, ...
    double(edf.event_table.latency(1)))
verifyEqual(testCase, report.edf_end_sample, ...
    result.edf_terminal_anchor_latency)
end

function testPartialStartPreservesNumberingAndContinuesIntoFullChunk(testCase)
[eeg, edf] = build_pair(complete_keep, partial_start_keep);
edf.validation_sequence_count = 0;
edf.validation_sequences = table( ...
    zeros(0, 1), zeros(0, 1), ...
    'VariableNames', {'start_latency', 'end_latency'});
edf.pre_authoritative_validation_event_count = 0;
[result, eeg, edf] = match_triggers(eeg, edf);

verifyEqual(testCase, result.startup_mode, "partial_cycle_late_start")
verifyEqual(testCase, result.comparison_history.outcome(1:2), ...
    ["partial_start_match"; "match"])
verifyEqual(testCase, result.comparison_history{1:2, 1:2}, ...
    [1 1; 2 2])
verifyEqual(testCase, result.initial_canonical_interval_number, 25)
verifyEqual(testCase, result.eeg_initial_anchor_event_index, 25)
verifyEqual(testCase, result.edf_initial_anchor_event_index, 1)

report = extract_segments(result, eeg, edf);

verifyEqual(testCase, height(report), 1)
verifyEqual(testCase, ...
    report{1, {'start_eeg_chunk', 'start_edf_chunk', 'start_interval', ...
    'end_eeg_chunk', 'end_edf_chunk', 'end_interval'}}, ...
    [1 1 25 5 5 49])
verifyEqual(testCase, report.eeg_start_sample, ...
    double(eeg.event_table.latency(25)))
verifyEqual(testCase, report.edf_start_sample, ...
    double(edf.event_table.latency(1)))
verifyEqual(testCase, report.eeg_end_sample, ...
    result.eeg_terminal_anchor_latency)
verifyEqual(testCase, report.edf_end_sample, ...
    result.edf_terminal_anchor_latency)
end

function testShortLossIsSplitAtCanonicalPrefixAndSuffix(testCase)
for affected_system = ["eeg", "edf"]
    if affected_system == "eeg"
        [eeg, edf] = build_pair(short_loss_keep, complete_keep);
    else
        [eeg, edf] = build_pair(complete_keep, short_loss_keep);
    end
    [result, eeg, edf] = match_triggers(eeg, edf);

    report = extract_segments(result, eeg, edf);

    verifyEqual(testCase, height(report), 2)
    verifyEqual(testCase, ...
        report{:, {'start_eeg_chunk', 'start_edf_chunk', ...
        'start_interval', 'end_eeg_chunk', 'end_edf_chunk', ...
        'end_interval'}}, ...
        [1 1 1 2 2 10; 2 2 14 5 5 49])

    loss = result.missing_periods(1, :);
    if affected_system == "eeg"
        verifyEqual(testCase, report.eeg_end_sample(1), ...
            loss.start_anchor_latency)
        verifyEqual(testCase, report.eeg_start_sample(2), ...
            loss.end_anchor_latency)
        verifyEqual(testCase, report.edf_end_sample(1), ...
            double(edf.event_table.latency(61)))
        verifyEqual(testCase, report.edf_start_sample(2), ...
            double(edf.event_table.latency(64)))
    else
        verifyEqual(testCase, report.edf_end_sample(1), ...
            loss.start_anchor_latency)
        verifyEqual(testCase, report.edf_start_sample(2), ...
            loss.end_anchor_latency)
        verifyEqual(testCase, report.eeg_end_sample(1), ...
            double(eeg.event_table.latency(61)))
        verifyEqual(testCase, report.eeg_start_sample(2), ...
            double(eeg.event_table.latency(64)))
    end
end
end

function testPenultimateLossSuffixContinuesThroughTerminal(testCase)
clean_keep = complete_keep;
clean_keep{5} = (1:30)';
loss_keep = clean_keep;
loss_keep{4} = setdiff((1:50)', (12:13)', 'stable');

for affected_system = ["eeg", "edf"]
    if affected_system == "eeg"
        [eeg, edf] = build_pair(loss_keep, clean_keep);
    else
        [eeg, edf] = build_pair(clean_keep, loss_keep);
    end
    [result, eeg, edf] = match_triggers(eeg, edf);

    report = extract_segments(result, eeg, edf);

    verifyEqual(testCase, height(report), 2)
    verifyEqual(testCase, ...
        report{:, {'start_eeg_chunk', 'start_edf_chunk', ...
        'start_interval', 'end_eeg_chunk', 'end_edf_chunk', ...
        'end_interval'}}, ...
        [1 1 1 4 4 10; 4 4 14 5 5 29])
    verifyEqual(testCase, report.eeg_end_sample(2), ...
        result.eeg_terminal_anchor_latency)
    verifyEqual(testCase, report.edf_end_sample(2), ...
        result.edf_terminal_anchor_latency)

    loss = result.missing_periods(1, :);
    if affected_system == "eeg"
        verifyEqual(testCase, report.eeg_end_sample(1), ...
            loss.start_anchor_latency)
        verifyEqual(testCase, report.eeg_start_sample(2), ...
            loss.end_anchor_latency)
    else
        verifyEqual(testCase, report.edf_end_sample(1), ...
            loss.start_anchor_latency)
        verifyEqual(testCase, report.edf_start_sample(2), ...
            loss.end_anchor_latency)
    end
end
end

function testMarkedExactCountShortenedIntervalSplitsAtAnchors(testCase)
[eeg, edf] = build_pair(complete_keep, complete_keep);
eeg_start_anchor = 70;
eeg_end_anchor = eeg_start_anchor + 1;
eeg.event_table.latency(eeg_end_anchor:end) = ...
    eeg.event_table.latency(eeg_end_anchor:end) - eeg.Fs;
eeg = attach_amplifier_interruptions( ...
    eeg, eeg_start_anchor, eeg_end_anchor, "paired");
[result, eeg, edf] = match_triggers(eeg, edf);

report = extract_segments(result, eeg, edf);

verifyEqual(testCase, height(report), 2)
verifyEqual(testCase, report.eeg_end_sample(1), ...
    double(eeg.event_table.latency(eeg_start_anchor)))
verifyEqual(testCase, report.eeg_start_sample(2), ...
    double(eeg.event_table.latency(eeg_end_anchor)))
verifyEqual(testCase, report.edf_end_sample(1), ...
    result.amplifier_interruptions.edf_start_anchor_latency)
verifyEqual(testCase, report.edf_start_sample(2), ...
    result.amplifier_interruptions.edf_end_anchor_latency)
end

function testMarkedMultiIntervalPerfectTimingStillSplits(testCase)
[eeg, edf] = build_pair(complete_keep, complete_keep);
eeg_start_anchor = 70;
eeg_end_anchor = eeg_start_anchor + 3;
eeg = attach_amplifier_interruptions( ...
    eeg, eeg_start_anchor, eeg_end_anchor, "paired");
[result, eeg, edf] = match_triggers(eeg, edf);

report = extract_segments(result, eeg, edf);

verifyEqual(testCase, height(report), 2)
verifyEqual(testCase, report.eeg_end_sample(1), ...
    double(eeg.event_table.latency(eeg_start_anchor)))
verifyEqual(testCase, report.eeg_start_sample(2), ...
    double(eeg.event_table.latency(eeg_end_anchor)))
end

function testMarkedCrossBoundaryBracketSplitsAtExplicitAnchors(testCase)
[eeg, edf] = build_pair(complete_keep, complete_keep);
eeg_start_anchor = 99;
eeg_end_anchor = 102;
eeg = attach_amplifier_interruptions( ...
    eeg, eeg_start_anchor, eeg_end_anchor, "paired");
[result, eeg, edf] = match_triggers(eeg, edf);

report = extract_segments(result, eeg, edf);

verifyEqual(testCase, height(report), 2)
verifyEqual(testCase, report.eeg_end_sample(1), ...
    double(eeg.event_table.latency(eeg_start_anchor)))
verifyEqual(testCase, report.eeg_start_sample(2), ...
    double(eeg.event_table.latency(eeg_end_anchor)))
end

function testMarkedSingleIntervalWithinToleranceStaysContinuous(testCase)
[eeg, edf] = build_pair(complete_keep, complete_keep);
eeg_start_anchor = 70;
eeg_end_anchor = eeg_start_anchor + 1;
tolerance_sec = 2 / min([eeg.Fs, edf.trigger_Fs]);
within_tolerance_samples = floor(tolerance_sec * eeg.Fs);
eeg.event_table.latency(eeg_end_anchor:end) = ...
    eeg.event_table.latency(eeg_end_anchor:end) - ...
    within_tolerance_samples;
eeg = attach_amplifier_interruptions( ...
    eeg, eeg_start_anchor, eeg_end_anchor, "paired");
[result, eeg, edf] = match_triggers(eeg, edf);

report = extract_segments(result, eeg, edf);

verifyEqual(testCase, height(report), 1)
verifyEqual(testCase, ...
    result.amplifier_interruptions.resolution_status, ...
    "retained_within_tolerance")
verifyFalse(testCase, ...
    result.amplifier_interruptions.exclude_from_segments)
end

function testTerminalMissingReconnectEndsAtDisconnectAnchor(testCase)
terminal_keep = complete_keep;
terminal_keep{end} = (1:30)';
[eeg, edf] = build_pair(terminal_keep, terminal_keep);
eeg_start_anchor = 211;
eeg = attach_amplifier_interruptions( ...
    eeg, eeg_start_anchor, NaN, "missing_reconnect");
[result, eeg, edf] = match_triggers(eeg, edf);

report = extract_segments(result, eeg, edf);

verifyEqual(testCase, height(report), 1)
verifyEqual(testCase, report.eeg_end_sample, ...
    double(eeg.event_table.latency(eeg_start_anchor)))
verifyEqual(testCase, report.edf_end_sample, ...
    result.amplifier_interruptions.edf_start_anchor_latency)
verifyEqual(testCase, ...
    result.amplifier_interruptions.resolution_status, ...
    "terminal_eeg_loss")
end

function testDuplicateSameIntervalInterruptionsExcludeOnlyOnce(testCase)
[eeg, edf] = build_pair(complete_keep, complete_keep);
eeg_start_anchor = 70;
eeg_end_anchor = eeg_start_anchor + 1;
eeg.event_table.latency(eeg_end_anchor:end) = ...
    eeg.event_table.latency(eeg_end_anchor:end) - eeg.Fs;
eeg = attach_amplifier_interruptions( ...
    eeg, repmat(eeg_start_anchor, 2, 1), ...
    repmat(eeg_end_anchor, 2, 1), repmat("paired", 2, 1));
[result, eeg, edf] = match_triggers(eeg, edf);

report = extract_segments(result, eeg, edf);

verifyEqual(testCase, height(result.amplifier_interruptions), 2)
verifyEqual(testCase, ...
    numel(unique(result.amplifier_interruptions.gap_group_index)), 1)
verifyEqual(testCase, height(result.missing_periods), 1)
verifyEqual(testCase, height(report), 2)
verifyEqual(testCase, report.eeg_end_sample(1), ...
    double(eeg.event_table.latency(eeg_start_anchor)))
verifyEqual(testCase, report.eeg_start_sample(2), ...
    double(eeg.event_table.latency(eeg_end_anchor)))
end

function testOverlappingLeadingInterruptionDoesNotLeakEarlyInterval(testCase)
[eeg, edf] = build_pair(complete_keep, complete_keep);
eeg.event_table.latency(1:2) = eeg.event_table.latency(1:2) + eeg.Fs;
eeg = attach_amplifier_interruptions( ...
    eeg, [NaN; 2], [3; 3], ["paired"; "paired"]);
[result, eeg, edf] = match_triggers(eeg, edf);

report = extract_segments(result, eeg, edf);

verifyEqual(testCase, report.eeg_start_sample(1), ...
    double(eeg.event_table.latency(3)))
verifyEqual(testCase, report.edf_start_sample(1), ...
    double(edf.event_table.latency(3)))
end

function testMissingBoundaryMapsSuffixToSkippedCleanChunk(testCase)
[eeg, edf] = build_pair(missing_boundary_keep, complete_keep);
[result, eeg, edf] = match_triggers(eeg, edf);

report = extract_segments(result, eeg, edf);

verifyEqual(testCase, height(report), 2)
verifyEqual(testCase, ...
    report{:, {'start_eeg_chunk', 'start_edf_chunk', 'start_interval', ...
    'end_eeg_chunk', 'end_edf_chunk', 'end_interval'}}, ...
    [1 1 1 2 2 49; 2 3 7 4 5 49])
verifyEqual(testCase, report.eeg_start_sample(2), ...
    result.missing_periods.end_anchor_latency)
verifyEqual(testCase, report.edf_start_sample(2), ...
    double(edf.event_table.latency(107)))
end

function testMissingBoundaryAndSkippedCounterpartTimingStepBothSplit(testCase)
for boundary_system = ["eeg", "edf"]
    for step_sec = [-1, 1]
        if boundary_system == "eeg"
            [eeg, edf] = build_pair(missing_boundary_keep, complete_keep);
            edf.event_table.latency(121:end) = ...
                edf.event_table.latency(121:end) + step_sec * edf.trigger_Fs;
            eeg_start_index = [1; 101; 115];
            edf_start_index = [1; 107; 121];
            eeg_end_index = [100; 114; height(eeg.event_table)];
            edf_end_index = [100; 120; height(edf.event_table)];
        else
            [eeg, edf] = build_pair(complete_keep, missing_boundary_keep);
            eeg.event_table.latency(121:end) = ...
                eeg.event_table.latency(121:end) + step_sec * eeg.Fs;
            eeg_start_index = [1; 107; 121];
            edf_start_index = [1; 101; 115];
            eeg_end_index = [100; 120; height(eeg.event_table)];
            edf_end_index = [100; 114; height(edf.event_table)];
        end
        [result, eeg, edf] = match_triggers(eeg, edf);

        report = extract_segments(result, eeg, edf);

        verifyEqual(testCase, height(report), 3)
        verifyEqual(testCase, report.eeg_start_sample, ...
            double(eeg.event_table.latency(eeg_start_index)))
        verifyEqual(testCase, report.edf_start_sample, ...
            double(edf.event_table.latency(edf_start_index)))
        verifyEqual(testCase, report.eeg_end_sample, ...
            double(eeg.event_table.latency(eeg_end_index)))
        verifyEqual(testCase, report.edf_end_sample, ...
            double(edf.event_table.latency(edf_end_index)))
        boundary_loss = result.missing_periods( ...
            result.missing_periods.reason == "suspected_missing_boundary", :);
        verifyEqual(testCase, boundary_loss.system, boundary_system)
        verifyEqual(testCase, boundary_loss.start_anchor_event_index, 100)
        verifyEqual(testCase, boundary_loss.end_anchor_event_index, 101)
        verify_segment_interval_timing(testCase, report, eeg, edf, ...
            result.interval_tolerance_sec)
    end
end
end

function testMissingBoundaryAndMultipleSkippedCounterpartStepsAllSplit(testCase)
for boundary_system = ["eeg", "edf"]
    if boundary_system == "eeg"
        [eeg, edf] = build_pair(missing_boundary_keep, complete_keep);
        edf.event_table.latency(121:end) = ...
            edf.event_table.latency(121:end) - edf.trigger_Fs;
        edf.event_table.latency(131:end) = ...
            edf.event_table.latency(131:end) + 2 * edf.trigger_Fs;
        eeg_start_index = [1; 101; 115; 125];
        edf_start_index = [1; 107; 121; 131];
        eeg_end_index = [100; 114; 124; height(eeg.event_table)];
        edf_end_index = [100; 120; 130; height(edf.event_table)];
    else
        [eeg, edf] = build_pair(complete_keep, missing_boundary_keep);
        eeg.event_table.latency(121:end) = ...
            eeg.event_table.latency(121:end) - eeg.Fs;
        eeg.event_table.latency(131:end) = ...
            eeg.event_table.latency(131:end) + 2 * eeg.Fs;
        eeg_start_index = [1; 107; 121; 131];
        edf_start_index = [1; 101; 115; 125];
        eeg_end_index = [100; 120; 130; height(eeg.event_table)];
        edf_end_index = [100; 114; 124; height(edf.event_table)];
    end
    [result, eeg, edf] = match_triggers(eeg, edf);

    report = extract_segments(result, eeg, edf);

    verifyEqual(testCase, height(report), 4)
    verifyEqual(testCase, report.eeg_start_sample, ...
        double(eeg.event_table.latency(eeg_start_index)))
    verifyEqual(testCase, report.edf_start_sample, ...
        double(edf.event_table.latency(edf_start_index)))
    verifyEqual(testCase, report.eeg_end_sample, ...
        double(eeg.event_table.latency(eeg_end_index)))
    verifyEqual(testCase, report.edf_end_sample, ...
        double(edf.event_table.latency(edf_end_index)))
    verify_segment_interval_timing(testCase, report, eeg, edf, ...
        result.interval_tolerance_sec)
end
end

function testMissingBoundaryAndSkippedAmplifierInterruptionBothSplit(testCase)
[eeg, edf] = build_pair(complete_keep, missing_boundary_keep);
eeg.event_table.latency(121:end) = ...
    eeg.event_table.latency(121:end) - eeg.Fs;
eeg = attach_amplifier_interruptions(eeg, 120, 121, "paired");
[result, eeg, edf] = match_triggers(eeg, edf);

report = extract_segments(result, eeg, edf);

verifyEqual(testCase, height(report), 3)
verifyEqual(testCase, report.eeg_end_sample(1:2), ...
    double(eeg.event_table.latency([100; 120])))
verifyEqual(testCase, report.eeg_start_sample(2:3), ...
    double(eeg.event_table.latency([107; 121])))
verifyEqual(testCase, report.edf_end_sample(1:2), ...
    double(edf.event_table.latency([100; 114])))
verifyEqual(testCase, report.edf_start_sample(2:3), ...
    double(edf.event_table.latency([101; 115])))
verifyEqual(testCase, result.amplifier_interruptions.resolution_status, ...
    "excluded_and_relocked")
verifyTrue(testCase, result.amplifier_interruptions.exclude_from_segments)
verifyEqual(testCase, sort(result.missing_periods.system), ["edf"; "eeg"])
verify_segment_interval_timing(testCase, report, eeg, edf, ...
    result.interval_tolerance_sec)
end

function testMissingBoundaryAndSkippedCounterpartMissingTriggerBothSplit(testCase)
counterpart_keep = complete_keep;
counterpart_keep{3} = setdiff((1:50)', 20, 'stable');
for boundary_system = ["eeg", "edf"]
    if boundary_system == "eeg"
        [eeg, edf] = build_pair(missing_boundary_keep, counterpart_keep);
        eeg_start_index = [1; 101; 115];
        edf_start_index = [1; 107; 120];
        eeg_end_index = [100; 113; height(eeg.event_table)];
        edf_end_index = [100; 119; height(edf.event_table)];
    else
        [eeg, edf] = build_pair(counterpart_keep, missing_boundary_keep);
        eeg_start_index = [1; 107; 120];
        edf_start_index = [1; 101; 115];
        eeg_end_index = [100; 119; height(eeg.event_table)];
        edf_end_index = [100; 113; height(edf.event_table)];
    end
    [result, eeg, edf] = match_triggers(eeg, edf);

    report = extract_segments(result, eeg, edf);

    verifyEqual(testCase, height(report), 3)
    verifyEqual(testCase, report.eeg_start_sample, ...
        double(eeg.event_table.latency(eeg_start_index)))
    verifyEqual(testCase, report.edf_start_sample, ...
        double(edf.event_table.latency(edf_start_index)))
    verifyEqual(testCase, report.eeg_end_sample, ...
        double(eeg.event_table.latency(eeg_end_index)))
    verifyEqual(testCase, report.edf_end_sample, ...
        double(edf.event_table.latency(edf_end_index)))
    verify_segment_interval_timing(testCase, report, eeg, edf, ...
        result.interval_tolerance_sec)
end
end

function testMissingBoundaryAndAdjacentSkippedMarkersExcludeWholeGroup(testCase)
[eeg, edf] = build_pair(complete_keep, missing_boundary_keep);
eeg = attach_amplifier_interruptions( ...
    eeg, [120; 121], [121; 122], "paired");
[result, eeg, edf] = match_triggers(eeg, edf);

report = extract_segments(result, eeg, edf);

verifyEqual(testCase, height(report), 3)
verifyEqual(testCase, report.eeg_end_sample(1:2), ...
    double(eeg.event_table.latency([100; 120])))
verifyEqual(testCase, report.eeg_start_sample(2:3), ...
    double(eeg.event_table.latency([107; 122])))
verifyEqual(testCase, report.edf_end_sample(1:2), ...
    double(edf.event_table.latency([100; 114])))
verifyEqual(testCase, report.edf_start_sample(2:3), ...
    double(edf.event_table.latency([101; 116])))
verifyEqual(testCase, result.amplifier_interruptions.resolution_status, ...
    repmat("excluded_and_relocked", 2, 1))
verifyTrue(testCase, all(result.amplifier_interruptions.exclude_from_segments))
verifyEqual(testCase, ...
    numel(unique(result.amplifier_interruptions.gap_group_index)), 1)
boundary_loss = result.missing_periods( ...
    result.missing_periods.reason == "suspected_missing_boundary", :);
verifyEqual(testCase, boundary_loss.system, "edf")
verifyEqual(testCase, boundary_loss.start_anchor_event_index, 100)
verifyEqual(testCase, boundary_loss.end_anchor_event_index, 101)
verify_segment_interval_timing(testCase, report, eeg, edf, ...
    result.interval_tolerance_sec)
end

function testMissingBoundaryAndSkippedMultiIntervalMarkersBothSplit(testCase)
for anchors = [120, 149; 123, 152]
    [eeg, edf] = build_pair(complete_keep, missing_boundary_keep);
    eeg.event_table.latency(anchors(2):end) = ...
        eeg.event_table.latency(anchors(2):end) - eeg.Fs;
    eeg = attach_amplifier_interruptions( ...
        eeg, anchors(1), anchors(2), "paired");
    [result, eeg, edf] = match_triggers(eeg, edf);

    report = extract_segments(result, eeg, edf);

    verifyEqual(testCase, height(report), 3)
    verifyEqual(testCase, report.eeg_end_sample(1:2), ...
        double(eeg.event_table.latency([100; anchors(1)])))
    verifyEqual(testCase, report.eeg_start_sample(2:3), ...
        double(eeg.event_table.latency([107; anchors(2)])))
    verifyEqual(testCase, report.edf_end_sample(1:2), ...
        double(edf.event_table.latency([100; anchors(1) - 6])))
    verifyEqual(testCase, report.edf_start_sample(2:3), ...
        double(edf.event_table.latency([101; anchors(2) - 6])))
    verifyEqual(testCase, result.amplifier_interruptions.resolution_status, ...
        "excluded_and_relocked")
    verifyTrue(testCase, result.amplifier_interruptions.exclude_from_segments)
    verify_segment_interval_timing(testCase, report, eeg, edf, ...
        result.interval_tolerance_sec)
end
end

function testSharedTriggerGlitchStaysContinuous(testCase)
[eeg, edf] = build_pair(short_loss_keep, short_loss_keep);
[result, eeg, edf] = match_triggers(eeg, edf);

report = extract_segments(result, eeg, edf);

verifyEqual(testCase, height(report), 1)
verifyEqual(testCase, report.start_interval, 1)
verifyEqual(testCase, report.end_interval, 49)
end

function testSharedGlitchDoesNotRelockAnOpenLoss(testCase)
[eeg_keep, edf_keep] = open_loss_then_shared_glitch_keep;
[eeg, edf] = build_pair(eeg_keep, edf_keep);
[result, eeg, edf] = match_triggers(eeg, edf);

verifyEqual(testCase, result.comparison_history.outcome, ...
    ["match"; "eeg_data_loss"; "shared_trigger_glitch"; ...
    "match"; "match"; "same_last_trigger"])
verifyEqual(testCase, result.missing_periods.start_anchor_time_sec, 3058)
verifyEqual(testCase, result.missing_periods.end_anchor_time_sec, 5175)
report = extract_segments(result, eeg, edf);

verifyEqual(testCase, height(report), 2)
verifyEqual(testCase, ...
    report{:, {'start_eeg_chunk', 'start_edf_chunk', 'start_interval', ...
    'end_eeg_chunk', 'end_edf_chunk', 'end_interval'}}, ...
    [1 1 1 2 2 43; 4 4 1 6 6 49])
end

function testSharedGlitchDoesNotRelockAnOpenEdfLoss(testCase)
[eeg_keep, edf_keep] = open_loss_then_shared_glitch_keep;
[eeg, edf] = build_pair(edf_keep, eeg_keep);
[result, eeg, edf] = match_triggers(eeg, edf);

verifyEqual(testCase, result.comparison_history.outcome, ...
    ["match"; "edf_data_loss"; "shared_trigger_glitch"; ...
    "match"; "match"; "same_last_trigger"])
verifyEqual(testCase, result.missing_periods.start_anchor_time_sec, 3158)
verifyEqual(testCase, result.missing_periods.end_anchor_time_sec, 5275)
report = extract_segments(result, eeg, edf);

verifyEqual(testCase, height(report), 2)
verifyEqual(testCase, ...
    report{:, {'start_eeg_chunk', 'start_edf_chunk', 'start_interval', ...
    'end_eeg_chunk', 'end_edf_chunk', 'end_interval'}}, ...
    [1 1 1 2 2 43; 4 4 1 6 6 49])
end

function testSharedTerminalGlitchKeepsCanonicalIntervalNumber(testCase)
[eeg, edf] = build_pair(shared_terminal_glitch_keep, ...
    shared_terminal_glitch_keep);
[result, eeg, edf] = match_triggers(eeg, edf);

report = extract_segments(result, eeg, edf);

verifyEqual(testCase, height(report), 1)
verifyEqual(testCase, report.end_interval, 49)
end

function testDurationsUseElapsedSampleSpans(testCase)
[eeg, edf] = build_pair(complete_keep, complete_keep);
edf.event_table.latency(end) = edf.event_table.latency(end) - 1;
[result, eeg, edf] = match_triggers(eeg, edf); %#ok<ASGLU>

display_text = evalc( ...
    'report = extract_segments(result, eeg, edf);');

eeg_span = report.eeg_end_sample - report.eeg_start_sample;
edf_span = report.edf_end_sample - report.edf_start_sample;
expected_eeg_duration = eeg_span ./ eeg.Fs;
expected_edf_duration = edf_span ./ edf.trigger_Fs;
verifyEqual(testCase, report.eeg_duration_sec, ...
    expected_eeg_duration, 'AbsTol', 1e-12)
verifyEqual(testCase, report.edf_duration_sec, ...
    expected_edf_duration, 'AbsTol', 1e-12)
verifyEqual(testCase, report.edf_minus_eeg_sec, ...
    expected_edf_duration - expected_eeg_duration, 'AbsTol', 1e-12)
expected_eeg_duration_text = compose('%.5f', report.eeg_duration_sec);
expected_edf_duration_text = compose('%.5f', report.edf_duration_sec);
verifyTrue(testCase, contains( ...
    string(display_text), expected_eeg_duration_text))
verifyTrue(testCase, contains( ...
    string(display_text), expected_edf_duration_text))
verifyClass(testCase, report.eeg_duration_sec, 'double')
verifyClass(testCase, report.edf_duration_sec, 'double')
end

function testClockDriftMetricsAndPlotDefault(testCase)
[eeg, edf] = build_pair(complete_keep, complete_keep);
[result, eeg, edf] = match_triggers(eeg, edf);

figure_count_before = numel(findall(groot, 'Type', 'figure'));
report = extract_segments(result, eeg, edf);
figure_count_after = numel(findall(groot, 'Type', 'figure'));

verifyEqual(testCase, figure_count_after, figure_count_before)
verifyEqual(testCase, report.clock_drift_r_squared, 1)
verifyEqual(testCase, report.clock_drift_max_abs_residual_sec, 0, ...
    'AbsTol', 1e-10)
end

function testLinearClockDriftIsReportedForEverySegment(testCase)
[eeg, edf] = build_pair(short_loss_keep, complete_keep);
[result, eeg, edf] = match_triggers(eeg, edf);
edf = add_clock_drift(edf, 50);

report = extract_segments(result, eeg, edf);

verifyEqual(testCase, height(report), 2)
verifyGreaterThan(testCase, report.clock_drift_r_squared, 0.98)
verifyLessThanOrEqual(testCase, ...
    report.clock_drift_max_abs_residual_sec, ...
    result.interval_tolerance_sec)
end

function testPlotOverlaysCumulativeDifferenceBySegment(testCase)
[eeg, edf] = build_pair(short_loss_keep, complete_keep);
[result, eeg, edf] = match_triggers(eeg, edf);
edf = add_clock_drift(edf, 50);

original_visibility = get(groot, 'defaultFigureVisible');
visibility_cleanup = onCleanup(@() set( ...
    groot, 'defaultFigureVisible', original_visibility));
set(groot, 'defaultFigureVisible', 'off')
figure_count_before = numel(findall(groot, 'Type', 'figure'));

report = extract_segments(result, eeg, edf, true);
figure_handle = gcf;
figure_cleanup = onCleanup(@() close(figure_handle));

verifyEqual(testCase, ...
    numel(findall(groot, 'Type', 'figure')), figure_count_before + 1)
verifyEmpty(testCase, figure_handle.Name)
verifyEmpty(testCase, figure_handle.Tag)
lines = findobj(figure_handle, 'Type', 'line');
verifyNumElements(testCase, lines, height(report))
verifyEqual(testCase, sort(string({lines.DisplayName})), ...
    ["segment index 1", "segment index 2"])
verifyTrue(testCase, all(string({lines.LineStyle}) == "-"))
verifyTrue(testCase, all(string({lines.Marker}) == "o"))
verifyTrue(testCase, all([lines.LineWidth] == 1))
axes_handle = findobj(figure_handle, 'Type', 'axes');
legend_handle = findobj(figure_handle, 'Type', 'legend');
verifyNumElements(testCase, axes_handle, 1)
verifyNumElements(testCase, legend_handle, 1)
verifyEqual(testCase, ...
    char(axes_handle.InteractionOptions.DatatipsSupported), 'off')
verifyEqual(testCase, ...
    char(axes_handle.InteractionOptions.ZoomSupported), 'on')
verifyEqual(testCase, axes_handle.FontSize, 16)
verifyEqual(testCase, axes_handle.Title.FontSize, 24)
verifyEqual(testCase, legend_handle.FontSize, 20)
verifyEqual(testCase, string(axes_handle.XLabel.String), ...
    "Time within continuous segment (hr)")

for line_i = 1:numel(lines)
    segment_i = sscanf(lines(line_i).DisplayName, 'segment index %d');
    verifyEqual(testCase, lines(line_i).XData(1), 0, ...
        'AbsTol', 1e-12)
    verifyEqual(testCase, lines(line_i).XData(end), ...
        report.eeg_duration_sec(segment_i) / 3600, 'AbsTol', 1e-12)
    verifyGreaterThan(testCase, diff(lines(line_i).XData), 0)
    verifyEqual(testCase, lines(line_i).YData(1), 0, ...
        'AbsTol', 1e-12)
    verifyEqual(testCase, lines(line_i).YData(end), ...
        report.edf_minus_eeg_sec(segment_i), 'AbsTol', 1e-12)
end
end

function testMutatedIntervalFailsBeforeDriftFit(testCase)
[eeg, edf] = build_pair(short_loss_keep, complete_keep);
[result, eeg, edf] = match_triggers(eeg, edf);
report = extract_segments(result, eeg, edf);

segment_event_index = find( ...
    double(edf.event_table.latency) >= report.edf_start_sample(2) & ...
    double(edf.event_table.latency) <= report.edf_end_sample(2));
step_event_index = segment_event_index(ceil(numel(segment_event_index) / 2));
step_samples = round(0.1 * edf.trigger_Fs);
edf.event_table.latency(step_event_index:end) = ...
    edf.event_table.latency(step_event_index:end) + step_samples;

verifyError(testCase, ...
    @() extract_segments(result, eeg, edf), ...
    'extract_segments:UnvalidatedIntervalPair')
end

function testCumulativeNonlinearityFailsWithLocallyMatchedIntervals(testCase)
[eeg, edf] = build_pair(short_loss_keep, complete_keep);
[result, eeg, edf] = match_triggers(eeg, edf);
report = extract_segments(result, eeg, edf);

eeg_event_index = find( ...
    double(eeg.event_table.latency) >= report.eeg_start_sample(2) & ...
    double(eeg.event_table.latency) <= report.eeg_end_sample(2));
edf_event_index = find( ...
    double(edf.event_table.latency) >= report.edf_start_sample(2) & ...
    double(edf.event_table.latency) <= report.edf_end_sample(2));
ramp_start = edf_event_index(ceil(numel(edf_event_index) / 2));
offset_samples = min(20, max(0, ...
    (1:height(edf.event_table))' - ramp_start + 1));
edf.event_table.latency = edf.event_table.latency + offset_samples;

interval_difference_sec = ...
    diff(double(edf.event_table.latency(edf_event_index))) / edf.trigger_Fs - ...
    diff(double(eeg.event_table.latency(eeg_event_index))) / eeg.Fs;
verifyLessThanOrEqual(testCase, abs(interval_difference_sec), ...
    result.interval_tolerance_sec)
verifyError(testCase, ...
    @() extract_segments(result, eeg, edf), ...
    'extract_segments:NonlinearClockDrift')
end

function [eeg, edf] = build_pair(eeg_keep, edf_keep)

eeg.Fs = 500;
eeg.event_table = build_event_table(eeg.Fs, 0, eeg_keep);
edf.Fs = 256;
edf.trigger_Fs = 128;
edf.event_table = build_event_table(edf.trigger_Fs, 100, edf_keep);

end

function verify_segment_interval_timing(testCase, report, eeg, edf, tolerance_sec)

for segment_i = 1:height(report)
    eeg_events = eeg.event_table( ...
        double(eeg.event_table.latency) >= report.eeg_start_sample(segment_i) & ...
        double(eeg.event_table.latency) <= report.eeg_end_sample(segment_i), :);
    edf_events = edf.event_table( ...
        double(edf.event_table.latency) >= report.edf_start_sample(segment_i) & ...
        double(edf.event_table.latency) <= report.edf_end_sample(segment_i), :);
    verifyEqual(testCase, eeg_events.type, edf_events.type)
    verifyLessThanOrEqual(testCase, ...
        abs(diff(double(eeg_events.latency)) / eeg.Fs - ...
        diff(double(edf_events.latency)) / edf.trigger_Fs), tolerance_sec)
end

end

function system = attach_amplifier_interruptions( ...
    system, eeg_start_anchor, eeg_end_anchor, pairing_status)

eeg_start_anchor = eeg_start_anchor(:);
eeg_end_anchor = eeg_end_anchor(:);
pairing_status = string(pairing_status(:));
n_interruptions = numel(eeg_start_anchor);
if isscalar(eeg_end_anchor) && n_interruptions > 1
    eeg_end_anchor = repmat(eeg_end_anchor, n_interruptions, 1);
end
if isscalar(pairing_status) && n_interruptions > 1
    pairing_status = repmat(pairing_status, n_interruptions, 1);
end
assert(numel(eeg_end_anchor) == n_interruptions && ...
    numel(pairing_status) == n_interruptions, ...
    'Synthetic interruption inputs must have equal lengths.')

interruption_id = (1:n_interruptions)';
disconnect_raw_event_index = ...
    height(system.event_table) + (1:2:(2 * n_interruptions))';
disconnect_latency = nan(n_interruptions, 1);
eeg_start_anchor_latency = nan(n_interruptions, 1);
finite_start = isfinite(eeg_start_anchor);
disconnect_latency(finite_start) = double( ...
    system.event_table.latency(eeg_start_anchor(finite_start)));
eeg_start_anchor_latency(finite_start) = ...
    disconnect_latency(finite_start);
disconnect_latency(~finite_start) = ...
    double(system.event_table.latency(1)) - 1;
disconnect_event_text = repmat( ...
    "9001, Amplifier disconnected", n_interruptions, 1);
reconnect_raw_event_index = nan(n_interruptions, 1);
reconnect_latency = nan(n_interruptions, 1);
reconnect_event_text = strings(n_interruptions, 1);
eeg_end_anchor_latency = nan(n_interruptions, 1);

for interruption_i = 1:n_interruptions
    if pairing_status(interruption_i) == "paired"
        reconnect_raw_event_index(interruption_i) = ...
            disconnect_raw_event_index(interruption_i) + 1;
        reconnect_latency(interruption_i) = double( ...
            system.event_table.latency(eeg_end_anchor(interruption_i)));
        reconnect_event_text(interruption_i) = ...
            "9002, Amplifier reconnected";
        eeg_end_anchor_latency(interruption_i) = ...
            reconnect_latency(interruption_i);
    else
        assert(pairing_status(interruption_i) == "missing_reconnect", ...
            'Unsupported synthetic amplifier interruption status.')
        eeg_end_anchor(interruption_i) = NaN;
    end
end

system.amplifier_interruptions = table( ...
    interruption_id, pairing_status, ...
    disconnect_raw_event_index, disconnect_latency, ...
    disconnect_event_text, ...
    reconnect_raw_event_index, reconnect_latency, reconnect_event_text, ...
    eeg_start_anchor, eeg_start_anchor_latency, ...
    eeg_end_anchor, eeg_end_anchor_latency, ...
    'VariableNames', { ...
        'interruption_id', 'pairing_status', ...
        'disconnect_raw_event_index', 'disconnect_latency', ...
        'disconnect_event_text', ...
        'reconnect_raw_event_index', 'reconnect_latency', ...
        'reconnect_event_text', ...
        'eeg_start_anchor_event_index', 'eeg_start_anchor_latency', ...
        'eeg_end_anchor_event_index', 'eeg_end_anchor_latency'});

end

function system = add_clock_drift(system, slope_ppm)

first_latency = double(system.event_table.latency(1));
relative_latency = double(system.event_table.latency) - first_latency;
system.event_table.latency = first_latency + round( ...
    relative_latency .* (1 + slope_ppm * 1e-6));

end

function event_table = build_event_table(Fs, start_time_sec, keep_by_cycle)

[canonical_type, canonical_interval_sec, canonical_time_sec] = ...
    canonical_template;
cycle_duration_sec = sum(canonical_interval_sec);
event_time_sec = zeros(0, 1);
event_type = strings(0, 1);

for cycle_i = 1:numel(keep_by_cycle)
    keep_index = keep_by_cycle{cycle_i};
    cycle_start_sec = start_time_sec + ...
        (cycle_i - 1) * cycle_duration_sec;
    event_time_sec = [event_time_sec; ...
        cycle_start_sec + canonical_time_sec(keep_index)]; %#ok<AGROW>
    event_type = [event_type; canonical_type(keep_index)]; %#ok<AGROW>
end

latency = round(event_time_sec .* Fs) + 1;
event_table = table(latency, event_type, ...
    'VariableNames', {'latency', 'type'});

end

function keep_by_cycle = complete_keep

keep_by_cycle = repmat({(1:50)'}, 5, 1);

end

function keep_by_cycle = partial_start_keep

keep_by_cycle = complete_keep;
keep_by_cycle{1} = (25:50)';

end

function keep_by_cycle = short_loss_keep

keep_by_cycle = complete_keep;
keep_by_cycle{2} = setdiff((1:50)', (12:13)', 'stable');

end

function keep_by_cycle = missing_boundary_keep

keep_by_cycle = complete_keep;
keep_by_cycle{3} = (7:50)';

end

function [eeg_keep, edf_keep] = open_loss_then_shared_glitch_keep

eeg_keep = repmat({(1:50)'}, 6, 1);
edf_keep = eeg_keep;
eeg_keep{2} = (1:44)';
eeg_keep{3} = setdiff((1:50)', (12:13)', 'stable');
edf_keep{3} = eeg_keep{3};

end

function keep_by_cycle = shared_terminal_glitch_keep

keep_by_cycle = complete_keep;
keep_by_cycle{end} = setdiff((1:50)', (12:13)', 'stable');

end

function [type, interval_to_next_sec, relative_time_sec] = ...
    canonical_template

type = [repmat("64", 4, 1); repmat("63", 46, 1)];
interval_to_next_sec = (10:59)';
relative_time_sec = [0; cumsum(interval_to_next_sec(1:end - 1))];

end
