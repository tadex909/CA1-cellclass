function [t_ds, X_ds_n, XSpeed] = behavior_vectors_from_traj(Traj)
%BEHAVIOR_VECTORS_FROM_TRAJ Build continuous behavior vectors from Traj.
%
% theta_s.m expects t_ds, X_ds_n, and XSpeed to already exist. Some
% TrajData.mat files only contain the trial struct, so this helper stitches
% the trial fields into session-level vectors.

allTime = [];
allPosition = [];
allSpeed = [];

assert(isfield(Traj, 'time'), ...
    'Traj does not contain the required field: time');
assert(isfield(Traj, 'VRtraj'), ...
    'Traj does not contain the required field: VRtraj');
assert(isfield(Traj, 'XSpeed') || isfield(Traj, 'Speed'), ...
    'Traj does not contain XSpeed or Speed.');

for k = 1:numel(Traj)
    trialTime = double(Traj(k).time(:));
    trialPosition = double(Traj(k).VRtraj(:));

    if isfield(Traj, 'XSpeed') && ~isempty(Traj(k).XSpeed)
        trialSpeed = double(Traj(k).XSpeed(:));
    else
        trialSpeed = double(Traj(k).Speed(:));
    end

    if isfield(Traj, 'tstart') && ~isempty(Traj(k).tstart)
        trialTime = trialTime + double(Traj(k).tstart);
    end

    n = min([numel(trialTime), numel(trialPosition), numel(trialSpeed)]);
    if n == 0
        continue;
    end

    allTime = [allTime; trialTime(1:n)]; %#ok<AGROW>
    allPosition = [allPosition; trialPosition(1:n)]; %#ok<AGROW>
    allSpeed = [allSpeed; trialSpeed(1:n)]; %#ok<AGROW>
end

assert(~isempty(allTime), 'Could not reconstruct behavior vectors from Traj.');

[t_ds, order] = sort(allTime);
X_ds_n = allPosition(order);
XSpeed = allSpeed(order);
end
