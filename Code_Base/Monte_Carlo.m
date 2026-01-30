%% Monte Carlo Simulation: Dynamic Mutant Landscape vs Wild Type
% This script compares protein folding dynamics between:
% 1. Control: Static wild-type energy landscape
% 2. Experiment: Dynamically switching mutant energy landscapes

clear; clc; close all;

%% PARAMETERS
% File and Simulation Settings
folderName = 'Free_Energies';
numSteps = 20*10000000;                        % Total Monte Carlo steps
k = 20;                                     %Number of switches
switchInterval = numSteps/k;                    % Steps between landscape switches
numResidues = 76;                          % Number of residues (mutants)
To_mask = [6,11,16,18,24,32,33,34,39,42,48,51,52,54,58,63,64,72,74];
% Physics Parameters
T = 310;                                   % Temperature (K)
R = 0.008314;                              % Gas constant (kJ/mol/K)
beta = 1 / (R * T);                        % Inverse thermal energy

% Output Settings
saveData = false;                          % Save trajectory data to file

%% 1. PRE-LOAD DATA
fprintf('Loading energy profiles...\n');

% Load Wild Type (WT) Profile
wtFile = fullfile(folderName, 'WT.txt');
if ~isfile(wtFile)
    error('Wild type file not found: %s', wtFile);
end

% Read WT data
wtData = load(wtFile);
if size(wtData, 2) >= 2
    reactionCoord = wtData(:, 1);
    wtEnergy = wtData(:, 2);
else
    wtEnergy = wtData(:, 1);
    reactionCoord = (0:length(wtEnergy)-1)';
    warning('WT file only had 1 column. Generated reaction coordinate.');
end

numStates = length(wtEnergy);
fprintf('  Loaded WT Profile: %d states\n', numStates);

% Pre-load All Mutant Profiles
mutantEnergies = zeros(numStates, numResidues);
fprintf('  Loading %d mutant profiles... ', numResidues);

for j = 1:length(To_mask)
    i = To_mask(j);
    mutantFile = fullfile(folderName, sprintf('Mutant_%d.txt', i));
    
    if isfile(mutantFile)
        mutData = load(mutantFile);
        
        % Extract energy column
        if size(mutData, 2) >= 2
            tmpEnergy = mutData(:, 2);
        else
            tmpEnergy = mutData(:, 1);
        end
        
        % Handle size mismatches
        if length(tmpEnergy) ~= numStates
            minLen = min(length(tmpEnergy), numStates);
            mutantEnergies(1:minLen, i) = tmpEnergy(1:minLen);
            if numStates > minLen
                % Pad with last value if mutant profile is shorter
                mutantEnergies(minLen+1:end, i) = tmpEnergy(end);
            end
            warning('Mutant_%d.txt size mismatch. Adjusted to match WT.', i);
        else
            mutantEnergies(:, i) = tmpEnergy;
        end
    else
        % If mutant file doesn't exist, use WT energy as fallback
        mutantEnergies(:, i) = wtEnergy;
        warning('Mutant_%d.txt not found. Using WT energy.', i);
    end
end
fprintf('Done.\n\n');

%% 2. INITIALIZATION
fprintf('Initializing Monte Carlo simulation...\n');

% Start both walkers at lowest energy state (most stable)
[~, startIdx] = min(wtEnergy);
idxWT = startIdx;
idxDynamic = startIdx;

% Preallocate trajectory arrays
trajWT = zeros(numSteps, 1);
trajDynamic = zeros(numSteps, 1);

% Log when landscape switches occur
numSwitches = floor(numSteps / switchInterval);
switchLog = zeros(numSwitches, 2);  % [step, mutantID]
switchCount = 0;

% Initialize dynamic walker's energy landscape
currentEnergy = wtEnergy;

% Pre-generate random numbers for efficiency
randMoveWT = rand(numSteps, 1);
randAcceptWT = rand(numSteps, 1);
randMoveDynamic = rand(numSteps, 1);
randAcceptDynamic = rand(numSteps, 1);

fprintf('  Starting position: state %d (lowest energy)\n', startIdx);
fprintf('  Temperature: %.1f K\n', T);
fprintf('  Landscape switches every %d steps\n', switchInterval);
fprintf('  Expected number of switches: %d\n\n', numSwitches);

%% 3. MONTE CARLO SIMULATION
fprintf('Running Monte Carlo simulation (%d steps)...\n', numSteps);
tic;

for t = 1:numSteps
    
    % --- DYNAMIC WALKER: Check for landscape switch ---
    if mod(t, switchInterval) == 0
        switchCount = switchCount + 1;
        
        
        randIdx = randi(length(To_mask)); % Pick a random index from the mask list
        mutantID = To_mask(randIdx); % Randomly select a mutant landscape
        currentEnergy = mutantEnergies(:, mutantID);
        
        % Log the switch
        switchLog(switchCount, :) = [t, mutantID];
        
        fprintf('  Step %d: Switched to mutant %d (%d/%d switches)\n', ...
                    t, mutantID, switchCount, numSwitches);
        
    end
    
    % --- WALKER 1: WILD TYPE (Control) ---
    % Propose a move: -1 (left) or +1 (right)
    if randMoveWT(t) >= 0.5
        step = 1;
    else
        step = -1;
    end
    proposedWT = idxWT + step;
    
    % Check boundaries
    if proposedWT >= 1 && proposedWT <= numStates
        % Calculate energy change
        deltaE = wtEnergy(proposedWT) - wtEnergy(idxWT);
        
        % Metropolis acceptance criterion
        if deltaE < 0 || randAcceptWT(t) <= exp(-deltaE * beta)
            idxWT = proposedWT;
        end
    end
    
    % --- WALKER 2: DYNAMIC (Experiment) ---
    % Propose a move: -1 (left) or +1 (right)
    if randMoveDynamic(t) >= 0.5
        step = 1;
    else
        step = -1;
    end
    proposedDynamic = idxDynamic + step;
    
    % Check boundaries
    if proposedDynamic >= 1 && proposedDynamic <= numStates
        % Calculate energy change on current landscape
        deltaE = currentEnergy(proposedDynamic) - currentEnergy(idxDynamic);
        
        % Metropolis acceptance criterion
        if deltaE < 0 || randAcceptDynamic(t) <= exp(-deltaE * beta)
            idxDynamic = proposedDynamic;
        end
    end
    
    % --- RECORD TRAJECTORIES ---
    trajWT(t) = reactionCoord(idxWT);
    trajDynamic(t) = reactionCoord(idxDynamic);
end

elapsedTime = toc;
fprintf('Simulation complete! (%.2f seconds)\n\n', elapsedTime);

%% 4. SAVE DATA (Optional)
if saveData
    fprintf('Saving trajectory data...\n');
    outputData = [(1:numSteps)', trajWT, trajDynamic];
    outputFile = 'MC_trajectories.txt';
    
    fid = fopen(outputFile, 'w');
    fprintf(fid, '%% Monte Carlo Trajectories\n');
    fprintf(fid, '%% Column 1: Step\n');
    fprintf(fid, '%% Column 2: WT Reaction Coordinate\n');
    fprintf(fid, '%% Column 3: Dynamic Reaction Coordinate\n');
    fprintf(fid, '%% Temperature: %.1f K\n', T);
    fprintf(fid, '%% Switch Interval: %d steps\n', switchInterval);
    fprintf(fid, '%10s %15s %15s\n', 'Step', 'WT_RC', 'Dynamic_RC');
    
    for i = 1:numSteps
        fprintf(fid, '%10d %15.6f %15.6f\n', outputData(i,:));
    end
    fclose(fid);
    
    fprintf('  Saved to: %s\n\n', outputFile);
end

%% 5. VISUALIZATION
fprintf('Generating plots...\n');

figure('Color', 'w', 'Position', [100, 100, 1200, 800]);

% Determine Y-axis limits
yMax = max(reactionCoord) + 5;
yMin = min(reactionCoord) - 5;

% --- Subplot 1: Wild Type (Top) ---
ax1 = subplot(2, 1, 1);
plot(1:numSteps, trajWT, 'k-', 'LineWidth', 0.5);
ylabel('Number of Folded Residues', 'FontSize', 12, 'FontWeight', 'bold');
title('Control: Wild Type Energy Landscape (Static)', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
ylim([yMin, yMax]);
set(gca, 'FontSize', 11);

% --- Subplot 2: Dynamic (Bottom) ---
ax2 = subplot(2, 1, 2);
plot(1:numSteps, trajDynamic, 'b-', 'LineWidth', 0.5);
hold on;

% Mark landscape switches with red dotted lines
for i = 1:size(switchLog, 1)
    if switchLog(i, 1) > 0
        xline(switchLog(i, 1), 'r:', 'LineWidth', 1, 'Alpha', 0.6);
    end
end

ylabel('Number of Folded Residues', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Monte Carlo Steps', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Experiment: Dynamic Mutant Landscape (Switch every %d steps)', switchInterval), ...
      'FontSize', 14, 'FontWeight', 'bold');
grid on;
ylim([yMin, yMax]);
set(gca, 'FontSize', 11);

% Link x-axes for synchronized zooming
linkaxes([ax1, ax2], 'x');

fprintf('Visualization complete!\n');

%% 6. SUMMARY STATISTICS
fprintf('\n=== SIMULATION SUMMARY ===\n');
fprintf('Total steps: %d\n', numSteps);
fprintf('Number of landscape switches: %d\n', switchCount);
fprintf('Unique mutants visited: %d\n', length(unique(switchLog(:,2))));
fprintf('\nWild Type Walker:\n');
fprintf('  Mean position: %.2f residues\n', mean(trajWT));
fprintf('  Std deviation: %.2f residues\n', std(trajWT));
fprintf('  Min/Max: %.0f / %.0f residues\n', min(trajWT), max(trajWT));
fprintf('\nDynamic Walker:\n');
fprintf('  Mean position: %.2f residues\n', mean(trajDynamic));
fprintf('  Std deviation: %.2f residues\n', std(trajDynamic));
fprintf('  Min/Max: %.0f / %.0f residues\n', min(trajDynamic), max(trajDynamic));
fprintf('==========================\n\n');