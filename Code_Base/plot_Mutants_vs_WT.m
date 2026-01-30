%% Script to plot Free Energy vs Folded Residues for all mutants vs WT
% This script reads all Mutant_X.txt files and plots each against WT.txt

clear; clc; close all;

%% Parameters
input_dir = 'Free_Energies'; % Directory containing the data files
output_dir = 'Plots'; % Directory to save plots

% Create output directory if it doesn't exist
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Load wildtype data
wt_file = fullfile(input_dir, 'WT.txt');
if ~exist(wt_file, 'file')
    error('Wildtype file not found: %s', wt_file);
end

% Read WT data (skip header lines starting with %)
wt_data = load(wt_file);
wt_nres = wt_data(:, 1);
wt_fe = wt_data(:, 2);

fprintf('Loaded wildtype data: %s\n', wt_file);
fprintf('  Number of data points: %d\n\n', length(wt_nres));

%% Find all mutant files
file_list = dir(fullfile(input_dir, 'Mutant_*.txt'));
num_mutants = length(file_list);

if num_mutants == 0
    error('No mutant files found in directory: %s', input_dir);
end

fprintf('Found %d mutant files\n\n', num_mutants);

%% Process each mutant
for i = 1:num_mutants
    mutant_file = fullfile(input_dir, file_list(i).name);
    
    % Extract residue number from filename
    tokens = regexp(file_list(i).name, 'Mutant_(\d+)\.txt', 'tokens');
    if isempty(tokens)
        warning('Could not parse filename: %s', file_list(i).name);
        continue;
    end
    res_num = str2double(tokens{1}{1});
    
    fprintf('Processing %s (Residue %d)...\n', file_list(i).name, res_num);
    
    % Load mutant data
    mutant_data = load(mutant_file);
    mutant_nres = mutant_data(:, 1);
    mutant_fe = mutant_data(:, 2);
    
    % Create figure
    fig = figure('Position', [100, 100, 800, 600]);
    hold on;
    
    % Plot WT
    plot(wt_nres, wt_fe, 'k-', 'LineWidth', 2, 'DisplayName', 'WT');
    
    % Plot Mutant
    plot(mutant_nres, mutant_fe, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('Mutant %d', res_num));
    
    % Formatting
    xlabel('Number of Folded Residues', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Free Energy (kJ/mol)', 'FontSize', 14, 'FontWeight', 'bold');
    title(sprintf('Free Energy Profile: Mutant %d vs WT', res_num), 'FontSize', 16, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 12);
    grid on;
    set(gca, 'FontSize', 12, 'LineWidth', 1.5);
    box on;
    
    % Save figure
    output_file_png = fullfile(output_dir, sprintf('FE_Mutant_%d_vs_WT.png', res_num));
    output_file_fig = fullfile(output_dir, sprintf('FE_Mutant_%d_vs_WT.fig', res_num));
    
    saveas(fig, output_file_png);
    saveas(fig, output_file_fig);
    
    fprintf('  Saved: %s\n', output_file_png);
    fprintf('  Saved: %s\n\n', output_file_fig);
    
    close(fig);
end

fprintf('All plots completed!\n');
fprintf('Plots saved in directory: %s\n', output_dir);