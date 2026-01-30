%% Script to compute free energies with different residue masking
% This script iterates through different residues to mask their charges
% and saves the resulting free energy profiles

clear; clc;

%% Parameters
T = 310; % Temperature in Kelvin

% Define which residues to mask (list of residue indices)
To_mask = [6,11,16,18,24,32,33,34,39,42,48,51,52,54,58,63,64,72,74];

% Create output directory if it doesn't exist
output_dir = 'Free_Energies';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Wild Type

fprintf('Processing Wild Type...\n');
    
% Call the function with single residue masked
tic;
fes = get_Free_Energy(T, []);
elapsed_time = toc;

fprintf('  Completed in %.2f seconds\n', elapsed_time);

% Save the free energy profile
filename = fullfile(output_dir, sprintf('WT.txt'));

% Write header and data
fid = fopen(filename, 'w');
fprintf(fid, '%% Free Energy Profile\n');
fprintf(fid, '%% Temperature: %d K\n', T);
fprintf(fid, '%% Wild Type');
fprintf(fid, '%% Column 1: Number of Structured Blocks\n');
fprintf(fid, '%% Column 2: Free Energy (kJ/mol)\n');
fprintf(fid, '%10s %15s\n', 'N_blocks', 'FE_kJ/mol');

for j = 1:size(fes, 1)
    fprintf(fid, '%10d %15.6f\n', fes(j, 1), fes(j, 2));
end

fclose(fid);

fprintf('  Saved to: %s\n\n', filename);

%% Iterate through each residue to mask
fprintf('Starting free energy calculations...\n');
fprintf('Temperature: %d K\n', T);
fprintf('Total residues to mask: %d\n\n', length(To_mask));

for i = 1:length(To_mask)
    res_to_mask = To_mask(i);
    
    fprintf('Processing residue %d (%d/%d)...\n', res_to_mask, i, length(To_mask));
    
    % Call the function with single residue masked
    tic;
    fes = get_Free_Energy(T, res_to_mask);
    elapsed_time = toc;
    
    fprintf('  Completed in %.2f seconds\n', elapsed_time);
    
    % Save the free energy profile
    filename = fullfile(output_dir, sprintf('Mutant_%d.txt', res_to_mask));
    
    % Write header and data
    fid = fopen(filename, 'w');
    fprintf(fid, '%% Free Energy Profile\n');
    fprintf(fid, '%% Temperature: %d K\n', T);
    fprintf(fid, '%% Masked Residue: %d\n', res_to_mask);
    fprintf(fid, '%% Column 1: Number of Structured Blocks\n');
    fprintf(fid, '%% Column 2: Free Energy (kJ/mol)\n');
    fprintf(fid, '%10s %15s\n', 'N_blocks', 'FE_kJ/mol');
    
    for j = 1:size(fes, 1)
        fprintf(fid, '%10d %15.6f\n', fes(j, 1), fes(j, 2));
    end
    
    fclose(fid);
    
    fprintf('  Saved to: %s\n\n', filename);
end

fprintf('All calculations completed!\n');
fprintf('Results saved in directory: %s\n', output_dir);