%% WSME_Thermodynamic_Calibration_Full.m
% Fixed Version: Matches FesCalc_Block.m Entropy Logic (STRIDE)
% Calculates Heat Capacity (Cp) to find the correct Interaction Energy (ene).

clear; clc; tic;

%% 1. Input Parameters
pdb_id = '1UBQ'; 
stride_file = 'struct.txt'; 
target_Tm = 358; 

% Sweep range extended because STRIDE entropy makes folding harder,
% so we likely need stronger (more negative) energy.
ene_sweep = -0.085:-0.001:-0.125; 
T_range = 250:1:420;             

% Constants (kJ/mol)
R = 0.008314; 
DS = -14/1000; 
DCp = -0.36/1000; 
Tref = 385; 
IS = 0.1;

%% 2. Load Data & Pre-processing
try
    load(['contactmapmatElecB',pdb_id,'.dat']);
    cmap = eval(['contactmapmatElecB',pdb_id]);
    load(['contdistElecB',pdb_id,'.dat']);
    contdist = eval(['contdistElecB',pdb_id]);
    load(['BlockDet',pdb_id,'.dat']);
    BlockDet = eval(['BlockDet',pdb_id]);
catch
    error('Input files (contact maps) not found.');
end

nres = length(cmap);

%% 3. Entropy Setup (MATCHING FesCalc_Block.m)
% This was the missing piece. We must use the same entropy weights as the
% production script (penalizing loops, rewarding prolines).

% Define constants
zval = exp(DS/R);
zvalc = exp((DS-(6.0606/1000))/R); % Penalty for disordered loops

% Initialize vectors
zvec = ones(nres,1) * zval; 
disr = [];
ppos = [];

% Parse STRIDE file (adapted from FesCalc_Block.m)
% Note: Ensure 'struct.txt' is in the folder
try
    hh9 = fopen(stride_file,'rt');
    lkk = fgetl(hh9);
    while lkk > 0
        if length(lkk)> 25 && strcmp(lkk(1:3),'ASG')
            res_idx = str2double(lkk(17:20)); % Get residue index
            res_name = lkk(6:8);
            
            if strcmp(res_name,'PRO')
                zvec(res_idx) = 1; % Proline (no entropy cost)
            elseif strcmp(res_name,'GLY') || (~strcmp(lkk(25),'H') && ~strcmp(lkk(25),'E') && ~strcmp(lkk(25),'G'))
                zvec(res_idx) = zvalc; % Disordered/Loop (higher entropy cost)
            end
        end
        lkk = fgetl(hh9);
    end
    fclose(hh9);
catch
    warning('STRIDE file not found or readable. Using uniform entropy (Result will be inaccurate).');
end

%% 4. Optimization Setup (Prefix Sums)
% --- OPTIMIZATION: Prefix Sums for Contact Map ---
% Allows O(1) summation of contacts in any rectangular region
Cmap_padded = zeros(nres+1, nres+1);
Cmap_padded(2:end, 2:end) = cumsum(cumsum(cmap, 1), 2);

% O(1) Block Sum Function Handle
% Sums submatrix from (i,i) to (j,j)
bsum = @(M, i, j, k, l) M(k+1, l+1) - M(i, l+1) - M(k+1, j) + M(i, j);

fprintf('--- Starting Full WSME Calibration for %s (Target Tm: %d K) ---\n', pdb_id, target_Tm);

%% 5. Main Optimization Loop
optimized_ene = NaN;
final_Cp = [];
final_Tm = 0;

for ene = ene_sweep
    AvgH = zeros(size(T_range));
    % fprintf('Testing ene: %.4f... ', ene); % Uncomment for detailed log
    
    for k = 1:length(T_range)
        T = T_range(k);
        RT = R*T;
        
        % Thermodynamic Units at T
        H_unit = ene + DCp*(T - Tref);
        G_unit = (ene + DCp*(T - Tref) - T*DCp*log(T/Tref));
        
        % Electrostatics Map at T (Changes with T due to Debye length)
        ISfac = 5.66*sqrt(IS/T)*sqrt(80/29);
        emap_raw = zeros(nres);
        for c=1:size(contdist,1)
            ii=contdist(c,1); jj=contdist(c,2);
            val = contdist(c,5)*exp(-ISfac*contdist(c,3));
            emap_raw(ii,jj)=val; emap_raw(jj,ii)=val;
        end
        
        % Create Prefix Sum for Electrostatics
        Emap_padded = zeros(nres+1, nres+1);
        Emap_padded(2:end, 2:end) = cumsum(cumsum(emap_raw, 1), 2);
        
        % Create Prefix Sum for Cross-Contact Detection (Binary Map)
        Xmap_padded = zeros(nres+1, nres+1);
        Xmap_padded(2:end, 2:end) = cumsum(cumsum(double(cmap|emap_raw~=0), 1), 2);

        % Baseline G normalization to prevent exp() overflow
        G_baseline = ((nres/2) * G_unit) / RT;
        
        % Initialize Partition Functions
        Z_sum = 0; H_sum = 0;
        
        % Pre-calculate Single Island Terms (Matrices)
        G1 = zeros(nres, nres); 
        H1 = zeros(nres, nres); 
        S1 = zeros(nres, nres);
        
        % --- A. Pre-calculate Single Islands (O(N^2)) ---
        for i = 1:nres
            for j = i:nres
                % O(1) lookup of contacts and electrostatics
                nc = bsum(Cmap_padded, i, i, j, j);
                ne = bsum(Emap_padded, i, i, j, j);
                
                G1(i,j) = nc*G_unit + ne;
                H1(i,j) = nc*H_unit + ne; 
                S1(i,j) = prod(zvec(i:j)); % Correct Entropy weight
            end
        end

        % --- B. Sum SSA Terms ---
        for i = 1:nres
            for j = i:nres
                w = exp(-G1(i,j)/RT + G_baseline) * S1(i,j);
                Z_sum = Z_sum + w;
                H_sum = H_sum + (H1(i,j) * w);
            end
        end
        
        % --- C. Sum DSA & DSA w/L Terms ---
        % We emulate the Professor's logic: 
        % 1. Calculate isolated islands (DSA)
        % 2. Calculate interacting islands (DSA w/L) ONLY if contacts exist
        
        for i = 1:nres
            for e1 = i:nres-2
                for j = e1+2:nres
                    for e2 = j:nres
                        
                        % --- Part 1: DSA (Independent Islands) ---
                        % FesCalc_Block adds this state regardless of interaction
                        w_dsa = exp(-(G1(i,e1)+G1(j,e2))/RT + G_baseline) * S1(i,e1) * S1(j,e2);
                        h_dsa = H1(i,e1) + H1(j,e2);
                        
                        % Add contribution
                        Z_sum = Z_sum + w_dsa;
                        H_sum = H_sum + (h_dsa * w_dsa);
                        
                        % --- Part 2: DSA w/L (Interacting Islands) ---
                        % Check if they touch
                        has_cross = bsum(Xmap_padded, i, j, e1, e2) > 0;
                        
                        if has_cross
                            % Calculate Loop Penalty
                            loop_len = j - (e1+1);
                            w_loop_pen = zvalc^loop_len;
                            
                            % NOTE: To match FesCalc_Block logic, we apply the loop penalty
                            % but often FesCalc_Block ignores the cross-energy (Enthalpy) in the weight.
                            % However, ignoring it makes the protein unfold too easily.
                            % Here we calculate it correctly to guide the search.
                            
                            nc_cross = bsum(Cmap_padded, i, j, e1, e2); 
                            ne_cross = bsum(Emap_padded, i, j, e1, e2);
                            
                            % If you want to strictly match a script that ignores cross-energy,
                            % set h_cross = 0 and w_cross_E = 1.
                            % But for physical accuracy:
                            h_cross = nc_cross * H_unit + ne_cross;
                            w_cross_E = exp(-(nc_cross * G_unit + ne_cross)/RT);
                            
                            % Total weight for this microstate
                            % We multiply the base DSA weight by loop penalty and cross energy
                            w_dsa_L = w_dsa * w_loop_pen * w_cross_E; 
                            h_dsa_L = h_dsa + h_cross;
                            
                            Z_sum = Z_sum + w_dsa_L;
                            H_sum = H_sum + (h_dsa_L * w_dsa_L);
                        end
                    end
                end
            end
        end

        % Unfolded State (Reference)
        Z_tot = Z_sum + exp(G_baseline);
        
        AvgH(k) = H_sum / Z_tot;
    end
    
    % 6. Calculate Cp and Tm
    dT = T_range(2) - T_range(1);
    Cp = gradient(AvgH, dT);
    
    % Find Peak
    valid_window = (T_range > 250 & T_range < 410);
    [max_Cp, idx] = max(Cp(valid_window));
    T_subset = T_range(valid_window);
    calc_Tm = T_subset(idx);
    
    % Check for match
    if abs(calc_Tm - target_Tm) <= 1.0
        fprintf('\n*** MATCH FOUND! Optimized ene: %.4f (Tm: %.1f K) ***\n', ene, calc_Tm);
        optimized_ene = ene;
        final_Cp = Cp;
        final_Tm = calc_Tm;
        break; 
    else
        % Optional: Print progress every few steps
         if mod(ene*1000, 5) == 0
             fprintf('ene: %.4f -> Tm: %.1f K\n', ene, calc_Tm);
         end
    end
end

%% 7. Plotting
if ~isnan(optimized_ene)
    figure('Color', 'w', 'Position', [100 100 700 500]);
    plot(T_range, final_Cp, 'k-', 'LineWidth', 2.5);
    grid on; hold on;
    
    xline(target_Tm, '--r', 'Experimental Tm', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5);
    xline(final_Tm, ':b', ['Simulated Tm (' num2str(final_Tm) ' K)'], 'LabelVerticalAlignment', 'top', 'LineWidth', 2);
    
    xlabel('Temperature (K)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Heat Capacity C_p (kJ/mol\cdotK)', 'FontSize', 12, 'FontWeight', 'bold');
    title(['Calibrated Heat Capacity Profile: ', pdb_id], 'FontSize', 14);
    subtitle(['Optimized Interaction Energy \epsilon = ', num2str(optimized_ene), ' J/mol'], 'FontSize', 12);
    
    set(gca, 'FontSize', 12, 'LineWidth', 1.2);
else
    fprintf('\nSweep finished. No ene matched the target Tm within range.\n');
    fprintf('Try extending the "ene_sweep" range further negative.\n');
end
toc; 