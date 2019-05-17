% RUN_MULTIPLE_EXPERIMENTS is a script for running Experiments 3-5
% repeatedly for different types of matrices to produce the plots in the
% supplementary material of our paper. The results are saved in mat files
% so that they can later be read and processed. Note that the appropriate
% lines need to be commented out in the scripts experiment3.m,
% experiment4.m and experiment5.m before running this code, since otherwise
% the random_seed and mat_type chosen in this script will just be
% overwritten by whatever is used in those scripts.

% Settings
random_seeds = [2 3 4 5 6];
mat_types = {'normal', 'uniform', 'adversarial_1', 'adversarial_2'};

% Determine number of experiments and number of matrix types
no_experiment = length(random_seeds);
no_mat_types = length(mat_types);

% Run Experiment 3
for mat = 1:no_mat_types
    mat_type = mat_types{mat};
    for expr = 1:no_experiment
        random_seed = random_seeds(expr);
        fprintf('Running Experiment 3. Matrix type: %s. Run: %d\n', mat_type, expr);
        experiment3;
        save(['experiment3-mat_type-', mat_type, '-run-', num2str(expr)]);
        clearvars -except random_seeds no_experiment mat_types no_mat_types mat mat_type expr random_seed
        fprintf('Finished Experiment 3. Matrix type: %s. Run: %d\n', mat_type, expr);
    end
end

% Run Experiment 4
for mat = 1:no_mat_types
    mat_type = mat_types{mat};
    for expr = 1:no_experiment
        fprintf('Running Experiment 4. Matrix type: %s. Run: %d\n', mat_type, expr);
        experiment4;
        save(['experiment4-mat_type-', mat_type, '-run-', num2str(expr)]);
        clearvars -except random_seeds no_experiment mat_types no_mat_types mat mat_type expr random_seed
        fprintf('Finished Experiment 4. Matrix type: %s. Run: %d\n', mat_type, expr);
    end
end

% Run Experiment 5
for mat = 1:no_mat_types
    mat_type = mat_types{mat};
    for expr = 1:no_experiment
        fprintf('Running Experiment 5. Matrix type: %s. Run: %d\n', mat_type, expr);
        experiment5;
        save(['experiment5-mat_type-', mat_type, '-run-', num2str(expr)]);
        clearvars -except random_seeds no_experiment mat_types no_mat_types mat mat_type expr random_seed
        fprintf('Finished Experiment 5. Matrix type: %s. Run: %d\n', mat_type, expr);
    end
end