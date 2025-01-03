%% For all shapes this demo implements matching such that M is smaler shape and N is larger shape. 

addpath(genpath('.')),
close all
clear all
clc


%% options
options = struct;
options.shot_num_bins = 10; % number of bins for shot
options.shot_radius = 5;

options.shuffle = false; % leave false as not implemented on different shapes for DIR
options.isometric = false;


options.option1.nb_iter_max = 30; % iteration number needed for fast marching %% this should be smallest
options.option2.nb_iter_max = 120; %% this largest
options.nb_iter_max = inf;
options.num_eigenvalues = 100;  % Specify the number of eigenvalues (adjust as needed)
options.power = 0.5;  % This is just an example, adjust based on your needs


%% loading... for all pairs this demo implements matching larger to smaller shape.
name ='';
% 
%% %%For Isometric uncoment any of the isometric sample pairs below as well as
%% the gt_in line and lines 40-42
% % i = 'mesh023';
% % j = 'mesh071';
% % gt_in = [1:12500]';
% % 
% i ='mesh03';
% j = 'mesh050';
% gt_in = [1:12500]';

% % i ='wolf1';
% % j = 'wolf2';
% % gt_in = [1:4344]';
% 
% options.isometric = true;
% gt = [gt_in, gt_in];
% disp(size(gt))

%% %%for non_isometric uncomment any of the pairs below and the gt_in line
%% as well as lines 58-61
% i = 'kid16'; 
% j = 'kid17';
% gt_in = [1:10988]';%11292

i = '3311_surface'; 
% j = '4022_surface';

% i = '5277_surface';
j = '5273_surface';

gt_in = [1:8515]';%11292

% i = 'kid19'; 
% j = 'kid20';
% gt_in = [1:8515]';%11292

filename1 = strcat(i, '.off');
filename2 = strcat(j, '.off');


check_off_file(filename1)
check_off_file(filename2)

gt_M_null = read_correspondence(strcat(name, i, '.mat'));% load
gt_N_null = read_correspondence(strcat(name, j, '.mat')); %load

disp(size(gt_M_null));
disp(size(gt_N_null));
gt = merge_ground_truth(gt_M_null, gt_N_null); % merge
disp(size(gt));

%% specs for GEM & DIR
if ~(options.isometric)
    options.maxIter = 60;
    options.spec_dim = 10;
    options.th = [100];
else
    options.maxIter = 10;
    options.spec_dim = 500;    
    low = 0.14; % local distortion lower bound
    gap = 5; % local distortion gaps
    options.th = 0.6-(0.5-low)/gap:-(0.5-low)/gap:low; % local distortion
    disp(options.th)
end



[N, M, n1, n2, diameters, corr_true] = load_and_preprocess(i, j,'', '', '', '', options);
M.surface.X = M.VERT(:,1);
M.surface.Y = M.VERT(:,2);
M.surface.Z = M.VERT(:,3);
M.surface.VERT = M.VERT;
M.surface.TRIV = M.TRIV;
M.nf = size(M.TRIV,1);
M.nv = size(M.VERT,1);

N.surface.X = N.VERT(:,1);
N.surface.Y = N.VERT(:,2);
N.surface.Z = N.VERT(:,3);
N.surface.VERT = N.VERT;
N.surface.TRIV = N.TRIV;
N.nf = size(N.TRIV,1);
N.nv = size(N.VERT,1);
options.constraint_map = ones(M.nv, 1); % Default to ones (or customize as needed)
%% GEM
% Get number of vertices
nverts = size(M.VERT, 1);  % or size(N.VERT, 1), depending on which mesh you are using

% Initialize options structure
options.option1.nb_iter_max = 30;
options.option2.nb_iter_max = 120;
options.nb_iter_max = inf;

% Set L (constraint map) if not already provided
if ~isfield(options, 'L')
    options.L = ones(nverts, 1);  % Initialize L with ones (or any default value)
end

% Check for NaN or Inf in L
if any(isnan(options.L(:))) || any(isinf(options.L(:)))
    disp('L (constraint map) contains NaN or Inf values.');
    options.L(isnan(options.L) | isinf(options.L)) = 0;  % Replace NaN/Inf with 0
end

% Define end_points explicitly if not provided in the options
if ~isfield(options.option1, 'end_points') || isempty(options.option1.end_points)
    options.option1.end_points = 1:size(M.VERT, 1);  % Use all vertices as end points
end

% Check M.VERT and M.TRIV for NaNs or Infs
if any(isnan(M.VERT(:))) || any(isinf(M.VERT(:)))
    error('M.VERT contains NaN or Inf values!');
end

if any(isnan(M.TRIV(:))) || any(isinf(M.TRIV(:)))
    error('M.TRIV contains NaN or Inf values!');
end

disp('M.VERT and M.TRIV are valid.');



%% Debugging
% Check for NaN or Inf in L 
if any(isnan(options.L(:))) || any(isinf(options.L(:))) 
    disp('L (constraint map) contains NaN or Inf values.'); 
    options.L(isnan(options.L) | isinf(options.L)) = 0; % Replace NaN/Inf with 0 
end 

% Define end_points explicitly if not provided in the options 
if ~isfield(options.option1, 'end_points') || isempty(options.option1.end_points) 
    options.option1.end_points = 1:size(M.VERT, 1); % Use all vertices as end points 
end 

% Check M.VERT and M.TRIV for NaNs or Infs 
assert(all(~isnan(M.VERT(:))), 'M.VERT contains NaN values'); 
assert(all(~isinf(M.VERT(:))), 'M.VERT contains Inf values'); 
assert(all(~isnan(M.TRIV(:))), 'M.TRIV contains NaN values'); 
assert(all(~isinf(M.TRIV(:))), 'M.TRIV contains Inf values'); 
assert(all(~isnan(N.VERT(:))), 'N.VERT contains NaN values'); 
assert(all(~isinf(N.VERT(:))), 'N.VERT contains Inf values'); 
assert(all(~isnan(N.TRIV(:))), 'N.TRIV contains NaN values'); 
assert(all(~isinf(N.TRIV(:))), 'N.TRIV contains Inf values'); 

disp('M.VERT, M.TRIV, N.VERT, and N.TRIV are valid.');

%% Parallel Processing 
% Profile the code execution
profile on;

% Divide the vertices of M and N into chunks for processing
chunk_size = 75;

% Define chunk size
num_chunks_M = ceil(n1 / chunk_size);
num_chunks_N = ceil(n2 / chunk_size);

% Initialize storage for distances
distances_M = cell(1, num_chunks_M);
distances_N = cell(1, num_chunks_N);

% Ensure any existing parallel pool is deleted
delete(gcp('nocreate'));

% Start parallel pool with 12 workers  
num_workers = 12; 
parpool('local', num_workers);

% Process chunks for M in parallel using parfor
parfor kk = 1:num_chunks_M
    % Define chunk range
    start_idx = (kk - 1) * chunk_size + 1;
    end_idx = min(kk * chunk_size, n1);

    % Extract chunk of start points
    start_points_chunk = start_idx:end_idx;

    % Debugging: Print progress
    fprintf('Processing chunk %d of M: start = %d, end = %d\n', kk, start_idx, end_idx);
    try
        % Call the fast marching function to calculate geodesic distances for M
        distances_chunk = perform_fast_marching_mesh(M.VERT', M.TRIV', start_points_chunk, options.option1);

        % Store distances
        distances_M{kk} = distances_chunk;
        fprintf('Finished chunk %d for M\n', kk);
    catch ME
        fprintf('Error occurred at chunk %d of M: %s\n', kk, ME.message);
        disp(ME.stack);
    end
end

% Combine all chunk results into a single matrix for M
M.distances = cell2mat(distances_M);

% Process chunks for N in parallel using parfor
parfor kk = 1:num_chunks_N
    % Define chunk range
    start_idx = (kk - 1) * chunk_size + 1;
    end_idx = min(kk * chunk_size, n2);

    % Extract chunk of start points
    start_points_chunk = start_idx:end_idx;

    % Debugging: Print progress
    fprintf('Processing chunk %d of N: start = %d, end = %d\n', kk, start_idx, end_idx);
    try
        % Call the fast marching function to calculate geodesic distances for N
        distances_chunk = perform_fast_marching_mesh(N.VERT', N.TRIV', start_points_chunk, options.option2);

        % Store distances
        distances_N{kk} = distances_chunk;
        fprintf('Finished chunk %d for N\n', kk);
    catch ME
        fprintf('Error occurred at chunk %d of N: %s\n', kk, ME.message);
        disp(ME.stack);
    end
end

% Combine all chunk results into a single matrix for N
N.distances = cell2mat(distances_N);

% Shut down the parallel pool
delete(gcp('nocreate'));

% Stop profiling and view results
profile off;
profile viewer;

disp('Geodesic processing for M and N done');


%% geodesic distances for plot
fprintf('geodesic processing Full N for plot... ');
tic
distances = geodesic_distance(N.TRIV,N.VERT); %Added by JX  
%     distances = perform_fast_marching_mesh(N.VERT', N.TRIV', 1:n2, options);
distances = sparse(distances);
fprintf('done \n')
toc

%% GEM
disp('------GEM------')
tic
% Assume M, N, corr_true, and options are already loaded
corr_GEM = GEM(M, N, corr_true, options);

toc
[corr_GEM, ~] = find(corr_GEM);
% corr = [gt_in, corr_GEM];
figure(1);
subplot(1,2,1); visualize_map_on_source(M, N, corr_GEM); title('Source');
subplot(1,2,2); visualize_map_on_target(M, N, corr_GEM); title('2D-GEM')

%% DIR

options.spec_dim = 500;
options.spec_dim_cut = 460;
if ~(options.isometric)
    options.th = [100];    
end
disp('------DIR------')
tic
corr_DIR = DIR(strcat(i, '.off'), strcat(j, '.off'), options, corr_true);
toc
figure(2);
subplot(1,2,1); visualize_map_on_source(M, N, corr_DIR); title('Source');
subplot(1,2,2); visualize_map_on_target(M, N, corr_DIR); title('DIR')

%% GRAMPA
disp('------GRAMPA------')
tic
corr_GRAMPA = GRAMPA(M, N, options, 1); %eta 1 as in original paper
[corr_GRAMPA, ~] = find(corr_GRAMPA);
toc
figure(3)
subplot(1,2,1); visualize_map_on_source(M, N, corr_GRAMPA); title('Source');
subplot(1,2,2); visualize_map_on_target(M, N, corr_GRAMPA); title('GRAMPA')

%% ZoomOut
options.spec_dim = 200;
disp('------ZoomOut------')
options.k_init = 20;
options.k_step = 1;
options.k_final = 200; % as in original paper from 4 or 20 to a max of 200 
tic
[corr_ZoomOut, ~, ~, ~] = zoomOut_refine(M, N, options);
toc
figure(4)
subplot(1,2,1); visualize_map_on_source(M, N, corr_ZoomOut); title('Source');
subplot(1,2,2); visualize_map_on_target(M, N, corr_ZoomOut); title('ZoomOut')

%% getting curves
all_corr = cell(4,1);

% Ensure dimensions match before concatenation
min_length = min([size(gt_in, 1), size(corr_GEM, 1), size(corr_GRAMPA, 1), size(corr_DIR, 1), size(corr_ZoomOut, 1)]);
gt_in = gt_in(1:min_length, :);
corr_GEM = corr_GEM(1:min_length, :);
corr_GRAMPA = corr_GRAMPA(1:min_length, :);
corr_DIR = corr_DIR(1:min_length, :);
corr_ZoomOut = corr_ZoomOut(1:min_length, :);

% Debugging check
disp('Adjusted size of gt_in:');
disp(size(gt_in));
disp('Adjusted size of corr_GEM:');
disp(size(corr_GEM));
disp('Adjusted size of corr_GRAMPA:');
disp(size(corr_GRAMPA));
disp('Adjusted size of corr_DIR:');
disp(size(corr_DIR));
disp('Adjusted size of corr_ZoomOut:');
disp(size(corr_ZoomOut));

% Perform the concatenation
all_corr{1} = [gt_in, corr_GEM];
all_corr{2} = [gt_in, corr_GRAMPA];
all_corr{3} = [gt_in, corr_DIR];
all_corr{4} = [gt_in, corr_ZoomOut];

all_curves = cell(4, 1);
for i = 1:4
    corr = cell2mat(all_corr(i));
    disp(['Size of corr ' num2str(i) ':']);
    disp(size(corr));
    
    errors = zeros(size(corr, 1), 1);
    for m = 1:size(corr, 1)
        gt_match = gt(gt(:, 1) == corr(m, 1), 2);
        match = corr(m, 2);

        if ~isempty(gt_match)
            % Using the geodesic distance of the second graph
            errors(m) = distances(gt_match, match); % TODO include your geodesics here
        else
            errors(m) = 200;
        end
    end
    thresholds = 0:0.01:0.25;
    errors = errors / diameters;
    
    % Debugging: Print error values
    disp(['Errors for corr ' num2str(i) ':']);
    disp(errors);

    curve = zeros(1, length(thresholds));
    for m = 1:length(thresholds)
        curve(m) = 100 * sum(errors <= thresholds(m)) / length(errors);
    end
    
    % Debugging: Print curve values
    disp(['Curve for corr ' num2str(i) ':']);
    disp(curve);

    all_curves{i} = curve;
end
figure(5);
subplot(1,6,1); visualize_map_on_source(M, N, corr_GEM); title('Source');
subplot(1,6,2); visualize_map_on_target(M, N, corr_GEM); title('GEM')
subplot(1,6,3); visualize_map_on_target(M, N, corr_GRAMPA); title('GRAMPA');
subplot(1,6,4); visualize_map_on_target(M, N, corr_DIR); title('DIR')
subplot(1,6,5); visualize_map_on_target(M, N, corr_ZoomOut); title('ZoomOut')

% Plot all curves
figure(7);
plot(thresholds', mean(all_curves{1}, 1)', thresholds', mean(all_curves{2}, 1)', thresholds', mean(all_curves{3}, 1)', thresholds', mean(all_curves{4}, 1)');
% ylim([0 100]);
line_width = 1.5;
hline = findobj(gcf, 'type', 'line');
set(hline, 'LineWidth', line_width);
legend({'2D-GEM', 'GRAMPA', 'DIR', 'ZoomOut'}, 'FontSize', 10);
