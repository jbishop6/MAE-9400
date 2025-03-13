%% For all shapes this demo implements matching such that M is smaler shape and N is larger shape. 

addpath(genpath('.')),
close all
clear all
clc


%% options
options = struct;
options.shot_num_bins = 15; % number of bins for shot increased from 10 to 15
options.shot_radius = 20; % Increased from 5 to 20

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
% 
i = '3311_surface'; 
j = '4022_surface';

% j = '5277_surface';
% j = '5273_surface';

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
disp(['M has ', num2str(size(M.VERT,1)), ' vertices']);
disp(['N has ', num2str(size(N.VERT,1)), ' vertices']);

disp(['Number of vertices in N: ', num2str(size(N.VERT,1))]);
disp(['Number of triangles in N: ', num2str(size(N.TRIV,1))]);




M.surface.X = M.VERT(:,1);
M.surface.Y = M.VERT(:,2);
M.surface.Z = M.VERT(:,3);
M.surface.VERT = M.VERT;
M.surface.TRIV = M.TRIV;
M.nf = size(M.TRIV,1);
M.nv = size(M.VERT,1);

[tform, N_aligned] = pcregistericp(pointCloud(N.VERT), pointCloud(M.VERT));
N.VERT = transformPointsForward(tform, N.VERT);


N.surface.X = N.VERT(:,1);
N.surface.Y = N.VERT(:,2);
N.surface.Z = N.VERT(:,3);
N.surface.VERT = N.VERT;
N.surface.TRIV = N.TRIV;
N.nf = size(N.TRIV,1);
N.nv = size(N.VERT,1);
options.constraint_map = ones(M.nv, 1); % Default to ones (or customize as needed)


% Compute shot descriptors

M_desc = calc_shot(M.VERT', M.TRIV', 1:M.n, options.shot_num_bins, options.shot_radius, 3)';
N_desc = calc_shot(N.VERT', N.TRIV', 1:N.n, options.shot_num_bins, options.shot_radius, 3)';


% Normalize descriptors
M_desc = M_desc ./ vecnorm(M_desc, 2, 2);
N_desc = N_desc ./ vecnorm(N_desc, 2, 2);


disp(['Size of M_desc: ', num2str(size(M_desc))]);
disp(['Size of N_desc: ', num2str(size(N_desc))]);

% Compute cosine similarity
similarity = M_desc * N_desc'; % Cosine similarity matrix
[~, idx] = max(similarity, [], 2); % Best match for each M vertex


k = 5; % Number of closest candidates per point
Mdl = KDTreeSearcher(N_desc); % Build KD-tree for fast searching
if size(M_desc,2) ~= size(N_desc,2) % Ensure feature dimensions match
    M_desc = M_desc';
    N_desc = N_desc';
end

[idx_all, D] = knnsearch(Mdl, M_desc, 'K', k); % Get distances too





% Compute Euclidean distance for each candidate and pick the best match
% Ensure idx_all is properly shaped before passing to pdist2
distances_matrix = zeros(size(M_desc,1), k); % Preallocate correct size
for l = 1:size(M_desc,1)
    distances_matrix(l, :) = pdist2(M_desc(l, :), N_desc(idx_all(l, :), :), 'euclidean');
end

distances_matrix = reshape(distances_matrix, size(M_desc,1), k); % Reshape to match M.n × k
[~, best_match] = min(distances_matrix, [], 2);


% Assign the best match from the K candidates
gt = [(1:M.n)', idx_all(sub2ind(size(idx_all), (1:M.n)', best_match))];
disp(['Size of gt: ', num2str(size(gt,1)), ' x ', num2str(size(gt,2))]);
disp('First 10 correspondences:');
disp(gt(1:10,:)); % Show first 10 mappings

num_repeats = size(gt,1) - length(unique(gt(:,2)));
disp(['Number of repeated matches in N before filtering: ', num2str(num_repeats)]);


% Compute distances between matched points
distances = sqrt(sum((M.VERT(gt(:,1),:) - N.VERT(gt(:,2),:)).^2, 2));

% Find median distance
median_dist = median(distances);

% Filter out matches that are too far from the median (e.g., 3x median)
valid_matches = distances < (2 * median_dist);
gt_filtered = gt(valid_matches, :);

num_repeats_filtered = size(gt_filtered,1) - length(unique(gt_filtered(:,2)));
disp(['Number of repeated matches in N after filtering: ', num2str(num_repeats_filtered)]);
disp('First 10 filtered correspondences:');
disp(gt_filtered(1:10,:));







% Checking the surfaces

figure;
subplot(1,2,1);
scatter3(M.VERT(gt(:,1),1), M.VERT(gt(:,1),2), M.VERT(gt(:,1),3), 30, 'r', 'filled');
title('Corresponding Points on M');
axis equal; view(3);

subplot(1,2,2);
scatter3(N.VERT(gt(:,2),1), N.VERT(gt(:,2),2), N.VERT(gt(:,2),3), 30, 'b', 'filled');
title('Corresponding Points on N');
axis equal; view(3);

figure;
patch('Faces', M.TRIV, 'Vertices', M.VERT, ...
    'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'black', 'EdgeAlpha', 0.3); % Light gray surface with black edges
title('Surface M with Wireframe');
axis equal; view(3); camlight; lighting gouraud;

figure;
patch('Faces', N.TRIV, 'Vertices', N.VERT, ...
    'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'black', 'EdgeAlpha', 0.3);
title('Surface N with Wireframe');
axis equal; view(3); camlight; lighting gouraud;

%% GEM
disp('---- GEM preprocess ----')
tic
fprintf('geodesic processing for M...'); 
M.distances = [];
vec = double(int32(linspace(n1/5,n1,5)));
%     disp(vec)
begining = 1;
for kk=1:length(vec)%1:n1
    ending = vec(kk);
    idx = begining:ending;
    begining = ending + 1;
    distances = perform_fast_marching_mesh(M.VERT', M.TRIV', idx, options.option1);
    M.distances = [M.distances, distances];
%         disp(size(M.distances))
end
fprintf('done \n'); 
fprintf('geodesic processing for N...'); 
N.distances = [];
vec = double(int32(linspace(n2/5,n2,5)));
begining = 1;
for kk=1:length(vec)%1:n2
   ending = vec(kk);
   idx = begining:ending;
   begining = ending + 1;
   distances = perform_fast_marching_mesh(N.VERT', N.TRIV', idx, options.option2);
   N.distances = [N.distances, distances];
end
fprintf('done \n');
toc


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
figure;
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
corr_DIR = DIR(strcat(i), strcat(j), options, corr_true);

toc
figure;
subplot(1,2,1); visualize_map_on_source(M, N, corr_DIR); title('Source');
subplot(1,2,2); visualize_map_on_target(M, N, corr_DIR); title('DIR')

%% GRAMPA
disp('------GRAMPA------')
tic
corr_GRAMPA = GRAMPA(M, N, options, 1); %eta 1 as in original paper
[corr_GRAMPA, ~] = find(corr_GRAMPA);
toc
figure;
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
figure;
subplot(1,2,1); visualize_map_on_source(M, N, corr_ZoomOut); title('Source');
subplot(1,2,2); visualize_map_on_target(M, N, corr_ZoomOut); title('ZoomOut')

%% getting curves
all_corr = cell(4,1);

% Ensure dimensions match before concatenation
gt_in = (1:min(size(M.VERT, 1), size(N.VERT, 1)))';


min_length = min([size(gt_in, 1), size(corr_GEM, 1), size(corr_GRAMPA, 1), size(corr_DIR, 1), size(corr_ZoomOut, 1)]);
gt_in = gt_in(1:min_length, :);
corr_GEM = corr_GEM(1:min_length, :);
corr_GRAMPA = corr_GRAMPA(1:min_length, :);
corr_DIR = corr_DIR(1:min_length, :);
corr_ZoomOut = corr_ZoomOut(1:min_length, :);




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

    curve = zeros(1, length(thresholds));
    for m = 1:length(thresholds)
        curve(m) = 100 * sum(errors <= thresholds(m)) / length(errors);
    end

    all_curves{i} = curve;
end
figure;
subplot(1,6,1); visualize_map_on_source(M, N, corr_GEM); title('Source');
subplot(1,6,2); visualize_map_on_target(M, N, corr_GEM); title('GEM')
subplot(1,6,3); visualize_map_on_target(M, N, corr_GRAMPA); title('GRAMPA');
subplot(1,6,4); visualize_map_on_target(M, N, corr_DIR); title('DIR')
subplot(1,6,5); visualize_map_on_target(M, N, corr_ZoomOut); title('ZoomOut')

% Plot all curves
figure;
plot(thresholds', mean(all_curves{1}, 1)', thresholds', mean(all_curves{2}, 1)', thresholds', mean(all_curves{3}, 1)', thresholds', mean(all_curves{4}, 1)');
line_width = 1.5;
hline = findobj(gcf, 'type', 'line');
set(hline, 'LineWidth', line_width);
legend({'2D-GEM', 'GRAMPA', 'DIR', 'ZoomOut'}, 'FontSize', 10);
