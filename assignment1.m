%%Assignment 1

%% Part 1: load & plot
load IOdownload.mat

swe = io.swe2000;
idn = io.idn2000;

n = length(swe);

G_swe = digraph(swe, name);
G_idn = digraph(idn, name);

% figure
% plot(G_swe)
% title('Graph Sweden')
% 
% figure
% plot(G_idn)
% title('Graph Indonesia')

%% 1a

%Sweden
in_deg_swe = swe'*ones(n,1);
out_deg_swe = swe*ones(n,1);

[sort_in, idx_in] = sort(abs(in_deg_swe), 'descend');
[sort_out, idx_out] = sort(abs(out_deg_swe), 'descend');


for i = 1:3
    fprintf('Sweden largest sector number %.2f:\n', i)
    fprintf('Factor in-centrality: %s\n', name{idx_in(i)})
    fprintf('Factor out-centrality: %s\n', name{idx_out(i)})
end

%Indonesia
in_deg_idn = idn'*ones(n,1);
out_deg_idn = idn*ones(n,1);

[sort_in, idx_in] = sort(abs(in_deg_idn), 'descend');
[sort_out, idx_out] = sort(abs(out_deg_idn), 'descend');


fprintf('\n')

for i = 1:3
    fprintf('Indonesia largest sector number %.2f:\n', i)
    fprintf('Factor in-centrality: %s\n', name{idx_in(i)})
    fprintf('Factor out-centrality: %s\n', name{idx_out(i)})
end

%% 1b


%Sweden 
%find nodes
[BINS, BINSIZE] = conncomp(G_swe);
[~, max_comp] = max(BINSIZE); %%max_comp is the largest component
nodes = find(BINS == max_comp); %find all indices of all nodes associated with largest component

%create component
comp_swe = swe(nodes, nodes); %create matrix for component
comp_name_swe = name(nodes);

%centrality 
[V, D] = eig(comp_swe');
[eig_max, idx_max]= max(D,[],'all', 'linear');
z_swe = V(:,idx_max);

z_swe = z_swe / sum(z_swe);
%sort
[z_swe_sort, idx_swe] = sort(z_swe, 'descend');

%print
for i = 1:3
    fprintf('Largest eigenvector sector number for %.2f:\n', i)
    fprintf('Factor Sweden: %s\n', comp_name_swe{idx_swe(i)})
end

%correct answer: vehicles, radio and other

%% 
% Indonesia
% Find nodes
[BINS, BINSIZE] = conncomp(G_idn);
[~, max_comp] = max(BINSIZE); % Max_comp is the largest component
nodes = find(BINS == max_comp); % Find all indices of all nodes associated with the largest component

% Create component
comp_idn = idn(nodes, nodes); % Create matrix for component
comp_name_idn = name(nodes);

% Centrality 
[V, D] = eig(comp_idn');
[eig_max, idx_max]= max(D,[],'all', 'linear');
z_idn = V(:,idx_max);

z_idn = z_idn / sum(z_idn);

% Sort
[z_idn_sort, idx_idn] = sort(z_idn, 'descend');

% Print
for i = 1:3
    fprintf('Largest eigenvector sector number for %.2f:\n', i)
    fprintf('Factor Indonesia: %s\n', comp_name_idn{idx_idn(i)})
end


%correct answer: Food, Hotels & Agriculture




%% plots
%plot graph
% subG_idn = digraph(comp_idn, comp_name);
% 
% figure
% plot(subG_idn)
% title('Graph of largest component Indonesia')


% %plot graph
% subG_swe = digraph(comp_swe, comp_name_swe);
% 
% figure
% plot(subG_swe)
% title('Graph of largest component Sweden')

%% 1c
beta = 0.15;
mu_1 = ones(n,1);
mu_2 = zeros(n,1);
mu_2(31) = 1; %31 = Wholesale & retail trade; repairs


%Sweden
[V, D] = eigs(swe);
lambda_w = D(1,1);
A = (eye(n) - (1 - beta)/lambda_w * swe');

z1_swe = A \ (beta*mu_1);
z2_swe = A \ (beta*mu_2);


[sort_in, idx_1] = sort(z1_swe, 'descend');
[sort_out, idx_2] = sort(z2_swe, 'descend');


for i = 1:3
    fprintf('Sweden largest sector number %.2f:\n', i)
    fprintf('Factor centrality mu_1: %s\n', name{idx_1(i)})
    fprintf('Factor centrality mu_2: %s\n', name{idx_2(i)})
end
fprintf('\n')

%Indonesia
[V, D] = eigs(idn);
lambda_w = D(1,1);
A = (eye(n) - (1 - beta)/lambda_w * idn');

z1_idn = A \ (beta*mu_1);
z2_idn = A \ (beta*mu_2);


[sort_in, idx_1] = sort(z1_idn, 'descend');
[sort_out, idx_2] = sort(z2_idn, 'descend');

for i = 1:3
    fprintf('Indonesia largest sector number %.2f:\n', i)
    fprintf('Factor centrality mu_1: %s\n', name{idx_1(i)})
    fprintf('Factor centrality mu_2: %s\n', name{idx_2(i)})
end


%% Part 2: load & matrix creation
load('twitter.mat', '-ascii');
load('users.mat', '-ascii')
W = spconvert(twitter);

[n, ~] = size(users); %number of users
W(1, n) = 0; %fixing size, assuming last nodes found are sources

%adding self-loops
%this part might be unnecessary since we know that it's always first node
%that has no out-neighbours, still might be good to have general method
%though

w = W*ones(n,1); %out-neighbours 
sinks = find(w == 0); %finding nodes without out neighbours

for k = sinks
    W(k, k) = 1; %adding self-loops to sinks
end

P = diag(sum(W,2)) \ W; %creation of stochastic matrix

%% 2a

%PageRank centrality
z = zeros(n,1);
Ps = sparse(eye(n));
% Parameters
beta = 0.15;

mu = ones(n,1);
k = 0;
dz = beta*mu;
z = z + dz;
while norm(dz) > 1e-6
    k = k + 1;
    Ps = P'*Ps;
    dz = beta*(1-beta)^k*Ps*mu;
    z = z + dz;
end

disp('Five largest PageRank centralities: ')
[sort, nodes] = sort(z, 'descend');

for i = 1:5
    nodes(i)
end

%nodes with highest page rank are 1, 2, 112, 9, 26

%% 2b-c: Simulation

% Stubborn nodes
% Number of iterations
niter = 3000;

% Stubborn and regular nodes
stubborn = [80 1];
regular = setdiff(1:n, stubborn);

% Input to the stubborn nodes
u = [0 1]';

% Submatrices
Q = P(regular, regular);
E = P(regular, stubborn);
x = zeros(n,niter);

%initial values
x(:,1) = 0.5*ones(n,1);
x(stubborn,1) = u;

for i = 2:niter
    x(regular, i) = Q*x(regular, i-1) + E*x(stubborn, i-1);
    x(stubborn, i) = x(stubborn, i-1);
end

%% 2b

nodes = [stubborn 77 99 40 112];
plot(x(nodes , :)')

title('Simulation of opinion dynamics')
xlabel('Iteration')
ylabel('Opinion')
legend(cellstr(num2str(nodes')))

%% 2c investigation of nodes
edges = 0:0.1:1;

histogram(x(:, niter), edges, 'Normalization', 'probability')
title_str = sprintf('Distribution of opinions with stubborn nodes [%s]'' with values u = [0 1].', num2str(stubborn));
title(title_str)
