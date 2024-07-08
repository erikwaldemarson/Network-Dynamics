%% 1.1

n = 100;
W = zeros(n);
W = W + diag(ones(n-1,1),1); % add ones on the +1 off-diagonal
W = W + diag(ones(n-1,1),-1); % add ones on the -1 off-diagonal
W = W + diag(ones(n-2,1),2); % add ones on the +2 off-diagonal
W = W + diag(ones(n-2,1),-2); % add ones on the -2 off-diagonal
W = W + diag(ones(1,1),n-1); % add ones on the +n-1 off-diagonal
W = W + diag(ones(1,1),1-n); % add ones on the -n+1 off-diagonal
W = W + diag(ones(2,1),n-2); % add ones on the +n-2 off-diagonal
W = W + diag(ones(2,1),2-n); % add ones on the -n+2 off-diagonal
W = sparse(W); % transform it into a sparse matrix
G = graph(W); % convert W into a graph (might not work on all MATLAB versions.
%plot(G) % probably best to do on a smaller graph...

S = 0;
I = 1;
R = 2;

beta = 0.3;
rho = 0.7;

T = 15; 
N = 100;

mean_infected = zeros(T+1,1);
mean_newly_infected = zeros(T+1,1);
mean_susceptible = zeros(T+1,1);
mean_recovered = zeros(T+1,1);


for k = 1:N
    %initial X
    X = zeros(n,1);
    nr_initial_inf = 10;
    X(randi(n, [nr_initial_inf,1])) = I; 

    mean_infected(1) = mean_infected(1) + nr_initial_inf;
    mean_susceptible(1) = mean_susceptible(1) + n - nr_initial_inf;

    for t = 1:T
        newly_infected = 0;
        i_infected = X == I;
        nr_inf_neighbours = W*i_infected;
        
        susceptible = find(X == S);
        for i = susceptible'
            m = nr_inf_neighbours(i);
            prob = 1 - (1-beta)^m;
            eps = rand();
            
            if prob >= eps
                X(i) = I;
                newly_infected = newly_infected + 1;
            end
        end
    
        infected = find(i_infected);
        for i = infected'
            prob = rho;
            eps = rand();
            
            if prob >= eps
                X(i) = R;
            end
        end

        mean_newly_infected(t+1) = mean_newly_infected(t+1) + newly_infected;
        mean_infected(t+1) = mean_infected(t+1) + sum(X==I);
        mean_susceptible(t+1) = mean_susceptible(t+1) + sum(X==S);
        mean_recovered(t+1) = mean_recovered(t+1) + sum(X==R);
    end
end

mean_newly_infected = mean_newly_infected / N;
mean_infected = mean_infected / N;
mean_susceptible = mean_susceptible / N;
mean_recovered = mean_recovered / N;



%% plot
figure
hold on
% Create a 2x2 grid of subplots
subplot(2, 2, 1); % First subplot
plot(0:T, mean_newly_infected);
title('Mean newly infected');
xlabel('t / Weeks')


subplot(2, 2, 2); % Third subplot
plot(0:T, mean_susceptible);
title('Mean susceptible');
xlabel('t / Weeks')


subplot(2, 2, 3); % Second subplot
plot(0:T, mean_infected);
title('Mean infected');
xlabel('t / Weeks')


subplot(2, 2, 4); % Fourth subplot
plot(0:T, mean_recovered);
title('Mean recovered');
xlabel('t / Weeks')


hold off

%% 1.2 

n = 900; % # of nodes
k = 2.3; % Average node degree
k0 = round(k+1); % Number of nodes in current graph
W = ones(k0,k0)-diag(ones(k0,1)); % Adjacency matrix
W = sparse(W); % Transform W into a sparse matrix
for i=(k0+1):n
    c = floor(k/2);
    if mod(k/2, 1) >= rand
        c = c + 1;
    end

    w = sum(W,2); % Degree vector of graph
    P = w./sum(w); % Probability vector for adding links

    for j=1:c % Choose c neighbours
        % Choose one neighbour with prob. prop. to node degrees
        neighbour = randsample(1:k0,1,true,full(P));
        % Assure not to choose the same neighbour again
        P(neighbour) = 0;
        W(k0+1,neighbour) = 1; % Add link (one direction)
        W(neighbour,k0+1) = 1; % Add link (other direction)
    end
    k0 = k0 + 1; % Number of nodes in current graph
end

disp("Average degree:" )
sum(sum(W)) / length(W)
%% Plot graph
G = graph(W);
plot(G)
title('Random graph from preferential attachment model')


%% 2 
n = 500; % # of nodes
k = 6; % Average node degree
k0 = round(k+1); % Number of nodes in current graph
W = ones(k0,k0)-diag(ones(k0,1)); % Adjacency matrix
W = sparse(W); % Transform W into a sparse matrix
for i=(k0+1):n
    c = floor(k/2);
    if mod(k/2, 1) >= rand
        c = c + 1;
    end

    w = sum(W,2); % Degree vector of graph
    P = w./sum(w); % Probability vector for adding links

    for j=1:c % Choose c neighbours
        % Choose one neighbour with prob. prop. to node degrees
        neighbour = randsample(1:k0,1,true,full(P));
        % Assure not to choose the same neighbour again
        P(neighbour) = 0;
        W(k0+1,neighbour) = 1; % Add link (one direction)
        W(neighbour,k0+1) = 1; % Add link (other direction)
    end
    k0 = k0 + 1; % Number of nodes in current graph
end

disp("Average degree:" )
sum(sum(W)) / length(W)
G = graph(W);
plot(G)

S = 0;
I = 1;
R = 2;

beta = 0.3;
rho = 0.7;


T = 15; 
N = 100;

mean_infected = zeros(T+1,1);
mean_newly_infected = zeros(T+1,1);
mean_susceptible = zeros(T+1,1);
mean_recovered = zeros(T+1,1);


for k = 1:N
    %initial X
    X = zeros(n,1);
    nr_initial_inf = 10;
    X(randi(n, [nr_initial_inf,1])) = I; 

    mean_infected(1) = mean_infected(1) + nr_initial_inf;
    mean_susceptible(1) = mean_susceptible(1) + n - nr_initial_inf;

    for t = 1:T
        newly_infected = 0;
        i_infected = X == I;
        nr_inf_neighbours = W*i_infected;
        
        susceptible = find(X == S);
        for i = susceptible'
            m = nr_inf_neighbours(i);
            prob = 1 - (1-beta)^m;
            eps = rand();
            
            if prob >= eps
                X(i) = I;
                newly_infected = newly_infected + 1;
            end
        end
    
        infected = find(i_infected);
        for i = infected'
            prob = rho;
            eps = rand();
            
            if prob >= eps
                X(i) = R;
            end
        end

        mean_newly_infected(t+1) = mean_newly_infected(t+1) + newly_infected;
        mean_infected(t+1) = mean_infected(t+1) + sum(X==I);
        mean_susceptible(t+1) = mean_susceptible(t+1) + sum(X==S);
        mean_recovered(t+1) = mean_recovered(t+1) + sum(X==R);
    end
end

mean_newly_infected = mean_newly_infected / N;
mean_infected = mean_infected / N;
mean_susceptible = mean_susceptible / N;
mean_recovered = mean_recovered / N;


%% plot
figure
hold on
% Create a 2x2 grid of subplots
subplot(2, 2, 1); % First subplot
plot(0:T, mean_newly_infected);
title('Mean newly infected');
xlabel('t / Weeks')


subplot(2, 2, 2); % Third subplot
plot(0:T, mean_susceptible);
title('Mean susceptible');
xlabel('t / Weeks')


subplot(2, 2, 3); % Second subplot
plot(0:T, mean_infected);
title('Mean infected');
xlabel('t / Weeks')


subplot(2, 2, 4); % Fourth subplot
plot(0:T, mean_recovered);
title('Mean recovered');
xlabel('t / Weeks')


hold off

%% 3

n = 500 ; % # of nodes
k = 6; % Average node degree
k0 = round(k+1); % Number of nodes in current graph
W = ones(k0,k0)-diag(ones(k0,1)); % Adjacency matrix
W = sparse(W); % Transform W into a sparse matrix
for i=(k0+1):n
    c = floor(k/2);
    if mod(k/2, 1) >= rand
        c = c + 1;
    end

    w = sum(W,2); % Degree vector of graph
    P = w./sum(w); % Probability vector for adding links

    for j=1:c % Choose c neighbours
        % Choose one neighbour with prob. prop. to node degrees
        neighbour = randsample(1:k0,1,true,full(P));
        % Assure not to choose the same neighbour again
        P(neighbour) = 0;
        W(k0+1,neighbour) = 1; % Add link (one direction)
        W(neighbour,k0+1) = 1; % Add link (other direction)
    end
    k0 = k0 + 1; % Number of nodes in current graph
end

disp("Average degree:" )
sum(sum(W)) / length(W);
G = graph(W);
figure
plot(G)

S = 0;
I = 1;
R = 2;

beta = 0.3;
rho = 0.7;
Vacc = [0, 0, 5, 15, 25, 35, 45, 55, 60, 60, 60, 60, 60, 60, 60, 60]; %adding extra 0 for t = 0

T = 15; 
N = 100;

mean_newly_infected = zeros(T+1,1);
mean_newly_vacc = zeros(T+1,1);
mean_infected = zeros(T+1,1);
mean_susceptible = zeros(T+1,1);
mean_recovered = zeros(T+1,1);
mean_vacc = zeros(T+1,1);


for k = 1:N
    %initial X
    X = zeros(n,1);
    nr_initial_inf = 10;
    X(randi(n, [nr_initial_inf,1])) = I; 
    %initial infected
    mean_infected(1) = mean_infected(1) + nr_initial_inf;
    mean_susceptible(1) = mean_susceptible(1) + n - nr_initial_inf;

    
    X_vacc = zeros(n,1); %binary vector 0 or 1

    %initial vaccinated
    nr_init_vacc = n*Vacc(1) / 100;
    not_vacc = find(X_vacc == 0);
    init_vacc = randsample(not_vacc, nr_init_vacc); %people who receives new vaccinations

    X_vacc(init_vacc) = 1;
    
    mean_vacc(1) = mean_vacc(1) + sum(X_vacc);

    for t = 2:T+1
        newly_infected = 0;
        nr_new_vacc = n*(Vacc(t) - Vacc(t-1)) / 100;
        not_vacc = find(X_vacc == 0);
        new_vacc = randsample(not_vacc, nr_new_vacc); %people who receives new vaccinations

        X_vacc(new_vacc) = 1;
        

        %S -> I
        i_infected = X == I;
        spread = i_infected;
        spread(X_vacc == 1) = 0; %vaccinated people can't spread

        nr_inf_neighbours = W*spread;
        susceptible = find(X == S);

        for i = susceptible'
            if(X_vacc(i) == 0)
                m = nr_inf_neighbours(i);
                prob = 1 - (1-beta)^m;
                eps = rand();
                
                if prob >= eps
                    X(i) = I;
                    newly_infected = newly_infected + 1;
                end
            end
        end
    
        %I -> R
        infected = find(i_infected);
        for i = infected'
            prob = rho;
            eps = rand();
            
            if prob >= eps
                X(i) = R;
            end
        end
        
        %Means
        mean_vacc(t) = mean_vacc(t) + sum(X_vacc);
        mean_newly_vacc(t) = mean_newly_vacc(t) + nr_new_vacc;

        mean_newly_infected(t) = mean_newly_infected(t) + newly_infected;
        mean_infected(t) = mean_infected(t) + sum(X==I);
        mean_susceptible(t) = mean_susceptible(t) + sum(X==S);
        mean_recovered(t) = mean_recovered(t) + sum(X==R);
    end
end

mean_newly_infected = mean_newly_infected / N;
mean_newly_vacc = mean_newly_vacc / N;
mean_infected = mean_infected / N;
mean_susceptible = mean_susceptible / N;
mean_recovered = mean_recovered / N;
mean_vacc = mean_vacc / N;


%% SIR test
beta = 0.5;
rho = 1;
Vacc = [0, 0, 5, 15, 25, 35, 45, 55, 60, 60, 60, 60, 60, 60, 60, 60]; %adding extra 0 for t = 0

T = 15; 
N = 100;

n = 500 ; % # of nodes
k = 6; % Average node degree

[mean_newly_infected, mean_susceptible, mean_infected, mean_recovered, mean_vacc, mean_newly_vacc] = SIR(n, k, beta, rho, Vacc, T, N);

%% plot
figure
hold on
% Create a 2x2 grid of subplots
subplot(3, 3, 1); % First subplot
plot(0:T, mean_newly_infected);
title('Mean newly infected');
xlabel('t / Weeks')


subplot(3, 3, 2); % Third subplot
plot(0:T, mean_susceptible);
title('Mean susceptible');
xlabel('t / Weeks')


subplot(3, 3, 3); % First subplot
plot(0:T, mean_newly_vacc);
title('Mean newly vaccinated');
xlabel('t / Weeks')


subplot(3, 3, 4); % Second subplot
plot(0:T, mean_infected);
title('Mean infected');
xlabel('t / Weeks')


subplot(3, 3, 5); % Fourth subplot
plot(0:T, mean_recovered);
title('Mean recovered');
xlabel('t / Weeks')


subplot(3, 3, 6); % Fourth subplot
plot(0:T, mean_vacc);
title('Mean vacc');
xlabel('t / Weeks')


hold off

%% 4: H1N1 Pandemic

Vacc = [5, 9, 16, 24, 32, 40, 47, 54, 59, 60, 60, 60, 60, 60, 60, 60];
I = [1, 1, 3, 5, 9, 17, 32, 32, 17, 5, 2, 1, 0, 0, 0, 0];

n = 934;
T = 15;
N = 10; 

k0 = 11;
beta0 = 0.18; 
rho0 = 0.6;

dk = 0.25;
dbeta = 0.025;
drho = 0.025;

run = true;
while run 
    k = [k0 - dk, k0, k0 + dk];
    beta = [beta0 - dbeta, beta0, beta0 + dbeta];
    rho = [rho0 - drho, rho0, rho0 + drho];

    if beta0 - dbeta < 0
        beta(1) = 0;
    elseif beta0 + dbeta > 1
        beta(3) = 1;
    end

    if rho0 - drho < 0
        rho(1) = 0;
    elseif rho0 + drho > 1
        rho(3) = 1;
    end
    
    RMSE = zeros(1, 3);
    RMSE(:,:, 3) = 0;

    for i = 1:3
        for j = 1:3
            for l = 1:3
                [mean_newly_infected, mean_susceptible, mean_infected, mean_recovered, mean_vacc] = SIR(n, k(i), beta(j), rho(l), Vacc, T, N);
                
                RMSE(i, j, l) = sqrt(mean( (mean_newly_infected(2:16)' - I(2:16)).^2 ));
            end
        end
    end

    [min_RMSE, ind] = min(RMSE, [], 'all');
    [i, j, l] = ind2sub(size(RMSE), ind);

    min_RMSE
    if (i == 2) && (j == 2) && (l == 2)
        run = false;
    else
        k0 = k(i)
        beta0 = beta(j)
        rho0 = rho(l)
    end

    [mean_newly_infected, mean_susceptible, mean_infected, mean_recovered, mean_vacc, mean_newly_vacc] = SIR(n, k0, beta0, k0, Vacc, T, N);
    % figure
    % hold on
    % plot(0:15, mean_newly_infected, 'b')
    % plot(0:15, I, 'r')
    % hold off
end

%%

figure
hold on
plot(0:15, mean_newly_infected, 'b')
plot(0:15, I, 'r')
hold off

%% test
Vacc = [5, 9, 16, 24, 32, 40, 47, 54, 59, 60, 60, 60, 60, 60, 60, 60];
I = [1, 1, 3, 5, 9, 17, 32, 32, 17, 5, 2, 1, 0, 0, 0, 0];

n = 934;
T = 15;
N = 10; 
[mean_newly_infected, mean_susceptible, mean_infected, mean_recovered, mean_vacc, mean_newly_vacc] = SIR(n, 11.5, 0.17, 0.55, Vacc, T, 100);
figure
title('Mean newly infected')
hold on
plot(0:15, mean_newly_infected, 'b')
plot(0:15, I, 'r')
xlabel('t / Weeks')
legend('Simulation','H1N1 data')
hold off


rmse = sqrt(mean( (mean_newly_infected(2:16)' - I(2:16)).^2 ))




%% plot
figure
hold on
% Create a 2x2 grid of subplots
subplot(3, 3, 1); % First subplot
plot(0:T, mean_newly_infected);
title('Mean newly infected');
xlabel('t / Weeks')


subplot(3, 3, 2); % Third subplot
plot(0:T, mean_susceptible);
title('Mean susceptible');
xlabel('t / Weeks')


subplot(3, 3, 3); % First subplot
plot(0:T, mean_newly_vacc);
title('Mean newly vaccinated');
xlabel('t / Weeks')


subplot(3, 3, 4); % Second subplot
plot(0:T, mean_infected);
title('Mean infected');
xlabel('t / Weeks')



subplot(3, 3, 5); % Fourth subplot
plot(0:T, mean_recovered);
title('Mean recovered');
xlabel('t / Weeks')


subplot(3, 3, 6); % Fourth subplot
plot(0:T, mean_vacc);
title('Mean vacc');
xlabel('t / Weeks')


hold off

%% 5: challenge

Vacc = [5, 9, 16, 24, 32, 40, 47, 54, 59, 60, 60, 60, 60, 60, 60, 60];
I = [1, 1, 3, 5, 9, 17, 32, 32, 17, 5, 2, 1, 0, 0, 0, 0];

n = 934;
T = 15;
N = 100; 
numClusters = 5;


p0 = 0.01;
q0 = 0.02;
beta0 = 0.15; 
rho0 = 0.6;

dp = 0.0005;
dq = 0.005;
dbeta = 0.01;
drho = 0.05;

run = true;
while run 
    p = [p0 - dp, p0, p0 + dp];
    q = [q0 - dq, q0, q0 + dq];
    beta = [beta0 - dbeta, beta0, beta0 + dbeta];
    rho = [rho0 - drho, rho0, rho0 + drho];

    if beta0 - dbeta < 0
        beta(1) = 0;
    elseif beta0 + dbeta > 1
        beta(3) = 1;
    end

    if rho0 - drho < 0
        rho(1) = 0;
    elseif rho0 + drho > 1
        rho(3) = 1;
    end
    
    RMSE = zeros(1, 3);
    RMSE(:,:, 3) = 0;
    RMSE(:,:,:,3) = 0;

    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    [mean_newly_infected, mean_susceptible, mean_infected, mean_recovered, mean_vacc] = SIRv2(n, numClusters, p(i), q(j), beta(k), rho(l), Vacc, T, N);
                
                    RMSE(i, j, k, l) = sqrt(mean( (mean_newly_infected(2:16)' - I(2:16)).^2 ));
                end
            end
        end
    end

    [min_RMSE, ind] = min(RMSE, [], 'all');
    [i, j, k, l] = ind2sub(size(RMSE), ind);

    min_RMSE
    if (i == 2) && (j == 2) && (k == 2) && (l == 2)
        run = false;
    else
        p0 = p(i)
        q0 = q(j)
        beta0 = beta(k)
        rho0 = rho(l)
    end

    %[mean_newly_infected, mean_susceptible, mean_infected, mean_recovered, mean_vacc, mean_newly_vacc] = SIR(n, k0, beta0, k0, Vacc, T, N);

end

%%
Vacc = [5, 9, 16, 24, 32, 40, 47, 54, 59, 60, 60, 60, 60, 60, 60, 60];
I = [1, 1, 3, 5, 9, 17, 32, 32, 17, 5, 2, 1, 0, 0, 0, 0];

n = 934;
T = 15;
N = 100; 
numClusters = 5;

[mean_newly_infected, mean_susceptible, mean_infected, mean_recovered, mean_vacc, mean_newly_vacc] = SIRv2(n, numClusters, 0.01, 0.02, 0.175, 0.95, Vacc, T, 200);

figure
title('Mean newly infected')
hold on
plot(0:15, mean_newly_infected, 'b')
plot(0:15, I, 'r')
xlabel('t / Weeks')
legend('Simulation','H1N1 data')
hold off



rmse = sqrt(mean( (mean_newly_infected(2:16)' - I(2:16)).^2 ))


%% Genetic Algorithm

options = optimoptions('ga', ...
    'PopulationSize', 25, ...
    'MaxGenerations', 15, ...
    'CrossoverFraction', 0.8, ...
    'EliteCount', 2, ...
    'MutationFcn', @mutationadaptfeasible)

nvars = 4; %p, q, beta, rho
lb = [0.005, 0.005, 0.1, 0.1];
ub = [0.025, 0.025, 0.5, 0.6];

[x, fval] = ga(@objectiveFunction, nvars, [], [], [], [], lb, ub, [], options);

disp(['Optimal solution: ', num2str(x)]);
disp(['Objective function value at optimal solution: ', num2str(fval)]);




%%

Vacc = [5, 9, 16, 24, 32, 40, 47, 54, 59, 60, 60, 60, 60, 60, 60, 60];
I = [1, 1, 3, 5, 9, 17, 32, 32, 17, 5, 2, 1, 0, 0, 0, 0];

n = 934;
T = 15;
N = 100; 
numClusters = 5;

[mean_newly_infected, mean_susceptible, mean_infected, mean_recovered, mean_vacc, mean_newly_vacc] = SIRv2(n, numClusters, 0.0061525, 0.0054639 , 0.49186, 0.5804, Vacc, T, 200);

figure
title('Mean newly infected')
hold on
plot(0:15, mean_newly_infected, 'b')
plot(0:15, I, 'r')
xlabel('t / Weeks')
legend('Simulation','H1N1 data')
hold off


rmse = sqrt(mean( (mean_newly_infected(2:16)' - I(2:16)).^2 ))

