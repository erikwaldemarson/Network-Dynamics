function [mean_newly_infected, mean_susceptible, mean_infected, mean_recovered, mean_vacc, mean_newly_vacc] = SIRv2(n, numClusters, p, q, beta, rho, Vacc, T, N)
%Simulate a pandemic using SIR model on a given graph described by W

W = clustered_graph(n, numClusters, p, q);

S = 0;
I = 1;
R = 2;

mean_newly_infected = zeros(T+1,1);
mean_newly_vacc = zeros(T+1,1);
mean_infected = zeros(T+1,1);
mean_susceptible = zeros(T+1,1);
mean_recovered = zeros(T+1,1);
mean_vacc = zeros(T+1,1);


for k = 1:N
    %initial X
    X = zeros(n,1);
    nr_initial_inf = 1;
    X(randi(n, [nr_initial_inf,1])) = I; 
    %initial infected
    mean_infected(1) = mean_infected(1) + nr_initial_inf;
    mean_susceptible(1) = mean_susceptible(1) + n - nr_initial_inf;
    mean_newly_infected(1) = mean_newly_infected(1) + nr_initial_inf;

    X_vacc = zeros(n,1); %binary vector 0 or 1

    %initial vaccinated
    nr_init_vacc = round(n*Vacc(1)/ 100);
    not_vacc = find(X_vacc == 0);
    init_vacc = randsample(not_vacc, nr_init_vacc); %people who receives new vaccinations

    X_vacc(init_vacc) = 1;
    
    mean_vacc(1) = mean_vacc(1) + sum(X_vacc);

    for t = 2:T+1
    
        newly_infected = 0;
        nr_new_vacc = round(n*(Vacc(t) - Vacc(t-1)) / 100);
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

end
