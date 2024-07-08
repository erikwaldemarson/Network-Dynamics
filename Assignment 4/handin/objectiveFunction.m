function RMSE = objectiveFunction(x)
%x1 = p, x2 = q, x3 = beta, x4 = rho
    Vacc = [5, 9, 16, 24, 32, 40, 47, 54, 59, 60, 60, 60, 60, 60, 60, 60];
    I = [1, 1, 3, 5, 9, 17, 32, 32, 17, 5, 2, 1, 0, 0, 0, 0];
    
    n = 934;
    T = 15;
    N = 100; 
    numClusters = 5;

    [mean_newly_infected, mean_susceptible, mean_infected, mean_recovered, mean_vacc] = SIRv2(n, numClusters, x(1), x(2), x(3), x(4), Vacc, T, N);

    RMSE = sqrt(mean( (mean_newly_infected(2:16)' - I(2:16)).^2 ));

end