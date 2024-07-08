function [W, avg_deg] = random_graph(n,k)
%Generates a random graph of n nodes with avg. degree k >= 2 using preferential
%attachment. Returns weight matrix W

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

avg_deg = sum(sum(W)) / length(W);

end

