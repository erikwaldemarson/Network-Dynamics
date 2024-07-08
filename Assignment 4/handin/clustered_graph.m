function [W, avg_deg] = clustered_graph(n, numClusters, p, q)

nodesPerCluster = round(n / numClusters); % nodes per cluster

% Initialize adjacency matrix
W = zeros(n);

% Intra-cluster edges
for c = 1:numClusters
    startIndex = (c-1)*nodesPerCluster + 1;
    endIndex = min(c*nodesPerCluster, n);
    clusterSize = endIndex - startIndex + 1;
    intraClusterEdges = rand(clusterSize) < p;
    intraClusterEdges = triu(intraClusterEdges, 1);
    W(startIndex:endIndex, startIndex:endIndex) = intraClusterEdges + intraClusterEdges';
end

% Inter-cluster edges
for c1 = 1:numClusters
    for c2 = c1+1:numClusters
        startIndex1 = (c1-1)*nodesPerCluster + 1;
        endIndex1 = min(c1*nodesPerCluster, n);
        startIndex2 = (c2-1)*nodesPerCluster + 1;
        endIndex2 = min(c2*nodesPerCluster, n);
        
        clusterSize1 = endIndex1 - startIndex1 + 1;
        clusterSize2 = endIndex2 - startIndex2 + 1;
        
        interClusterEdges = rand(clusterSize1, clusterSize2) < q;
        W(startIndex1:endIndex1, startIndex2:endIndex2) = interClusterEdges;
        W(startIndex2:endIndex2, startIndex1:endIndex1) = interClusterEdges';
    end
end

avg_deg = sum(sum(W)) / length(W);

end

