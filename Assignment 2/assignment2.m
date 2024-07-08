%% load
load traffic.mat
load capacities.mat -ascii
load traveltime.mat -ascii
load flow.mat -ascii

%% 1a

[num_nodes, num_edges] = size(traffic);

head = zeros(1, num_edges);
tail = zeros(1, num_edges);

for i = 1:num_edges
    edge = traffic(:,i);
    
    head(i) = find(edge==1);
    tail(i) = find(edge==-1);
end

G = graph(head, tail, traveltime);

figure
p = plot(G,'EdgeLabel',G.Edges.Weight);
title('Shortest path (red) from Santa Monica (1) to Santa Ana (17)')
path = shortestpath(G,1,17);
highlight(p, path,'EdgeColor','red')

%% 1b

G_flow = graph(head, tail, capacities);
max_flow = maxflow(G_flow, 1, 17)

%% 1c

nu_ext = traffic*flow;
net_inflow = nu_ext(1);

%% 1d

B = traffic;
nu = zeros(num_nodes, 1);

nu(1) = net_inflow;
nu(num_nodes) = -net_inflow;

cvx_begin
    variables f(num_edges)
    
    minimize sum(capacities.*traveltime.*inv_pos(1 - f./capacities) - capacities.*traveltime)
    
    subject to 
        B*f == nu;
        0 <= f <= capacities;
cvx_end




%% 1e 

cvx_begin
    variables f_0(num_edges)
    
    minimize sum(-capacities.*traveltime.*log(1 - f_0./capacities))
    
    subject to 
        B*f_0 == nu;
        0 <= f_0 <= capacities;
cvx_end

%% cost
cost_social_opt = sum(capacities.*traveltime.*inv_pos(1 - f./capacities) - capacities.*traveltime);
cost_wardrop = sum(capacities.*traveltime.*inv_pos(1 - f_0./capacities) - capacities.*traveltime);
PoA_0 = cost_wardrop / cost_social_opt;

%% 1f


omega = f.*(traveltime./capacities)./ (1 - f./capacities).^2;

cvx_begin
    variables f_omega(num_edges)
    
    minimize sum(-capacities.*traveltime.*log(1 - f_omega./capacities) + omega.*f_omega)
    
    subject to 
        B*f_omega == nu;
        0 <= f_omega <= capacities;
cvx_end

%% cost
cost_wardrop_tolls = sum(capacities.*traveltime.*inv_pos(1 - f_omega./capacities) - capacities.*traveltime);
PoA_omega = cost_wardrop_tolls / cost_social_opt;


%% 1g

cvx_begin
    variables f(num_edges)
    
    minimize sum(capacities.*traveltime.*inv_pos(1 - f./capacities) - capacities.*traveltime - f.*traveltime)
    
    subject to 
        B*f == nu;
        0 <= f <= capacities;
cvx_end



omega = f.*(traveltime./capacities)./ (1 - f./capacities).^2;

cvx_begin
    variables f_omega(num_edges)
    
    minimize sum(-capacities.*traveltime.*log(1 - f_omega./capacities) -f_omega.*traveltime + omega.*f_omega)
    
    subject to 
        B*f_omega == nu;
        0 <= f_omega <= capacities;
cvx_end


cost_social_opt = sum(capacities.*traveltime.*inv_pos(1 - f./capacities) - capacities.*traveltime - f.*traveltime);
cost_wardrop_tolls = sum(capacities.*traveltime.*inv_pos(1 - f_omega./capacities) - capacities.*traveltime - f_omega.*traveltime);

PoA_omega = cost_wardrop_tolls / cost_social_opt;


