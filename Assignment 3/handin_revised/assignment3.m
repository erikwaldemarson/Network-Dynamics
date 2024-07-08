%% matrix

Lambda = [0, 2/5, 1/5, 0, 0;
         0, 0, 3/4, 1/4, 0;
         1/2, 0, 0, 1/2, 0;
         0, 0, 1/3, 0, 2/3;
         0, 1/3, 0, 1/3, 0];
     
w = sum(Lambda, 2);
wstar = max(w);
Pbar = Lambda./wstar;
Pbar = Pbar + diag(ones(length(w),1) - sum(Pbar, 2));

%% 1a

a = 2; %which node to start in, in this case a=2

niterations = 10^5;
return_time = 0;
cumpbar = cumsum(Pbar, 2);

for i = 1:niterations
    pos = a; %starts in node a
    tnext = -log(rand())/wstar;

    while pos == a %initialize, must leave a
        pos = find(cumpbar(pos, :) > rand(),1);
        return_time = return_time + tnext;
        tnext = -log(rand())/wstar;
    end

    while pos ~= a  %how long time to return
        pos = find(cumpbar(pos, :) > rand(),1);
        return_time = return_time + tnext;
        tnext = -log(rand())/wstar;
    end
end

avg_return_time = return_time / niterations;

%% 1b
[V,D] = eig(Pbar');
lambda = diag(D);
[~,i] = sort(lambda,1,'descend');
pibar = V(:,i(1));
pibar = pibar/sum(pibar);

1 / (pibar(2)*w(2))

%% 1c

o = 1; %which node to start in, o =1
d = 5; %which node to end in, d = 5

niterations = 10^5;

hitting_time = 0;
cumpbar = cumsum(Pbar, 2);


for i = 1:niterations
    pos = o; %starts in node o
    tnext = -log(rand())/wstar;

    while pos ~= d  %how long to return
        pos = find(cumpbar(pos, :) > rand(),1);
        hitting_time = hitting_time + tnext;
        tnext = -log(rand())/wstar;
    end
end

avg_hitting_time = hitting_time / niterations

%% 1d
[n,~] = size(Lambda);
P = diag(w) \ Lambda;
o = 1; %which node to start in, o =1
d = 5; %which node to end in, d = 5

S = d;
R = setdiff(1:n, S);
Phat = P(R,R);
what = w(R);
xhat = (eye(length(R)) - Phat) \ (what).^(-1);
Ts = zeros(n, 1);
Ts(R) = xhat;


%% 2a
n = 10;
W = diag(ones(n-1,1),1) + diag(ones(n-1,1),-1); %simple line graph
states = [1, 2]; %red = 1, green = 2

cost = @(s, x) s == x; %cost function
eta = @(t) t / 100; %inverse noise

niterations = 500; %iterations
X = ones(n, niterations); %intialize as red

U = zeros(1, niterations); %potential function

for i = 1:n 
     for j = 1:n
            U(1) = U(1) + 0.5*W(i,j)*c(X(i,1), X(j, 1)); %initial U
     end
end

for t = 1:niterations-1
    i = randi(n); %particle / node chosen

    z_t = zeros(1,length(states));
    for s = states

        cost_tot = 0;
        for j = 1:n
            cost_tot = cost_tot  + W(i,j)*cost(s,X(j,t));
        end

        z_t(s) = exp(-eta(t)*cost_tot); %proportional to prob of state s

    end

    probs = z_t ./ sum(z_t);
    cumprobs = cumsum(probs);

    X(:,t+1) = X(:,t);
    X(i,t+1) = find(cumprobs >= rand(), 1);

    for i = 1:n 
        for j = 1:n
            U(t+1) = U(t+1) + 0.5*W(i,j)*cost(X(i,t+1), X(j, t+1));
        end
    end
end
figure
plot(U,'-o')
title("Potential function \Phi(x) with \eta = {\it t} / 100 for t \in [1, 500].")
xlabel('\it t')
ylabel('\Phi(x)')

%% 2b
load('coord.mat', '-ascii');
load('wifi.mat', '-ascii');

W = sparse(wifi);
[n, ~] = size(W);
states = 1:8;

eta = @(t) t/10; %inverse noise

niterations = 500; %iterations
X = zeros(n, niterations); 
X(:,1) = randi(8,[n,1]);%intialize everyone as random state

U = zeros(1, niterations); %potential function


for i = 1:n 
     for j = 1:n
            U(1) = U(1) + 0.5*W(i,j)*c(X(i,1), X(j, 1)); %initial U
     end
end


for t = 1:niterations-1
    i = randi(n); %particle / node chosen

    z_t = zeros(1,length(states));
    for s = states

        cost_tot = 0;
        for j = 1:n
            cost_tot = cost_tot  + W(i,j)*c(s,X(j,t));
        end

        z_t(s) = exp(-eta(t)*cost_tot); %proportional to prob of state s
    end

    probs = z_t ./ sum(z_t);
    cumprobs = cumsum(probs);

    X(:,t+1) = X(:,t);
    X(i,t+1) = find(cumprobs >= rand(), 1);
    
    for i = 1:n 
        for j = 1:n
            U(t+1) = U(t+1) + 0.5*W(i,j)*c(X(i,t+1), X(j, t+1));
        end
    end

end
figure
plot(U,'-o')
title("Potential function U({\itt}) with \eta = {\it t} / 10 for t \in [1, 500].")
xlabel('\it t')
ylabel('U({\itt})')


%%
figure
hold on
gplot(W, coord ,'-k')
colors = {'r', 'g', 'b', 'y', 'm', 'c', 'w', 'k'};
x = X(:,end);
for i = 1:n
    scatter(coord(i,1),coord(i,2),200,'markeredgecolor','k','markerfacecolor', cell2mat(colors(x(i))));
end
title('Best response for \eta = {\it t} / 10 for t = 500.')
hold off



%% 

x = ones(100, 1)';

c(1, x)


function cost = c(s,x)
    if s == x
        cost = 2;
    elseif abs(s-x) == 1
        cost = 1;
    else 
        cost = 0;
    end
end









