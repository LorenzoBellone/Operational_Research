%% Implementation of a greedy heuristic algorithm for the resolution of the
% LTD problem in case of unbalanced traffic matrix.
clc; clear all; close all; 

% Definition of the parameters
Delta = 4; % Number of transmitters and receivers for each node.
% Ranges for the distribution of the traffic.
% High traffic 
Ua_HT = 5;  
Ub_HT = 15;
% Low traffic
Ua_LT = 0;
Ub_LT = 3;

% Probability that 2 nodes exchange a high or low traffic.
HT = 0.1;
LT = 1-HT;
rng(25)
Nodes = 40;

% Creation of the unbalanced traffic matrix.
Tsd_HT = (Ub_HT-Ua_HT).*rand(Nodes,Nodes)+Ua_HT;
Tsd_HT = Tsd_HT - diag(diag(Tsd_HT));
Tsd_LT = (Ub_LT-Ua_LT).*rand(Nodes,Nodes)+Ua_LT;
Tsd_LT = Tsd_LT - diag(diag(Tsd_LT));

% Two different matrices, one with a low traffic and the another one with
% high traffic.
Tsd=zeros(Nodes, Nodes);
% Each entry of the traffic matrix will have with probability 0.1 a high
% traffic, with probability 0.9 a low traffic.
for row = 1:Nodes
    for col = 1:Nodes
        random_decision = rand;
        if random_decision < LT
            Tsd(row,col) = Tsd_LT(row,col);
        else
            Tsd(row,col) = Tsd_HT(row,col);
        end     
    end
end

b = zeros(Nodes,Nodes); % Adjacent matrix
fsd = zeros(Nodes,Nodes,Nodes,Nodes); % Traffic on link ij belonging to nodes s-d
fij = zeros(Nodes,Nodes); % Total traffic on link ij
DeltaTx = zeros(Nodes,1); % Out-degree of node i
DeltaRx = zeros(Nodes,1); % In-degree of node j
trans = [];
rec = [];

i = 1;
% Implementation of the algorithm:
while (max(Tsd(:)) >  0) % While there is still a positive traffic in the traffic matrix:
    [M, I] = max(Tsd(:));
    [s, d] = ind2sub(size(Tsd),I); % Pick the maximum traffic between two nodes s and d
    if (DeltaTx(s, 1) < Delta) && (DeltaRx(d, 1) < Delta) % If the degree constraints are satisfied:
        DeltaTx(s, 1) = DeltaTx(s, 1) + 1;
        DeltaRx(d, 1) = DeltaRx(d, 1) + 1;
        b(s, d) = 1; % Add a link between node s and node d
        trans(i) = s;
        rec(i) = d;
        i = i+1;
    end
    Tsd(s, d) = -Tsd(s, d);  % Make that traffic entry negative in order to check all the traffic in the matrix.
end
Tsd(:, :) = -Tsd(:, :); % Reset the traffic matrix to the original values
disp('End While')

G1 = digraph(trans, rec);
figure(1);
plot(G1);
% Implementation of a routing strategy based on the shortest path algorithm:
for s = 1:Nodes
    for d = 1:Nodes
        path = shortestpath(G1, s, d); % For any pair of nodes, evaluate the shortest path given the topology created before.
        % For each link in the shortest path, add the amount of traffic
        % belonging to nodes s-d.
        for i=1:(length(path)-1)
            fij(path(i), path(i+1)) = fij(path(i), path(i+1)) + Tsd(s, d);
            fsd(s, d, path(i), path(i+1)) = Tsd(s, d);
        end
    end
end
result = max(fij(:)); % The most congested link is the higher fij.

