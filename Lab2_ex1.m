%% Implementation of a greedy heuristic algorithm for the resolution of the 
%  LTD problem.

clc; clear all; close all; 
% Definition of the parameters
Delta = 4; % Number of transmitters and receivers for each node.
Ua = 0;
Ub = 4; % Ranges for the distribution of the traffic

rng(25)
Nodes = 40; % Number of nodes in the network
% Creation of the uniform traffic matrix:
Tsd = (Ub-Ua).*rand(Nodes,Nodes)+Ua;
Tsd = Tsd - diag(diag(Tsd));
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
    Tsd(s, d) = -Tsd(s, d); % Make that traffic entry negative in order to check all the traffic in the matrix.
end
Tsd(:, :) = -Tsd(:, :); % Reset the traffic matrix to the original values
disp('End While');
G = digraph(trans, rec);
figure(1);
plot(G);
title(['Nodes =  ', num2str(Nodes),',  \Delta =  ', num2str(Delta)]);
% Implementation of a routing strategy based on the shortest path algorithm:
for s = 1:Nodes
    for d = 1:Nodes 
        path = shortestpath(G, s, d); % For any pair of nodes, evaluate the shortest path given the topology created before.
        if isempty(path)
            disp("Not optimal"); % If there is no path, the graph is disconnected and the solution cannot be evaluated.
        end
        % For each link in the shortest path, add the amount of traffic
        % belonging to nodes s-d.
        for i=1:(length(path)-1)
            fij(path(i), path(i+1)) = fij(path(i), path(i+1)) + Tsd(s, d);
            fsd(s, d, path(i), path(i+1)) = Tsd(s, d);
        end
    end
end
result = max(fij(:)); % The most congested link is the higher fij.
%saveas(figure(1), "Topology_Lab2_ex1_Delta4.png");