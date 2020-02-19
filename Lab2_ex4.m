%% Implementation of a local search algorithm for the LTD problem in which the 
% topology is a bidirectional Manhattan.
clc; close all; clear all;

% Definition of the parameters
Nodes = 40;
Delta = 4; % Number of transmitters and receivers per node
Ua = 0;
Ub = 4;
Row = 5; %Manhattan rows
Col = 8; %Manhattan columns
% Table in which each node of the Manhattan topology will be associated
% with a node of the considered network
names = table('Size', [Nodes, 1], 'VariableTypes', ["int8"], 'VariableNames', ["Names"]);

rng(25)
% Creation of the traffic matrix
Tsd = (Ub-Ua).*rand(Nodes,Nodes)+Ua;
Tsd = Tsd - diag(diag(Tsd));

b = zeros(Nodes,Nodes); % Adjacent matrix
G1 = digraph; 
iter = 100; % Number of iterations for the local search algorithm.
opt_solution = sum(sum(Tsd)); % The initial solution is the worst possible case, in which the entire traffic is routed on one single link.
%% Creation of the Manhattan Topology
%Horizontal Links of Manhattan Topology
for i=1:Nodes
    if (i == 1) || (mod(i-1, Col)==0)
        G1 = addedge(G1, i, i+1);
        b(i, i+1) = 1;
        G1 = addedge(G1, i, i+Col-1);
        b(i, i+Col-1) = 1;
    else
        if mod(i, Col) == 0
            G1 = addedge(G1, i, i-1);
            G1 = addedge(G1, i, i-Col+1);
            b(i, i-1) = 1;
            b(i, i-Col+1) = 1;
        else
            G1 = addedge(G1, i, i+1);
            G1 = addedge(G1, i, i-1);
            b(i, i+1) = 1;
            b(i, i-1) = 1;
        end
    end
end
%Vertical Links
for i=1:Nodes
    if (i<=Col)
        G1 = addedge(G1, i, i+Col);
        G1 = addedge(G1, i, Nodes-Col+i);
        b(i, i+Col) = 1;
        b(i, Nodes-Col+i) = 1;
    else
        if (i > Nodes-Col)
            G1 = addedge(G1, i, i-Col);
            G1 = addedge(G1, i, i-Nodes+Col);
            b(i, i-Col) = 1;
            b(i, i-Nodes+Col) = 1;
        else
            G1 = addedge(G1, i, i+Col);
            G1 = addedge(G1, i, i-Col);
            b(i, i+Col) = 1;
            b(i, i-Col) = 1;
        end
    end
end
figure(1);
plot(G1);
%% Smart Positioning of the nodes in the Manhattan
while (max(Tsd(:)) >  0) % While there is still a positive traffic in the traffic matrix:
    [M, I] = max(Tsd(:));
    [s, d] = ind2sub(size(Tsd),I); % Pick the maximum traffic between two nodes s and d
    %% First Scenario: Both the nodes are already included in the Topology.
    if (ismember([s,d], names.Names) == 1)
        Tsd(s, d) = -Tsd(s, d); % Do nothing, simply mark that traffic as already selected.
    end
    %% Second Scenario: Both the nodes are not in the Topology
    if (ismember(s, names.Names) == 0) && (ismember(d, names.Names) == 0)
        % Pick a random node that has not been associated with a node of
        % the network
        ok = 0;
        while ok == 0
            source = randi([1, Nodes]);
            if names.Names(source) == 0
                ok = 1;
            end
        end
        ok = 0;
        names.Names(source) = s; % Associate that node with the current node s.
        % Try to place the other node d as close as possible to s.
        clos_neigh = nearest(G1, source, Nodes);
        i = 1;
        while ok == 0
            dest = clos_neigh(i);
            if names.Names(dest) == 0
                ok = 1;
            else
                i = i+1;
            end
        end
        names.Names(dest) = d; % Associate the closest available node to node d.
        Tsd(s, d) = -Tsd(s, d);
    end
    %% Third scenario: One of the two nodes have already been placed in the topology.
    if (ismember(s, names.Names) == 1) && (ismember(d, names.Names) == 0)
    % If the placed node is s, try to place node d as close as possible
    % to node s.
        clos_neigh = nearest(G1, find(names.Names == s), Nodes);
        ok = 0;
        i = 1;
        while ok == 0
            dest = clos_neigh(i);
            if names.Names(dest) == 0
                ok = 1;
            else
                i = i+1;
            end
        end
        names.Names(dest) = d;
        Tsd(s, d) = -Tsd(s, d);
    end
    if (ismember(s, names.Names) == 0) && (ismember(d, names.Names) == 1)
    % If the placed node is d, try to place node s as close as possible to
    % node d.
        clos_neigh = nearest(G1, find(names.Names == d), Nodes);
        ok = 0;
        i = 1;
        while ok == 0
            source = clos_neigh(i);
            if names.Names(source) == 0
                ok = 1;
            else
                i = i+1;
            end
        end
        names.Names(source) = s;
        Tsd(s, d) = -Tsd(s, d);
    end
end
Tsd = -Tsd; % Reset the traffic matrix to the original values
G1.Nodes = names;
figure(2);
plot(G1, 'NodeLabel', names.Names);
% Implementation of the local search.
for N=1:iter
    %% Routing Strategy, route the network according to the current node placement.
    fij = zeros(Nodes, Nodes);
    fsd = zeros(Nodes, Nodes, Nodes, Nodes);
    for s = 1:Nodes
        for d = 1:Nodes
            % Find the shortest path between node s and node d
            path = shortestpath(G1, find(names.Names == s), find(names.Names == d));
            path = names.Names(path);
            % For each link in the shortest path, add the amount of traffic
            % belonging to nodes s-d.
            for i=1:(length(path)-1)
                fij(path(i), path(i+1)) = fij(path(i), path(i+1)) + Tsd(s, d); %Per ogni hop tra s e d incrementa il traffico sul link i j
                fsd(s, d, path(i), path(i+1)) = Tsd(s, d);
            end
        end
    end
    result = max(fij(:)); % Solution of the LTD according to the current node placement
    % If the solution is better than the current one, update the best
    % solution
    if result < opt_solution 
        opt_solution = result;
    % If the solution is not the optimal one, put the two nodes previously
    % swapped in their original position. At the first iteration this step
    % cannot be performed, since the initial solution will be always better
    % than the current one.
    else
        names.Names([swap2 swap1]) =  names.Names([swap1 swap2]);
    end
    % In order to generate neighbor solutions, pick to random nodes and
    % swap them. Then route again the traffic in order to see if there are improvements 
    % in the solution.
    swap1 = randi([1, Nodes]);
    swap2 = randi([1, Nodes]);
    names.Names([swap1 swap2]) =  names.Names([swap2 swap1]);
    % Repeat these steps for a number of iterations, and then take the last
    % current best solution as the optimal one. 
end
G1.Nodes = names;
figure(3);
plot(G1, 'NodeLabel', names.Names);