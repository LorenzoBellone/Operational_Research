clc; clear all; close all; 
% roba in più ^_^
Ua_HT = 5;  % boundarys for the traffic tsd
Ub_HT = 15;
Ua_LT = 0;
Ub_LT = 3;
HT = 0.1;
LT = 1-HT;

%rng('default')
Nodes = 20;
Tsd_HT = (Ub_HT-Ua_HT).*rand(Nodes,Nodes)+Ua_HT;
Tsd_HT = Tsd_HT - diag(diag(Tsd_HT));
Tsd_LT = (Ub_LT-Ua_LT).*rand(Nodes,Nodes)+Ua_LT;
Tsd_LT = Tsd_LT - diag(diag(Tsd_LT));

Tsd=zeros(Nodes, Nodes);

%generate a random number such that 90% of case Tsd = Tsd_LT
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

% parameters
Delta = 4;
Row = 4; %Manhattan rows
Col = 5; %Manhattan columns
names = table('Size', [Nodes, 1], 'VariableTypes', ["int8"], 'VariableNames', ["Names"]);
b = zeros(Nodes,Nodes);
fsd = zeros(Nodes,Nodes,Nodes,Nodes);
fij = zeros(Nodes,Nodes);
G1 = digraph;
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
while (max(Tsd(:)) >  0)
    [M, I] = max(Tsd(:));
    [s, d] = ind2sub(size(Tsd),I);
    if (ismember([s,d], names.Names) == 1)
        Tsd(s, d) = -Tsd(s, d);
    end
    if (ismember(s, names.Names) == 0) && (ismember(d, names.Names) == 0)
        ok = 0;
        while ok == 0
            source = randi([1, Nodes]);
            if names.Names(source) == 0
                ok = 1;
            end
        end
        ok = 0;
        names.Names(source) = s;
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
        names.Names(dest) = d;
        Tsd(s, d) = -Tsd(s, d);
    end
    if (ismember(s, names.Names) == 1) && (ismember(d, names.Names) == 0)
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
Tsd = -Tsd;
G1.Nodes = names;
figure(2);
plot(G1, 'NodeLabel', names.Names);
%% Routing Strategy
for s = 1:Nodes
    for d = 1:Nodes
        path = shortestpath(G1, find(names.Names == s), find(names.Names == d)); %Calcola il percorso più breve tra s e d
        path = names.Names(path);
        for i=1:(length(path)-1)
            fij(path(i), path(i+1)) = fij(path(i), path(i+1)) + Tsd(s, d); %Per ogni hop tra s e d incrementa il traffico sul link i j
            fsd(s, d, path(i), path(i+1)) = Tsd(s, d);
        end
    end
end
result = max(fij(:));
Row = 4; %Manhattan rows
Col = 5; %Manhattan columns
names = table('Size', [Nodes, 1], 'VariableTypes', ["int8"], 'VariableNames', ["Names"]);

rng('default')
b = zeros(Nodes,Nodes);
fsd = zeros(Nodes,Nodes,Nodes,Nodes);
fij = zeros(Nodes,Nodes);
G1 = digraph;
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
while (max(Tsd(:)) >  0)
    [M, I] = max(Tsd(:));
    [s, d] = ind2sub(size(Tsd),I);
    if (ismember([s,d], names.Names) == 1)
        Tsd(s, d) = -Tsd(s, d);
    end
    if (ismember(s, names.Names) == 0) && (ismember(d, names.Names) == 0)
        ok = 0;
        while ok == 0
            source = randi([1, Nodes]);
            if names.Names(source) == 0
                ok = 1;
            end
        end
        ok = 0;
        names.Names(source) = s;
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
        names.Names(dest) = d;
        Tsd(s, d) = -Tsd(s, d);
    end
    if (ismember(s, names.Names) == 1) && (ismember(d, names.Names) == 0)
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
Tsd = -Tsd;
G1.Nodes = names;
figure(2);
plot(G1, 'NodeLabel', names.Names);
%% Routing Strategy
for s = 1:Nodes
    for d = 1:Nodes
        path = shortestpath(G1, find(names.Names == s), find(names.Names == d)); %Calcola il percorso più breve tra s e d
        path = names.Names(path);
        for i=1:(length(path)-1)
            fij(path(i), path(i+1)) = fij(path(i), path(i+1)) + Tsd(s, d); %Per ogni hop tra s e d incrementa il traffico sul link i j
            fsd(s, d, path(i), path(i+1)) = Tsd(s, d);
        end
    end
end
result = max(fij(:));