clc; clear all; close all; 
% parameters
Ua = 0;
Ub = 4;
Nodes = 40;

rng(25)
Tsd = (Ub-Ua).*rand(Nodes,Nodes)+Ua;
Tsd = Tsd - diag(diag(Tsd));
b = zeros(Nodes,Nodes);
fsd = zeros(Nodes,Nodes,Nodes,Nodes);
fij = zeros(Nodes,Nodes);
G1 = digraph;

%% Creation of the Manhattan Topology
%Horizontal Links of Manhattan Topology
Row = 5; %Manhattan rows
Col = 8; %Manhattan columns
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
for s = 1:Nodes
    for d = 1:Nodes
        path = shortestpath(G1, s, d); %Calcola il percorso più breve tra s e d
        for i=1:(length(path)-1)
            %Per ogni hop tra s e d incrementa il traffico sul link i j
            fij(path(i), path(i+1)) = fij(path(i), path(i+1)) + Tsd(s, d); 
            fsd(s, d, path(i), path(i+1)) = Tsd(s, d);
        end
    end
end
result = max(fij(:)); % Il link più caricato

%local search a priori a posteriori
