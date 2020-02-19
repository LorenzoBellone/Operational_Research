clc; clear all; close all; 
% parameters
Delta = 4;

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

Tsd=zeros(Nodes)

%generate a random number such that 90% of case Tsd = Tsd_LT
for row = 1:Nodes
    for col = 1:Nodes
        random_decision = rand;
        if random_decision < LT
            Tsd(row,col) = Tsd_LT(row,col)
        else
            Tsd(row,col) = Tsd_HT(row,col)
        end     
    end
end

% roba in più ^_^

b = zeros(Nodes,Nodes);
fsd = zeros(Nodes,Nodes,Nodes,Nodes);
fij = zeros(Nodes,Nodes);
DeltaTx = zeros(Nodes,1);
DeltaRx = zeros(Nodes,1);

%mahattan 4x5
% 1   2   3   4   5 
% 6   7   8   9   10
% 11  12  13  14  15
% 16  17  18  19  20
for i = 1:Nodes
    r = mod(i,5);
    %horizontal
    if r ~= 0
        b(i,i+1) = 1;
    else
        b(i,i-4) = 1;
    end
    if r == 1
        b(i,i+4) = 1;
    else
        b(i,i-1) = 1;
    end
    %vertical
    if i < 16 
        b(i, i+5) = 1;
    else
        b(i,i-15) = 1;
    end
    if i > 5 
        b(i, i-5) = 1;
    else
        b(i,i+15) = 1;
    end
end

%check sum(b,1) sum(b,2)

disp('End For')

G = digraph(b);
plot(G);
for s = 1:Nodes
    for d = 1:Nodes
        path = shortestpath(G, s, d); %Calcola il percorso più breve tra s e d
        for i=1:(length(path)-1)
            fij(path(i), path(i+1)) = fij(path(i), path(i+1)) + Tsd(s, d); %Per ogni hop tra s e d incrementa il traffico sul link i j
            fsd(s, d, path(i), path(i+1)) = Tsd(s, d);
        end
    end
end
result = max(fij(:)) % Il link più caricato