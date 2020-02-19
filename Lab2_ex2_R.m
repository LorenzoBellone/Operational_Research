%% Random generated topology
clc; clear all; close all;
Delta = 4;

% roba in più ^_^
Ua_HT = 5;  % boundarys for the traffic tsd
Ub_HT = 15;
Ua_LT = 0;
Ub_LT = 3;
HT = 0.1;
LT = 1-HT;

rng(25)
Nodes = 40;
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
% roba in più ^_^

b = zeros(Nodes,Nodes);
fsd = zeros(Nodes,Nodes,Nodes,Nodes);
fij = zeros(Nodes,Nodes);
DeltaTx = zeros(Nodes,1);
DeltaRx = zeros(Nodes,1);
G2 = digraph;

for i=1:Nodes
    for j=1:Nodes
        ok = 0;
        while ok == 0
            % Prendi due nodi a caso nel grafo (tra 1 e Nodes)
            s = randi([1, Nodes]);
            d = randi([1, Nodes]);
            if (s ~= d) && (b(s, d) == 0) % Se i nodi che hai selezionato lo erano già stati o hai preso due numeri ugali riprova
                ok = 1;
                if (DeltaTx(s, 1) < Delta) && (DeltaRx(d, 1) < Delta) % Condizione sul delta
                    G2 = addedge(G2, s, d); % Se tutto va bene aggiungi il link
                    b(s, d) = 1;
                    DeltaTx(s, 1) = DeltaTx(s, 1) +1;
                    DeltaRx(d, 1) = DeltaRx(d, 1) +1;
                end
            end
        end
    end
end
figure(2);
plot(G2);
% Stessa risoluzione di prima..
for s = 1:Nodes
    for d = 1:Nodes
        path = shortestpath(G2, s, d);
        for i=1:(length(path)-1)
            fij(path(i), path(i+1)) = fij(path(i), path(i+1)) + Tsd(s, d);
            fsd(s, d, path(i), path(i+1)) = Tsd(s, d);
        end
    end
end
result = max(fij(:));