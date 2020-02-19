%% Random generated topology
clc; clear all; close all;
Delta = 4;
Ua = 0;
Ub = 4;

rng(25)
Nodes = 40;
Tsd = (Ub-Ua).*rand(Nodes,Nodes)+Ua;
Tsd = Tsd - diag(diag(Tsd));
b = zeros(Nodes,Nodes);
fsd = zeros(Nodes,Nodes,Nodes,Nodes);
fij = zeros(Nodes,Nodes);
DeltaTx = zeros(Nodes,1);
DeltaRx = zeros(Nodes,1);
G = digraph;

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
                    G = addedge(G, s, d); % Se tutto va bene aggiungi il link
                    b(s, d) = 1;
                    DeltaTx(s, 1) = DeltaTx(s, 1) +1;
                    DeltaRx(d, 1) = DeltaRx(d, 1) +1;
                    
                end
            end
            end
    end
end
plot(G);
% Stessa risoluzione di prima..
for s = 1:Nodes
    for d = 1:Nodes
        path = shortestpath(G, s, d);
        for i=1:(length(path)-1)
            fij(path(i), path(i+1)) = fij(path(i), path(i+1)) + Tsd(s, d);
            fsd(s, d, path(i), path(i+1)) = Tsd(s, d);
        end
    end
end
result = max(fij(:));