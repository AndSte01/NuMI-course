function quadplot(nodes,elements,sol)
arguments
    nodes (:,2) double
    elements (:,4) double
    sol (:,1) double
end
%quadplot  Plots a 3D surface
%
% Inputs:
%   nodes:    Nodes
%                 [x_1, y_1; ...; x_n, y_n] (Index is row in matrix)
%   elements: Elements (ElementId, [local_node_id])
%   sol:      Solution vector (Row is node_id)
%
% Exercise 1.3
%
% © 2024, Andreas Steger

%% umrechnen der Vierecke in Dreiecke
% zählen im Uhrzeigersinn von unten links aus (ul: 1, ur: 2, or: 3, ol: 4)
% ein Viereck zerfällt folglich in zwei Dreiecke: (1, 2, 4) und (2, 3, 4)

% abfragen wie viele Vierecke vorhanden sind
n_elements = height(elements); % Hilfsvariable ist hier leicht schneller (tradeoff Speicher gegen Berechnungszeit gerechtfertigt)

% erzeugen der leeren Matrix für die neuen Dreiecke, wichtig weil die
% dynamische Änderung der Größe einer Matrix schlecht für die Performance
% ist
triangles = zeros([n_elements*2 3]);

% Möglichkeit 1: berechnet die Dreiecke in einer sortierten Reihenfolge
% (jene zwei Dreiecke welche zuvor ein Viereck bildeten sind hintereinander
% angeordnet in der `triangles` Matrix)
%for i=1:n_elements
%    triangles(2*i,:) = [elements(i,1), elements(i,2), elements(i,4)];
%    triangles(2*i-1,:) = [elements(i,2), elements(i,3), elements(i,4)];
%end

% Möglichkeit 2: Dreiecke werden fragmentiert abgespeichert (in der ersten
% Hälfte der `triangles` Matrix befinden sich die _ersten_ Dreiecke der
% Vierecke, in der zweiten Hälfte der Matrix analog die _zweiten_) dadurch
% erleichtert sich aber die Berechnung der Zeilenindizes, sie kondensiert
% sich auf die Berechnung eines einfachen Offsets für die zweiten Dreiecke
for i=1:n_elements
    triangles(i,:) = [elements(i,1), elements(i,2), elements(i,4)];
    triangles(n_elements + i,:) = [elements(i,2), elements(i,3), elements(i,4)];
end

%% plotten der generierten Dreiecke
trisurf(triangles, nodes(:,1), nodes(:,2), sol);
end