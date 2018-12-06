% Written by Henry Hu for CSCI 1430 @ Brown and CS 4495/6476 @ Georgia Tech

% Find the best fundamental matrix using RANSAC on potentially matching
% points

% 'matches_a' and 'matches_b' are the Nx2 coordinates of the possibly
% matching points from pic_a and pic_b. Each row is a correspondence (e.g.
% row 42 of matches_a is a point that corresponds to row 42 of matches_b.

% 'Best_Fmatrix' is the 3x3 fundamental matrix
% 'inliers_a' and 'inliers_b' are the Mx2 corresponding points (some subset
% of 'matches_a' and 'matches_b') that are inliers with respect to
% Best_Fmatrix.

% For this section, use RANSAC to find the best fundamental matrix by
% randomly sample interest points. You would reuse
% estimate_fundamental_matrix() from part 2 of this assignment.

% If you are trying to produce an uncluttered visualization of epipolar
% lines, you may want to return no more than 30 points for either left or
% right images.

function [ Best_Fmatrix, inliers_a, inliers_b] = ransac_fundamental_matrix(matches_a, matches_b)

ptsPerItr = 8;          #pontos selecionados em cada iteracao  
maxInliers = 0;         #melhores pontos inLayers 0 para iniciar
errThreshold = 0.05;    #erro
numIterations = 1000;   #numero de iteracoes



Best_Fmatrix = zeros(3,3);                    #alocando melhor matriz com 0
xa = [matches_a ones(size(matches_a,1),1)];   
xb = [matches_b ones(size(matches_b,1),1)];   

for i = 1:numIterations                                 #iniciando iteracoes
    ind = randi(size(matches_a,1), [ptsPerItr,1]);      #selecionando pontos aleatorios 
    FmatrixEstimate = estimate_fundamental_matrix(matches_a(ind,:), matches_b(ind,:)); #estimando a matriz fundamental  
    err = sum((xb .* (FmatrixEstimate * xa')'),2);      #calculando o erro da matriz fundamental
    currentInliers = size( find(abs(err) <= errThreshold) , 1);   #numeros de pontos inLayer na atual matriz fundamental 
                                                                  #se o valor absoluto do erro e menor que o  
                                                                  
    if (currentInliers > maxInliers)    #se houver mais ponstos corretos entao a nova matriz e atualizada
       Best_Fmatrix = FmatrixEstimate;    
       maxInliers = currentInliers;
    end    
end

maxInliers
err = sum((xb .* (Best_Fmatrix * xa')'),2);
[Y,I]  = sort(abs(err),'ascend');
inliers_a = matches_a(I(1:30),:);
inliers_b = matches_b(I(1:30),:);

end

