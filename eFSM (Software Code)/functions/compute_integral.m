function integral_value = compute_integral(MSH,weight, ArchBool, n, fL, Hip_Flexor)
% compute integral of f on triangular mesh MSH via midpoint rule
% MSH.p : n x 2: list of coordinates of nodes
% MSH.TR : n x 3: list of elements
% f: n x 1: list of nodal values of function to be integrated

% unit test:
%     load('MSH_test.mat') % circle of radius 1
%     f = @(x) ones(size(x));
%     f_nodalvalues = f(MSH_test.p);
%     ComputeIntegral(MSH_test,f_nodalvalues) % should be the area of unit
%                                             % circle (roughly  pi)
if ArchBool || Hip_Flexor
    weights = weight(MSH.p(:,1), MSH.p(:,2));
    weights = n * weights;
end
% plot3(MSH.p(:,1), MSH.p(:,2), weights,'r.','MarkerSize', 10); hold on
integral_value = 0;
for i = 1:size(MSH.TR,1)
    % add local mass matrices
    coordinates = MSH.p(MSH.TR(i,:),:);
    eAB = coordinates(2,:) - coordinates(1,:);
    eAC = coordinates(3,:) - coordinates(1,:);
    % det/2 = area of triangle
    det_D = abs(eAB(1)*eAC(2) - eAB(2)*eAC(1));
    % f evaluated at midpoint ( by taking mean value )
    %         f_midpoint = 1 / 3 * sum(f(MSH.p(MSH.TR(i,:),1), MSH.p(MSH.TR(i,:),2)));
    % additionally take weights in account !
% % % % % % % % % % % %     vals = f(MSH.p(MSH.TR(i,:),1), MSH.p(MSH.TR(i,:),2));
    vals = fL(MSH.TR(i,:));
%     norm(vals - valsN)
    if ArchBool || Hip_Flexor
        f_midpoint = 1 / 3 * sum(weights(MSH.TR(i,:)) .* vals);
    else
        f_midpoint = 1 / 3 * sum(vals);
    end
    % midpoint rule
    integral_value = integral_value + det_D / 2 * f_midpoint;
end

end