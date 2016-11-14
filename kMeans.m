%First K-Means algorithm is implemented in Matlab
%Then few matrices APIs to evaluate the algorithm are used. 
k = 4;    
dataSetLength = 150;    %Total number of data, numbers of rows if each row is one observation.
numberOfFeature = 4;    %Number of Features or Dimensions  
numberOfIterations = 10;    %Stop recalculating new centres after this many iterations. Efficiency sake. 

l = load('data.mat');   
data = l.L.data;

%Initializing the memberShipMatrix and centers
memberShipMatrix = zeros(dataSetLength,k);
centroidMatrix = zeros(k,numberOfFeature);
objectiveFunctionSum = zeros(1,numberOfIterations);
iterations = 1;

% Randomly assinging centeroid 
%%{
for i = 1:k
    %TODO: try to randomnly change i in data matrix for randomization
    %As of below it will use intial K values from data matrix as centers.
    centroidMatrix(i,:) = data(i,:);
end
%%}
%Use below intialization if you explicitly want to initalize centers/seeds.
%{
centroidMatrix(1,:) = data( 1, :);
centroidMatrix(2,:) = data( 68, :);
centroidMatrix(3,:) = data( 38, :);
centroidMatrix(4,:) = data( 100, :);
 %}

%Loop until iterations < numberOfIterations
while 1
    memberShipMatrix = zeros(dataSetLength,k);
    % Loop through each data point to find distance between current
    % center(centroidMatrix) and data point. Assing data point to center it
    % is nearest to(min).
    %calculating MemberShipMatrix(which shows which cluster a point belongs to)  and objective function

    for i = 1:size(data,1)
        currentDistance = zeros(1, k);
        for j = 1:size(centroidMatrix,1)
            currentDistance(j) = euclideanDistance(data(i,:), centroidMatrix(j,:));
        end
        [minimumDistance, I] = min(currentDistance);
        memberShipMatrix(i,I) = 1;
        objectiveFunctionSum(iterations) = objectiveFunctionSum(iterations) + minimumDistance ;
    end
    clusterSum = zeros(k,numberOfFeature);
    numberOfElementsPerCluster = zeros(1,k);
    
    %Finding Number of data points in each cluster 
    %Also calulating sum of each point in each cluster across dimensions( to calculate new
    %center)
    for y = 1:size(memberShipMatrix,1)
        for z = 1:size(memberShipMatrix,2)
            if(memberShipMatrix(y, z)==1)
                clusterNumber = z;
                numberOfElementsPerCluster(z) = numberOfElementsPerCluster(z) + 1;
            end        
        end
        clusterSum(clusterNumber,:) = clusterSum(clusterNumber,:) ... 
            + data(y,:);    
    end
    
    newCentroid = zeros(k,numberOfFeature);
    %Calculating new Centroid    
    for u = 1:size(newCentroid,1)
        newCentroid(u,:) = clusterSum(u,:)/numberOfElementsPerCluster(u);
    end
    
    centroidMatrix = newCentroid;    
    iterations = iterations + 1; 
    
    if iterations > numberOfIterations
        break;
    end
    
end

objectiveFunctionSum
numberOfElementsPerCluster

figure();
plot(2:numberOfIterations, objectiveFunctionSum(2:numberOfIterations) )
titleString =sprintf('Objective Function vs number of iterations for K = %d', k);
title(titleString)
xlabel('Number of Iterations') % x-axis label
ylabel('Objective Function') % y-axis label




%l = load('data.mat');
%Calculating Dunns Index for different values of K 
for k=2:10
    [X,C] = kmeans(l.L.data,k);
    % Generating memebership matrix 
    A = zeros(k,150);

    for i=1:150
        x = X(i);
        A(x,i) = 1;
    end
    distM=squareform(pdist(l.L.data));
    v(k) = dunnsM(k, distM, X);
end
figure();
plot(2:10,v(2:10))
title('Dunn Index vs Number of cluster');
xlabel('Number of Clusters');
ylabel('Dunn Index ');


%DavisBouldin Index using Kmeans
eva = evalclusters((l.L.data),'kmeans','DaviesBouldin','klist',[2:10])

figure();
plot(2:10,eva.CriterionValues)
title('Davies&Bouldin Index vs Number of cluster');
xlabel('Number of Clusters');
ylabel('Davies&Bouldin Index');

tree = linkage(l.L.data,'average');
figure()
dendrogram(tree)





