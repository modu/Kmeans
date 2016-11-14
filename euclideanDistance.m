function distance = euclideanDistance(firstPoint, secondPoint)
    distance  = sqrt(sum((firstPoint - secondPoint) .^ 2));
end