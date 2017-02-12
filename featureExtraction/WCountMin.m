function Y = WCountMin(X, thresh)
    Y = sum(X<thresh,2);
end