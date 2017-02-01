function Y = Wpca1(X)
    [~,Y] = pca(X,'NumComponents',1);
end