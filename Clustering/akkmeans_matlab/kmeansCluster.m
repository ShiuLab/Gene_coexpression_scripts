function kmeansCluster(data, max_iter, k, ntry)
% kmeans cluster
N = size(data, 1);
m = round(N/10);
lambda = 1/100;

%Sample m data points
perm = randperm(N);
indices = perm(1:m);

%Compute m x N RBF kernel
Krect = exp(-lambda*(ones(m,1)*sum(data.^2,2)' - 2 * data(indices,:)*data' + sum(data(indices,:).^2,2)*ones(1,N)));
labels=approx_kkmeans(Krect,k,max_iter,indices);


fid = fopen([num2str(k), '_', num2str(ntry), '_labels.txt'], 'w');
num = length(labels);
for i = 1 : num
    fprintf(fid, '%d\t', labels(i));
end
fclose(fid);






        


