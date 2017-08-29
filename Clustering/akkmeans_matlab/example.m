% An example to show how approximate  kernel kmeans  is invoked
% Generates two-dimensional Gaussian mixture data of size 500 and finds the 2
% clusters using 50 samples

max_iter = 100;
k =2;
N = 500;
m = 50;
lambda = 1/100;

data = randn(N/2,2);
data = [data; randn(N/2,2) + 5];

plot(data(:,1),data(:,2),'.b');

%Sample m data points
perm = randperm(N);
indices = perm(1:m);

%Compute m x N RBF kernel
Krect = exp(-lambda*(ones(m,1)*sum(data.^2,2)' - 2 * data(indices,:)*data' + sum(data(indices,:).^2,2)*ones(1,N)));
labels=approx_kkmeans(Krect,k,max_iter,indices);
figure
plot(data(labels==1,1),data(labels==1,2),'.r',data(labels==2,1),data(labels==2,2),'.b')