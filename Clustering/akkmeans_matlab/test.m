

% read the data into .mat file, please specify the name of the file

ntry = 10;
max_iter = 100;
k = [5, 10, 25, 50, 100, 200, 300, 400, 500, 1000, 2000];

filename = '/mnt/home/uygunsah/1_Expression_Database/expression_matrix/5_Stress_FC';
data = readInputFile(filename);

for try0 = 1 : ntry
    for i = 1 : length(k)
        tic
        k0 = k(i);
        kmeansCluster(data, max_iter, k0, try0)
        toc
    end
end



        


