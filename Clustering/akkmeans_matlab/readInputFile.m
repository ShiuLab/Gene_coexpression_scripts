function data = readInputFile(filename)
    
fid = fopen(filename);

tline = fgetl(fid);
tline = fgetl(fid);

tabchar = sprintf('\t');
row = 0;
while ischar(tline)
    row = row + 1;
    index = find(tline == tabchar);
    
    num = length(index);
    for i = 1 : num-1
        data(row, i) = str2double(tline(index(i)+1:index(i+1)-1));
    end
    data(row, num) = str2double(tline(index(num)+1:end));
    tline = fgetl(fid);
    if mod(row, 100)==1
        disp(num2str(row))
    end
end
fclose(fid);
