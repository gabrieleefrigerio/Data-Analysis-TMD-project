function [DATA, HAMMER] = clearFaults(data, hammer, folderName) 

faulty_acc = readmatrix([folderName,'\faulty_acc.txt']);
data(:,faulty_acc,:) = [];

faulty_tests = readmatrix([folderName,'\faulty_tests.txt']);
data(:,:,faulty_tests) = [];
hammer(:,faulty_tests) = [];

DATA = data;
HAMMER = hammer;

end