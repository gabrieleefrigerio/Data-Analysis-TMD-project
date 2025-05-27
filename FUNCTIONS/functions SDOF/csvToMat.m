function [] = csvToMat(folderName)

baseFolder = pwd;
csvPath = [folderName,'\csv'];
cd(csvPath)
namecsv = uigetfile('*.csv',"MultiSelect","on");


for ii = 1:length(namecsv)
    temp = readmatrix(namecsv{ii});
    hammer(:,ii) = temp(2:end,end);
    temp = temp(2:end,2:end-1);

    for jj = 1:size(temp,1)*size(temp,2)     % remove outiliers
        if abs(temp(jj))>1e5
            if jj ~= 1
                temp(jj) = temp(jj-1);
            else 
                temp(jj) = 0;
            end
        end
    end

    data(:,:,ii) = temp;
end

cd('..')
temp = readmatrix('setup.txt');
fsamp = temp(end,end);

dt = 1/fsamp;
time = 0:dt:(size(hammer,1)-1)*dt;
time = time';

save('data.mat','fsamp','time','hammer','data')
cd(baseFolder)