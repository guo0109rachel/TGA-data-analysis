clear
clc
close('all')

%Read the data to what we want
[FileName,PathName,~]=uigetfile('*.mol2*','select the file directly from TGA software');
%readingdata=readtable([PathName,FileName],'FileType','text','Delimiter',{'\s+','\t'},'HeaderLines',15);
readingdata=fileread([PathName,FileName]);
rawdata=string(readingdata);
rawdata=splitlines(rawdata);
%Cut the set to what we want
startline=find(strcmp(rawdata,'@<TRIPOS>ATOM'));
endline=find(strcmp(rawdata,'@<TRIPOS>BOND'));
head=rawdata(1:startline);
last=rawdata(endline:end);
rawdata=rawdata((startline+1):(endline-1));
rawdata=strtrim(rawdata);
p={'           ','          ','         ','        ','       ','      ','     ','    ','   ','  ',' '};
rawdata=replace(rawdata,p,',');
dataset=split(rawdata,',');

%get dimension values
Dimensions=str2double(dataset(:,3:5));

%Calculation for the appropriate charge value
for i=1:length(dataset)
    %find appropriate O atoms
    if (dataset(i,2)=='O')
        for j=1:length(dataset)
            distance(j,1)=(Dimensions(j,1)-Dimensions(i,1))^2+(Dimensions(j,2)-Dimensions(i,2))^2+(Dimensions(j,3)-Dimensions(i,3))^2;
        end
        [~,minIndex]=sort(distance);
        a=minIndex(2);
        if (dataset(a,2)=='H')
           dataset(i,9)='-1.0500';
        else
           dataset(i,9)='-0.9500';
        end
    %find Al atoms
        else if (dataset(i,2)=='AL')
            dataset(i,9)='1.5750';
            %find Si atoms
            else if (dataset(i,2)=='SI')
                dataset(i,9)='2.1000';
                %find H atoms
                else if(dataset(i,2)=='H')
                    dataset(i,9)='0.4250';
                    end
                end
            end
    end
end

newFileName=strrep(FileName,FileName((end-4):end),'_handle.mol2');
fileID=fopen([PathName,newFileName],'w');
FormatSpec='%7s %2s %15s %9s %9s %3s %5s %7s %11s\r\n';
fprintf(fileID,'%s\r\n',head);
fprintf(fileID,FormatSpec,dataset');
fprintf(fileID,'%s\r\n',last);
fclose(fileID);

