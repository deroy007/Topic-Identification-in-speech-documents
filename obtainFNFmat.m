function [FNF,S_exp,S_obt,B2,B23] = obtainFNFmat(phnFilename,T1,winS) 

phnFile = importdata(phnFilename);
split1 = cellfun(@(x) strsplit(x,' '), phnFile( ),'uni',0);%online solution
 command=sprintf('cat %s | cut -d " " -f 1-2 > col12.txt',phnFilename);
system (command);
B2 = importdata('col12.txt');
command2=sprintf('cat %s | cut -d " " -f 3- > col3.txt',phnFilename);
system (command2);
FDC = T1; %fricative discriminative coefficient
% sFDC = zeros(winS*length(FDC),1);
s11 = repmat(FDC,[1,winS]);
s12 = s11';
sFDC = s12(:);
for v=1:length(B2)
    a2 = B2(v,:);
    if (v==length(B2)-1 ) && (a2(2)>length(sFDC))
        FNF(v,:) = mean(sFDC(a2(1)+1:end));
    else
        if (v==length(B2)) && ((a2(1)>length(sFDC)))
                FNF(v,:) = 0;
            
        else
            FNF(v,:) = mean(sFDC(a2(1)+1:end));
        end
        
    end

   clear a2; 
%         FNF(v,:) = mean(sFDC(a2(1)+1:a2(2)));
end
    
% end
fricvsAll = zeros(length(FNF),1);
fricInx = find(FNF>0.2);
fricvsAll(fricInx) = 1;
% figure,stem(fricvsAll)
B23 = importdata('col3.txt');
B23 = [B23;' '];
fafList = {'sh','s','th','f','ch','zh','z','dh','v','jh'};
SS_inx = [];
for lis = 1:length(fafList)
strS = strfind(B23,fafList{1,lis});
S_inx = find(~cellfun(@isempty,strS));
SS_inx = [SS_inx;S_inx(:)]; %fricatives + affricates index
S_exp(lis) = length(S_inx);
S_obt(lis) = sum(fricvsAll(S_inx));
end
NF_inx = setdiff([1:length(B23)-1],SS_inx); %finds index other than fricatives and affricates
S_exp(lis+1) = length(NF_inx);
S_obt(lis+1) = sum(fricvsAll(NF_inx));