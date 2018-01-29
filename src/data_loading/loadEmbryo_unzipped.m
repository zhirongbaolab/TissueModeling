function [embinfo,errors ]...
    = loadEmbryo_unzipped(file,starttime,endtime)

% ----------------------------------------------
% Given the first file and an int representing the final timepoint,
% iterate over all files (shared naming paradigm) and process the nuclei
% to load an entire embryo
% ----------------------------------------------

errors=[];

%templocation='temp_unzip\';
%{
%unzip zipfile to temp file


if ~exist(zipfile,'file')
    errors.zipfilemissing=1;
    return
end
try
    unzip(zipfile,templocation);
catch exception
    errors.zipfilecorrupted=1;
    return
end
%}


% examine all divisions and see if A cell is more anterior based on long
% axis of variance and axis label in auxinfo file


embinfo=[];
allpoints=[];

%load entire embryo
templocation=file{:};
for t=starttime:endtime % iterate over all timepoints
    
    % load t current, first expecting no file extension
    nuclei=[templocation,num2str(t,'t%d'),'-nuclei'];
    if ~exist(nuclei, 'file')
        errors.nucleiFiles=['expected nuclei file missing time: ',nuclei];
        return
    end


    % process all nuclei at this timepoint retrieving cell data, cell names, and forced name (if any)
    [celldata,cellnames,forcednames]=readnuclei_allfields(nuclei);

    % save embryo info
    embinfo(t).celldata=celldata;
    embinfo(t).cellnames=cellnames;
    p1_sucessors=celldata(:,9:10);
    celllocations=celldata(:,4:6);%pull nuclei from labeled data

    % find 4 cell stage cells ABa, P2
    ABaindex=find(strcmp('ABa',cellnames));
    P2index=find(strcmp('P2',cellnames));

    % WHAT IS GOING ON HERE
    if(~isempty(ABaindex)&&~isempty(P2index))
        ABapos=celllocations(ABaindex,:);
        P2pos=celllocations(P2index,:);
    end
    
    
    % check if this is before the last timepoint
    if(t<endtime)

        %load t+1, the time point after this
        nuclei=[templocation,num2str(t+1,'t%d'),'-nuclei'];
        if ~exist(nuclei)
            errors.nucleiFiles=['expected nuclei file missing time: ',nuclei];
            return
        end

        % load data from the time point after the current time. cell data, cell names, and their corresponding indices
        [celldata_c2,cellnames2]=readnuclei(nuclei);
        indiciesp2=celldata_c2(:,2);
        
        % translate acetree sucessor indicies to array indicies
        s=size(p1_sucessors); % the number of ancestors
        p1_sucessors_t=[];
        
        % iterate from 1 to the the number of items in p1_successors (s is the size of this data array, and the second value represents WHAT)
        for j=1:s(1);

            % assume the successors are INVALID
            suc1=-1;
            suc2=-1;

            % compute an AND operation on [the remaining successors (specifically the VALID/INVALID 1/-1 index] AND [-1 i.e. invalid])
            % if there are any valid cells in the rest of the array of p1 successors, this operation will be TRUE, otherwise FALSE
            if(p1_sucessors(j,1)~=-1)
                suc1=find(indiciesp2==p1_sucessors(j,1));
                if isempty(suc1)
                    suc1=-1;
                    'invalid nucleus pointed to by valid nucleus, should not happen in new files'
                end
            end
            if(p1_sucessors(j,2)~=-1)
                suc2=find(indiciesp2==p1_sucessors(j,2));
                if isempty(suc2)
                    suc2=-1;
                    'invalid nucleus pointed to by valid nucleus, should not happen in new files'
                end
            end
            p1_sucessors_t=[p1_sucessors_t;suc1,suc2];
        end
    else
        p1_sucessors_t=p1_sucessors;% dont try and translate if endtime
    end
    embinfo(t).finalpoints=celllocations;
    embinfo(t).suc=p1_sucessors_t;
    embinfo(t).names=cellnames;
    embinfo(t).forcednames=forcednames;
end

end

