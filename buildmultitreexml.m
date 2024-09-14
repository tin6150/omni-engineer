%%% 1. Change the function signature:
%%%```
%%%/edit line 1
function buildmultitreexml(fastafiles, outputfile)
%%%```

%%%2. Add input validation for the new parameters:
%%%```
%%%/insert after line 1
if ~iscell(fastafiles) || isempty(fastafiles)
    error('Input fastafiles must be a non-empty cell array of file names');
end
if ~ischar(outputfile) || isempty(outputfile)
    error('Output file name must be a non-empty string');
end
%%%```


clear
rng(187989)


%## 3. Replace the single FASTA reading with a loop to read multiple files:
%## ```
%## /edit line starting with "fasta = fastaread" until line starting with "compregion ="
%## all_fasta = cell(1, length(fastafiles));
%## for i = 1:length(fastafiles)
    %## all_fasta{i} = fastaread(fastafiles{i});
%## end

% Determine the common region length (assuming all files have the same length)
%## compregion = 1:length(all_fasta{1}(1).Sequence);
%## ```



%// fasta = fastaread('data/coregenome_snp_aln.fasta');
%%

%// compregion = 1:length(fasta(1).Sequence);

all_fasta = cell(1, length(fastafiles));
for i = 1:length(fastafiles)
    all_fasta{i} = fastaread(fastafiles{i});
end

% Determine the common region length (assuming all files have the same length)
compregion = 1:length(all_fasta{1}(1).Sequence);
```


%//for a = 1 : length(fasta)
    %//disp(a)
    %//for b = a+1:length(fasta)
        %//indices = find(fasta(a).Sequence(compregion)~=fasta(b).Sequence(compregion));        
        %//red_indices = find(fasta(a).Sequence(indices)=='-' | fasta(b).Sequence(indices)=='-');
        %//pdist(a,b) = length(indices)-length(red_indices);
        %//pdist(b,a) = pdist(a,b);
    %//end
%//end
%%

%##4. Update the distance calculation to work with multiple FASTA files:
%##```
%## /edit line starting with "for a = 1 : length(fasta)" until line containing "end" after "pdist(b,a) = pdist(a,b);"
pdist = cell(1, length(all_fasta));
for file_idx = 1:length(all_fasta)
    fasta = all_fasta{file_idx};
    pdist{file_idx} = zeros(length(fasta));
    for a = 1 : length(fasta)
        disp(['File ' num2str(file_idx) ', Sequence ' num2str(a)])
        for b = a+1:length(fasta)
            indices = find(fasta(a).Sequence(compregion)~=fasta(b).Sequence(compregion));        
            red_indices = find(fasta(a).Sequence(indices)=='-' | fasta(b).Sequence(indices)=='-');
            pdist{file_idx}(a,b) = length(indices)-length(red_indices);
            pdist{file_idx}(b,a) = pdist{file_idx}(a,b);
        end
    end
end
%##```


%## 5. Update the clustering process to work with multiple FASTA files:
%## ```
%## /edit line starting with "members = cell(0,0);" until line containing "end" after "c=c+1;"
members = cell(1, length(all_fasta));
for file_idx = 1:length(all_fasta)
    fasta = all_fasta{file_idx};
    members{file_idx} = cell(0,0);
    already_clustered = [];
    
    min_dist=200;
    
    c = 1;
    cl_size = 0;
    for i = 1:length(fasta)
        if ~ismember(i, already_clustered)
            newmembers = find(pdist{file_idx}(i, :)<min_dist);
            members{file_idx}{c} = newmembers;
     
            if sum(ismember(newmembers, already_clustered))>0
                error('alal');
            end
            while ~isempty(newmembers)
                tmp = [];
                for j = 1 : length(newmembers)
                    tmp = [tmp,find(pdist{file_idx}(newmembers(j), :)<min_dist)];
                end
                tmp = unique(tmp);
                newmembers = tmp(~ismember(tmp,members{file_idx}{c}));
                
                if sum(ismember(newmembers, already_clustered))>0
                    error('alal');
                end
    
                members{file_idx}{c} = [members{file_idx}{c} newmembers];
            end
            cl_size = cl_size+length(members{file_idx}{c});
            already_clustered = [already_clustered, members{file_idx}{c}];
            c=c+1;    
        end
    end
end
%## ```




%%

f = fopen('data/SF_metadata.csv');fgets(f);

id = cell(0,0);
while ~feof(f)
    line = strsplit(fgets(f), ',');
    id{end+1,1} = [line{1} '_' line{3}];
    date_number = datenum(line{4}, 'mm/dd/yy');
    id{end,2} = datestr(date_number, 'yyyy-mm-dd');   
    id{end,3} = line{3};   

end
fclose(f);

%%


%## 6. Update the XML generation loop to handle multiple FASTA files:
%## ```
%## /edit line starting with "for rep = 0:2" until end of file
for rep = 0:2
    f = fopen('TemplateMultiCluster.xml');
    g = fopen(sprintf('%s_rep%d.xml', outputfile, rep), 'w');
    while ~feof(f)
        line = fgets(f);
        if contains(line, 'insert_sequences')
            for file_idx = 1:length(all_fasta)
                fasta = all_fasta{file_idx};
                for j = 1 : length(members{file_idx})
                    fprintf(g,'\t\t<data id="cluster%d_%d" spec="Alignment" name="alignment">\n', file_idx, j);
                    for i = members{file_idx}{j}
                        name = strrep(fasta(i).Header, '.','-');
                        fprintf(g,'\t\t\t<sequence id="seq_%s" spec="Sequence" taxon="%s" totalcount="4" value="%s"/>\n', name, name, fasta(i).Sequence(compregion));
                    end
                    fprintf(g,'\t\t</data>\n');
                end
            end
        elseif contains(line, 'insert_trees')
            for file_idx = 1:length(all_fasta)
                fasta = all_fasta{file_idx};
                for j = 1 : length(members{file_idx})
                    vals = '';
                    for i = members{file_idx}{j}
                        ind = find(ismember(id(:,1),fasta(i).Header));
                        vals = [vals ',' id{ind,1} '=' id{ind,2}];
                    end
                    
                    fprintf(g,'\t\t<tree id="Tree.t:%d_%d" spec="beast.evolution.tree.Tree" name="stateNode">\n', file_idx, j);
                    fprintf(g,'\t\t\t<trait id="dateTrait.t:%d_%d" spec="beast.evolution.tree.TraitSet" dateFormat="yyyy-M-dd" traitname="date" value="%s">\n', file_idx, j, vals(2:end));
                    fprintf(g,'\t\t\t\t<taxa id="TaxonSet.%d_%d" spec="TaxonSet">\n', file_idx, j);
                    fprintf(g,'\t\t\t\t\t<alignment idref="cluster%d_%d"/>\n', file_idx, j);
                    fprintf(g,'\t\t\t\t</taxa>\n');
                    fprintf(g,'\t\t\t</trait>\n');
                    fprintf(g,'\t\t\t<taxonset idref="TaxonSet.%d_%d"/>\n', file_idx, j);
                    fprintf(g,'\t\t</tree>\n');
                end
            end
        % Continue updating the rest of the XML generation process similarly,
        % replacing single-file references with nested loops for all files
        % ...
        else
            fprintf(g, line);
        end
    end
    fclose(f);
    fclose(g);   
end
%## ```

