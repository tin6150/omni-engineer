% buildmultitreexml_annot.m is a matlab script that read in fasta file and convert to xml per a template file.  explain what the various block of codes are doing.

%% >> %%

% 🤖 Assistant:
// Certainly! I'll explain the main blocks of code in the buildmultitreexml_annot.m script:


// This script essentially takes sequence data and metadata, performs some clustering and distance calculations, and then generates XML files suitable for phylogenetic analysis, likely for use with BEAST (Bayesian Evolutionary Analysis Sampl



%% << %%

clear
rng(187989)


fasta = fastaread('data/coregenome_snp_aln.fasta');
%%

// 1. Initial setup and data loading:
//    - Clears workspace and sets random number generator seed
//    - Reads a FASTA file containing core genome SNP alignment

// 2. Pairwise distance calculation:
//    - Computes pairwise distances between sequences
//    - Stores distances in the 'pdist' matrix



compregion = 1:length(fasta(1).Sequence);

for a = 1 : length(fasta)
    disp(a)
    for b = a+1:length(fasta)
        indices = find(fasta(a).Sequence(compregion)~=fasta(b).Sequence(compregion));        
        red_indices = find(fasta(a).Sequence(indices)=='-' | fasta(b).Sequence(indices)=='-');
        pdist(a,b) = length(indices)-length(red_indices);
        pdist(b,a) = pdist(a,b);
    end
end
%%


// 3. Clustering sequences:
//    - Groups sequences into clusters based on a minimum distance threshold
//    - Stores cluster members in the 'members' cell array


members = cell(0,0);
already_clustered = [];

min_dist=200;

c = 1;
cl_size = 0;
for i = 1:length(fasta)
    if ~ismember(i, already_clustered)
        newmembers = find(pdist(i, :)<min_dist);
        members{c} = newmembers;
 
        if sum(ismember(newmembers, already_clustered))>0
            error('alal');
        end
        while ~isempty(newmembers)
            tmp = [];
            for j = 1 : length(newmembers)
                tmp = [tmp,find(pdist(newmembers(j), :)<min_dist)];
            end
            tmp = unique(tmp);
            newmembers = tmp(~ismember(tmp,members{c}));
            
            if sum(ismember(newmembers, already_clustered))>0
                error('alal');
            end

            
            members{c} = [members{c} newmembers];
        end
        cl_size = cl_size+length(members{c});
        already_clustered = [already_clustered, members{c}];
        c=c+1;    
    end
end

%%


// 4. Metadata loading:
//    - Reads metadata from a CSV file
//    - Stores sequence IDs, dates, and other information in the 'id' cell array


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






// 5. XML generation loop:
//    - Iterates 3 times (0 to 2) to create multiple XML files
//    - Opens a template XML file and a new output XML file

// 6. XML content generation:
//    - Processes the template file line by line
//    - Replaces placeholder text with generated content:
//      a. insert_sequences: Adds sequence data for each cluster
//      b. insert_trees: Creates tree structures for each cluster
//      c. insert_rootlengths: Adds root length parameters
//      d. insert_init: Initializes random trees for clusters with multiple sequences
//      e. insert_types: Adds sequence type information
//      f. insert_taxa: Adds taxon information
//      g. insert_rootstrees: References trees and root lengths
//      h. insert_clock: Sets the clock rate
//      i. insert_treelikelihood: Adds tree likelihood distributions
//      j. insert_tree_operators: Adds various tree operators for MCMC analysis









%%
for rep = 0:2
    f = fopen('TemplateMultiCluster.xml');
    g = fopen(sprintf('xmls/Ecoli_sf_rep%d.xml', rep), 'w');
    while ~feof(f)
        line = fgets(f);
        if contains(line, 'insert_sequences')
            for j = 1 : length(members)
                fprintf(g,'\t\t<data id="cluster%d" spec="Alignment" name="alignment">\n',j);
                for i = members{j}
                    name = strrep(fasta(i).Header, '.','-');
                    fprintf(g,'\t\t\t<sequence id="seq_%s" spec="Sequence" taxon="%s" totalcount="4" value="%s"/>\n', name, name, fasta(i).Sequence(compregion));
                end
                fprintf(g,'\t\t</data>\n');

            end
            
        elseif contains(line, 'insert_trees')
            for j = 1 : length(members)
                vals = '';
                for i = members{j}
                    ind = find(ismember(id(:,1),fasta(i).Header));
                    vals = [vals ',' id{ind,1} '=' id{ind,2}];
                end
                
                fprintf(g,'\t\t<tree id="Tree.t:%d" spec="beast.evolution.tree.Tree" name="stateNode">\n', j);
                fprintf(g,'\t\t\t<trait id="dateTrait.t:%d" spec="beast.evolution.tree.TraitSet" dateFormat="yyyy-M-dd" traitname="date" value="%s">\n', j, vals(2:end));
                fprintf(g,'\t\t\t\t<taxa id="TaxonSet.%d" spec="TaxonSet">\n', j);
                fprintf(g,'\t\t\t\t\t<alignment idref="cluster%d"/>\n', j);
                fprintf(g,'\t\t\t\t</taxa>\n');
                fprintf(g,'\t\t\t</trait>\n');
                fprintf(g,'\t\t\t<taxonset idref="TaxonSet.%d"/>\n', j);
                fprintf(g,'\t\t</tree>\n');
            end
        elseif contains(line, 'insert_rootlengths')
            for j = 1 : length(members)
                fprintf(g,'\t\t\t<parameter id="rootLength:%d" name="stateNode" dimension="1">1</parameter>\n',j);
            end
        elseif contains(line, 'insert_init')
            for j = 1 : length(members)
                if length(members{j})>1
                    fprintf(g,'\t\t\t<init id="RandomTree.t:%d" spec="beast.evolution.tree.RandomTree" estimate="false" initial="@Tree.t:%d" taxa="@cluster%d">\n',j,j,j);
                    fprintf(g,'\t\t\t\t<populationModel id="ConstantPopulation0.t:%d" spec="ConstantPopulation">\n',j);
                    fprintf(g,'\t\t\t\t\t<parameter id="randomPopSize.t:%d" spec="parameter.RealParameter" name="popSize">0.1</parameter>\n',j);
                    fprintf(g,'\t\t\t\t</populationModel>\n');
                    fprintf(g,'\t\t\t</init>\n');
                end
            end
        elseif contains(line, 'insert_types')
            vals = '';
            for i = 1 :size(id,1)
                vals = [vals ',' id{i,1} '=' id{i,3}];
            end
            fprintf(g, strrep(line, 'insert_types', vals(2:end)));
        elseif contains(line, 'insert_taxa')
            for i = 1 : size(id,1)
                fprintf(g,'\t\t\t<taxon id="%s" spec="beast.evolution.alignment.Taxon"/>\n',id{i,1});
            end
        elseif contains(line, 'insert_rootstrees')
            for j = 1 : length(members)
                fprintf(g,'\t\t\t\t<tree idref="Tree.t:%d"/>\n', j);
                fprintf(g,'\t\t\t\t<rootLength idref="rootLength:%d"/>\n', j);
            end       
        elseif contains(line, 'insert_clock')            
            fprintf(g, strrep(line, 'insert_clock', num2str(10^-6*4800000/length(fasta(1).Sequence)) ));
        elseif contains(line, 'insert_treelikelihood')
            for j = 1 : length(members)
                if length(members{j})>1
                    fprintf(g,'\t\t\t\t<distribution id="treeLikelihood.%d" spec="ThreadedTreeLikelihood" data="@cluster%d" tree="@Tree.t:%d" siteModel="@SiteModel.s:coregenome_snp_aln_masked" branchRateModel="@StrictClock.c:coregenome_snp_aln_masked"/>\n',j,j,j);              
                end
            end
        elseif contains(line, 'insert_tree_operators')
            for j = 1 : length(members)
                fprintf(g,'\t\t\t<operator id="TreeRootScaler.t:%d" spec="ScaleOperator" parameter="@rootLength:%d" scaleFactor="0.5" weight="1"/>\n',j,j);

                if length(members{j})>1            
                    fprintf(g,'\t\t\t<operator id="MascotTreeScaler.t:%d" spec="ScaleOperator" scaleFactor="0.5" tree="@Tree.t:%d" weight="0.30"/>\n',j,j);
                    fprintf(g,'\t\t\t<operator id="MascotTreeRootScaler.t:%d" spec="ScaleOperator" rootOnly="true" scaleFactor="0.5" tree="@Tree.t:%d" weight="0.30"/>\n',j,j);
                    if length(members{j})>2
                        fprintf(g,'\t\t\t<operator id="MascotUniformOperator.t:%d" spec="Uniform" tree="@Tree.t:%d" weight="3.0"/>\n',j,j);
                        fprintf(g,'\t\t\t<operator id="MascotSubtreeSlide.t:%d" spec="SubtreeSlide" tree="@Tree.t:%d" weight="1.50"/>\n',j,j);
                        fprintf(g,'\t\t\t<operator id="MascotNarrow.t:%d" spec="Exchange" tree="@Tree.t:%d" weight="1.50"/>\n',j,j);
                        fprintf(g,'\t\t\t<operator id="MascotWide.t:%d" spec="Exchange" isNarrow="false" tree="@Tree.t:%d" weight=".30"/>\n',j,j);
                        fprintf(g,'\t\t\t<operator id="MascotWilsonBalding.t:%d" spec="WilsonBalding" tree="@Tree.t:%d" weight=".30"/>\n',j,j);
                    end
                end
            end
        else
            fprintf(g, line);
        end
    end
    fclose(f);
    fclose(g);   
end

// 7. File closing:
//    - Closes the template and output XML files
