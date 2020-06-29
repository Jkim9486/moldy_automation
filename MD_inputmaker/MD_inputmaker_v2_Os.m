clear all
load constant
Kcal2KJ=4.184;% 1Kcal/mol = 4.184 KJ/mol
rm2sigma=1/(2^(1/6));% sigma=rm/(2^(1/6))
KbK2KJ=8.314472e-3;% Kb*K=R*K/Na=8.314472 J/mol;
% Moldy LJ input unit is 4 * KJ/mol and Angstrom (sigma)
% [Kcal2KJ,rm2sigma,KbK2KJ,LJ]=allocate;
%%% User options
numofsolv=512;
solvent_name={
    'Cyclohexane'
     };
solvent_list={
    'input/Cyclohexane_md.xyz'
    };

coord_list={
  'input/Os3CO12_DFT_md.xyz'
  'input/Os3CO11axi_DFT_md.xyz'
  'input/Os3CO11bri_DFT_md.xyz'
  'input/Os3CO11equ_DFT_md.xyz'
  'input/Os3CO10equequ_md.xyz'
  'input/Os3CO10axiequ_md.xyz'
  'input/Os3CO10axiaxi_md.xyz'
  'input/CO_DFT_md.xyz'
  };

solute_names={
    'Os3CO12'
    'Os3CO11_axi'
    'Os3CO11_bri'
    'Os3CO11_equ'
    'Os3CO10_equequ'
    'Os3CO11_axiequ'
    'Os3CO11_axiaxi'
    'CO'
    };

charge{1}=[
    -1.33468
    -1.33435
    -1.33435
    0.75757
    0.75757
    0.78058
    0.78058
    0.75752
    0.75752
    0.78056
    0.7806
    0.75752
    0.75752
    0.7806
    0.78056
    -0.44192
    -0.44192
    -0.42889
    -0.42889
    -0.44198
    -0.44198
    -0.42892
    -0.42897
    -0.44198
    -0.44198
    -0.42897
    -0.42892
    ];

charge{2}=[
    -1.44087
-1.44087
-0.65365
0.75582
0.75452
0.78095
0.77528
0.75582
0.75452
0.77528
0.78095
0.73707
0.70631
0.70631
-0.43756
-0.43897
-0.43249
-0.42664
-0.43756
-0.43897
-0.42664
-0.43249
-0.42077
-0.42767
-0.42767
    ];

charge{3}=[
    -1.23402
-1.01061
-1.00328
0.76492
0.76492
0.78556
0.78041
0.69767
0.69767
0.73182
0.70289
0.70289
0.53912
0.72558
-0.43359
-0.43359
-0.42537
-0.42634
-0.42071
-0.42071
-0.41807
-0.42073
-0.42073
-0.40839
-0.41729
];

charge{4}=[
    -0.6868
-1.23663
-1.42219
0.66066
0.66066
0.66981
0.76785
0.76785
0.78581
0.78661
0.74076
0.74076
0.76769
0.77212
-0.42391
-0.42391
-0.46146
-0.42275
-0.42275
-0.41351
-0.42168
-0.45488
-0.45488
-0.43611
-0.43914
    ];

charge{5}=[
    -0.2917
-1.3032
-1.3032
0.7727
0.77328
0.79496
0.78573
0.77328
0.77269
0.78573
0.79497
-0.46485
-0.45655
-0.44692
-0.44959
-0.45655
-0.46486
-0.44959
-0.44692
0.58576
0.58576
-0.44546
-0.44546
    ];

charge{6}=[
    -0.77806
-1.33154
-0.68212
0.70268
0.6828
0.68762
0.75078
0.75882
0.77076
0.79313
0.69992
0.65704
0.57884
-0.40438
-0.40395
-0.43689
-0.43036
-0.42495
-0.41609
-0.42734
-0.4309
-0.44241
-0.47342
];

charge{7}=[
    -0.12972
-1.30019
-1.54942
0.77626
0.77173
0.78891
0.78262
0.71031
0.72252
0.75828
0.7768
-0.41694
-0.41498
-0.41311
-0.41875
-0.45152
-0.45372
-0.4323
-0.43959
0.61806
0.63065
-0.45951
-0.45639
];

charge{8}=[
    0.5161
    -0.5161
    ];

charge_solvent{1}=[
    -0.12
    -0.12
    -0.12
    -0.12
    -0.12
    -0.12
    0.06
    0.06
    0.06
    0.06
    0.06
    0.06
    0.06
    0.06
    0.06
    0.06
    0.06
    0.06
    ];

label_change{1}=[
    "Hsolvent"
    "Csolvent"
    ];

label_change{2}=[
    "H_99"
    "C_99"
];

%% main func.
input_checker_batch(charge,coord_list,solute_names);
fprintf('               Solute coordinations process.\n')
coord=read_coord(coord_list,constant);
coord=LJ_weight_setter_batch(coord,charge,constant,coord_list);
fprintf('               Solvent coordination process.\n')
coord_solvent=read_coord(solvent_list,constant);
coord_solvent=LJ_weight_setter_batch(coord_solvent,charge_solvent,constant,solvent_list);



% coord_solvent_rearrange=coord_rearranger(coord_solvent,LJ,LJ_elements);
% coord_rearrange=coord_rearranger(coord,LJ,LJ_elements);
%%% Pair generator
coord=pair_generator_batch(coord,coord_solvent);

results_maker(coord,coord_solvent,solute_names,solvent_name,numofsolv,label_change);

% [Kcal2KJ,rm2sigma,KbK2KJ,LJ]=allocate;
%% Specifications %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LJ(1,:)=[33.2*4*KbK2KJ;3.5]; %% C at cyclohexane
% LJ(2,:)=[15.09*4*KbK2KJ;2.5]; %% H at cyclohexane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% File write
% dlmwrite('LJ.txt',pair,'\t');
% function [Kcal2KJ,rm2sigma,KbK2KJ,LJ]=allocate
% Kcal2KJ=4.184;% 1Kcal/mol = 4.184 KJ/mol
% rm2sigma=1/(2^(1/6));% sigma=rm/(2^(1/6))
% KbK2KJ=8.314472e-3;% Kb*K=R*K/Na=8.314472 J/mol;
% LJ=nan(1,2); %% Fix
% end
function coord=LJ_weight_setter_batch(coord,charge,constant,list)
numofcoord=length(coord);
for coord_index=1:numofcoord
    disp('*=============================================================*')
    if coord_index==1
        fprintf('               The 1st coordination: %s\n',list{coord_index})
    elseif coord_index==2
        fprintf('               The 2nd coordination: %s\n',list{coord_index})
    elseif coord_index==3
        fprintf('               The 3rd coordination: %s\n',list{coord_index})
    else
        fprintf('               The %dth coordination: %s\n',coord_index,list{coord_index})
    end
    coord{coord_index}.charge=charge{coord_index};
    coord{coord_index}=LJ_weight_setter(coord{coord_index},constant);
end
% disp('*=============================================================*')
% disp('               Input parameter check and process completed.')
disp('*=============================================================*')
    function coord_temp=LJ_weight_setter(coord_temp,constant)
        numofatom=length(coord_temp.num);
        coord_temp.weight=nan(numofatom,1);
        %         coord_temp.LJ=nan(numofatom,2);
        for atom_index=1:numofatom
            coord_temp.weight(atom_index)=constant.masses(coord_temp.num(atom_index));
            %             coord_temp{coord_index}.charge=charge_temp{coord_index};
            %             coord_temp.LJ(atom_index,:)=constant.LJ(coord_temp.num(atom_index),:);
        end
        atom_all=string(coord_temp.textdata);
        atom_unique=unique(atom_all);
        numofLJ=length(atom_unique);
        num_atom_unique=sym2num_batch(atom_unique,constant);
        if numofatom~=length(coord_temp.charge)
            error('*=============================================================*\nError !!\nThe number of atom is not match.\nCheck coordination and charge.\n*=============================================================*',0);
        end
        LJ=nan(numofLJ,2);
        for LJ_index=1:numofLJ
%             posi_temp=find(constant.num==num_atom_unique(LJ_index));
            %             posi_temp=find(constant.==atom_unique(LJ_index));
%             posi=posi_temp(1);
            LJ(LJ_index,:)=constant.LJ(num_atom_unique(LJ_index),:);
        end
        coord_temp.LJ=LJ;
        coord_temp.LJ_atom=atom_unique;
        fprintf('               The number of LJ potential: %d\n',numofLJ)
        fprintf('               The number of atom in coordination: %d\n',numofatom)
        function num=sym2num_batch(textdata,constant)
            numofatom_temp=length(textdata);
            num=nan(numofatom_temp,1);
            for atom_index_temp=1:numofatom_temp
                num(atom_index_temp)=sym2num(textdata{atom_index_temp},constant);
            end
            function num=sym2num(element,constant)
                element_cut=sym_cutter(element);
                num=find(strcmp(element_cut,constant.elements));
                
                function sym_cut=sym_cutter(sym)
                    posi_cut=strfind(sym,'_');
                    if isempty(posi_cut)
                        sym_cut=sym;
                    else
                        sym_cut=sym(1:posi_cut-1);
                    end
                end
            end
        end
    end
end


function input_checker_batch(charge,coord_list,solute_names)
    numofcoord=length(coord_list);
    if ~all([numofcoord==length(solute_names),numofcoord==length(charge)])
        error('*=============================================================*\nError !!\nThe number of coordination is not match.\nCheck coord list and solute names.\n*=============================================================*',0);
    end
    disp('*=============================================================*')
    fprintf('               The number of coordination: %d\n',numofcoord)
    disp('               Input parameter check completed.')
    disp('*=============================================================*')
end

% function [numofLJs,numofcoord]=input_checker_batch(LJ_elements,four_epsilon,sigma,charge,atomic_weight,coord_list,solute_names)
%     numofcoord=length(coord_list);
%     numofLJs=cell(numofcoord,1);
%     if ~(numofcoord==length(solute_names))
%         error('*=============================================================*\nError !!\nThe number of coordination is not match.\nCheck coord list and solute names.\n*=============================================================*',0);
%     end
% for coord_index=1:numofcoord
%     numofLJs{coord_index}=input_checker(LJ_elements{coord_index},four_epsilon{coord_index},sigma{coord_index},charge{coord_index},atomic_weight{coord_index},coord_index);
% end
%      fprintf('               The number of coordination: %d\n',numofcoord)
%     disp('*=============================================================*')
% 
% function numofLJ=input_checker(LJ_elements,four_epsilon,sigma,charge,atomic_weight,coord_index)
%     numofLJ=length(LJ_elements);
%     if ~all([numofLJ==length(four_epsilon),numofLJ==length(sigma),numofLJ==length(charge),numofLJ==length(atomic_weight)])
%         error('*=============================================================*\nError !!\nThe number of LJ is not match.\nCheck LJ_elements,four epsilon, sigma, charge and atomic weight.\n*=============================================================*',0);
%     end
%     disp('*=============================================================*')
%     if coord_index==1
%     fprintf('               The 1st coordination\n')
%     elseif coord_index==2
%     fprintf('               The 2nd coordination\n')
%     elseif coord_index==3
%     fprintf('               The 3rd coordination\n')
%     else
%     fprintf('               The %dth coordination\n',coord_index)
%     end
%     fprintf('               The number of LJ potential: %d\n',numofLJ)
% 
%     disp('               Input parameter check completed.')
%     disp('*=============================================================*')
% end
% end

function coord_modify=coord_rearranger(coord,LJ,LJ_element)
numofcoord=length(coord);
coord_modify=cell(1,numofcoord);
for coord_index=1:numofcoord
    numofLJ=size(LJ,1);
%     numofatom=size(coord{coord_index}.data,1);
    
    coord_modify{coord_index}.data=[];
    coord_modify{coord_index}.textdata={};
    
    for LJ_index=1:numofLJ
        LJ_bool=strcmp(LJ_element{LJ_index},coord{coord_index}.textdata);
        coord_modify{coord_index}.data=[coord_modify{coord_index}.data;coord{coord_index}.data(LJ_bool,:)];
        coord_modify{coord_index}.textdata=[coord_modify{coord_index}.textdata;coord{coord_index}.textdata(LJ_bool)];
    end
end
end

function results_maker(coord,coord_solvent,solute_names,solvent_name,numofsolv,label_change)
numofcoord=length(coord);
[~,~]=mkdir('results');
for coord_index=1:numofcoord
    
    LJ_elements=[coord_solvent{1}.LJ_atom;coord{coord_index}.LJ_atom];
    numofLJ=length(LJ_elements);
    LJ_checker=zeros(1,numofLJ);
    
    temp_dir=sprintf('results/%s',solute_names{coord_index});
    [~,~]=mkdir(temp_dir);
    temp=sprintf('results/%s/%s.in',solute_names{coord_index},solute_names{coord_index});
    fp=fopen(temp,'w');
    
    fprintf(fp,"#MD input file for simulation of %s in %s\n",solute_names{coord_index},solvent_name{1});
    fprintf(fp,"#Structure\n");
    fprintf(fp,"%s %d\n",solvent_name{1},numofsolv);
    fprintf(fp,"#Pair_number\tX\tY\tZ\tMass(a.u)\tCharge(C)\tAtom_name\n");
    
    numofatom=length(coord_solvent{1}.textdata);
    for index_atom=1:numofatom
        bool_LJ=strcmp(coord_solvent{1}.textdata(index_atom),LJ_elements);
        num_LJ=find(bool_LJ);
        if isempty(num_LJ)
            error('*=============================================================*\nError !!\nThe name of LJ is not match.\nCheck the name of LJ and coordination files.\n*=============================================================*',0);
        end
        posi=find(strcmp(coord_solvent{1}.textdata{index_atom},label_change{1}), 1);
        if isempty(posi)
            text_temp=coord_solvent{1}.textdata(index_atom,1);
        else
            text_temp=label_change{2}{posi};
        end
        fprintf(fp,"%d\t%f\t%f\t%f",num_LJ,coord_solvent{1}.data(index_atom,1), ...
            coord_solvent{1}.data(index_atom,2), ...
            coord_solvent{1}.data(index_atom,3));
        if LJ_checker(num_LJ)==0
        fprintf(fp,"\t%f\t%f\t%s",coord_solvent{1}.weight(index_atom),coord_solvent{1}.charge(index_atom),text_temp);
        LJ_checker(num_LJ)=1;
        end
        fprintf(fp,'\n');
    end
    
    fprintf(fp,"%s %d\n",solute_names{coord_index},1);
    
    numofatom=length(coord{coord_index}.textdata);
    for index_atom=1:numofatom
        bool_LJ=strcmp(coord{coord_index}.textdata(index_atom),LJ_elements);
        num_LJ=find(bool_LJ);
        if isempty(num_LJ)
            error('*=============================================================*\nError !!\nThe name of LJ is not match.\nCheck the name of LJ and coordination files.\n*=============================================================*',0);
        end
        fprintf(fp,"%d\t%f\t%f\t%f",num_LJ,coord{coord_index}.data(index_atom,1), ...
            coord{coord_index}.data(index_atom,2), ...
            coord{coord_index}.data(index_atom,3));
        if LJ_checker(num_LJ)==0
            fprintf(fp,"\t%f\t%f\t%s",coord{coord_index}.weight(index_atom),coord{coord_index}.charge(index_atom),coord{coord_index}.textdata{index_atom});
            LJ_checker(num_LJ)=1;
        end
        fprintf(fp,'\n');
    end
    
%     numofatom=length(coord{coord_index}.textdata);
%     for index_atom=1:numofatom
%         bool_LJ=strcmp(coord{coord_index}.textdata(index_atom),LJ_elements);
%         num_LJ=find(bool_LJ);
%         fprintf(fp,"%d\t%f\t%f\t%f\t%f\t%f\t%s\n",num_LJ,coord{coord_index}.data(index_atom,1), ...
%             coord{coord_index}.data(index_atom,2), ...
%             coord{coord_index}.data(index_atom,3), ...
%             atomic_weight{bool_LJ},charge{bool_LJ},coord{coord_index}.textdata{index_atom});
%     end
    
    fprintf(fp,"end\n");
    fprintf(fp,"Lennard-Jones\n");
    
    pair=coord{coord_index}.pair;
    numofpair=size(pair,1);
    
    for pair_index=1:numofpair
        fprintf(fp,'%d \t %d \t %f \t %f \n',pair(pair_index,1),pair(pair_index,2),pair(pair_index,3),pair(pair_index,4));
    end
    
    fprintf(fp,"end\n");
    fclose(fp);
    fprintf('               The results file made.\n')
    disp('*=============================================================*')
end
end

% function LJ=LJ_setter(four_epsilon,sigma)
% % LJ=cell(1,numofcoord);
% % for coord_index=1:numofcoord
%     numofLJ=length(four_epsilon);
%     LJ=nan(numofLJ,2);
%     for LJ_index=1:numofLJ
%         LJ(LJ_index,:)=[four_epsilon{LJ_index};sigma{LJ_index}];
%     end
% %     LJ_temp=nan(numofatom,2);
% %     for atom_index=1:numofatom
% %         LJ_index=strcmp(coord{coord_index}.textdata(atom_index),LJ_elements);
% %         LJ_temp(atom_index,:)=[four_epsilon{LJ_index};sigma{LJ_index}];
% %     end
% %     LJ{coord_index}=LJ_temp;
% % end
% end

function coord=pair_generator_batch(coord,coord_solvent)
numofcoord=length(coord);
for coord_index=1:numofcoord
    LJ=[coord_solvent{1}.LJ;coord{coord_index}.LJ];
    pair=pair_generator(LJ);
    coord{coord_index}.pair=pair;
end
    function pair=pair_generator(LJ)
        numofatom=size(LJ,1);
        numofpair=sum([1:1:numofatom]);
        pair=nan(numofpair,4);
        pair_count=0;
        
        for ii=1:numofatom
            pair_count=pair_count+1;
            pair(pair_count,1)=ii;
            pair(pair_count,2)=ii;
            pair(pair_count,3)=LJ(ii,1);
            pair(pair_count,4)=LJ(ii,2);
        end
        
        for ii=1:numofatom-1
            for jj=ii+1:numofatom
                pair_count=pair_count+1;
                pair(pair_count,1)=ii;
                pair(pair_count,2)=jj;
                pair(pair_count,3)=sqrt(LJ(ii,1)*LJ(jj,1));
                pair(pair_count,4)=(LJ(ii,2)+LJ(jj,2))/2;
            end
        end
        fprintf('               The generation of LJ pair completed.\n')
        fprintf('               The number of pair: %d\n',size(pair,1))
        disp('*=============================================================*')
    end
end

function coord=read_coord(coord_list,constant)
numofcoord=length(coord_list);
coord=cell(1,numofcoord);
for ii=1:length(coord_list)
    coord{ii}=importdata(coord_list{ii});
    %     coord{ii}.num=sym2num(coord{ii}.textdata,constant);
    coord{ii}.num=sym2num_batch(coord{ii}.textdata,constant);
end
    function num=sym2num_batch(textdata,constant)
        numofatom=length(textdata);
        num=nan(numofatom,1);
        for atom_index=1:numofatom
            num(atom_index)=sym2num(textdata{atom_index},constant);
        end
        function num=sym2num(element,constant)
            element_cut=sym_cutter(element);
            num=find(strcmp(element_cut,constant.elements));
            
            function sym_cut=sym_cutter(sym)
                posi=strfind(sym,'_');
                if isempty(posi)
                    sym_cut=sym;
                else
                    sym_cut=sym(1:posi-1);
                end
            end
        end
    end
end
% Ref
% #Y. M. Munoz-Munoz, G. Guevara-Carrion, M. Llano-Restrepo and J. Vrabec, Fluid Phase Equilibria, 2015, 404, 150-160.