clear all
load constant
Kcal2KJ=4.184;% 1Kcal/mol = 4.184 KJ/mol
rm2sigma=1/(2^(1/6));% sigma=rm/(2^(1/6))
KbK2KJ=8.314472e-3;% Kb*K=R*K/Na=8.314472 J/mol;
% Moldy LJ input unit is 4 * KJ/mol and Angstrom (sigma)
% [Kcal2KJ,rm2sigma,KbK2KJ,LJ]=allocate;
%%% User options
numofsolv=1024;
solvent_name={
    'Cyclohexane'
    };
solvent_list={
    'cyclohexane_DFT.xyz'
    };

coord_list={
    'coord/Azo_t_S0.xyz'
    'coord/Azo_t_T1.xyz'
    'coord/Azo_t_S1.xyz'
    'coord/Azo_c_S0.xyz'
    'coord/Azo_c_T1.xyz'
    'coord/Azo_c_S1.xyz'
    
    'coord/DiBrAzo_t_S0.xyz'
    'coord/DiBrAzo_t_T1.xyz'
    'coord/DiBrAzo_t_S1.xyz'
    'coord/DiBrAzo_c_S0.xyz'
    'coord/DiBrAzo_c_T1.xyz'
    'coord/DiBrAzo_c_S1.xyz'
    
    };

solute_names={
    'Azo_c_S0'
    'Azo_c_S1'
    'Azo_c_T1'
    'Azo_t_S0'
    'Azo_t_S1'
    'Azo_t_T1'
    'DiBrAzo_c_S0'
    'DiBrAzo_c_S1'
    'DiBrAzo_c_T1'
    'DiBrAzo_t_S0'
    'DiBrAzo_t_S1'
    'DiBrAzo_t_T1'
    };

charge{1}=[
    -0.17718
    -0.17709
    0.0744
    0.07437
    -0.17109
    -0.17116
    -0.19307
    -0.19307
    -0.20385
    -0.20381
    -0.19595
    -0.19594
    -0.18308
    -0.18311
    0.21352
    0.21351
    0.22405
    0.22408
    0.2051
    0.20508
    0.20384
    0.20383
    0.20332
    0.20331
    ];
charge{2}=[
    -0.14273
    -0.14259
    0.07485
    0.0747
    -0.17875
    -0.1787
    -0.20856
    -0.2086
    -0.20155
    -0.20156
    -0.19405
    -0.19402
    -0.19135
    -0.19135
    0.2132
    0.21319
    0.21469
    0.21473
    0.20566
    0.20567
    0.2053
    0.20527
    0.20327
    0.20326
    ];
charge{3}=[
    -0.1818
    -0.18246
    0.0985
    0.09855
    -0.21522
    -0.21533
    -0.23557
    -0.23569
    -0.245
    -0.24495
    -0.23874
    -0.23869
    -0.22935
    -0.22952
    0.25409
    0.25402
    0.25773
    0.25753
    0.24624
    0.24616
    0.2459
    0.24586
    0.24389
    0.24386
    ];
charge{4}=[
    0.05511
    -0.19959
    -0.19209
    -0.19957
    -0.18945
    -0.20879
    -0.11182
    -0.11187
    0.05454
    -0.19811
    -0.19249
    -0.19879
    -0.18962
    -0.20841
    0.21459
    0.20549
    0.20448
    0.20557
    0.21506
    0.21479
    0.20558
    0.2045
    0.20554
    0.21536
    ];
charge{5}=[
    0.07438
    -0.179
    -0.20148
    -0.19159
    -0.19395
    -0.20879
    -0.14213
    -0.14243
    0.07462
    -0.17888
    -0.20146
    -0.19152
    -0.19398
    -0.20876
    0.21329
    0.20572
    0.20334
    0.20533
    0.21482
    0.21331
    0.20568
    0.20335
    0.20533
    0.21479
    ];
charge{6}=[
    0.09311
    -0.2276
    -0.23898
    -0.23935
    -0.23349
    -0.24682
    -0.144
    -0.14514
    0.09389
    -0.22916
    -0.23853
    -0.24087
    -0.23302
    -0.24885
    0.25409
    0.24641
    0.24386
    0.24655
    0.24904
    0.25374
    0.24622
    0.24369
    0.24645
    0.24876
    ];
charge{7}=[
    -0.17578
    -0.17575
    0.07014
    0.07016
    -0.15403
    -0.15413
    -0.17584
    -0.17582
    -0.2422
    -0.2422
    -0.23361
    -0.23364
    -0.059
    -0.05896
    0.21998
    0.21997
    0.23006
    0.23006
    0.22367
    0.22366
    0.22232
    0.22232
    0.07428
    0.07434
    ];
charge{8}=[
    -0.14158
    -0.14177
    0.0721
    0.07208
    -0.16224
    -0.16217
    -0.19304
    -0.19316
    -0.23831
    -0.23833
    -0.23087
    -0.23088
    -0.0697
    -0.06953
    0.21967
    0.21965
    0.22007
    0.22002
    0.22383
    0.22382
    0.22345
    0.22342
    0.0766
    0.07685
    ];
charge{9}=[
    -0.18465
    -0.18066
    0.09682
    0.09743
    -0.20206
    -0.20212
    -0.22173
    -0.22196
    -0.25437
    -0.25434
    -0.24829
    -0.24817
    -0.11514
    -0.11531
    0.26126
    0.26112
    0.26361
    0.26359
    0.26399
    0.26396
    0.26345
    0.26339
    0.07526
    0.07492
    ];
charge{10}=[
    0.04919
    -0.18162
    -0.23018
    -0.07271
    -0.22806
    -0.1925
    -0.10892
    -0.10936
    0.04843
    -0.17944
    -0.2307
    -0.07178
    -0.22845
    -0.19183
    0.22162
    0.2243
    0.22461
    0.22067
    0.22192
    0.22439
    0.22462
    0.22105
    0.07268
    0.07207
    ];
charge{11}=[
    0.07166
    -0.1625
    -0.23822
    -0.06953
    -0.23079
    -0.19339
    -0.14122
    -0.14131
    0.07193
    -0.16235
    -0.23829
    -0.06942
    -0.2309
    -0.19317
    0.21978
    0.22386
    0.22349
    0.22019
    0.21976
    0.22387
    0.2235
    0.22014
    0.07652
    0.0764
    ];
charge{12}=[
    0.09235
    -0.21466
    -0.24835
    -0.12482
    -0.2427
    -0.23496
    -0.14227
    -0.14241
    0.09357
    -0.21731
    -0.24771
    -0.12645
    -0.24191
    -0.23876
    0.26128
    0.26428
    0.26438
    0.2549
    0.26081
    0.26404
    0.2644
    0.25477
    0.07303
    0.07451
    ];

charge_solvent{1}=[
    0.108985
    0.108985
    0.108985
    0.108985
    0.108985
    0.108985
    -0.036082
    -0.036082
    -0.036082
    -0.036082
    -0.036082
    -0.036082
    -0.072903
    -0.072903
    -0.072903
    -0.072903
    -0.072903
    -0.072903
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

results_maker(coord,coord_solvent,solute_names,solvent_name,numofsolv);

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

function results_maker(coord,coord_solvent,solute_names,solvent_name,numofsolv)
numofcoord=length(coord);
for coord_index=1:numofcoord
    
    LJ_elements=[coord_solvent{1}.LJ_atom;coord{coord_index}.LJ_atom];
    numofLJ=length(LJ_elements);
    LJ_checker=zeros(1,numofLJ);
    [~,~]=mkdir('results');
    temp=sprintf('%s/%s.in','results',solute_names{coord_index});
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
        fprintf(fp,"%d\t%f\t%f\t%f",num_LJ,coord_solvent{1}.data(index_atom,1), ...
            coord_solvent{1}.data(index_atom,2), ...
            coord_solvent{1}.data(index_atom,3));
        if LJ_checker(num_LJ)==0
            fprintf(fp,"\t%f\t%f\t%s",coord_solvent{1}.weight(index_atom),coord_solvent{1}.charge(index_atom),coord_solvent{1}.textdata{index_atom});
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