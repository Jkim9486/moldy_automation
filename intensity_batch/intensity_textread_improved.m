%   program intensity.c
%
%       version 0.1 28/02/2004 01h45
%
%       written by Rodolphe Vuilleumier
%       email: vuilleum@lptl.jussieu.fr
%
%
%   calculates x-ray elastic scattering intensities and partials from
%   radial distribution functions
%
%
%   the atomic form factors are defined in file atoms.h generated from
%   /users/opid09/lib/sfac.dat
%   see header of atoms.h for more information
%
%   NOTE: a fake atom M was added in the list of atom types to be used as
%         a ghost atom in moldy

% intensity.c --> intensity.cpp Rewriten by Dr. Hosung Ki, Nano-Bio Structural dynamics Lab.
% intensity.cpp --> intensity.m Rewriten by Jungmin Kim, Nano-Bio Structural dynamics Lab.
% Email : jkim9486@gmail.com
clear all
global atom
global constant
% Get constant
% temp2=sprintf('%s/constant.mat',temp);
% load(temp2);
load('constant.mat')

plot_option=0;

% User Input Arguments
% input_filename='rdf_results/Average_Os3CO12.rdf';

q=[0.005:0.01:10]'; %#ok<NBRAK>
Box_length=35.9113; %% Check in moldy '~.out' file.
solvent={'C_100','H_100','H_200'}; %% For off cage calculate;
numofatom_solv=[256*6,256*6,256*6]; 
coord_list={
    'MD_input/coord/Azo_t_S0.xyz'
    'MD_input/coord/Azo_t_T1.xyz'
    'MD_input/coord/Azo_t_S1.xyz'
    'MD_input/coord/Azo_c_S0.xyz'
    'MD_input/coord/Azo_c_T1.xyz'
    'MD_input/coord/Azo_c_S1.xyz'
    
    'MD_input/coord/DiBrAzo_t_S0.xyz'
    'MD_input/coord/DiBrAzo_t_T1.xyz'
    'MD_input/coord/DiBrAzo_t_S1.xyz'
    'MD_input/coord/DiBrAzo_c_S0.xyz'
    'MD_input/coord/DiBrAzo_c_T1.xyz'
    'MD_input/coord/DiBrAzo_c_S1.xyz'
    };

input_filenames={
    'input/Azo_c_S0/rdf_results/Average_Azo_c_S0.rdf'
    'input/Azo_c_S1/rdf_results/Average_Azo_c_S1.rdf'
    'input/Azo_c_T1/rdf_results/Average_Azo_c_T1.rdf'
    'input/Azo_t_S0/rdf_results/Average_Azo_t_S0.rdf'
    'input/Azo_t_S1/rdf_results/Average_Azo_t_S1.rdf'
    'input/Azo_t_T1/rdf_results/Average_Azo_t_T1.rdf'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    'input/DiBrAzo_c_S0/rdf_results/Average_DiBrAzo_c_S0.rdf'
    'input/DiBrAzo_c_S1/rdf_results/Average_DiBrAzo_c_S1.rdf'
    'input/DiBrAzo_c_T1/rdf_results/Average_DiBrAzo_c_T1.rdf'
    'input/DiBrAzo_t_S0/rdf_results/Average_DiBrAzo_t_S0.rdf'
    'input/DiBrAzo_t_S1/rdf_results/Average_DiBrAzo_t_S1.rdf'
    'input/DiBrAzo_t_T1/rdf_results/Average_DiBrAzo_t_T1.rdf'
    };
% input_filename='Azo_c_S0rdf_results/Average_Os3CO12.rdf';

% atoms{1}=[
%     'Os_1',1;
%     ];

% atom=({
%     'Os_1',1;
%     'Os_2',1;
%     'Os_3',1;
%     %%%%%%%%%
%     'C_1',1;
%     'C_2',1;
%     'C_3',1;
%     'C_4',1;
%     'C_5',1;
%     'C_6',1;
%     'C_7',1;
%     'C_8',1;
%     'C_9',1;
%     'C_10',1;
%     'C_11',1;
%     'C_12',1;
%     %%%%%%%%%
%     'O_1',1;
%     'O_2',1;
%     'O_3',1;
%     'O_4',1;
%     'O_5',1;
%     'O_6',1;
%     'O_7',1;
%     'O_8',1;
%     'O_9',1;
%     'O_10',1;
%     'O_11',1;
%     'O_12',1;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     'C_100',1024*6;
%     'H_100',1024*6;
%     'H_200',1024*6;
%     });

atoms=atoms_maker(coord_list,solvent,constant,numofatom_solv);
numofcoord=length(coord_list);

for coord_index=1:numofcoord
    atom=atoms{coord_index};
    input_filename=input_filenames{coord_index};
    intensity_main(input_filename,q,Box_length,solvent,plot_option);
    fprintf('              Run %d.\n',coord_index);
end
fprintf('              End.\n')

%% Data treat_arguments
function intensity_main(input_filename,q,Box_length,solvent,plot_option)
global atom
global constant

cage_option=isempty(solvent);
numofq=length(q);
[dir_origianl,name,~]=fileparts(input_filename);
path_current=pwd;
temp=sprintf('%s/codes',path_current);
addpath(temp)
V=Box_length^3;
%% File import-importdata
data_temp=importdata(input_filename);

r=data_temp.data(:,1); %% r
gr=data_temp.data(:,2:end); %% Hg-I_1
[numofr,numofpair]=size(gr);
dr=diff(r); %% Get dr
%% File import-textscan
% fid=fopen(input_filename,'r');
data_temp=extractAfter(string(cell2mat(data_temp.textdata)), 'r');
pair_cell=cell(1,numofpair); %% (1,:)=pair, (2,:) First atom symbol (3,:) Second atom symbol (4,:) First atom (5,:) Second atom (6,:) First atom number (7,:) Second atom number
pairs = regexp(data_temp,'[A-Z][a-z_0-9]*-[A-Z][a-z_0-9]*','match');
for ii=1:length(pairs)
    pair_cell{ii}=char(pairs(ii));
end
[pair,Z,N]=pair_analyzer(pair_cell);
%% File import-textscan
% fid=fopen(input_filename,'r');
% data_temp=textscan(fid,'%s');
% 
% pair_cell=cell(1,numofpair); %% (1,:)=pair, (2,:) First atom symbol (3,:) Second atom symbol (4,:) First atom (5,:) Second atom (6,:) First atom number (7,:) Second atom number
% for ii=1:numofpair
%     pair_cell{ii}=data_temp{1}{2+ii,1};
% end
% [pair,Z,N]=pair_analyzer(pair_cell);
%% Formfactor
f=nan(2,numofpair,numofq);
for ii=1:numofpair
    for jj=1:2
        for kk=1:numofq
            f(jj,ii,kk)=form(Z(jj,ii),q(kk));
        end
    end
end
%% RDF 2 S(q)
signal=zeros(numofq,numofpair);
for ii=1:numofq
    for jj=1:numofpair
        fact=2;
        if 1==strcmp(pair{1,jj},pair{2,jj})
            signal(ii,jj)=signal(ii,jj)+f(1,jj,ii)^2*N(1,jj);
            fact=1;
        else
        end
        if gr(end,jj)>0.25
            substration=1;
        else
            substration=0;
        end
        for kk=1:numofr-1
            signal(ii,jj)=signal(ii,jj)+N(1,jj)*N(2,jj)/V*f(1,jj,ii)*f(2,jj,ii)*4*pi*r(kk)*sin(q(ii)*r(kk))/q(ii)*(gr(kk,jj)-substration)*fact*dr(kk);
        end
    end
end
results=[q,sum(signal,2),signal];
%% Solvent checker
if 1~=cage_option
    count=0;
    count2=0;
    for ii=1:numofpair
        checker=0;
        count=count+1;
        for jj=1:2
            for kk=1:length(solvent)
            if 1==strcmp(pair{jj,ii},solvent{kk})
                checker=checker+1;
            end
            end
        end
        if checker==1
            count2=count2+1;
        count_array(count2)=count; %#ok<SAGROW>
        end
    end
    pair_cell_cage=cell(1,length(count_array));
    signal_cage=nan(numofq,length(count_array));
    for ii=1:length(count_array)
        pair_cell_cage{ii}=pair_cell{count_array(ii)};
        signal_cage(:,ii)=signal(:,count_array(ii));
    end
    results_cage_element=[q,signal_cage];
    results_cage=[q,sum(signal_cage,2)];
end
%% Making Output file
dir=sprintf('%s/%s_results_intensity',dir_origianl,name);
[~,~]=mkdir(dir);

% temp=cell2mat(pair_cell);
filename=sprintf('%s/Sq_%s.dat',dir,name);
fid=fopen(filename,'w');
fprintf(fid, '%-6s','q');
fprintf(fid,' \t%-12s','sum');
fprintf(fid,' \t%-12s',pair_cell{:});
      for ii=1:numofq
             fprintf(fid,'\n');
             fprintf(fid,'%-4g ',q(ii));
             fprintf(fid,'\t%-12g',results(ii,2:end));
      end
fclose(fid);

if 1~=cage_option
    filename=sprintf('%s/Sq_%s_cage_elements.dat',dir,name);
    fid=fopen(filename,'w');
    fprintf(fid, '%-6s','q');
    fprintf(fid,' \t%-12s',pair_cell_cage{:});
          for ii=1:numofq
                 fprintf(fid,'\n');
                 fprintf(fid,'%-4g ',q(ii));
                 fprintf(fid,'\t%-12g',results_cage_element(ii,2:end));
          end
    fclose(fid);

    filename=sprintf('%s/Sq_%s_cage.dat',dir,name);
    fid=fopen(filename,'w');
    fprintf(fid, '%-6s','q');
    fprintf(fid,' \t%-12s','cage');
          for ii=1:numofq
                 fprintf(fid,'\n');
                 fprintf(fid,'%-4g ',q(ii));
                 fprintf(fid,'\t%-12g',results_cage(ii,2:end));
          end
    fclose(fid);
end
%% Ploting
if plot_option
    figure
    plot(q,results(:,2:end))
    legend(['Static',pair_cell])
    figure_title=sprintf('Sq_%s',name);
    set(gcf,'numbertitle','off','name', figure_title);
    if 1~=cage_option
        
        figure
        plot(q,results_cage(:,2:end))
        figure_title=sprintf('Sq_%s_cage',name);
        set(gcf,'numbertitle','off','name', figure_title);
        
        figure
        plot(q,results_cage_element(:,2:end))
        legend(pair_cell_cage)
        figure_title=sprintf('Sq_%s_cage_element',name);
        set(gcf,'numbertitle','off','name', figure_title);
        
    end
end
end
function atoms=atoms_maker(coord_list,solvent,constant,numofatom_solv)
coord=read_coord(coord_list,constant);
numofcoord=length(coord_list);
numofkind_solv=length(solvent);
atoms=cell(1,numofcoord);
for coord_index=1:numofcoord
    numofatom=length(coord{coord_index}.textdata);
    atom_temp=cell(numofatom+numofkind_solv,2);
    for atom_index=1:numofatom
        atom_temp{atom_index,1}=coord{coord_index}.textdata{atom_index};
        atom_temp{atom_index,2}=1;
    end
    for atom_index=numofatom+1:numofatom+numofkind_solv
        atom_temp{atom_index,1}=solvent{atom_index-numofatom};
        atom_temp{atom_index,2}=numofatom_solv(atom_index-numofatom);
    end
    atoms{coord_index}=atom_temp;
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
end
