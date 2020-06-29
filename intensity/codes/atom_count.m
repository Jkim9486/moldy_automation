function number=atom_count(text)
global atom
[numofatom,~]=size(atom);
for ii=1:numofatom
    if 1==strcmp(atom{ii,1},text)
        number=atom{ii,2};
    else
    end
end
end