function [pair,Z,N]=pair_analyzer(pair_cell)
    numofpair=length(pair_cell);
    pair=cell(2,numofpair);
    Z=nan(2,numofpair);
    N=nan(2,numofpair);
    for ii=1:numofpair
        posi1=strfind(pair_cell{1,ii},'-');
        pair{1,ii}=pair_cell{1,ii}(1:posi1-1);
        pair{2,ii}=pair_cell{1,ii}(posi1+1:end);
         for jj=1:2
            Z(jj,ii)=sym2num(pair{jj,ii});
            N(jj,ii)=atom_count(pair{jj,ii});
         end
    end
end