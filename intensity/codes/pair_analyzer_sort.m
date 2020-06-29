function [pair_results,pair_atom_results]=pair_analyzer_sort(pair_cell)
    numofpair=length(pair_cell);
    pair_results=cell(3,numofpair);
    pair_atom_results=nan(2,numofpair);
    for ii=1:numofpair
        posi1=strfind(pair_cell{1,ii},'-');
        pair_results{1,ii}=pair_cell{1,ii}(1:posi1-1);
        pair_results{2,ii}=pair_cell{1,ii}(posi1+1:end);
        temp_string=string({pair_cell{1,ii}(1:posi1-1),pair_cell{1,ii}(posi1+1:end)});
        sort_string=sort(temp_string);
        pair_results{1,ii}=pair_cell{1,ii};
         for jj=1:2
            pair_results{1+jj,ii}=convertStringsToChars(sort_string(jj));
            pair_atom_results(jj,ii)=sym2num(convertStringsToChars(sort_string(jj)));
         end
    end
end