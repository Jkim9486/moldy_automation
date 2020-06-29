function num=sym2num(text)
global constant
posi=strfind(text,'_');
if 1==isempty(posi)
    textdata=text;
else
    textdata=text(1:posi-1);
end
for i=1:104
    if strcmp(textdata,constant.elements{i})==1
        num=i;
    else
    end
end
end