function Results = debye_solver_single(coord,q)
numofq=length(q);
[numofatom,~]=size(coord.data);
x=coord.data(:,1);
y=coord.data(:,2);
z=coord.data(:,3);%원자의 X,Y,Z정보
atom=coord.num(:);

% Termfrist=zeros(numofq,1); %form_factor file 요구됨
%     for n=1:numofq
%         for i=1:numofatom
%         Termfrist(n)=Termfrist(n)+form_factor(atom(i),q(n))^2;
%         end
%     end
    
Termsecond=zeros(numofq,1);
for n=1:numofq
    for i=1:numofatom-1
        for j=i+1:numofatom
            r=sqrt((x(i)-x(j))^2+(y(i)-y(j))^2+(z(i)-z(j))^2);
            Termsecond(n)=Termsecond(n)+form_factor(atom(i),q(n))*form_factor(atom(j),q(n))*sin(q(n)*r)/(q(n)*r);
        end
    end
end
    Results=2*Termsecond;
end