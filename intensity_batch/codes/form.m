function[results]=form(atomic_number,q)
global constant
s=q/(4*pi);
%% Simple calculate for An,Bn,Cn term (D. Waasmaier and A. Kirfel)
temp=0;
for k=1:5
        temp=temp+constant.Xfactor(atomic_number,k)*exp(-constant.Xfactor(atomic_number,k+6)*s*s);
end
results=temp+constant.Xfactor(atomic_number,6);