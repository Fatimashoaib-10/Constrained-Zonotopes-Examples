%% Mohammad Khajenejad
% Given f(.), cretaes a sign-stable g(x)=f(x)+Bx and then finds an interval enclosing the image of f(.) 

%clear; close all; clc;
function [Fmax,Fmin]=decomp_signstable(f,xup,xlow,dfup,dflow)
%clear; close all; clc;
n=size(xlow,1);m=size(dfup,1);%lambda=ppermute(n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
JBpp=zeros(size(dfup));Bp=zeros(size(dfup));
%% sign-stability conditions
for j=1:n
    for i=1:m
        Bpp(i,j)=-eps-dfup(i,j);%negative sign-stable
        Bp(i,j)=eps-dflow(i,j);%positive sign-stable
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% generating all the binary m by n matrices

MATRIX_SIZE_1 = m;MATRIX_SIZE_2 = n;

NumBits = MATRIX_SIZE_1*MATRIX_SIZE_2;
MaxNumericRepr = (2^NumBits)-1;

BinaryMatrices = cell(MaxNumericRepr + 1, 1);

for i = 0:MaxNumericRepr
    MatrixDigits = dec2bin(i, NumBits);
    MatrixDigits = reshape(MatrixDigits, [NumBits 1]);
    MatrixDigits = str2num(MatrixDigits);

    BinaryMatrices{i+1} = reshape(MatrixDigits, [MATRIX_SIZE_1 MATRIX_SIZE_2]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fmax=Inf(1,m);Fmin=-Inf(1,m);
%% enumerating all the possible Bs
for s=1:2^(m*n)
    B=BinaryMatrices{s}.*Bp+(ones(m,n)-BinaryMatrices{s}).*Bpp;
    Bpl=max(B,zeros(size(B)));Bne=Bpl-B;
    xcup=zeros(m,n);xclow=zeros(m,n);
    for i=1:m
        for j=1:n
            xcup(i,j)=0;xclow(i,j)=0;
            if BinaryMatrices{s}(i,j)==1
                xcup(i,j)=xup(j);
                xclow(i,j)=xlow(j);
            else
                xcup(i,j)=xlow(j);
                xclow(i,j)=xup(j);
            end
          %  if Bpl(i,j)<0
          %      Bpl(i,j)=0;
           % end
        end
    end
   % Bne=Bpl-B;
   fup=zeros(size(Fmax));flow=zeros(size(Fmin));
    for i=1:m
        fcup=f(xcup(i,:)');fclow=f(xclow(i,:)');
        fup(i)=fcup(i)+B(i,:)*xcup(i,:)'+Bne(i,:)*xup-Bpl(i,:)*xlow;
        flow(i)=fclow(i)+B(i,:)*xclow(i,:)'+Bne(i,:)*xlow-Bpl(i,:)*xup;
    end
    Fmax=min(Fmax,fup);
    Fmin=max(Fmin,flow);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for i=1:m
 %   Fmax(i)=Inf;Fmin(i)=-Inf;
  %  for k=1:size(lambda,1)
 
   % B(i,:)=lambda(k,:).*Bp(i,:)+(ones(1,n)-lambda(i,:)).*Bpp(i,:);
   % Bpl(i,:)=B(i,:);
   % for s=1:size(B(i,:),2)
    %    if Bpl(i,s)<0
     %       Bpl(i,s)=0;
      %  end
  %  end
   % Bne(i,:)=Bpl(i,:)-B(i,:);
  %  xcup=(lambda(i,:).*xup'+(ones(1,n)-lambda(i,:)).*xlow')';
   % xclow=(lambda(i,:).*xlow'+(ones(1,n)-lambda(i,:)).*xup')';
   % fxcup=f(xcup);fxclow=f(xclow);
   % fup=fxcup(i)+B(i,:)*xcup+Bne(i,:)*max(xcup,xlow)-Bpl(i,:)*min(xcup,xlow);
  %  flow=fxclow(i)+B(i,:)*xclow+Bne(i,:)*min(xcup,xlow)-Bpl(i,:)*max(xcup,xlow);
   % Fmax(i)=min(fup,Fmax(i));Fmin(i)=max(flow,Fmin(i));
    %,f(xclow));flow=min(f(xcup),f(xclow));
    %gup=fup(i)
   % end
%end
%end
%const=[const;dg+(Bp+Bpp)*(xup-xlow)<=theta*ones(m,1)];
%% corner points
%xc=xlow;xcc=xup;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% comparison with Liren
 %ours

end

