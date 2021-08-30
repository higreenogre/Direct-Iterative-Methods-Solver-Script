%                            BY HIGREENOGRE
%                            FUNCTION  FILE
%%                     DIRECT ELIMINATION method for:
%%
clc; clear all;
%% 1.GAUSS ELIMINATION method for, eg.- 2x+3z=5; z=1; 3x+2z=5;

disp(' The solution by Gauss Elimination method for 2x+3z=5; z=1; 3x+2z=5 is:')
A=[ 2 0 3 ;0 0 1;  3 0 2 ];
B=[5;1;5];
solution1=gauss_el(A,B)
  
%% 2.GAUSS-JORDAN ELIMINATION method for,eg: 1x+2y+3z=6; 3x+6y+1z=10; x+3y+2z=6;

disp('The solution by the Gauss Jordan Elimination method for 1x+2y+3z=6; 3x+6y+1z=10; x+3y+2z=6 is:')
A=[1 2 3 ;3 6 1; 1 3 2 ];
B=[6;10;6];
solution2=gj(A,B)

%% 3.THE MATRIX INVERSE method for,eg: 1x+2y+3z=6; 3x+6y+1z=10; x+3y+2z=6;

disp('The solution by the matrix inverse method for 1x+2y+3z=6; 3x+6y+1z=10; x+3y+2z=6 is:')
A=[1 2 3 ;3 6 1; 1 3 2 ];
B=[6;10;6];
solution3=inversemethod(A,B)

%% 4. DOOLITTLE LU FACTORIZATION, eg: 2x-1y-2z=-1; -4x+6y+5z=7; -4x-2y+8z=2;

disp('The solution by doolittle LU factorization method for 2x-1y-2z=-1; -4x+6y+5z=7; -4x-2y+8z=2 is:')
A=[2 -1 -2 ;-4 6 5 ; -4 -2 8 ];
B=[-1;7;2];
UpperTriangle=lu(A,2)
LowerTriangle=lu(A,1)
 
UpperSolution=lu(A,2)*B
LowerSolution=lu(A,1)*B

%%                       ITERATIVE method for
%% 1.JACOBIAN ITERATION  eg: 5x-y+2z=12; 3x+8y-2z=-25; x+y+4z=6

disp('The solution by Jacobian Iteration method for 5x-y+2z=12; 3x+8y-2z=-25; x+y+4z=6 for 15 iterations is:')
A=[5 -1 2 ;3 8 -2; 1 1 4 ];
B=[12;-25;6];
max=15;
solution5=jacobiiteration(A,B,max)

%% 2.GAUSS SIEDAL ITERATION eg: 2x+2y+3z=7; 3x+6y+1z=10; x+3y+2z=6

disp('The solution by Gauss Siedel Iteration method for 2x+2y+3z=7; 3x+6y+1z=10; x+3y+2z=6 for 10 iterations is:')
A=[2 2 3 ;3 6 1; 1 3 2 ];
B=[7;10;6];
max=10
solution6=gs(A,B,max)

%% 3.SUCCESSIVE-OVER-RELAXATION eg: 2x+2y+3z=7; 3x+6y+1z=10; x+3y+2z=6

disp('The solution by Successive-over-relaxation method for 2x+2y+3z=7; 3x+6y+1z=10; x+3y+2z=6 for 100 iterations is:')
A=[2 2 3 ;3 6 1; 1 3 2 ];
B=[7;10;6];
weight=1.1;
max=100;
solution7=sor(A,B,weight,max)

%%                         FUNCTIONS
%GAUSS ELIMINATION method for

function final=gauss_el(A,B)
    trimmed_mat=trim1(A,B);
    C=[A,B];
    dim=size(C);
    dimofmat=size(trimmed_mat{1});
    if dimofmat==0
        final=zeros(dim(1),dim(2));
        return
    end
    final=convert1(A,B,trimmed_mat{1},trimmed_mat{2},trimmed_mat{3});
end

function final=trim1(A,B)
    C=[A,B];
    dim=size(C);
    removecol=[];
    for i=2:dim(1)
        C(1,:)=C(1,:)+C(i,:);
    end
    for i=1:dim(1)
        n=i;
        while C(1,i)==0 && n<=dim(1)
            if n==i
                n=n+1;
                continue
            end
            C(i,:)=C(i,:)+C(n,:);
            n=n+1;
        end
        if n==dim(1)+1
            removecol=[removecol,i];
        end
    end
    C; % by the end of this statement the inactive variables are known
    removecol ;% the columns that dont contain variables
    n=dim(1);
    removerow=[];
    for i=removecol
        C(:,i-dim(1)+n)=[];
        C(n,:)=[];
        removerow=[removerow,n];
        n=n-1;
    end
    removerow;
    mat=C;
    dimofmat=size(mat);
    if dimofmat(1) == 0
        final={mat,removerow,removecol};
        return %if there are no variables then a zero matriz will be returned and the program ends here
    end
    mat(1,:)=mat(1,:)/mat(1,1);
    for i=2:dimofmat(1)
        for j=i:dimofmat(1)
            mat(j,:)=mat(j,:)-mat(j,i-1)*mat(i-1,:);
        end
        n=i+1;
        while mat(i,i)==0 & n<=dimofmat
            mat(i,:)=mat(i,:)+mat(n,:);
            n=n+1;
        end
        mat(i,:)=mat(i,:)/mat(i,i);
    end
    final={mat,removerow,removecol};
end

function final=convert1(A,B,mat,removerow,removecol)
    C=[A,B];
    dim=size(C);
    final=zeros(dim(1),dim(2));
    globalrow=1:dim(1);
    globalcol=1:dim(2);
    for i=removerow
        globalrow(i)=[];
    end
    count=0;
    for j=removecol
        globalcol(j-count)=[];
        count=count+1;
    end
    for i=1:length(globalrow)
        for j=1:length(globalcol)
            final(globalrow(i),globalcol(j))=mat(i,j);
        end
    end
end

% GAUSS JORDAN  ELIMINATION method for

function final=gj(A,B)
    trimmed_mat=trim2(A,B);
    C=[A,B];
    dim=size(C);
    dimofmat=size(trimmed_mat{1});
    if dimofmat==0
        final=zeros(dim(1),dim(2));
        return
    end
    final=convert2(A,B,trimmed_mat{1},trimmed_mat{2},trimmed_mat{3});
end

function final=trim2(A,B)
    C=[A,B];
    dim=size(C);
    removecol=[];
    for i=1:dim(1)
        n=i;
        while C(1,i)==0 && n<=dim(1)
            if n==i
                n=n+1;
                continue
            end
            C(i,:)=C(i,:)+C(n,:);
            n=n+1;
        end
        if n==dim(1)+1
            removecol=[removecol,i];
        end
    end
    C; % by the end of this statement the inactive variables are known
    removecol ;% the columns that dont contain variables
    n=dim(1);
    removerow=[];
    for i=removecol
        C(:,i-dim(1)+n)=[];
        C(n,:)=[];
        removerow=[removerow,n];
        n=n-1;
    end
    removerow;
    mat=C;
    dimofmat=size(mat);
    if dimofmat(1) == 0
        final={mat,removerow,removecol};
        return %if there are no variables then a zero matriz will be returned and the program ends here
    end
    mat(1,:)=mat(1,:)/mat(1,1);
    for i=2:dimofmat(1)
        for j=i:dimofmat(1)
            mat(j,:)=mat(j,:)-mat(j,i-1)*mat(i-1,:);
        end
        n=i+1;
        while mat(i,i)==0 & n<=dimofmat
            mat(i,:)=mat(i,:)+mat(n,:);
            n=n+1;
        end
        mat(i,:)=mat(i,:)/mat(i,i);
    end
    for i=dimofmat(1)-1:-1:1
        for j=i:-1:1
            mat(j,:)=mat(j,:)-mat(i+1,:)*mat(j,i+1);
        end
    end
    
    final={mat,removerow,removecol};
end

function final=convert2(A,B,mat,removerow,removecol)
    C=[A,B];
    dim=size(C);
    final=zeros(dim(1),dim(2));
    globalrow=1:dim(1);
    globalcol=1:dim(2);
    for i=removerow
        globalrow(i)=[];
    end
    count=0;
    for j=removecol
        globalcol(j-count)=[];
        count=count+1;
    end
    for i=1:length(globalrow)
        for j=1:length(globalcol)
            final(globalrow(i),globalcol(j))=mat(i,j);
        end
    end
end


% THE MATRIX INVERSE method for

function final=inversemethod(A,B)
    final=Inversematrix(A)*B;
end

function final=Inversematrix(mat)
    if determinant(mat)~=0
       dim=size(mat);
        C=eye(dim);
        d=mat;
        mat(1,:)=mat(1,:)/mat(1,1);
        C(1,:)=C(1,:)/d(1,1);
        for i=2:dim
        for j=i:dim
            d=mat;
            mat(j,:)=mat(j,:)-mat(j,i-1)*mat(i-1,:);
            C(j,:)=C(j,:)-d(j,i-1)*C(i-1,:);
        end
        n=i+1;
        while mat(i,i)==0 & n<=dim
            mat(i,:)=mat(i,:)+mat(n,:);
            C(i,:)=C(i,:)+C(n,:);
            n=n+1;
        end
        d=mat;
        mat(i,:)=mat(i,:)/mat(i,i);
        C(i,:)=C(i,:)/d(i,i);
    end
    for i=dim-1:-1:1
        for j=i:-1:1
            d=mat;
            mat(j,:)=mat(j,:)-mat(i+1,:)*mat(j,i+1);
            C(j,:)=C(j,:)-C(i+1,:)*d(j,i+1);
        end
    end
    final=C;
    end
end

% THE DETERMINANT

function final=determinant(A)
    final=0;
    if size(A,1)==size(A,2) && size(A,1)==1
        final=A;
        return
    else 
        for i=1:size(A,1)
            removed=A;
            removed(1,:)=[];
            removed(:,i)=[];
            final=final-(-1)^(i)*A(1,i)*determinant(removed);
        end
    end
end

% DOOLITTLE LU FACTORIZATION

function final = lu(mat,n)
    dim = size(mat,1);
    lower = zeros(dim);
    upper = zeros(dim);
    for i = 1:dim
        % to build the upper triangle matrix
        for k = i:dim
            sum = 0;
            for j = 1:i
                sum = sum+ (lower(i,j)*upper(j,k));
            end
            upper(i,k) = mat(i,k) - sum;
        end
        
        % to build the lower triangle matrix
        for k = i:dim
            if i==k
                lower(i,k) = 1;
            else
                sum = 0;
                for j = 1:i
                    sum = sum+ (lower(k,j)*upper(j,i));
                end
                lower(k,i) = ( mat(k,i) - sum) / upper(i,i);
            end
        end
    end
    if n==1
        final = lower;
    else
        final = upper;
    end
end

% JACOBI ITERATION

function final=jacobiiteration(A,B,times)
    dim=size(A,1);
    x=ones(dim,1);
    y=zeros(dim,1);
    for  i  = 1:times
        for j=1:dim
            remove=1:dim;
            remove(j)=[];
            value=0;
            for k=remove
                value=value+A(j,k)*x(k);
            end
            y(j)=(B(j)-value)/A(j,j);
        end
        x=y;
    end
    final=x;
end

% GAUSS SIEDEL ITERATION

function final=gs(A,B,times)
    dim=size(A,1);
    x=zeros(dim,1);
    for  i  = 1:times
        for j=1:dim
            remove=1:dim;
            remove(j)=[];
            value=0;
            for k=remove
                value=value+A(j,k)*x(k);
            end
            x(j)=(B(j)-value)/A(j,j);
        end
    end
    final=x;
end

% SUCCESSIVE OVER RELAXATION

function final=sor(A,B,w,times)
    dim=size(A,1);
    x=zeros(dim,1);
    for  i  = 1:times
        for j=1:dim
            remove=1:dim;
            remove(j)=[];
            value=0;
            for k=remove
                value=value+A(j,k)*x(k);
            end
            x(j)=w*(B(j)-value)/A(j,j)+(1-w)*x(j);
        end
    end
    final=x;
end