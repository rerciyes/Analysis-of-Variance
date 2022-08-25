soil_a=transpose([8 7 7 6 7 13 15 14 15 14]);%
soil_b=transpose([5 6 6 4 5 12 14 13 14 13]);
%toxin_3=transpose([18 21 20 22 24 0]);
control=transpose([11 14 11 16 0 0]);
data_set=[soil_a soil_b]
ANOVA1_partition_TSS(data_set)
ANOVA1_test_equality(data_set,0.1)
%c=[1 -4 3];
%ANOVA1_is_contrast(c)
n=[4 6 5 4];
c1=[1 -4 3 3];c2=[1 -4 3 0];
ANOVA1_is_orthagonal(n,c1,c2)
%data set for two way anova
X=[];
X(:,1,:)=[8 7 7 6 7;13 15 14 15 14];
X(:,2,:)=[5 6 6 4 5;12 14 13 14 13];
%X(:,3,:)=[10 12 11 19;12 13 10 13];
ANOVA2_partition_TSS(X)
ANOVA2_test_equality(X,0.1,"AB")
C=[1 0 0 -1;-1 0 0 1];
d=[4 5 2 3];

%ANOVA1_CI_linear_combs(data_set,C,0.1,"Best")
ANOVA1_test_linear_combs(data_set,C,0.1,d,"Best")
ANOVA2_MLE(X)

%% 
% 

function [values] = ANOVA1_partition_TSS(data_set)
    %overall_sample_mean
    x_bar=0;
    counter=0;
    for i=1:size(data_set,1)
        for j=1:size(data_set,2)
          x_bar=x_bar+data_set(i,j);
          if data_set(i,j)~=0
              counter=counter+1;
          end
        end
    end
    overall_mean=x_bar/counter;
    
    %sample_mean_of_ith_sample
    sample_mean_of_ith=[];
    column_sizes=[];
    for i=1:size(data_set,2)
        x_ith_bar=0;
        counter=0;
        for j=1:size(data_set,1)
          x_ith_bar=x_ith_bar+data_set(j,i);
          if data_set(j,i)~=0
              counter=counter+1;
          end
        end
        sample_mean_of_ith(i)=x_ith_bar/counter;
        column_sizes(i)=counter;
    end
    
    %ss_total
    ss_total=0;
    for i=1:size(data_set,1)
        for j=1:size(data_set,2)
          if i-1<column_sizes(j)
             ss_total=ss_total+((data_set(i,j)- overall_mean)^2);
          end
        end
    end
    
    %ss_w
    ss_w=0;
    for i=1:size(data_set,1)
        for j=1:size(data_set,2)
          if i-1<column_sizes(j)
             ss_w=ss_w+((data_set(i,j)- sample_mean_of_ith(j))^2);
          end
        end
    end
    
    %ss_b
    ss_b=0;
    for j=1:size(data_set,2)
       ss_b=ss_b+(column_sizes(j)*(sample_mean_of_ith(j)-overall_mean)^2);
    end
    values=[ss_total ss_w ss_b];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [anova_table] = ANOVA1_test_equality(data_set,alpha)
x=ANOVA1_partition_TSS(data_set);
%finding ni's
n=0;
for i=1:size(data_set,1)
    for j=1:size(data_set,2)
        if(data_set(i,j)~=0)
            n=n+1;
        end
    end
end

ss_total=x(1);ss_w=x(2);ss_b=x(3);
df_of_b=size(data_set,2)-1;
df_of_w=n-size(data_set,2); 
df_total=n-1;
ms_b=ss_b/df_of_b;
ms_w=ss_w/df_of_w;
F=ms_b/ms_w;
f=finv(1-alpha,df_of_b,df_of_w);
decision=1;
if F>f
    decision=0;
end
pvalue = 1-fcdf(F,df_of_b,df_of_w);
anova_table=[df_of_b ss_b ms_b F;df_of_w ss_w ms_w NaN;df_total ss_total pvalue NaN];

if decision==0
    disp("Hypothesis Rejected")
else
    disp("Hypothesis not Rejected")
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function new_alpha=Bonferroni_correction(alpha,m)
new_alpha=alpha/m;
end

function new_alpha=Sidak_correction(alpha,m)
new_alpha=1-((1-alpha)^(1/m));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bool = ANOVA1_is_contrast(c)
total_sum=sum(c);
if total_sum==0
    disp("It's contrast")
    bool=true;
else
    disp("Not contrast")
    bool=false;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bool=ANOVA1_is_orthagonal(n,c1,c2)
%should return a warning if any of the linear combinations is not a contrast.
%n c1 and c2 must be same size
bool1=ANOVA1_is_contrast(c1);bool2=ANOVA1_is_contrast(c2);
if(bool1&&bool2)
    total_sum=c1*transpose(c2)
    if total_sum==0
        disp("It's orthagonal")
        bool=true;
    else
        disp("Not orthagonal")
        bool=false;
    end
else
    bool=false;
    disp("One of c is not contrast")
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [values]=ANOVA2_partition_TSS(X)

x_bar_ij=0;
for k=1:size(X,3)
    x_bar_ij=x_bar_ij+X(:,:,k);
end
x_bar_ij=x_bar_ij/size(X,3);

x_bar_i=0;
for j=1:size(X,2)
    for k=1:size(X,3)
    x_bar_i=x_bar_i+X(:,j,k);
    end
end
x_bar_i=x_bar_i/(size(X,3)*size(X,2));

x_bar_j=0;
for i=1:size(X,1)
    for k=1:size(X,3)
    x_bar_j=x_bar_j+X(i,:,k);
    end
end
x_bar_j=x_bar_j/(size(X,3)*size(X,1));

x_bar=0;
for i=1:size(X,1)
    for j=1:size(X,2)
        for k=1:size(X,3)
        x_bar=x_bar+X(i,j,k);
        end
    end
end
x_bar=x_bar/(size(X,3)*size(X,1)*size(X,2));

ssa=0;
for i=1:size(X,1)
    ssa=ssa+(x_bar_i(i)-x_bar)^2;
end
ssa=size(X,2)*size(X,3)*ssa;

ssb=0;
for j=1:size(X,2)
    ssb=ssb+(x_bar_j(j)-x_bar)^2;
end
ssb=size(X,1)*size(X,3)*ssb;

ssab=0;
for i=1:size(X,1)
    for j=1:size(X,2)
    ssab=ssab+(x_bar_ij(i,j)-x_bar_j(j)-x_bar_i(i)+x_bar)^2;
    end
end
ssab=ssab*size(X,3);

sse=0;
for i=1:size(X,1)
    for j=1:size(X,2)
        for k=1:size(X,3)
        sse=sse+(X(i,j,k)-x_bar_ij(i,j))^2;
        end
    end
end

sstotal=0;
for i=1:size(X,1)
    for j=1:size(X,2)
        for k=1:size(X,3)
        sstotal=sstotal+(X(i,j,k)-x_bar)^2;
        end
    end
end
values=[sstotal ssa ssb ssab sse];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [twoway_anovatable_desired_part]=ANOVA2_test_equality(X,alpha,choice)
values=ANOVA2_partition_TSS(X);
%prompt = "Choice for hypothesis: ";
%choice = input(prompt)
msa=values(2)/(size(X,1)-1);
msb=values(3)/(size(X,2)-1);
msab=values(4)/((size(X,1)-1)*(size(X,2)-1));
mse=values(5)/((size(X,1))*(size(X,2))*(size(X,3)-1));
f_1=msa/mse;
f_2=msb/mse;
f_3=msab/mse;
values=[size(X,1)-1 values(2) msa f_1;size(X,2)-1 values(3) msb f_2;(size(X,1)-1)*(size(X,2)-1) values(4) msab f_3;(size(X,1))*(size(X,2))*(size(X,3)-1) values(5) mse NaN;(size(X,1))*(size(X,2))*(size(X,3))-1 values(1) NaN NaN];
if(choice=="A")
    twoway_anovatable_desired_part=values(1,:);
end
if(choice=="B")
    twoway_anovatable_desired_part=values(2,:);
end
if(choice=="AB")
    twoway_anovatable_desired_part=values(3,:);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ANOVA2_MLE(X)
x_bar_ij=0;
for k=1:size(X,3)
    x_bar_ij=x_bar_ij+X(:,:,k);
end
x_bar_ij=x_bar_ij/size(X,3);

x_bar_i=0;
for j=1:size(X,2)
    for k=1:size(X,3)
    x_bar_i=x_bar_i+X(:,j,k);
    end
end
x_bar_i=x_bar_i/(size(X,3)*size(X,2));

x_bar_j=0;
for i=1:size(X,1)
    for k=1:size(X,3)
    x_bar_j=x_bar_j+X(i,:,k);
    end
end
x_bar_j=x_bar_j/(size(X,3)*size(X,1));

x_bar=0;
for i=1:size(X,1)
    for j=1:size(X,2)
        for k=1:size(X,3)
        x_bar=x_bar+X(i,j,k);
        end
    end
end
x_bar=x_bar/(size(X,3)*size(X,1)*size(X,2));

mean=x_bar
a_bar=x_bar_i-x_bar
b_bar=x_bar_j-x_bar
sigma=x_bar_ij-x_bar_i-x_bar_j+x_bar
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [bounds]=ANOVA1_CI_linear_combs(dataset,C,alpha,choice)
X=transpose(dataset);
onewayanova=ANOVA1_partition_TSS(dataset);
%deriving the contrasts
contrastbool=1;
for i=1:size(C,1)
    bool1=ANOVA1_is_contrast(C(i,:));
    if(~bool1)
        contrastbool=0;
    end
end

%finding ni's
for i=1:size(X,1)
    n(i)=0;
    for j=1:size(X,2)
        if(X(i,j)~=0)
            n(i)=n(i)+1;
        end
    end
end

%deriving the orthagonality
orthagonalbool=1;
for i=1:size(C,1)
    for j=1:size(C,2)
        bool2=ANOVA1_is_orthagonal(n,C(i,:),C(j,:));
        if(~bool2)
            orthagonalbool=0;
            break
        end
    end
    if(~bool2)
        break
    end
end 

%pairwise difference check
pairwisebool=0;
for i=1:size(C,1)
    nzeros=numel(C(i,:))-nnz(C(i,:));
    x=sum(C(i,:));
    if(x==0&&nzeros==size(C,2)-2)
        pairwisebool=1;
    end
end
for i=1:size(n,2)-1
    for j=i+1:size(n,2)
        if(n(i)~=n(j))
            pairwisebool=0;
            break
        end
    end
    if(pairwisebool==0)
        break
    end
end
matrix_for_the_best=zeros(2,8);
%Scheffe
if(contrastbool)
    fcontrast=finv(1-alpha,size(C,2)-1,sum(n)-size(C,2));
    sc1c=sqrt(onewayanova(2)*sum((C.^2)/n)/(sum(n)-size(C,2)))*sqrt(size(C,2)*fcontrast);
    sc2c=-sqrt(onewayanova(2)*sum((C.^2)/n)/(sum(n)-size(C,2)))*sqrt(size(C,2)*fcontrast);
    matrix_for_the_best(1,1)=sc1c;matrix_for_the_best(2,1)=sc2c;
    
    fcontrastb=finv(1-(alpha/size(C,1)),size(C,2)-1,sum(n)-size(C,2));
    sc1cb=sqrt(onewayanova(2)*sum((C.^2)/n)/(sum(n)-size(C,2)))*sqrt(size(C,2)*fcontrastb);
    sc2cb=-sqrt(onewayanova(2)*sum((C.^2)/n)/(sum(n)-size(C,2)))*sqrt(size(C,2)*fcontrastb);
    matrix_for_the_best(1,2)=sc1cb;matrix_for_the_best(2,2)=sc2cb;
else
    f=finv(1-alpha,size(C,2),sum(n)-size(C,2));
    sc1=sqrt(onewayanova(2)*sum((C.^2)/n)/(sum(n)-size(C,2)))*sqrt(size(C,2)*f);
    sc2=-sqrt(onewayanova(2)*sum((C.^2)/n)/(sum(n)-size(C,2)))*sqrt(size(C,2)*f);
    matrix_for_the_best(1,3)=sc1;matrix_for_the_best(2,3)=sc2;
    
    fb=finv(1-(alpha/size(C,1)),size(C,2),sum(n)-size(C,2));
    sc1b=sqrt(onewayanova(2)*sum((C.^2)/n)/(sum(n)-size(C,2)))*sqrt(size(C,2)*fb);
    sc2b=-sqrt(onewayanova(2)*sum((C.^2)/n)/(sum(n)-size(C,2)))*sqrt(size(C,2)*fb);
    matrix_for_the_best(1,4)=sc1b;matrix_for_the_best(2,4)=sc2b;
end


%Sidak
if(orthagonalbool)
    newalpha=Sidak_correction(alpha,size(C,1));
    fs=finv(1-newalpha,size(C,2),sum(n)-size(C,2));
    sc1s=sqrt(onewayanova(2)*sum((C.^2)/n)/(sum(n)-size(C,2)))*sqrt(size(C,2)*fs);
    sc2s=-sqrt(onewayanova(2)*sum((C.^2)/n)/(sum(n)-size(C,2)))*sqrt(size(C,2)*fs);
    matrix_for_the_best(1,5)=sc1s;matrix_for_the_best(2,5)=sc2s;
end

%Tukey
if(pairwisebool)
    q=qdist(alpha,size(C,2),sum(n)-size(C,2));
    tk1=q*sqrt(onewayanova(2))/sqrt((sum(n)-size(C,2))*n(1));
    tk2=-q*sqrt(onewayanova(2))/sqrt((sum(n)-size(C,2))*n(1));
    matrix_for_the_best(1,6)=tk1;matrix_for_the_best(2,6)=tk2;
    
    qb=qdist(alpha/size(C,1),size(C,2),sum(n)-size(C,2));
    tk1b=qb*sqrt(onewayanova(2))/sqrt((sum(n)-size(C,2))*n(1));
    tk2b=-qb*sqrt(onewayanova(2))/sqrt((sum(n)-size(C,2))*n(1));
    matrix_for_the_best(1,7)=tk1b;matrix_for_the_best(2,7)=tk2b;
    if(orthagonalbool)
        newalpha=Sidak_correction(alpha,size(C,1));
        qb=qdist(newalpha/size(C,1),size(C,2),sum(n)-size(C,2));
        tk1s=qb*sqrt(onewayanova(2))/sqrt((sum(n)-size(C,2))*n(1));
        tk2s=-qb*sqrt(onewayanova(2))/sqrt((sum(n)-size(C,2))*n(1));
        matrix_for_the_best(1,8)=tk1s;matrix_for_the_best(2,8)=tk2s;
    end
end
%prompt = "Choose methods: ";
%choice = input(prompt,'s')
choice="Best";

if(choice=="Scheffe")
    if(contrastbool)
        bounds=[sc1c sc2c];
    else
        bounds=[sc1 sc2];
    end
end

if(choice=="Bonferroni")
    if(contrastbool)
        bounds=[sc1cb sc2cb];
    elseif(pairwisebool)
        bounds=[tk1b tk2b];
    else
        bounds=[sc1b sc2b]
    end
end

if(choice=="Tukey")
    if(pairwisebool)
        bounds=[tk1 tk2];
    end
end

if(choice=="Sidak")
    if(orthagonalbool)
        bounds=[sc1s sc2s];
    elseif(pairwisebool)
        bounds=[tk1s tk2s];
    end
end

if(choice=="Best")
    matrix_for_the_best
    a=matrix_for_the_best(1,:);b=matrix_for_the_best(2,:);
    c=abs(a-b);nonzerominimum=min(c(c > 0));
    k=find(c==nonzerominimum);
    bound1=matrix_for_the_best(1,k);bound2=matrix_for_the_best(2,k);
    bounds=[bound1 bound2];
end

end

function [bounds]=ANOVA1_test_linear_combs(dataset,C,alpha,d,choice)
X=transpose(dataset);
onewayanova=ANOVA1_partition_TSS(dataset);
%deriving the contrasts
contrastbool=1;
for i=1:size(C,1)
    bool1=ANOVA1_is_contrast(C(i,:));
    if(~bool1)
        contrastbool=0;
    end
end
%determine ci x x bar i in sum
sample_mean_of_ith=[];
column_sizes=[];
for i=1:size(dataset,2)
    x_ith_bar=0;
    counter=0;
    for j=1:size(dataset,1)
       x_ith_bar=x_ith_bar+dataset(j,i);
        if dataset(j,i)~=0
           counter=counter+1;
        end
    end
    sample_mean_of_ith(i)=x_ith_bar/counter
    column_sizes(i)=counter;
end
a=sum(C*transpose(sample_mean_of_ith))-sum(d);

%finding ni's
for i=1:size(X,1)
    n(i)=0;
    for j=1:size(X,2)
        if(X(i,j)~=0)
            n(i)=n(i)+1;
        end
    end
end

%deriving the orthagonality
orthagonalbool=1;
for i=1:size(C,1)
    for j=1:size(C,2)
        bool2=ANOVA1_is_orthagonal(n,C(i,:),C(j,:));
        if(~bool2)
            orthagonalbool=0;
            break
        end
    end
    if(~bool2)
        break
    end
end 

%pairwise difference check
pairwisebool=0;
for i=1:size(C,1)
    nzeros=numel(C(i,:))-nnz(C(i,:));
    x=sum(C(i,:));
    if(x==0&&nzeros==size(C,2)-2)
        pairwisebool=1;
    end
end
for i=1:size(n,2)-1
    for j=i+1:size(n,2)
        if(n(i)~=n(j))
            pairwisebool=0;
            break
        end
    end
    if(pairwisebool==0)
        break
    end
end
matrix_for_the_best=zeros(2,8);
%Scheffe
if(contrastbool)
    fcontrast=finv(1-alpha,size(C,2)-1,sum(n)-size(C,2));
    sc1c=sqrt(onewayanova(2)*sum((C.^2)/n)/(sum(n)-size(C,2)))*sqrt(size(C,2)*fcontrast);
    sc2c=-sqrt(onewayanova(2)*sum((C.^2)/n)/(sum(n)-size(C,2)))*sqrt(size(C,2)*fcontrast);
    ratio=a/(sc1c/sqrt(fcontrast));pvalue=1-tcdf(ratio,size(C,2)-1);
    matrix_for_the_best(1,1)=sc1c;matrix_for_the_best(2,1)=sc2c;matrix_for_the_best(3,1)=pvalue;
    
    fcontrastb=finv(1-(alpha/size(C,1)),size(C,2)-1,sum(n)-size(C,2));
    sc1cb=sqrt(onewayanova(2)*sum((C.^2)/n)/(sum(n)-size(C,2)))*sqrt(size(C,2)*fcontrastb);
    sc2cb=-sqrt(onewayanova(2)*sum((C.^2)/n)/(sum(n)-size(C,2)))*sqrt(size(C,2)*fcontrastb);
    ratio=a/(sc1cb/sqrt(fcontrastb));pvalue=1-tcdf(ratio,size(C,2)-1);
    matrix_for_the_best(1,2)=sc1cb;matrix_for_the_best(2,2)=sc2cb;matrix_for_the_best(3,2)=pvalue;
else
    f=finv(1-alpha,size(C,2),sum(n)-size(C,2));
    sc1=sqrt(onewayanova(2)*sum((C.^2)/n)/(sum(n)-size(C,2)))*sqrt(size(C,2)*f);
    sc2=-sqrt(onewayanova(2)*sum((C.^2)/n)/(sum(n)-size(C,2)))*sqrt(size(C,2)*f);
    ratio=a/(sc1/sqrt(f));pvalue=1-tcdf(ratio,size(C,2));
    matrix_for_the_best(1,3)=sc1;matrix_for_the_best(2,3)=sc2;matrix_for_the_best(3,3)=pvalue;
    
    fb=finv(1-(alpha/size(C,1)),size(C,2),sum(n)-size(C,2));
    sc1b=sqrt(onewayanova(2)*sum((C.^2)/n)/(sum(n)-size(C,2)))*sqrt(size(C,2)*fb);
    sc2b=-sqrt(onewayanova(2)*sum((C.^2)/n)/(sum(n)-size(C,2)))*sqrt(size(C,2)*fb);
    ratio=a/(sc1b/sqrt(fb));pvalue=1-tcdf(ratio,size(C,2));
    matrix_for_the_best(1,4)=sc1b;matrix_for_the_best(2,4)=sc2b;matrix_for_the_best(3,4)=pvalue;
end


%Sidak
if(orthagonalbool)
    newalpha=Sidak_correction(alpha,size(C,1));
    fs=finv(1-newalpha,size(C,2),sum(n)-size(C,2));
    sc1s=sqrt(onewayanova(2)*sum((C.^2)/n)/(sum(n)-size(C,2)))*sqrt(size(C,2)*fs);
    sc2s=-sqrt(onewayanova(2)*sum((C.^2)/n)/(sum(n)-size(C,2)))*sqrt(size(C,2)*fs);
    ratio=a/(sc1s/sqrt(fs));pvalue=1-tcdf(ratio,size(C,2));
    matrix_for_the_best(1,5)=sc1s;matrix_for_the_best(2,5)=sc2s;matrix_for_the_best(3,5)=pvalue;
end

%Tukey
if(pairwisebool)
    q=qinv(1-alpha,size(C,2),sum(n)-size(C,2));
    tk1=q*sqrt(onewayanova(2))/sqrt((sum(n)-size(C,2))*n(1));
    tk2=-q*sqrt(onewayanova(2))/sqrt((sum(n)-size(C,2))*n(1));
    ratio=a/(tk1/q);pvalue=1-tcdf(ratio,size(C,2));
    matrix_for_the_best(1,6)=tk1;matrix_for_the_best(2,6)=tk2;matrix_for_the_best(3,6)=pvalue;
    
    qb=qinv(1-(alpha/size(C,1)),size(C,2),sum(n)-size(C,2));
    tk1b=qb*sqrt(onewayanova(2))/sqrt((sum(n)-size(C,2))*n(1));
    tk2b=-qb*sqrt(onewayanova(2))/sqrt((sum(n)-size(C,2))*n(1));
    ratio=a/(tk1b/qb);pvalue=1-tcdf(ratio,size(C,2));
    matrix_for_the_best(1,7)=tk1b;matrix_for_the_best(2,7)=tk2b;matrix_for_the_best(3,7)=pvalue;
    if(orthagonalbool)
        newalpha=Sidak_correction(alpha,size(C,1));
        qs=qinv(1-newalpha,size(C,2),sum(n)-size(C,2));
        tk1s=qs*sqrt(onewayanova(2))/sqrt((sum(n)-size(C,2))*n(1));
        tk2s=-qs*sqrt(onewayanova(2))/sqrt((sum(n)-size(C,2))*n(1));
        ratio=a/(tk1s/qs);pvalue=1-tcdf(ratio,size(C,2));
        matrix_for_the_best(1,8)=tk1s;matrix_for_the_best(2,8)=tk2s;matrix_for_the_best(3,8)=pvalue;
    end
end
prompt = "Choose methods: ";
%choice = input(prompt,'s')
%choice="Best";

if(choice=="Scheffe")
    if(contrastbool)
        bounds=matrix_for_the_best(:,1);
    else
        bounds=matrix_for_the_best(:,3);
    end
end

if(choice=="Bonferroni")
    if(contrastbool)
        bounds=matrix_for_the_best(:,2);
    elseif(pairwisebool)
        bounds=matrix_for_the_best(:,7);
    else
        bounds=matrix_for_the_best(:,4);
    end
end

if(choice=="Tukey")
    if(pairwisebool)
        bounds=matrix_for_the_best(:,6);
    end
end

if(choice=="Sidak")
    if(orthagonalbool)
        bounds=matrix_for_the_best(:,5);
    elseif(pairwisebool)
        bounds=matrix_for_the_best(:,8);
    end
end

if(choice=="Best")
    matrix_for_the_best
    a=matrix_for_the_best(1,:);b=matrix_for_the_best(2,:);
    c=abs(a-b);nonzerominimum=min(c(c > 0));
    k=find(c==nonzerominimum);
    bound1=matrix_for_the_best(1,k);bound2=matrix_for_the_best(2,k);pvalue=matrix_for_the_best(3,k);
    bounds=[bound1 bound2 pvalue];
end


end