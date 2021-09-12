clear all;
close all;

%N=������ �������� �����������
N=25;
error_true=0; %������ ��� ������������� � ������� ���������

%��������� �(�) ���� ���� �������� dp ����� fp
%��������� �(�) ���� ���� �������� dc ����� fs
dp=1;
dc=0;
fp=0.4*pi;
fs=0.5*pi;

%������ �����������
rot=0;
%����������������� �������� ����������
f=0:(pi/(100*N)):pi;

%������� ��� ���������� ��� ����� ���������
%��� ������������ ��� ������ ���������� ��� fd
k=1;
for i=1:1:length(f)
    if f(i)<=fp || f(i)>=fs
        fd(k)=f(i);
        k=k+1;
    end        
end

pi=0;
o=0;
%�������� 0, �
for i=1:1:length(fd)
    if fd(i)==pi
        p=1; %������ �������
    end
    if fd(i)==0
        o=1; %������ �������
    end
end
if p==0
   fd(length(fd)+1)=pi;
end
if o==0
   fd(length(fd)+1)=0;
end
fd=sort(fd);

%��� ��� ��������� ������
w=zeros(1,length(fd));
for i=1:1:length(fd)
    if fd(i)<=fp
        w(i)=1;
    end
    if fd(i)>=fs
        w(i)=2;
    end
end

%����������� ��������
d=ones(1,length(fd));
for i=1:1:length(fd)
    if fd(i)<=fp
        d(i)=dp;
    else
        d(i)=dc;        
    end
end

%i-> frequency
%j-> number of functions

for j=1:1:N+1
    if j==1
        for i=1:length(fd)
            matrix(i,j)=1;
        end
    else
        for i=1:length(fd)
            matrix(i,j)=2*cos((j-1)*fd(i));
        end
    end
end

[c,de]=linear_prog(N,d,w,matrix,fd);

%��������� �������
for i=1:1:length(fd)
    filter(i)=matrix(i,:)*c;
end
%������������ ��������� ���������
for i=1:1:length(fd)
        E(i)=w(i)*(d(i)-filter(i));
end
error=abs(max(E)-abs(min(E)));
plot(fd,filter); title('Filter'); xlabel('Frequency'); grid on;

%-----------------------------------------------------------------------------------------
function [c,de]=linear_prog(N,d,w,matrix,fd)
%������� F
f=zeros(1,N+2);
f(N+2)=1;


linear_matrix=zeros(2*length(fd),N+2);
delta=zeros(2*length(fd),1);
k=1;
j=1;
while j<=(2* length(fd)-1)
    linear_matrix(j,1:N+1)=matrix(k,:);
    linear_matrix(j,N+2)=-1/w(k);
    linear_matrix(j+1,1:N+1)=-matrix(k,:);
    linear_matrix(j+1,N+2)=-1/w(k);
    delta(j)=d(k);
    delta(j+1)=-d(k);
    j=j+2;
    k=k+1;
end

c=linprog(f,linear_matrix,delta);
de=c(length(c));
c=c(1:length(c)-1);
end



