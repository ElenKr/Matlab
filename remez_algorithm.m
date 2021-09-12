clear all;
close all;

%N=������ �������� �����������
N=25;
%��������� ��������� ��� ������������� � ������� ���������
error_true=0;
%������ �����������
rot=0;
%��������� �(�) ���� ���� �������� dp ����� fp
%��������� �(�) ���� ���� �������� dc ����� fs
dp=1;
dc=0;
fp=0.4*pi;
fs=0.5*pi;

%����������������� �������� ����������
f=0:(pi/(100*N)):pi;

%������� ��� ���������� ��� ����� ��������
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
%����������� �� �������� ���������� ���� ��� �������� ��� �,�
fd=sort(fd);

%����������� ��������� ������
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
            array(i,j)=1;
        end
    else
        for i=1:length(fd)
            array(i,j)=2*cos((j-1)*fd(i));
        end
    end
end

[E,c,error,rot,de]=remezalg(w,d,array,N,fd,error_true,rot);
for i=1:1:length(fd)
    
    filter(i)=array(i,:)*c;
end
%plot(fd,filter); xlabel('Frequency'); ylabel('������'); title('Filter'); grid on;

%------------------------------------------------------------------------

function [E,c,error,rot,de]=remezalg(w,d,array,N,fd,error_true,rot)

%���� 1
%����������� ����������� �� ������ ����� ����������
A=zeros(N+1,1);
for i=1:1:length(fd)
    A=A+d(i)*(w(i).^2)*(array(i,:).');
end

B=zeros(N+1,N+1);
for i=1:1:length(fd)
    x=array(i,:).';
    B=B+(w(i).^2)*(x*array(i,:));
end

%c=inv(B)*A;
c=linsolve(B,A);

%-----------------------------------------------------------------
%���� 2
%������ ��������� ��� ���������� 
%-----------------------------------------------------------------
%���� 3
while error_true==0
    
    for i=1:1:length(fd)
        pros(i)=array(i,:)*c;
    end

    %������������ ��������� ���������
    for i=1:1:length(fd)
        E(i)=w(i)*(d(i)-pros(i));
    end
    plot(fd,E); title('������������ ��������� ���������'); xlabel('Frequency'); grid on;
    axis([0 pi min(E)-0.001 max(E)+0.001]);
    %������ ������� ���������
    [megista,position11]=findpeaks(E); %���������� �� ������ �������
    e=-E;
    [elaxista,position22]=findpeaks(e); %���������� ��� ���� ��� ������� ���������
    elaxista=-elaxista;
    position1=position11;
    position2=position22;
    position11(length(position11)+1)=1;
    position11(length(position11)+1)=2252;
    %������������ ��� ������ ��� ������ �� ����������
    for i=1:1:length(position1)
        position1(i)=fd(position1(i));
    end

    for i=1:1:length(position2)
        position2(i)=fd(position2(i));
    end

    %���������� ���������� ���� ��� ������ ��� ���������
    position=zeros(length(position1)+length(position2),2);
    for i=1:1:length(position1)
        position(i,:)=[megista(i) position1(i)];
    end

    j=1;

    for i=length(position1)+1:1:length(position)
        position(i,:)=[elaxista(j) position2(j)];
        j=j+1;
    end
    position(length(position)+1,1)= -0.00014;
    position(length(position),2)= 0;
    position(length(position)+1,1)= -0.00011;
    position(length(position),2)= 3.1416;
    
    position=sortrows(position,2);

    %-----------------------------------------------------------------
    %���� 4
    %������� ��� ������ ��� ����������
    while length(position)~=N+2

        if mod(length(position)-N-1,2)==1

            if abs(position(1,1))<abs(position(length(position),1))
                position=position(2:length(position),:);
            else
                position=position(1:length(position)-1,:);

            end
        else     
            [m,index]=min(abs(position(:,1)));

            if index==length(position) 
                m_right=abs(position(1,1));
                m_left=abs(position(index-1,1));
                if m_right<m_left
                    position=position(2:length(position)-1,:);
                else
                   position=position(1:length(position)-2,:);       
                end
            elseif index==1
                m_right=abs(position(2,1));
                m_left=abs(position(length(position),1));  
                if m_right<m_left
                    position=position(3:length(position),:);
                else
                   position=position(2:length(position)-1,:);       
                end
            else               
                m_right=abs(position(index+1,1));
                m_left=abs(position(index-1,1));
            end
            if m_right>m_left
                temp=zeros(length(position)-2,2);
                for i=1:1:index-2
                    temp(i,:)=position(i,:);
                end
                for i=index+1:1:length(temp)
                    temp(i-1,:)=position(i,:);
                end
                position=temp;
            end
        end
    end



    %------------------------------------------------------------------
    %���� 5
    error=abs(max(position(:,1))-abs(min(position(:,1))));
    %��������� �� ������������� �� ��� ������������ ��� ��� ��� ���� ��� �����������
    if error>(abs(min(position(:,1))/100))
        error_true=0;
        remez_matrix=zeros(length(position),N+2);        
        new=zeros(1,length(position11)+length(position22));

        for i=1:1:length(position11)
            new(i)=position11(i);
        end
        j=1;

        for i=length(position11)+1:1:length(position22)+length(position11)
            new(i)=position22(j);
            j=j+1;
        end

        new=sort(new);
        delta=ones(1,length(position11)+length(position22));

        for i=1:1:length(new)
            delta(i)=d(new(i));
        end
        
        for i=1:1:length(new)
            w2(i)=w(new(i));
        end
        
        for i=1:1:N+2
            for j=1:1:length(position)
                if j==1
                    remez_matrix(i,j)=1;
                end
                if j==(N+2)
                    remez_matrix(i,j)=((-1)^(i+1))/(w2(i));
                end
                if j~=1 && j~=N+2
                    remez_matrix(i,j)=2*cos((j-1)*position(i,2));    
                end
             end
        end
        remez_matrix=inv(remez_matrix);
        delta=delta.';
        %����������� ��� ���� �����������
        c=remez_matrix*delta;
        de=c(length(c));
        c=c(1:length(c)-1);
        rot=rot+1;
    else
        error_true=1;
    end    
end
end






