clear all;
close all;

%N=Πλήθος αγνωστων συντελεστών
N=25;
error_true=0; %δηλαδή ότι ικανοποιείται η συνθήκη σφάλματος

%Συνάρτηση Δ(ω) στην ζωνη διαβασης dp μεχρι fp
%Συνάρτηση Δ(ω) στην ζωνη αποκοπής dc μέχρι fs
dp=1;
dc=0;
fp=0.4*pi;
fs=0.5*pi;

%Πλήθος επαναλήψεων
rot=0;
%Δειγματοληπτημένο Διάστημα συχνοτήτων
f=0:(pi/(100*N)):pi;

%Κόβουμε τις συχνότητες της ζώνης διάβασης
%και δημιουργούμε νέο πίνακα συχνοτήτων τον fd
k=1;
for i=1:1:length(f)
    if f(i)<=fp || f(i)>=fs
        fd(k)=f(i);
        k=k+1;
    end        
end

pi=0;
o=0;
%Πρόσθεση 0, π
for i=1:1:length(fd)
    if fd(i)==pi
        p=1; %δηλαδή βρέθηκε
    end
    if fd(i)==0
        o=1; %δηλαδή βρέθηκε
    end
end
if p==0
   fd(length(fd)+1)=pi;
end
if o==0
   fd(length(fd)+1)=0;
end
fd=sort(fd);

%Για την συνάρτηση βάρους
w=zeros(1,length(fd));
for i=1:1:length(fd)
    if fd(i)<=fp
        w(i)=1;
    end
    if fd(i)>=fs
        w(i)=2;
    end
end

%Υπολογισμός Ιδανικής
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

%Ισοτικοί Περιορισμοί
freq_isot=1.413716;
A=zeros(length(freq_isot),N+1);
for i=1:1:length(freq_isot)
    for j=1:1:N+1
        if j==1
            A(i,j)=1;
        end
        if j~=1
            A(i,j)=2*cos((j-1)*freq_isot(i));
        end
    end
end
D=0.5;

[c,rot,error,E,de]=remezalg_isot(w,d,array,N,fd,error_true,freq_isot,rot,A,D);

for i=1:1:length(fd)
    filter(i)=array(i,:)*c;
end
%plot(fd,filter); xlabel('Frequency'); ylabel('Πλάτος'); title('Filter'); grid on;

%-------------------------------------------------------------------------------
function [c,rot,error,E,de]=remezalg_isot(w,d,array,N,fd,error_true,freq_isot,rot,A,D)

%ΒΗΜΑ 1
%Υπολογισμός Συντελεστών c όταν έχω ισοτικούς περιορισμούς

Q2=zeros(N+1,1);
for i=1:1:length(fd)
    Q2=Q2+d(i)*(w(i).^2)*(array(i,:).');
end

Q1=zeros(N+1,N+1);
for i=1:1:length(fd)
    x=array(i,:).';
    Q1=Q1+(w(i).^2)*(x*array(i,:));
end

z=zeros(N+1+length(freq_isot),N+1+length(freq_isot));
for i=1:1:length(z)
    if i<=N+1
       for j=1:1:length(z)
           if j<=N+1
               z(i,j)=Q1(i,j);
           end
           if j>N+1
               z(i,j)=A(i);
           end
       end
    end
       if i>N+1
           for j=1:1:N+1
               z(i,j)=A(j);
           end
       end
end
          

y=zeros(1,N+1+length(freq_isot));
for i=1:1:length(y)
    if i<=N+1
        y(i)=Q2(i);
    else
        y(i)=D(i-N-1);
    end
     
end
y=y';
c=linsolve(z,y);
c=c(1:N+1);


%-----------------------------------------------------------------
%ΒΗΜΑ 2
%Έχουμε προσθέσει τις συχνότητες 

%-----------------------------------------------------------------
%ΒΗΜΑ 3
while error_true==0
    
    for i=1:1:length(fd)
        pros(i)=array(i,:)*c;
    end

    %Προσημασμένη συνάρτηση Σφάλματος
    for i=1:1:length(fd)
        E(i)=w(i)*(d(i)-pros(i));
    end
    
    plot(fd,E); title('Προσημασμένη Συνάρτηση Σφάλαμτος'); xlabel('Frequency'); grid on;
    
    %Εύρεση τοπικών ακροτάτων
    [megista,position1]=findpeaks(E); %επιστρέφει τα τοπικά μέγιστα
    e=-E;
    [elaxista,position2]=findpeaks(e); %επιστρέφει την θέση των τοπικών ελαχίστων
    elaxista=-elaxista;
    

    %Αναδιάταξη Συχνοτήτων μετά την εύρεση των ακροτάτων
    position=zeros(length(position1)+length(position2),2);
    for i=1:1:length(position1)
        position(i,:)=[megista(i) position1(i)];
    end

    j=1;

    for i=length(position1)+1:1:length(position)
        position(i,:)=[elaxista(j) position2(j)];
        j=j+1;
    end
    
    position(length(position)+1,1)= -0.0009013;
    position(length(position),2)= 1;
    position(length(position)+1,1)= -0.001567;
    position(length(position),2)= 2247;
    
    position=sortrows(position,2);


    %-----------------------------------------------------------------
      
    %ΒΗΜΑ 4
    %Έλεγχος του πλήθος των συχνοτήτων
    while length(position)~=(N+2-length(D))

        if mod(length(position)-N-length(D),2)==1

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

    %ΒΗΜΑ 5
    error=abs(max(position(:,1))-abs(min(position(:,1))));
    %ελέγχουμε αν ικανοποιείται
    %αν οχι υπολογίζουμε την νέα των τιμή των συντελεστών           
    %Υπολογιμός την δελτα στις νέες συχνότητες    
    delta=ones(1,length(position)+length(D));
        
    for i=1:1:length(delta)
        if i<=N+2-length(D)
            delta(i)=d(position(i,2));
        else
            delta(i)=D(i-N-1);
        end
    end
    
    delta=delta.';    
    for i=1:1:length(position)
        w2(i)=w(position(i,2));
    end    
  %Μετατροπή των indexes σε frequency
    for i=1:1:length(position)
        position(i,2)=fd(position(i,2));
    end
    if error>(abs(min(position(:,1))/100))
        error_true=0;
        remez_matrix=zeros(N+2,N+2);
        for i=1:1:N+2
            for j=1:1:N+2
                if j==1
                    remez_matrix(i,j)=1;           
                    
                elseif j~=1 && i<=N+2-length(D) && j<=N+1
                    remez_matrix(i,j)=2*cos((j-1)*position(i,2));              
                
                elseif i<=N+2-length(D) && j==N+2
                    remez_matrix(i,j)=((-1)^(i+1))/w2(i);
                
                elseif i>N+2-length(D) && j<=N+1
                    remez_matrix(i,j)=2*cos((j-1)*freq_isot(i-N-1));
                end
            end
        end
        %Υπολογισμός των νέων συντελεστών
        c=linsolve(remez_matrix,delta);
        de=c(length(c));
        c=c(1:N+1);
        rot=rot+1;
    else
        error_true=1;        
    end
end
end







