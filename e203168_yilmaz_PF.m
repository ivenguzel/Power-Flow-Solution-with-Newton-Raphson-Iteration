function [Total_time,V,del,Tol,Iter] = e203168_yilmaz_PF(path)


%% Variables
tic

N_bus = [];               % bus number matrix
G = [];                   % shunt conductance
B = [];                   % shunt susceptance
Y_bus = [];               % Y Bus Matrix
bus = 1;                  
row = 0;                  
I_transpose=[];
Z_bus = [];
unknown = [];
PQ_bus=[];
PV_bus=[];
teta = [];
volt = [];
q=1;
qq=1;
qqqq=1;
number_bus=0;
NNZ=0;
BMva=100;
Iter=0;

%% Reading CDF FÝLE

content = fopen(path);   
line = fgetl(content);
file=line;

while ischar(line)
line = fgetl(content);
file =char(file,line);
end

fclose(content);


%% Determining BUS1 row number

for initial=1:1:10               
    if str2double(file(initial,(1:4))) == 1
        break;
    end
end

row=initial;

%% Classify bus parameters(G,B)


while(1)                             
 
 N_bus(bus,1) = str2double(file(row,(1:4)));
 G(bus,1) = str2double(file(row,(107:114)));
 B(bus,1) = str2double(file(row,(115:122)));
      
 Y_bus(bus,bus) = G(bus,:) + i*B(bus,:);        %diagonal entries are entered to Y_bus
 bus = bus+1;
 row = row+1;
    
   
    if file(row,(1:4)) == '-999'
        break;
    end    
       
end

number_bus = row-initial;    % Number of busses are calculated



row = row+2;   %passed to branch data
bus = 1;


%% Classify branch parameters(R,X,Type,B)


while(1)
       
    N_sending(bus,1) = str2double(file(row,(1:4)));
    N_recieving(bus,1) = str2double(file(row,(6:9)));
    type(bus,1) = str2double(file(row,(19)));
    R(bus,1) = str2double(file(row,(20:29)));
    X(bus,1) = str2double(file(row,(30:40)));
    B_line(bus,1) = 1i*str2double(file(row,(41:50)))/2;
    Turn_ratio(bus,1) = str2double(file(row,(77:82)));
            
bus = bus+1;
row = row+1;

 if file(row,(1:4)) == '-999'
        break;
 end
    
end


[row,~] = size(N_sending);

for k = 1:1:row    
     
    
    if file(row,(1:4)) == '-999'
        break;
    end
    
    
    Y_1 = find(N_bus == N_sending(k,1));           % sending matrix'indeki elemanlarýn N_bus da hangi row da olduðu 
    Y_2 = find(N_bus == N_recieving(k,1));       % recieving matrix'indeki elemanlarýn N_bus da hangi row da olduðu 
    
   if(type(k,1) == 0)         %transmission line type
        Y_bus(Y_1,Y_1) = Y_bus(Y_1,Y_1) + (1/(R(k,1) + i* X(k,1)) + B_line(k,1)) ;
        Y_bus(Y_2,Y_2) = Y_bus(Y_2,Y_2)+ ((1/(R(k,1) + i* X(k,1))) + B_line(k,1));
       
        Y_bus(Y_1,Y_2) = Y_bus(Y_1,Y_2) - (1/(R(k,1) + i* X(k,1)));
        Y_bus(Y_2,Y_1) = Y_bus(Y_2,Y_1)- (1/(R(k,1) + i* X(k,1)));       
   end
   
   if(type(k,1) == 1||type(k,1) == 2||type(k,1) == 3||type(k,1) == 4)
      
         Y_bus(Y_1,Y_1) = Y_bus(Y_1,Y_1) + ((1/(R(k,1) + i* X(k,1))) / (Turn_ratio(k,1)^2)) ;
         Y_bus(Y_2,Y_2) = Y_bus(Y_2,Y_2)+ (1/(R(k,1) + i* X(k,1)));
       
         Y_bus(Y_1,Y_2) = Y_bus(Y_1,Y_2) - ((1/(R(k,1) + i* X(k,1))) / Turn_ratio(k,1));
         Y_bus(Y_2,Y_1) = Y_bus(Y_2,Y_1) - ((1/(R(k,1) + i* X(k,1))) /Turn_ratio(k,1)) ;     
      
   end      

end


spy(Y_bus,5,'r')

G=real(Y_bus);B=imag(Y_bus);

Tol=22222;
        counter=1;

%% Classify bus parameters(V,P,Q,Theta)
        
        
        
for r=initial:1:1000
    
        if file(r,(1:4)) == '-999'
        break;
        end
        
         V(counter,1) = str2double(file(r,(28:33)));
         del(counter,1) = str2double(file(r,(34:40)));
         Pg(counter,1) = str2double(file(r,(60:67)));
         Qg(counter,1) = str2double(file(r,(68:75)));
         Pl(counter,1) = str2double(file(r,(41:49)));
         Ql(counter,1) = str2double(file(r,(50:59)));

        counter=counter+1;
end

Pg=Pg/BMva;         % per unit conversion
Pl=Pl/BMva;
Qg=Qg/BMva;
Ql=Ql/BMva;




P=zeros(number_bus,1);
Q=zeros(number_bus,1);


 for r=1:number_bus
     del(r)=pi*del(r)/180; % Conversion radian to degree 
 end

for k = 1:1:row    

     if file(k,(1:4)) == '-999'
        break;
     end
    
    if number_bus==300
       
    if (str2double(file(k,(25:26))) == 0 || str2double(file(k,(25:26))) == 1 )
      PQ_bus(qq,1)= str2double(file(k,(128:132)));
      qq=qq+1;
    end
    
    if (str2double(file(k,(25:26))) == 2 )
      PV_bus(qqqq,1)= str2double(file(k,(128:132)));
      qqqq=qqqq+1;
    end
    
      if (str2double(file(k,(25:26))) == 3 )
     Slack= str2double(file(k,(128:132)));
      end
    end   
    
    
    if number_bus~=300
       
    if (str2double(file(k,(25:26))) == 0 || str2double(file(k,(25:26))) == 1 )
      PQ_bus(qq,1)= str2double(file(k,(1:4)));
      qq=qq+1;
    end
    
    if (str2double(file(k,(25:26))) == 2 )
      PV_bus(qqqq,1)= str2double(file(k,(1:4)));
      qqqq=qqqq+1;
    end
    
      if (str2double(file(k,(25:26))) == 3 )
     Slack= str2double(file(k,(1:4)));
      end
    end   
end   

[number_PQ,~] = size(PQ_bus);
[number_PV,~] = size(PV_bus);
number_unknown = number_PQ+number_PV ; 


while Tol>10


  
for r=1:number_bus
    for j=1:number_bus
 P(r)=P(r)+V(r)*V(j)*(G(r,j)*cos(del(r)-del(j))+B(r,j)*sin(del(r)-del(j)));
 Q(r)=Q(r)+V(r)*V(j)*(G(r,j)*sin(del(r)-del(j))-B(r,j)*cos(del(r)-del(j)));
        end
end


Active_Power_specified=Pg-Pl; Reactive_Power_specified=Qg-Ql;

dPa=Active_Power_specified-P;
dQa=Reactive_Power_specified-Q;
dP=dPa([1:Slack-1,Slack+1:end],1);

dQ=zeros(number_PQ,1);
 
 
    for r=1:number_bus
        if str2double(file(k,(25:26)))==0
            dQ(k,1)=dQa(r);
            k=k+1;
        end
    end
    M=[dP;dQ]; % delta Matrix 


%% Formation Of J1 

   
J1=zeros(number_bus,number_bus);

 for r=1:number_bus
        m=r;
        for j=1:number_bus;
            n=j;
            if m==n
                for n=1:number_bus 
                J1(r,j)=J1(r,j)+V(m)*V(n)*(-G(m,n)*sin(del(m)-del(n))+B(m,n)*cos(del(m)-del(n)));
                end
                J1(r,j)=J1(r,j)-V(m)^2*B(m,m);
            else
                J1(r,j)=V(m)*V(n)*(G(m,n)*sin(del(m)-del(n))-B(m,n)*cos(del(m)-del(n)));
            end
        end
 end
%% Formation Of J2
    J2=zeros(number_bus,number_PQ);
    for r=1:number_bus
        m=r;
        for j=1:number_PQ
            n=PQ_bus(j);
            if m==n
                for n=1:number_bus
                    J2(r,j)=J2(r,j)+V(n)*(G(m,n)*cos(del(m)-del(n))+B(m,n)*sin(del(m)-del(n)));
                end
                J2(r,j)=J2(r,j)+V(m)*G(m,m);
            else
                J2(r,j)=V(m)*(G(m,n)*cos(del(m)-del(n))+B(m,n)*sin(del(m)-del(n)));
            end
        end
    end
	%% Formation Of J3
    J3=zeros(number_PQ,number_bus);
    for r=1:number_PQ
        m=PQ_bus(r);
        for j=1:number_bus
            n=j;
            if m==n
                for n=1:number_bus
                    J3(r,j)=J3(r,j)+V(m)*V(n)*(G(m,n)*cos(del(m)-del(n))+B(m,n)*sin(del(m)-del(n)));
                end
                J3(r,j)=J3(r,j)-V(m)^2*G(m,m);
            else
                J3(r,j)=V(m)*V(n)*(-G(m,n)*cos(del(m)-del(n))-B(m,n)*sin(del(m)-del(n)));
            end
        end
    end
	%% Formation Of J4
    J4=zeros(number_PQ,number_PQ);
    for r=1:number_PQ
        m=PQ_bus(r);
        for j=1:number_PQ
            n=PQ_bus(j);
            if m==n
                for n=1:number_bus
                J4(r,j)=J4(r,j)+V(n)*(G(m,n)*sin(del(m)-del(n))-B(m,n)*cos(del(m)-del(n)));
                end
                J4(r,j)=J4(r,j)-V(m)*B(m,m);
            else
                J4(r,j)=V(m)*(G(m,n)*sin(del(m)-del(n))-B(m,n)*cos(del(m)-del(n)));
            end
        end              
    end
   J=[J1 J2;J3 J4]; % Jacobian Matrix 
    
   J=J([1:Slack-1,Slack+1:end],:);
   J=J(:,[1:Slack-1,Slack+1:end]);

   X=inv(J)*M;

   
   
    dTh=X([1:Slack-1,Slack+1:number_bus],:); % Change in angle 
    dV=X(number_bus:end);	% change in volatge mag 
    del([1:Slack-1,Slack+1:end],:)=del([1:Slack-1,Slack+1:end],:)+dTh;% Voltage angle update 
	% voltage mag update 
    k=1;
    for n=1:number_bus
        if str2double(file(k,(25:26))) == 0
            V(n)=V(n)+dV(k);
            k=k+1;
        end
    end
   
   
    
    
     Iter=Iter+1;
    Tol=max(abs(M))
    
    end
    
    
    Q=zeros(number_bus,1);
 for r=1:number_bus
        for j=1:number_bus
            P(r)=P(r)+V(r)*V(j)*(G(r,j)*cos(del(r)-del(j))+B(r,j)*sin(del(r)-del(j)));
            Q(r)=Q(r)+V(r)*V(j)*(G(r,j)*sin(del(r)-del(j))-B(r,j)*cos(del(r)-del(j)));
        end
 end
 for r=1:number_bus
     del(r)=180*del(r)/pi; % Converion radian to degree 
 end
    

 
 disp('----------------------------------------');
disp('  Newton Raphson Loadflow Solution    ');
disp('----------------------------------------');
disp(' |Bus |   |Voltage|    |Angle |');
disp(' | No.|   |pu     |    |Degree|');
disp('----------------------------------------');
for m=1:number_bus 
    fprintf(' %3g   ' ,m);
    fprintf(' %8.3f    ' ,V(m));
    fprintf(' %8.3f  ' ,del(m));
     fprintf(' %8.3f  ',Pg(m)*BMva);
    if str2double(file(k,(25:26))) == 2
    fprintf(' %8.3f  ',(Q(m)+Ql(m))*BMva);
     
    end
    fprintf('\n');
end
disp('----------------------------------------');
 fprintf( 'Number Of Ieration %3g \n',Iter)
toc  
Total_time=toc;



end
