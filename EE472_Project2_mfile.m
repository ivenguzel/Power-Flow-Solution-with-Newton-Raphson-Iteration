
%% PROJECT1 
%READING DATA FROM THE TEXT FILE
clear all;
path_ieee_cdf = 'C:\Users\Iven\Documents\git\Power-Flow-Solution-with-Newton-Raphson-Iteration\Test Cases for Ybus\test7.dat';

fileID = fopen(path_ieee_cdf); %opens the file for binary read access by giving file identifier
tline = fgetl(fileID); %returns the next line of fileID (removes nexline char)
file = tline;

while ischar(tline)
    tline = fgetl(fileID);
    file = char(file,tline);
end

fclose(fileID);

%% FIND the NUMBER of BUSES
first_row = 1;

while(1)   %Find where BUS info starts
    if file(first_row,1:3) == 'BUS'
        break;
    else
        first_row = first_row + 1;
    end
end

first_row = first_row + 1; %Because bus data starts in the next line from 'BUS'
row_num = first_row;

while(1)  %Count the rows until you see '-999'
    if file(row_num,(1:4)) == '-999'
        break;
    end
row_num = row_num + 1;
end

Nbus_first = first_row;
Nbus_last = row_num -1 ;  %row_num is the line with -999
Nbus = Nbus_last - Nbus_first +1; %Number of buses in the system

%% FIND the NUMBER of BRANCHES
 Nbranch_start = row_num + 2;
 row_num = Nbranch_start;
 
 while(1)     %Look for -999
    if file(row_num,(1:4)) == '-999'
        break;
    end
    row_num=row_num+1;
 end

Nbranch_last = row_num - 1; %row_num is where -999 is.
Nbranch = Nbranch_last - Nbranch_start +1; %Number of branches

%% GETTING BRANCH INFO and UPDATING YBUS
% We will use these from branch info:
% Sending bus: Column 1-4
% Receiving bus: Column 6-9
% Type of the line: Column 19
% R: Column 20-29
% X: Column 30-40
% Total B: Column 41-50
% 
% We will use these from bus info:
% G: Column 107-114
% B: Column 115-122

%I dont want to store any of these values 

BusIndex = zeros(Nbus,1);
Ybus = sparse(Nbus,Nbus);

for m=1:Nbus
    G = str2double(file(Nbus_first + (m-1) ,(107:114)));
    B = str2double(file(Nbus_first + (m-1),(115:122)));
    Ybus(m,m)= G + 1i*B ;
    BusIndex(m,1) =str2double(file(Nbus_first+ (m-1) ,(1:4))) ; %Bus names might be different than our index, so store these.
end
    

for m=1:Nbranch
    send_bus = str2double(file(Nbranch_start+m-1,(1:4)));
    rec_bus = str2double(file(Nbranch_start+m-1,(6:9)));
    R = str2double(file(Nbranch_start+m-1,(20:29)));
    X = str2double(file(Nbranch_start+m-1,(30:40)));
    B = str2double(file(Nbranch_start+m-1,(41:50)));
    Turn_ratio = str2double(file(Nbranch_start+m-1,(77:82)));
    
    TR_type = str2double(file(Nbranch_start+m-1,(19)));
    
    if TR_type == 4
        Phase_angle = str2double(file(Nbranch_start+m-1,(84:90)));
        Turn_ratio = (Turn_ratio*cos(Phase_angle*0.0174532925)+ Turn_ratio*sin(Phase_angle*0.0174532925) ); %For Phase shifters 1 deg=0.0174532925 radians
    end
    
	if Turn_ratio == 0 %they put zero sometimes to indicate no transformers
       Turn_ratio = 1;
    end
    
    Y1 = find(BusIndex == send_bus);           % To see which index sending bus name corresponds
    Y2 = find(BusIndex == rec_bus);       % To see which index receiving bus name corresponds 
    
    Ybus(Y1,Y1) = Ybus(Y1,Y1) + (1/(R + 1i* X)) / (Turn_ratio^2) + 1i*B/2; % When the type is TL it does not hange anything since n=1
    Ybus(Y2,Y2) = Ybus(Y2,Y2)+ (1/(R + 1i* X)) + 1i*B/2;
    Ybus(Y1,Y2) = Ybus(Y1,Y2) - ((1/(R + 1i* X)) / (conj(Turn_ratio)));
    Ybus(Y2,Y1) = Ybus(Y2,Y1) - ((1/(R + 1i* X)) /Turn_ratio) ;   
  
end

Yb=Ybus;

%% PROJECT 2:





    
