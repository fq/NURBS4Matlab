function nurbs2iges(nurbs, igesfile, units)
%Convert NURBS curves or surfaces to an IGES file.
%
% Parameters:
%   nurbs    - NURBS curves or surfaces, objects created by FQ's NURBS code
%   igesfile - Name of the IGES file to which the NURBS objects are being saved
%   units    - Units, 'in','inch' or 'mm', default is 'inch'
%
% Example:
%
%   Save a NURBS curve object 'crv' to IGES file 'my_iges.igs'
%       nurbs2iges(crv,'my_iges.igs','mm');
%
%   Save a NURBS surface object 'srf' to IGES file 'my_iges.igs'
%       nurbs2iges(srf,'my_iges.igs','mm');
%
%   Save several NURBS objects (curves or surfaces) to IGES file 'my_iges.igs'
%       nurbs2iges([srf1 srf2 crv1 crv2],'my_iges.igs','mm');
% 

%  Fu Qiang
%  2006-10-02

if ~isjava(nurbs)
    error('First input: Not a object create by NURBS code!')
end

if ~ischar(igesfile)
    error('Second input: Invalid filename.')
end

if nargin < 3
    units = 'IN';
    fprintf(' Nurbs2IGES: Using default units "Inches"\n');
end
units = upper(units);
switch units
   case {'IN','INCH','MM'}
      ;
   otherwise
      error('Unknown units.')
end

dt = datestr(now,'yyyymmdd.HHMMSS');

%IGES file structure
%Start Section, single line
S = 'IGES from Matlab NURBS code. Written by Fu,Qiang.';
%Global Section
G{1} = '1H,';      % Parameter Deliminator Character
G{2} = '1H;';      % Record Delimiter Character
G{3} = HString('Matlab Nurbs2IGES');      % Product ID from Sender
G{4} = HString(igesfile);   % File Name
G{5} = HString('Matlab NURBS'); % System ID
G{6} = HString('Nurbs2IGES V1.0'); % Pre-processor Version
G{7} = '32';       % Number of Bits for Integers
G{8} = '75';        % Single Precision Magnitude
G{9} = '6';       % Single Precision Significance
G{10}= '75';       % Double Precision Magnitude
G{11}= '15';       % Double Precision Significance
G{12}= HString('Nurbs from Matlab');  % Product ID for Receiver
G{13}= '1.0';  % Model Space Scale
G{14}= '3';    % Unit Flag (1 = inches, 3 = look to index 15 name)
G{15}= HString(units);     % Units  (Inches = "IN")
G{16}= '1000';        % Maximum Number of Line Weights
G{17}= '1.0'; % Size of Maximum Line Width
G{18}= HString(dt);   % Date and Time Stamp
G{19}= '0.001'; % Minimum User-intended Resolution
G{20}= '10000.0';   % Approximate Maximum Coordinate
G{21}= HString('Fu Qiang');     % Name of Author
G{22}= HString('CADCAM Group - Concordia U'); % Author's Organization
G{23}= '11'; % IGES Version Number
G{24}= '0';  % Drafting Standard Code
G{25}= HString(dt); % Model Creation/Change Date
% Convert section array to lines
SectionG = make_section(G,72);

%Directory Entry Section
D = [];
for i=1:length(nurbs)
    D(i).type = nurbs(i).EntityType;
    D(i).id = i*2-1;
    D(i).p_start = 0;
    D(i).p_count = 0;
end

%Parameter Data Section.
SectionP = {};
for i =1:length(nurbs)
    P = make_section_array(nurbs(i)); %finish one entity
    % Convert section array to lines
    SP = make_section(P,64);
    D(i).p_count = length(SP);
    if i == 1
        D(i).p_start = 1;
    else
        D(i).p_start = D(i-1).p_start + D(i-1).p_count;
    end
    SectionP{i} = SP;
end

% save sections to file
fid = fopen(igesfile,'w');
fprintf(fid,'%-72sS%7d\n',S,1);
for i=1:length(SectionG)
    fprintf(fid,'%-72sG%7d\n',SectionG{i},i);
end

for i=1:length(D)
    fprintf(fid,'%8d%8d%8d%8d%8d%8d%8d%8d%8dD%7d\n', ...
            D(i).type,D(i).p_start,0,0,0,0,0,0,0,i*2-1);
    fprintf(fid,'%8d%8d%8d%8d%8d%8s%8s%8s%8dD%7d\n', ...
            D(i).type,0,0,D(i).p_count,0,' ',' ',' ',0,i*2);
end

lines_p = 0;
for j = 1:length(D)
    sec = SectionP{j};
    for i=1:length(sec)
        lines_p = lines_p + 1;
        fprintf(fid,'%-64s %7dP%7d\n',sec{i},D(j).id,lines_p);
    end
end
sec_t = sprintf('S%7dG%7dD%7dP%7d', ...
                1,length(SectionG),2*length(D),lines_p);
fprintf(fid,'%-72sT%7d\n',sec_t,1);
fclose(fid);



function hs = HString(str)
hs = sprintf('%dH%s',length(str),str);

function P = make_section_array(nurbs)
P = {};
switch nurbs.EntityType
    case 126
        % Rational B-Spline Curve Entity
        cp = nurbs.getControlPoints;
        deg = nurbs.getDegree;
        knots = nurbs.getKnots;
        uspan = nurbs.getParamExtents;
        isplanar = ~any(cp(:,3));
        P{1} = '126';
        P{2} = int2str(size(cp,1)-1);
        P{3} = int2str(deg); 
        P{4} = int2str(isplanar);
        P{5} = '0';
        P{6} = '0';
        P{7} = '0';
        index = 8;
        for i = 1:length(knots)
            P{index} = sprintf('%f',knots(i));
            index = index + 1;
        end
        for i = 1:size(cp,1)
            P{index} = sprintf('%f',cp(i,4));
            index = index + 1;
        end
        for i = 1:size(cp,1)
            P{index} = sprintf('%f',cp(i,1));
            index = index + 1;
            P{index} = sprintf('%f',cp(i,2));
            index = index + 1;
            P{index} = sprintf('%f',cp(i,3));
            index = index + 1;
        end
        P{index} = sprintf('%f',uspan(1));index = index +1;
        P{index} = sprintf('%f',uspan(2));index = index +1;
        P{index} = '0.0';index = index +1;
        P{index} = '0.0';index = index +1;
        if isplanar
            P{index} = '1.0';
        else
            P{index} = '0.0';
        end
        index = index + 1;
        P{index} = '0';index = index + 1;
        P{index} = '0';
    case 128
        % Rational B-Spline Surface Entity
        cp = nurbs.getControlPoints;
        degU = nurbs.getDegreeU;
        degV = nurbs.getDegreeV;
        knotsU = nurbs.getKnotsU;
        knotsV = nurbs.getKnotsV;
        uspan = nurbs.getParamExtentsU;
        vspan = nurbs.getParamExtentsV;
        P{1} = '128';
        P{2} = int2str(size(cp,1)-1);
        P{3} = int2str(size(cp,2)-1);
        P{4} = int2str(degU); 
        P{5} = int2str(degV); 
        P{6} = '0';
        P{7} = '0';
        P{8} = '0';
        P{9} = '0';
        P{10} = '0';
        index = 11;
        for i = 1:length(knotsU)
            P{index} = sprintf('%f',knotsU(i));
            index = index + 1;
        end
        for i = 1:length(knotsV)
            P{index} = sprintf('%f',knotsV(i));
            index = index + 1;
        end
        for j = 1:size(cp,2)
            for i = 1:size(cp,1)
                P{index} = sprintf('%f',cp(i,j,4));
                index = index + 1;
            end
        end
        for j = 1:size(cp,2)
            for i = 1:size(cp,1)
                P{index} = sprintf('%f',cp(i,j,1));
                index = index + 1;
                P{index} = sprintf('%f',cp(i,j,2));
                index = index + 1;
                P{index} = sprintf('%f',cp(i,j,3));
                index = index + 1;
            end
        end
        P{index} = sprintf('%f',uspan(1));index = index +1;
        P{index} = sprintf('%f',uspan(2));index = index +1;
        P{index} = sprintf('%f',vspan(1));index = index +1;
        P{index} = sprintf('%f',vspan(2));index = index +1;
        P{index} = '0';index = index + 1;
        P{index} = '0';
    otherwise
        error('First input: Unknown type!')
end


function sec = make_section(fields,linewidth)
sec = {};
index = 1;
line = '';
num = length(fields);
for i = 1:num
    if i<num
        newitem = [fields{i} ','];
    else
        newitem = [fields{i} ';'];
    end
    len = length(line) + length(newitem);
    if ( len>linewidth )
        % new line
        sec{index} = line;
        index = index + 1;
        line = '';
    end
    line = [line newitem];
end
sec{index} = line;
