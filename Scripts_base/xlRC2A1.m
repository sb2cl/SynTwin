function CELL = xlRC2A1(ROW,COL)
%% Returns the column characters of Excel given a certain column number
% Input ROW : row number
%       COL : number of column
% Output CHAR : Character combination in Excel
% Obtained from: https://es.mathworks.com/matlabcentral/answers/54153-dynamic-ranges-using-xlswrite

    if COL <= 26                        % [A..Z]
        CHAR = char(mod(COL-1,26)+1+64);
        CELL = [CHAR num2str(ROW)];
    elseif COL <= 702                   % [AA..ZZ]
        COL = COL-26;    
        CHAR1 = char(floor((COL-1)/26)+1+64);
        CHAR0 = char(mod(COL-1,26)+1+64);
        CHAR = [CHAR1 CHAR0];
        CELL = [CHAR num2str(ROW)];
    elseif COL <= 16384                 % [AAA..XFD]
        COL = COL-702; 
        CHAR2 = char(floor((COL-1)/676)+1+64);
        COL=COL-(floor((COL-1)/676))*676;
        CHAR1 = char(floor((COL-1)/26)+1+64);
        CHAR0 = char(mod(COL-1,26)+1+64);
        CHAR = [CHAR2 CHAR1 CHAR0];
        CELL = [CHAR num2str(ROW)];
    else
        disp('Column does not exist in Excel!');
    end
end