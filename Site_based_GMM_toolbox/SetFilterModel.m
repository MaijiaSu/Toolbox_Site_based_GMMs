function [GMFilter,FilterNofPar,par_type] = SetFilterModel(type)
% par_type = 1,2,3,nan
% = 1, frequency type with support [0,inf]
% = 2, damper ratio with support [0,1]
% = 3, convex combination factor with support [0,1]

if strcmp(type, 'IIorderFilter')
    %  second order filter
    GMFilter = @(b,x)((b(1).^4)./...
        ((b(1).^2.-x.^2).^2.+4.*b(2).^2.*b(1).^2.*x.^2)).*b(3);
    FilterNofPar = 3;
    par_type = [1,2,nan];

elseif strcmp(type, 'KanaiTajimi')
    %  Kanai Tajimi filter
    GMFilter = @(b,x)( (b(1).^4+4.*b(2).^2.*b(1).^2.*x.^2)./...
        ((b(1).^2.-x.^2).^2.+4.*b(2).^2.*b(1).^2.*x.^2)).*b(3);
    FilterNofPar = 3;
    par_type = [1,2,nan];
    
elseif strcmp(type, 'CloughPenzien')
    %  Clough Penzien filter
     GMFilter = @(b,x)( (b(1).^4+4.*b(2).^2.*b(1).^2.*x.^2)./((b(1).^2.-x.^2).^2.+4.*b(2).^2.*b(1).^2.*x.^2).*... 
                         (x.^4)./((b(3).^2.-x.^2).^2.+4.*b(4).^2.*b(3).^2.*x.^2)).*b(5);
     FilterNofPar = 5;
     par_type = [1,2,1,2,nan];

elseif strcmp(type, 'CloughPenzien_hf')
    %  Clough Penzien filter
     GMFilter = @(b,x)( (b(1).^4+4.*b(2).^2.*b(1).^2.*x.^2)./((b(1).^2.-x.^2).^2.+4.*b(2).^2.*b(1).^2.*x.^2).*... 
                         (x.^4)./((b(3).^2.-x.^2).^2.+4.*1.^2.*b(3).^2.*x.^2)).*b(4);
     FilterNofPar = 4;
      par_type = [1,2,1,nan];

elseif strcmp(type, 'ConvexII-2')
    %  Convex combination of second order filter
    GMFilter = @(b,x)( b(5)*((b(1).^4)./((b(1).^2.-x.^2).^2.+4.*b(2).^2.*b(1).^2.*x.^2))+...
        ( (1-b(5))*(b(3).^4)./((b(3).^2.-x.^2).^2.+4.*b(4).^2.*b(3).^2.*x.^2)) ) .*b(6);
    FilterNofPar = 6;
     par_type = [1,2,1,2,3,nan];
    
elseif strcmp(type, 'CascadeII-2')
    %  Cascade combination of second order filter
    GMFilter = @(b,x)( ((b(1).^4)./((b(1).^2.-x.^2).^2.+4.*b(2).^2.*b(1).^2.*x.^2)).*...
        ( (b(3).^4)./((b(3).^2.-x.^2).^2.+4.*b(4).^2.*b(3).^2.*x.^2)) ) .*b(5);
    FilterNofPar = 5;
    par_type = [1,2,1,2,3,nan];

elseif strcmp(type, 'Convex_2-KT')
    %  Kanai Tajimi filter
    GMFilter = @(b,x)( b(5)* (b(1).^4+4.*b(2).^2.*b(1).^2.*x.^2)./((b(1).^2.-x.^2).^2.+4.*b(2).^2.*b(1).^2.*x.^2)+...
        ((1-b(5))* (b(3).^4+4.*b(4).^2.*b(3).^2.*x.^2)./((b(3).^2.-x.^2).^2.+4.*b(3).^2.*b(4).^2.*x.^2)) ).*b(6);
    FilterNofPar = 6;
    par_type = [1,2,1,2,3,nan];

elseif strcmp(type, 'Cascade_2-KT')
    %  Kanai Tajimi filter
    GMFilter = @(b,x)( (b(1).^4+4.*b(2).^2.*b(1).^2.*x.^2)./((b(1).^2.-x.^2).^2.+4.*b(2).^2.*b(1).^2.*x.^2).*...
        (b(3).^4+4.*b(4).^2.*b(3).^2.*x.^2)./((b(3).^2.-x.^2).^2.+4.*b(3).^2.*b(4).^2.*x.^2) ) .*b(5);
    FilterNofPar = 5;   
    par_type = [1,2,1,2,nan];

elseif strcmp(type, 'Convex-CP_II')
    %  Convex combination of one CP filter and one II-order filter
%     GMFilter = @(b,x)( b(6)*((b(1).^4)./((b(1).^2.-x.^2).^2.+4.*b(2).^2.*b(1).^2.*x.^2).*... 
%                          (x.^4)./((b(3).^2.-x.^2).^2.+4.*1.^2.*b(3).^2.*x.^2))+...  % CP-filter, 3 Pars; the bandwidth of second filter is 1
%         ( (1-b(6))*(b(4).^4)./((b(4).^2.-x.^2).^2.+4.*b(5).^2.*b(4).^2.*x.^2)) )...    % II-order filter
%                                                                               .*b(7);  % Normalization factor
%     FilterNofPar = 7;
    % version 2
    GMFilter = @(b,x)( b(7)*((b(1).^4)./((b(1).^2.-x.^2).^2.+4.*b(2).^2.*b(1).^2.*x.^2).*... 
                         (x.^4)./((b(3).^2.-x.^2).^2.+4.*b(4).^2.*b(3).^2.*x.^2))+...  % CP-filter, 4 Pars;
        ( (1-b(7))*(b(5).^4)./((b(5).^2.-x.^2).^2.+4.*b(6).^2.*b(5).^2.*x.^2)))...    % II-order filter
                                                                              .*b(8);  % Normalization factor
    FilterNofPar = 8;
    par_type = [1,2,1,2,1,2,3,nan];

end

end