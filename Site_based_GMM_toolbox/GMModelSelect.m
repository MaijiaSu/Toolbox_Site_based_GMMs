function FitOption = GMModelSelect(GMModelNO)
% Baseline model
FitOption.Tmodulate.type = 'spline';                  % 'spline', 'Gamma'
FitOption.Fmodulate.type = 'IIorderFilter';           % 'IIorderFilter','KanaiTajimi','ConvexII-2','CascadeII-2','CloughPenzien'
% FitOption.Fmodulate.Paramtericmodel = 'ReducedLinearmodel';  % 'Polymodel','Linearmodel','Constantmodel','ReducedLinearmodel'
FitOption.HighPassFilter.Type = 'Oscillator';         % 'Oscillator','Butterworth'
FitOption.GenMethod = 'spectral';                     % spectral or temporal
FitOption.Tmodulate.s_time = 0;

switch GMModelNO
    case 1
        FitOption.Fmodulate.Paramtericmodel = [2,1];
    case 2
        FitOption.Fmodulate.Paramtericmodel = [2,2];
    case 3
        FitOption.Fmodulate.type = 'IIorderFilter'; 
        FitOption.Fmodulate.Paramtericmodel = [3,3];
    case 4
        FitOption.Fmodulate.type = 'KanaiTajimi';
        FitOption.Fmodulate.Paramtericmodel = [2,1];
    case 5
        FitOption.Fmodulate.type = 'ConvexII-2';
        FitOption.Fmodulate.Paramtericmodel = [2,2,2,2,2];
    case 6
        FitOption.Fmodulate.type = 'ConvexII-2';
        FitOption.Fmodulate.Paramtericmodel = [3,3,3,3,3];
    case 7
        FitOption.Fmodulate.type = 'CascadeII-2';
        FitOption.Fmodulate.Paramtericmodel = [2,2,2,2];
    case 8
        FitOption.Fmodulate.type = 'CascadeII-2';
        FitOption.Fmodulate.Paramtericmodel = [3,3,3,3];
    case 9
        FitOption.Tmodulate.type = 'Gamma'; 
        FitOption.Fmodulate.Paramtericmodel = [2,1];
end
