function [Twb,Teq,epott]=WetBulb(TemperatureC,Pressure,Humidity,HumidityMode,ConvergenceMode)

% [Twb,Teq,epott]=WetBulb(TemperatureC,Pressure,Humidity,[HumidityMode])
%
% Calculate wet-bulb temperature, equivalent temperature, and equivalent
% potential temperature from temperature, pressure, and relative or
% specific humidity.
%
%
% INPUTS:
%
%    TemperatureC	2-m air temperature (degrees Celsius)
%    Pressure	        Atmospheric Pressure (Pa)
%    Humidity           Humidity -- meaning depends on HumidityMode
%    HumidityMode
%               0 (Default): Humidity is specific humidity (kg/kg)
%               1: Humidity is relative humidity (%, max = 100)
%
% TemperatureC, Pressure, and Humidity should either be scalars or arrays of
% identical dimension.
%
% OUTPUTS:
%
%    Twb	wet bulb temperature (C)
%    Teq	Equivalent Temperature (K)
%    epott 	Equivalent Potential Temperature (K)
%
%
% METHODOLOGY:
%
% Calculates Wet Bulb Temperature, Theta_wb, Theta_e, Moist Pot Temp, 
% Lifting Cond Temp, and Equiv Temp using Davies-Jones 2008 Method.
% 1st calculates the lifting cond temperature (Bolton 1980 eqn 22).  
% Then calculates the moist pot temp (Bolton 1980 eqn 24). Then 
% calculates Equivalent Potential Temperature (Bolton 1980 eqn 39).  
% From equivalent pot temp, equiv temp and Theta_w (Davies-Jones 
% 2008 eqn 3.5-3.8).  An accurate 'first guess' of wet bulb temperature
% is determined (Davies-Jones 2008 eqn 4.8-4.11). Newton-Raphson
% is used for 2 iterations, determining final wet bulb temperature 
% (Davies-Jones 2008 eqn 2.6).
%
% Reference:  Bolton: The computation of equivalent potential temperature. 
% 	      Monthly weather review (1980) vol. 108 (7) pp. 1046-1053
%	      Davies-Jones: An efficient and accurate method for computing the 
%	      wet-bulb temperature along pseudoadiabats. Monthly Weather Review 
%	      (2008) vol. 136 (7) pp. 2764-2785
% 	      Flatau et al: Polynomial fits to saturation vapor pressure. 
%	      Journal of Applied Meteorology (1992) vol. 31 pp. 1507-1513
% Note: Pressure needs to be in mb, mixing ratio needs to be in 
% 	kg/kg in some equations, and in g/kg in others.  
% Calculates Iteration via Newton-Raphson Method.  Only 2 iterations.
% Reference:  Davies-Jones: An efficient and accurate method for computing the 
%	      wet-bulb temperature along pseudoadiabats. Monthly Weather Review 
%	      (2008) vol. 136 (7) pp. 2764-2785
% 	      Flatau et al: Polynomial fits to saturation vapor pressure. 
%	      Journal of Applied Meteorology (1992) vol. 31 pp. 1507-1513
% Note: Pressure needs to be in mb, mixing ratio needs to be in 
% 	kg/kg in some equations. 
%
% Ported from HumanIndexMod 04-08-16 by Jonathan R Buzan
% MATLAB port by Robert Kopp
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Wed Jun 08 18:03:02 EDT 2016

    SHR_CONST_TKFRZ = 273.15;
    TemperatureK = TemperatureC + SHR_CONST_TKFRZ;

    constA = 2675; 	% Constant used for extreme cold temparatures (K)
    grms = 1000; 	% Gram per Kilogram (g/kg)
    p0 = 1000;   	% surface pressure (mb)

    kappad = 0.2854;	% Heat Capacity

    C = SHR_CONST_TKFRZ;		% Freezing Temperature
    pmb = Pressure*0.01;		% pa to mb
    T1 = TemperatureK;			% Use holder for T

    if ~exist('ConvergenceMode','var'); ConvergenceMode=1; end
    if ~exist('HumidityMode','var'); HumidityMode=0; end

    [es_mb,rs]=QSat_2(TemperatureK, Pressure);
   
    if HumidityMode==0
        qin=Humidity; % specific humidity
        relhum = 100*qin./rs; % relative humidity (%)
        vapemb = es_mb .* relhum * 0.01; % vapor pressure (mb) 
    elseif HumidityMode==1
        relhum=Humidity;  % relative humidity (%)
        qin = rs .* relhum * 0.01; % specific humidity
        vapemb = es_mb .* relhum * 0.01; % vapor pressure (mb) 
    end

    mixr = qin * grms; % change specific humidity to mixing ratio (g/kg)

    %    real(r8) :: k1;		        % Quadratic Parameter (C)
    %    real(r8) :: k2;		 	% Quadratic Parameter scaled by X (C) 
    %    real(r8) :: pmb;		 	% Atmospheric Surface Pressure (mb)
    %    real(r8) :: D;		 	% Linear Interpolation of X

    %   real(r8) :: hot	% Dimensionless Quantity used for changing temperature regimes
    %    real(r8) :: cold	% Dimensionless Quantity used for changing temperature regimes    

    %   real(r8) :: T1     	 		% Temperature (K)
    %   real(r8) :: vapemb        		% Vapour Pressure (mb)
    %   real(r8) :: mixr        		% Mixing Ratio (g/kg)

    %   real(r8) :: es_mb_teq		% saturated vapour pressure for wrt TEQ (mb)
    %   real(r8) :: de_mbdTeq		% Derivative of Saturated Vapour pressure wrt TEQ (mb/K)
    %   real(r8) :: dlnes_mbdTeq		% Log derivative of the sat. vap pressure wrt TEQ (mb/K)
    %   real(r8) :: rs_teq			% Mixing Ratio wrt TEQ (kg/kg)
    %   real(r8) :: rsdTeq			% Derivative of Mixing Ratio wrt TEQ (kg/kg/K)
    %   real(r8) :: foftk_teq		% Function of EPT wrt TEQ 
    %   real(r8) :: fdTeq			% Derivative of Function of EPT wrt TEQ 

    %   real(r8) :: wb_temp			% Wet Bulb Temperature First Guess (C)
    %   real(r8) :: es_mb_wb_temp		% Vapour Pressure wrt Wet Bulb Temp (mb)
    %   real(r8) :: de_mbdwb_temp		% Derivative of Sat. Vapour Pressure wrt WB Temp (mb/K)
    %   real(r8) :: dlnes_mbdwb_temp	% Log Derivative of sat. vap. pressure wrt WB Temp (mb/K)
    %   real(r8) :: rs_wb_temp		% Mixing Ratio wrt WB Temp (kg/kg)
    %   real(r8) :: rsdwb_temp		% Derivative of Mixing Ratio wrt WB Temp (kg/kg/K)
    %   real(r8) :: foftk_wb_temp		% Function of EPT wrt WB Temp
    %   real(r8) :: fdwb_temp		% Derivative of function of EPT wrt WB Temp

    %   real(r8) :: tl		 	% Lifting Condensation Temperature (K)
    %   real(r8) :: theta_dl	 	% Moist Potential Temperature (K)
    %   real(r8) :: pnd		 	% Non dimensional Pressure
    %   real(r8) :: X		 	% Ratio of equivalent temperature to freezing scaled by Heat Capacity

    %
    %-----------------------------------------------------------------------

    % Calculate Equivalent Pot. Temp (pmb, T, mixing ratio (g/kg), pott, epott)	
    % Calculate Parameters for Wet Bulb Temp (epott, pmb)
    pnd = (pmb./p0).^(kappad);
    D = 1./(0.1859.*pmb./p0 + 0.6512);
    k1 = -38.5.*pnd.*pnd +137.81.*pnd -53.737;
    k2 = -4.392.*pnd.*pnd +56.831.*pnd -0.384;

    % Calculate lifting condensation level.  first eqn 
    % uses vapor pressure (mb)
    % 2nd eqn uses relative humidity.  
    % first equation: Bolton 1980 Eqn 21.
    %   tl = (2840/(3.5*log(T1) - log(vapemb) - 4.805)) + 55;
    % second equation: Bolton 1980 Eqn 22.  relhum = relative humidity
    tl = (1./((1./((T1 - 55))) - (log(relhum/100)/2840))) + 55;

    % Theta_DL: Bolton 1980 Eqn 24.
    theta_dl = T1.*((p0./(pmb-vapemb)).^kappad).*((T1./tl).^(mixr*0.00028));
    % EPT: Bolton 1980 Eqn 39.  
    epott = theta_dl.*exp(((3.036./tl)-0.00178).*mixr.*(1 + 0.000448*mixr));
    Teq = epott.*pnd;			% Equivalent Temperature at pressure
    X = (C./Teq).^3.504;

    % Calculates the regime requirements of wet bulb equations.
    invalid = (Teq > 600) + (Teq < 200);
    hot = (Teq > 355.15);
    cold = ((X>=1).*(X<=D));
    X(invalid==1)=NaN; Teq(invalid==1)=NaN;


    % Calculate Wet Bulb Temperature, initial guess
    % Extremely cold regime if X.gt.D then need to 
    % calculate dlnesTeqdTeq 

    [es_mb_teq,rs_teq,de_mbdTeq, dlnes_mbdTeq, rsdTeq, foftk_teq, fdTeq]=QSat_2(Teq, Pressure);
    wb_temp = Teq - C - ((constA*rs_teq)./(1 + (constA.*rs_teq.*dlnes_mbdTeq)));
    sub=find(X<=D);
    wb_temp(sub) = (k1(sub) - 1.21 * cold(sub) - 1.45 * hot(sub) - (k2(sub) - 1.21 * cold(sub)) .* X(sub) + (0.58 ./ X(sub)) .* hot(sub));
    wb_temp(invalid==1) = NaN;

    % Newton-Raphson Method

    maxiter=3;
    iter=0;
    delta=1e6;
    while (max(delta(:))>.01)&&(iter<=maxiter)
        [es_mb_wb_temp,rs_wb_temp,de_mbdwb_temp, dlnes_mbdwb_temp, rsdwb_temp, foftk_wb_temp, fdwb_temp]=QSat_2(wb_temp+C, Pressure);
        delta=real((foftk_wb_temp - X)./fdwb_temp);
        delta=min(10,delta); delta=max(-10,delta);
        wb_temp = wb_temp - delta;
        wb_temp(invalid==1) = NaN;
        Twb = wb_temp;
        iter=iter+1;
    end
    
    % ! 04-06-16: Adding iteration constraint.  Commenting out original code.
    % but in the MATLAB code, for sake of speed, we only do this for the values
    % that didn't converge

    if ConvergenceMode
        
        convergence = 0.00001; maxiter=20000;

        [es_mb_wb_temp,rs_wb_temp,de_mbdwb_temp, dlnes_mbdwb_temp, rsdwb_temp, foftk_wb_temp, fdwb_temp]=QSat_2(wb_temp+C, Pressure);
        delta=real((foftk_wb_temp - X)./fdwb_temp);
        subdo=find(abs(delta)>convergence);

        iter = 0;
        while (length(subdo)>0)&&(iter<=maxiter)
            iter = iter + 1;
            
            wb_temp(subdo) = wb_temp(subdo) - 0.1*delta(subdo);

            [es_mb_wb_temp,rs_wb_temp,de_mbdwb_temp, dlnes_mbdwb_temp, rsdwb_temp, foftk_wb_temp, fdwb_temp]=QSat_2(wb_temp(subdo)+C, Pressure(subdo));
            delta=0*wb_temp;
            delta(subdo)=real((foftk_wb_temp - X(subdo))./fdwb_temp);
            subdo=find(abs(delta)>convergence);
        end

        Twb=wb_temp;
        if length(subdo)>0
            Twb(subdo) = TemperatureK(subdo)-C;
            for www=subdo(:)'
                disp(sprintf('WARNING-Wet_Bulb failed to converge. Setting to T: WB, P, T, RH, Delta: %0.2f, %0.0f, %0.2f, %0.1f, %0.2g, %0.1f, %0.2g', [Twb(www), Pressure(www), TemperatureK(www), relhum(www), ...
                                    delta(www)]));
            end
        end

    end
    
    Twb=real(Twb);
    
end

function [es_mb,rs,de_mbdT,dlnes_mbdT,rsdT,foftk,fdT]=QSat_2(T_k, p_t)

% [es_mb,rs,de_mbdT,dlnes_mbdT,rsdT,foftk,fdt]=QSat_2(T_k, p_t)
%
% DESCRIPTION:
%  Computes saturation mixing ratio and the change in saturation
%  mixing ratio with respect to temperature.  Uses Bolton eqn 10, 39.
%  Davies-Jones eqns 2.3,A.1-A.10
%  Reference:  Bolton: The computation of equivalent potential temperature. 
%  	      Monthly weather review (1980) vol. 108 (7) pp. 1046-1053
% 	      Davies-Jones: An efficient and accurate method for computing the 
% 	      wet-bulb temperature along pseudoadiabats. Monthly Weather Review 
% 	      (2008) vol. 136 (7) pp. 2764-2785
% 
% INPUTS:
%
%     T_k        temperature (K)
%     p_t        surface atmospheric pressure (pa)
%
% T_k and p_t should be arrays of identical dimensions.
%
% OUTPUTS:
%
%     es_mb      vapor pressure (pa)
%     rs       	 humidity (kg/kg)
%     de_mbdT    d(es)/d(T)
%     dlnes_mbdT dln(es)/d(T)
%     rsdT     	 d(qs)/d(T)
%     foftk      Davies-Jones eqn 2.3
%     fdT     	 d(f)/d(T)
%
% Ported from HumanIndexMod by Jonathan R Buzan 08/08/13
% MATLAB port by Robert Kopp
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Wed Sep 02 22:22:25 EDT 2015
%

    SHR_CONST_TKFRZ = 273.15;

    lambd_a = 3.504;	% Inverse of Heat Capacity
    alpha = 17.67;	% Constant to calculate vapour pressure
    beta = 243.5;		% Constant to calculate vapour pressure
    epsilon = 0.6220;	% Conversion between pressure/mixing ratio
    es_C = 6.112;		% Vapour Pressure at Freezing STD (mb)
    vkp = 0.2854;		% Heat Capacity
    y0 = 3036;		% constant
    y1 = 1.78;		% constant
    y2 = 0.448;		% constant
    Cf = SHR_CONST_TKFRZ;	% Freezing Temp (K)
    refpres = 1000;	% Reference Pressure (mb)

% $$$  p_tmb			% Pressure (mb)
% $$$  ndimpress		% Non-dimensional Pressure
% $$$  prersdt			% Place Holder for derivative humidity
% $$$  pminuse			% Vapor Pressure Difference (mb)
% $$$  tcfbdiff		% Temp diff ref (C)
% $$$  p0ndplam		% dimensionless pressure modified by ref pressure
% $$$     
% $$$  rsy2rs2			% Constant function of humidity
% $$$  oty2rs			% Constant function of humidity
% $$$  y0tky1			% Constant function of Temp
% $$$ 
% $$$  d2e_mbdT2		% d2(es)/d(T)2
% $$$  d2rsdT2			% d2(r)/d(T)2
% $$$  goftk			% g(T) exponential in f(T)
% $$$  gdT			% d(g)/d(T)
% $$$  d2gdT2			% d2(g)/d(T)2
% $$$ 
% $$$  d2fdT2			% d2(f)/d(T)2  (K)
%
%-----------------------------------------------------------------------
% Constants used to calculate es(T)
% Clausius-Clapeyron
    p_tmb = p_t*0.01;
    tcfbdiff = T_k - Cf + beta;
    es_mb = es_C.*exp(alpha.*(T_k - Cf)./(tcfbdiff));
    dlnes_mbdT = alpha.*beta./((tcfbdiff).*(tcfbdiff));
    pminuse = p_tmb - es_mb;
    de_mbdT = es_mb.*dlnes_mbdT;
    d2e_mbdT2 = dlnes_mbdT.*(de_mbdT - 2.*es_mb./(tcfbdiff));

    % Constants used to calculate rs(T)
    ndimpress = (p_tmb./refpres).^vkp;
    p0ndplam = refpres.*ndimpress.^lambd_a;
    rs = epsilon.*es_mb./(p0ndplam - es_mb + eps);
    prersdt = epsilon.*p_tmb./((pminuse).*(pminuse));
    rsdT = prersdt.*de_mbdT;
    d2rsdT2 = prersdt.*(d2e_mbdT2 -de_mbdT.*de_mbdT.*(2./(pminuse)));

    % Constants used to calculate g(T)
    rsy2rs2 = rs + y2.*rs.*rs;
    oty2rs = 1 + 2.*y2.*rs;
    y0tky1 = y0./T_k - y1;
    goftk = y0tky1.*(rs + y2.*rs.*rs);
    gdT = - y0.*(rsy2rs2)./(T_k.*T_k) + (y0tky1).*(oty2rs).*rsdT;
    d2gdT2 = 2.*y0.*rsy2rs2./(T_k.*T_k.*T_k) - 2.*y0.*rsy2rs2.*(oty2rs).*rsdT + ...
             y0tky1.*2.*y2.*rsdT.*rsdT + y0tky1.*oty2rs.*d2rsdT2;

    % Calculations for used to calculate f(T,ndimpress)
    foftk = ((Cf./T_k).^lambd_a).*(1 - es_mb./p0ndplam).^(vkp.*lambd_a).* ...
            exp(-lambd_a.*goftk);
    fdT = -lambd_a.*(1./T_k + vkp.*de_mbdT./pminuse + gdT);
    d2fdT2 = lambd_a.*(1./(T_k.*T_k) - vkp.*de_mbdT.*de_mbdT./(pminuse.*pminuse) - ...
                       vkp.*d2e_mbdT2./pminuse - d2gdT2);

    % avoid bad numbers
    rs(rs>1)=NaN;
    rs(rs<0)=NaN;

end
