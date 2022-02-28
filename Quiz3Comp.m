function [t,VAR,Output] = Quiz3Comp
%===========================================================================
% File: Quiz3Comp.m created Nov 07 2021 by MotionGenesis 6.1.
% Portions copyright (c) 2009-2021 Motion Genesis LLC.  Rights reserved.
% MotionGenesis Student Licensee: Jacob Sandler. (until August 2024).
% Paid-up MotionGenesis Student licensees are granted the right
% to distribute this code for legal student-academic (non-professional) purposes only,
% provided this copyright notice appears in all copies and distributions.
%===========================================================================
% The software is provided "as is", without warranty of any kind, express or    
% implied, including but not limited to the warranties of merchantability or    
% fitness for a particular purpose. In no event shall the authors, contributors,
% or copyright holders be liable for any claim, damages or other liability,     
% whether in an action of contract, tort, or otherwise, arising from, out of, or
% in connection with the software or the use or other dealings in the software. 
%===========================================================================
eventDetectedByIntegratorTerminate1OrContinue0 = [];
IzzA=0; IzzB=0; kd=0; kp=0; qcDt=0; qaDDt=0; qbDDt=0; yDDt=0; ysDDt=0; pen=0; TA=0; penDt=0; FHC=0; yAcm=0;


%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
bs                              =  1.5e4;                  % N*s/m               Constant
g                               =  9.80655;                % m/s^2               Constant
khc                             =  1e7;                    % UNITS               Constant
ks                              =  1e5;                    % N/m                 Constant
L                               =  0.8;                    % UNITS               Constant
LB                              =  .05;                    % m                   Constant
LC                              =  .1;                     % m                   Constant
lhc                             =  5e6;                    % UNITS               Constant
LN                              =  .05;                    % m                   Constant
mA                              =  10;                     % kg                  Constant
mB                              =  4;                      % kg                  Constant
mC                              =  4;                      % kg                  Constant
mD                              =  2;                      % kg                  Constant
mP                              =  1;                      % kg                  Constant
nhc                             =  1.5;                    % UNITS               Constant
qades                           =  0;                      % rad                 Constant
qbdesD                          =  20;                     % rad/s               Constant
rA                              =  .025;                   % m                   Constant
wn                              =  20;                     % UNITS               Constant

qa                              =  0;                      % deg                 Initial Value
qb                              =  90;                     % deg                 Initial Value
qc                              =  90;                     % deg                 Initial Value
y                               =  0.0;                    % m                   Initial Value
ys                              =  0.05;                   % m                   Initial Value
qaDt                            =  0;                      % deg                 Initial Value
qbDt                            =  0;                      % UNITS               Initial Value
yDt                             =  0;                      % UNITS               Initial Value
ysDt                            =  0;                      % UNITS               Initial Value

tInitial                        =  0;                      % s                   Initial Time
tFinal                          =  4;                      % s                   Final Time
tStep                           =  .02;                    % s                   Integration Step
printIntScreen                  =  1;                      % 0 or +integer       0 is NO screen output
printIntFile                    =  1;                      % 0 or +integer       0 is NO file   output
absError                        =  1.0E-05;                %                     Absolute Error
relError                        =  1.0E-08;                %                     Relative Error
%-------------------------------+--------------------------+-------------------+-----------------

% Unit conversions
DEGtoRAD = pi / 180.0;
qa = qa * DEGtoRAD;
qb = qb * DEGtoRAD;
qc = qc * DEGtoRAD;
qaDt = qaDt * DEGtoRAD;

% Evaluate constants
IzzA = 0.5*mA*rA^2;
IzzB = 0.08333333333333333*mB*LB^2;
kp = wn^2;
kd = 2*L*wn;
qcDt = 0;


VAR = SetMatrixFromNamedQuantities;
[t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, tStep, absError, relError, VAR, printIntScreen, printIntFile );
OutputToScreenOrFile( [], 0, 0 );   % Close output files.
if( printIntFile ~= 0 ),  PlotOutputFiles;  end


%===========================================================================
function sys = mdlDerivatives( t, VAR, uSimulink )
%===========================================================================
SetNamedQuantitiesFromMatrix( VAR );
qbDDt = -(2*(mA+mB+mC+mD)*(g*LB*mB*cos(qb)+2*g*LB*mA*cos(qb)-2*kp*(qbdesD-qbDt)-LB*(mB+2*mA)*(LC*sin(qb+qc)*qcDt^2+cos(qb+qc)*(LB*  ...
cos(qb)*qbDt^2+LC*cos(qc)*qcDt^2)/sin(qc)))+LB*(mB+2*mA)*cos(qb)*(2*ks*(LN-ys)+mC*(LC*sin(qc)*qcDt^2+(LB*cos(qb)*qbDt^2+LC*cos(qc)*qcDt^2)/  ...
tan(qc))-2*g*mA-2*g*mB-2*g*mC-2*g*mD-2*bs*ysDt-2*mA*(LB*sin(qb)*qbDt^2-LC*sin(qc)*qcDt^2-(LB*cos(qb)*qbDt^2+LC*cos(qc)*qcDt^2)/tan(  ...
qc))-mB*(LB*sin(qb)*qbDt^2-2*LC*sin(qc)*qcDt^2-2*(LB*cos(qb)*qbDt^2+LC*cos(qc)*qcDt^2)/tan(qc))))/(LB^2*(mB+2*mA)*cos(qb)*(mC*sin(  ...
qb)/tan(qc)+2*mA*(cos(qb)+sin(qb)/tan(qc))+mB*(cos(qb)+2*sin(qb)/tan(qc)))+(mA+mB+mC+mD)*(4*kd-4*IzzB-4*mA*LB^2*(1+sin(qb)*cos(qb+  ...
qc)/sin(qc))-mB*LB^2*(1+2*sin(qb)*cos(qb+qc)/sin(qc))));

% Quantities previously specified in MotionGenesis.
TA = kp*(qades-qa) - kd*qaDt;
pen = -y;
penDt = -yDt;

qaDDt = (TA+kp*(qbdesD-qbDt)+kd*(2*(mA+mB+mC+mD)*(g*LB*mB*cos(qb)+2*g*LB*mA*cos(qb)-2*kp*(qbdesD-qbDt)-LB*(mB+2*mA)*(LC*sin(qb+qc)*qcDt^2+  ...
cos(qb+qc)*(LB*cos(qb)*qbDt^2+LC*cos(qc)*qcDt^2)/sin(qc)))+LB*(mB+2*mA)*cos(qb)*(2*ks*(LN-ys)+mC*(LC*sin(qc)*qcDt^2+(LB*cos(qb)*qbDt^2+  ...
LC*cos(qc)*qcDt^2)/tan(qc))-2*g*mA-2*g*mB-2*g*mC-2*g*mD-2*bs*ysDt-2*mA*(LB*sin(qb)*qbDt^2-LC*sin(qc)*qcDt^2-(LB*cos(qb)*qbDt^2+LC*  ...
cos(qc)*qcDt^2)/tan(qc))-mB*(LB*sin(qb)*qbDt^2-2*LC*sin(qc)*qcDt^2-2*(LB*cos(qb)*qbDt^2+LC*cos(qc)*qcDt^2)/tan(qc))))/(LB^2*(mB+2*  ...
mA)*cos(qb)*(mC*sin(qb)/tan(qc)+2*mA*(cos(qb)+sin(qb)/tan(qc))+mB*(cos(qb)+2*sin(qb)/tan(qc)))+(mA+mB+mC+mD)*(4*kd-4*IzzB-4*mA*LB^2*(  ...
1+sin(qb)*cos(qb+qc)/sin(qc))-mB*LB^2*(1+2*sin(qb)*cos(qb+qc)/sin(qc)))))/IzzA;
if (pen > 0 )
    FHC = khc*pen^nhc + lhc*pen^nhc*penDt;
else
    FHC = 0
end
yDDt = -g - (ks*(LN-ys)-FHC-bs*ysDt)/mP;
ysDDt = g + 0.5*(LB^2*(mB+2*mA)*cos(qb)*(mC*sin(qb)/tan(qc)+2*mA*(cos(qb)+sin(qb)/tan(qc))+mB*(cos(qb)+2*sin(qb)/tan(qc)))+(mA+mB+  ...
mC+mD+mP)*(4*kd-4*IzzB-4*mA*LB^2*(1+sin(qb)*cos(qb+qc)/sin(qc))-mB*LB^2*(1+2*sin(qb)*cos(qb+qc)/sin(qc))))*(2*ks*(LN-ys)+mC*(LC*sin(  ...
qc)*qcDt^2+(LB*cos(qb)*qbDt^2+LC*cos(qc)*qcDt^2)/tan(qc))-2*g*mA-2*g*mB-2*g*mC-2*g*mD-2*bs*ysDt-2*mA*(LB*sin(qb)*qbDt^2-LC*sin(qc)*qcDt^2-(  ...
LB*cos(qb)*qbDt^2+LC*cos(qc)*qcDt^2)/tan(qc))-mB*(LB*sin(qb)*qbDt^2-2*LC*sin(qc)*qcDt^2-2*(LB*cos(qb)*qbDt^2+LC*cos(qc)*qcDt^2)/tan(  ...
qc)))/(mP*(LB^2*(mB+2*mA)*cos(qb)*(mC*sin(qb)/tan(qc)+2*mA*(cos(qb)+sin(qb)/tan(qc))+mB*(cos(qb)+2*sin(qb)/tan(qc)))+(mA+mB+mC+mD)*(  ...
4*kd-4*IzzB-4*mA*LB^2*(1+sin(qb)*cos(qb+qc)/sin(qc))-mB*LB^2*(1+2*sin(qb)*cos(qb+qc)/sin(qc))))) - 0.5*(2*FHC+mC*(LC*sin(qc)*qcDt^2+(  ...
LB*cos(qb)*qbDt^2+LC*cos(qc)*qcDt^2)/tan(qc))-2*g*mA-2*g*mB-2*g*mC-2*g*mD-2*mA*(LB*sin(qb)*qbDt^2-LC*sin(qc)*qcDt^2-(LB*cos(qb)*qbDt^2+  ...
LC*cos(qc)*qcDt^2)/tan(qc))-mB*(LB*sin(qb)*qbDt^2-2*LC*sin(qc)*qcDt^2-2*(LB*cos(qb)*qbDt^2+LC*cos(qc)*qcDt^2)/tan(qc)))/mP  ...
- LB*(mC*sin(qb)/tan(qc)+2*mA*(cos(qb)+sin(qb)/tan(qc))+mB*(cos(qb)+2*sin(qb)/tan(qc)))*(g*LB*mB*cos(qb)+2*g*LB*mA*cos(qb)-2*kp*(  ...
qbdesD-qbDt)-LB*(mB+2*mA)*(LC*sin(qb+qc)*qcDt^2+cos(qb+qc)*(LB*cos(qb)*qbDt^2+LC*cos(qc)*qcDt^2)/sin(qc)))/(LB^2*(mB+2*mA)*cos(qb)*(  ...
mC*sin(qb)/tan(qc)+2*mA*(cos(qb)+sin(qb)/tan(qc))+mB*(cos(qb)+2*sin(qb)/tan(qc)))+(mA+mB+mC+mD)*(4*kd-4*IzzB-4*mA*LB^2*(1+sin(qb)*  ...
cos(qb+qc)/sin(qc))-mB*LB^2*(1+2*sin(qb)*cos(qb+qc)/sin(qc))));

sys = transpose( SetMatrixOfDerivativesPriorToIntegrationStep );
end



%===========================================================================
function VAR = SetMatrixFromNamedQuantities
%===========================================================================
VAR = zeros( 1, 9 );
VAR(1) = qa;
VAR(2) = qb;
VAR(3) = qc;
VAR(4) = y;
VAR(5) = ys;
VAR(6) = qaDt;
VAR(7) = qbDt;
VAR(8) = yDt;
VAR(9) = ysDt;
end


%===========================================================================
function SetNamedQuantitiesFromMatrix( VAR )
%===========================================================================
qa = VAR(1);
qb = VAR(2);
qc = VAR(3);
y = VAR(4);
ys = VAR(5);
qaDt = VAR(6);
qbDt = VAR(7);
yDt = VAR(8);
ysDt = VAR(9);
end


%===========================================================================
function VARp = SetMatrixOfDerivativesPriorToIntegrationStep
%===========================================================================
VARp = zeros( 1, 9 );
VARp(1) = qaDt;
VARp(2) = qbDt;
VARp(3) = qcDt;
VARp(4) = yDt;
VARp(5) = ysDt;
VARp(6) = qaDDt;
VARp(7) = qbDDt;
VARp(8) = yDDt;
VARp(9) = ysDDt;
end



%===========================================================================
function Output = mdlOutputs( t, VAR, uSimulink )
%===========================================================================
yAcm = y + ys + LC*sin(qc) - LB*sin(qb);

Output = zeros( 1, 4 );
Output(1) = t;
Output(2) = y;
Output(3) = ys;
Output(4) = yAcm;
end


%===========================================================================
function OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile )
%===========================================================================
persistent FileIdentifier hasHeaderInformationBeenWritten;

if( isempty(Output) ),
   if( ~isempty(FileIdentifier) ),
      fclose( FileIdentifier(1) );
      clear FileIdentifier;
      fprintf( 1, '\n Output is in the file Quiz3Comp.1\n\n' );
   end
   clear hasHeaderInformationBeenWritten;
   return;
end

if( isempty(hasHeaderInformationBeenWritten) ),
   if( shouldPrintToScreen ),
      fprintf( 1,                '%%       t              y             ys            yacm\n' );
      fprintf( 1,                '%%     (sec)           (m)            (m)            (m)\n\n' );
   end
   if( shouldPrintToFile && isempty(FileIdentifier) ),
      FileIdentifier(1) = fopen('Quiz3Comp.1', 'wt');   if( FileIdentifier(1) == -1 ), error('Error: unable to open file Quiz3Comp.1'); end
      fprintf(FileIdentifier(1), '%% FILE: Quiz3Comp.1\n%%\n' );
      fprintf(FileIdentifier(1), '%%       t              y             ys            yacm\n' );
      fprintf(FileIdentifier(1), '%%     (sec)           (m)            (m)            (m)\n\n' );
   end
   hasHeaderInformationBeenWritten = 1;
end

if( shouldPrintToScreen ), WriteNumericalData( 1,                 Output(1:4) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(1), Output(1:4) );  end
end


%===========================================================================
function WriteNumericalData( fileIdentifier, Output )
%===========================================================================
numberOfOutputQuantities = length( Output );
if( numberOfOutputQuantities > 0 ),
   for( i = 1 : numberOfOutputQuantities ),
      fprintf( fileIdentifier, ' %- 14.6E', Output(i) );
   end
   fprintf( fileIdentifier, '\n' );
end
end



%===========================================================================
function PlotOutputFiles
%===========================================================================
figure;
data = load( 'Quiz3Comp.1' ); 
plot( data(:,1),data(:,2),'-b', data(:,1),data(:,3),'-.g', data(:,1),data(:,4),'--r', 'LineWidth',3 );
legend( 'y (m)', 'ys (m)', 'yacm (m)' );
xlabel('t (sec)');   % ylabel('Some y-axis label');   title('Some plot title');
clear data;
end



%===========================================================================
function [functionsToEvaluateForEvent, eventTerminatesIntegration1Otherwise0ToContinue, eventDirection_AscendingIs1_CrossingIs0_DescendingIsNegative1] = EventDetection( t, VAR, uSimulink )
%===========================================================================
% Detects when designated functions are zero or cross zero with positive or negative slope.
% Step 1: Uncomment call to mdlDerivatives and mdlOutputs.
% Step 2: Change functionsToEvaluateForEvent,                      e.g., change  []  to  [t - 5.67]  to stop at t = 5.67.
% Step 3: Change eventTerminatesIntegration1Otherwise0ToContinue,  e.g., change  []  to  [1]  to stop integrating.
% Step 4: Change eventDirection_AscendingIs1_CrossingIs0_DescendingIsNegative1,  e.g., change  []  to  [1].
% Step 5: Possibly modify function EventDetectedByIntegrator (if eventTerminatesIntegration1Otherwise0ToContinue is 0).
%---------------------------------------------------------------------------
% mdlDerivatives( t, VAR, uSimulink );        % UNCOMMENT FOR EVENT HANDLING
% mdlOutputs(     t, VAR, uSimulink );        % UNCOMMENT FOR EVENT HANDLING
functionsToEvaluateForEvent = [];
eventTerminatesIntegration1Otherwise0ToContinue = [];
eventDirection_AscendingIs1_CrossingIs0_DescendingIsNegative1 = [];
eventDetectedByIntegratorTerminate1OrContinue0 = eventTerminatesIntegration1Otherwise0ToContinue;
end


%===========================================================================
function [isIntegrationFinished, VAR] = EventDetectedByIntegrator( t, VAR, nIndexOfEvents )
%===========================================================================
isIntegrationFinished = eventDetectedByIntegratorTerminate1OrContinue0( nIndexOfEvents );
if( ~isIntegrationFinished ),
   SetNamedQuantitiesFromMatrix( VAR );
%  Put code here to modify how integration continues.
   VAR = SetMatrixFromNamedQuantities;
end
end



%===========================================================================
function [t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, tStep, absError, relError, VAR, printIntScreen, printIntFile )
%===========================================================================
OdeMatlabOptions = odeset( 'RelTol',relError, 'AbsTol',absError, 'MaxStep',tStep, 'Events',@EventDetection );
t = tInitial;                 epsilonT = 0.001*tStep;                   tFinalMinusEpsilonT = tFinal - epsilonT;
printCounterScreen = 0;       integrateForward = tFinal >= tInitial;    tAtEndOfIntegrationStep = t + tStep;
printCounterFile   = 0;       isIntegrationFinished = 0;
mdlDerivatives( t, VAR, 0 );
while 1,
   if( (integrateForward && t >= tFinalMinusEpsilonT) || (~integrateForward && t <= tFinalMinusEpsilonT) ), isIntegrationFinished = 1;  end
   shouldPrintToScreen = printIntScreen && ( isIntegrationFinished || printCounterScreen <= 0.01 );
   shouldPrintToFile   = printIntFile   && ( isIntegrationFinished || printCounterFile   <= 0.01 );
   if( isIntegrationFinished || shouldPrintToScreen || shouldPrintToFile ),
      Output = mdlOutputs( t, VAR, 0 );
      OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile );
      if( isIntegrationFinished ), break;  end
      if( shouldPrintToScreen ), printCounterScreen = printIntScreen;  end
      if( shouldPrintToFile ),   printCounterFile   = printIntFile;    end
   end
   [TimeOdeArray, VarOdeArray, timeEventOccurredInIntegrationStep, nStatesArraysAtEvent, nIndexOfEvents] = ode45( @mdlDerivatives, [t tAtEndOfIntegrationStep], VAR, OdeMatlabOptions, 0 );
   if( isempty(timeEventOccurredInIntegrationStep) ),
      lastIndex = length( TimeOdeArray );
      t = TimeOdeArray( lastIndex );
      VAR = VarOdeArray( lastIndex, : );
      printCounterScreen = printCounterScreen - 1;
      printCounterFile   = printCounterFile   - 1;
      if( abs(tAtEndOfIntegrationStep - t) >= abs(epsilonT) ), warning('numerical integration failed'); break;  end
      tAtEndOfIntegrationStep = t + tStep;
      if( (integrateForward && tAtEndOfIntegrationStep > tFinal) || (~integrateForward && tAtEndOfIntegrationStep < tFinal) ) tAtEndOfIntegrationStep = tFinal;  end
   else
      t = timeEventOccurredInIntegrationStep( 1 );    % time  at firstEvent = 1 during this integration step.
      VAR = nStatesArraysAtEvent( 1, : );             % state at firstEvent = 1 during this integration step.
      printCounterScreen = 0;
      printCounterFile   = 0;
      [isIntegrationFinished, VAR] = EventDetectedByIntegrator( t, VAR, nIndexOfEvents(1) );
   end
end
end


%==================================
end    % End of function Quiz3Comp
%==================================
