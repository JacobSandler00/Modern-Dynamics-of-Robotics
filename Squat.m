function [t,VAR,Output] = Squat
%===========================================================================
% File: Squat.m created Nov 30 2021 by MotionGenesis 6.1.
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
TA=0; TB=0; TC=0; qADt=0; qBDt=0; qADDt=0; qBDDt=0; qC=0; qCDt=0; yDt=0; qCDDt=0; yDDt=0; xQ=0; yQ=0;


%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
g                               =  9.8;                    % m/s^2               Constant
IzzA                            =  .3;                     % kg                  Constant
IzzB                            =  .3;                     % kg                  Constant
IzzC                            =  .3;                     % kg                  Constant
LA                              =  0.5;                    % m                   Constant
LB                              =  0.5;                    % m                   Constant
LC                              =  0.5;                    % m                   Constant
LX                              =  0.1;                    % m                   Constant
LY                              =  0.01;                   % m                   Constant
mA                              =  .3;                     % kg                  Constant
mB                              =  .4;                     % kg                  Constant
mC                              =  .37;                    % kg                  Constant
mQ                              =  .1;                     % kg                  Constant
yFinal                          =  1.1;                    % m                   Constant

qA                              =  0.5235987755982988;     % rad                 Initial Value
qB                              =  0;                      % deg                 Initial Value

tInitial                        =  0;                      % s                   Initial Time
tFinal                          =  10;                     % s                   Final Time
tStep                           =  0.02;                   % s                   Integration Step
printIntScreen                  =  1;                      % 0 or +integer       0 is NO screen output
printIntFile                    =  1;                      % 0 or +integer       0 is NO file   output
absError                        =  1.0E-05;                %                     Absolute Error
relError                        =  1.0E-08;                %                     Relative Error
%-------------------------------+--------------------------+-------------------+-----------------

% Unit conversions
DEGtoRAD = pi / 180.0;
qB = qB * DEGtoRAD;

VAR = SetMatrixFromNamedQuantities;
[t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, tStep, absError, relError, VAR, printIntScreen, printIntFile );
OutputToScreenOrFile( [], 0, 0 );   % Close output files.
if( printIntFile ~= 0 ),  PlotOutputFiles;  end


%===========================================================================
function sys = mdlDerivatives( t, VAR, uSimulink )
%===========================================================================
SetNamedQuantitiesFromMatrix( VAR );
% Quantities previously specified in MotionGenesis.
yDt = -0.04413397459621557 + 0.1*yFinal + 0.04413397459621557*(1-2.265828104423092*yFinal)*cos(0.6283185307179586*t);
yDDt = -0.02773019407303788*(1-2.265828104423092*yFinal)*sin(0.6283185307179586*t);
qC = 2.617993877991494;
qCDt = 0;
qCDDt = 0;

qADt = -((LY*cos(qC)+(LC-LX)*sin(qC))*qCDt*cos(qB)-(yDt+LY*sin(qC)*qCDt-(LC-LX)*cos(qC)*qCDt)*sin(qB))/(LA*sin(qA+qB));
qBDt = ((LY*cos(qC)+(LC-LX)*sin(qC))*qCDt*cos(qA)+(yDt+LY*sin(qC)*qCDt-(LC-LX)*cos(qC)*qCDt)*sin(qA))/(LB*sin(qA+qB));

sys = transpose( SetMatrixOfDerivativesPriorToIntegrationStep );
end



%===========================================================================
function VAR = SetMatrixFromNamedQuantities
%===========================================================================
VAR = zeros( 1, 2 );
VAR(1) = qA;
VAR(2) = qB;
end


%===========================================================================
function SetNamedQuantitiesFromMatrix( VAR )
%===========================================================================
qA = VAR(1);
qB = VAR(2);
end


%===========================================================================
function VARp = SetMatrixOfDerivativesPriorToIntegrationStep
%===========================================================================
VARp = zeros( 1, 2 );
VARp(1) = qADt;
VARp(2) = qBDt;
end



%===========================================================================
function Output = mdlOutputs( t, VAR, uSimulink )
%===========================================================================
qADDt = (LB*qBDt^2+cos(qB)*((LY*sin(qC)-(LC-LX)*cos(qC))*qCDt^2-(LY*cos(qC)+(LC-LX)*sin(qC))*qCDDt-LA*cos(qA)*qADt^2)-sin(qB)*((LC-  ...
LX)*cos(qC)*qCDDt-LY*cos(qC)*qCDt^2-(LC-LX)*sin(qC)*qCDt^2-yDDt-LY*sin(qC)*qCDDt-LA*sin(qA)*qADt^2))/(LA*sin(qA+qB));
qBDDt = (LA*qADt^2+cos(qA)*((LY*cos(qC)+(LC-LX)*sin(qC))*qCDDt-(LY*sin(qC)-(LC-LX)*cos(qC))*qCDt^2-LB*cos(qB)*qBDt^2)-sin(qA)*((LC-  ...
LX)*cos(qC)*qCDDt-LY*cos(qC)*qCDt^2-(LC-LX)*sin(qC)*qCDt^2-yDDt-LY*sin(qC)*qCDDt-LB*sin(qB)*qBDt^2))/(LB*sin(qA+qB));
TA = 0.5*g*LC*mC*cos(qC) + g*LA*mB*cos(qA) + g*LA*mC*cos(qA) + g*LA*mQ*cos(qA) + 0.5*g*LA*mA*cos(qA) + 0.5*LB*(LA*mB*sin(qA+qB)*qADt^2+  ...
mC*(LC*sin(qC+qB)*qCDt^2+2*LA*sin(qA+qB)*qADt^2)+2*mQ*(LY*cos(qC+qB)*qCDt^2+(LC-LX)*sin(qC+qB)*qCDt^2+LA*sin(qA+qB)*qADt^2))  ...
+ 0.25*(4*IzzC+mC*LC^2+4*mQ*(LY^2+(LC-LX)^2))*qCDDt + 0.5*LA*(LC*mC*cos(qC-qA)-2*mQ*(LY*sin(qC-qA)-(LC-LX)*cos(qC-qA)))*qCDDt  ...
+ 0.5*LA*LB*(mB+2*mC+2*mQ)*cos(qA+qB)*qBDDt + 0.25*(4*IzzA+mA*LA^2+4*mB*LA^2+4*mC*LA^2+4*mQ*LA^2)*qADDt + 0.5*LB*(LC*mC*cos(qC+qB)-  ...
2*mQ*(LY*sin(qC+qB)-(LC-LX)*cos(qC+qB)))*qBDDt + 0.5*LA*(LC*mC*cos(qC-qA)-2*mQ*(LY*sin(qC-qA)-(LC-LX)*cos(qC-qA)))*qADDt  ...
- g*mQ*(LY*sin(qC)-(LC-LX)*cos(qC)) - g*LB*mC*cos(qB) - g*LB*mQ*cos(qB) - 0.5*g*LB*mB*cos(qB) - 0.5*LC*mC*(LB*sin(qC+qB)*qBDt^2-LA*  ...
sin(qC-qA)*qADt^2) - mQ*(LB*LY*cos(qC+qB)*qBDt^2+LB*(LC-LX)*sin(qC+qB)*qBDt^2-LA*LY*cos(qC-qA)*qADt^2-LA*(LC-LX)*sin(qC-qA)*qADt^2)  ...
- 0.5*LA*(LB*mB*sin(qA+qB)*qBDt^2+mC*(LC*sin(qC-qA)*qCDt^2+2*LB*sin(qA+qB)*qBDt^2)+2*mQ*(LY*cos(qC-qA)*qCDt^2+(LC-LX)*sin(qC-qA)*qCDt^2+  ...
LB*sin(qA+qB)*qBDt^2)) - 0.5*LB*(LC*mC*cos(qC+qB)-2*mQ*(LY*sin(qC+qB)-(LC-LX)*cos(qC+qB)))*qCDDt - 0.5*LA*LB*(mB+2*mC+2*mQ)*cos(qA+  ...
qB)*qADDt - 0.25*(4*IzzB+mB*LB^2+4*mC*LB^2+4*mQ*LB^2)*qBDDt;
TB = g*mQ*(LY*sin(qC)-(LC-LX)*cos(qC)) + g*LB*mC*cos(qB) + g*LB*mQ*cos(qB) + 0.5*g*LB*mB*cos(qB) + 0.5*LC*mC*(LB*sin(qC+qB)*qBDt^2-  ...
LA*sin(qC-qA)*qADt^2) + mQ*(LB*LY*cos(qC+qB)*qBDt^2+LB*(LC-LX)*sin(qC+qB)*qBDt^2-LA*LY*cos(qC-qA)*qADt^2-LA*(LC-LX)*sin(qC-qA)*qADt^2)  ...
+ 0.5*LB*(LC*mC*cos(qC+qB)-2*mQ*(LY*sin(qC+qB)-(LC-LX)*cos(qC+qB)))*qCDDt + 0.5*LA*LB*(mB+2*mC+2*mQ)*cos(qA+qB)*qADDt + 0.25*(4*IzzB+  ...
mB*LB^2+4*mC*LB^2+4*mQ*LB^2)*qBDDt - 0.5*g*LC*mC*cos(qC) - 0.5*LB*(LA*mB*sin(qA+qB)*qADt^2+mC*(LC*sin(qC+qB)*qCDt^2+2*LA*sin(qA+qB)*qADt^2)+  ...
2*mQ*(LY*cos(qC+qB)*qCDt^2+(LC-LX)*sin(qC+qB)*qCDt^2+LA*sin(qA+qB)*qADt^2)) - 0.25*(4*IzzC+mC*LC^2+4*mQ*(LY^2+(LC-LX)^2))*qCDDt  ...
- 0.5*LB*(LC*mC*cos(qC+qB)-2*mQ*(LY*sin(qC+qB)-(LC-LX)*cos(qC+qB)))*qBDDt - 0.5*LA*(LC*mC*cos(qC-qA)-2*mQ*(LY*sin(qC-qA)-(LC-LX)*  ...
cos(qC-qA)))*qADDt;
TC = 0.5*g*LC*mC*cos(qC) + 0.25*(4*IzzC+mC*LC^2+4*mQ*(LY^2+(LC-LX)^2))*qCDDt + 0.5*LB*(LC*mC*cos(qC+qB)-2*mQ*(LY*sin(qC+qB)-(LC-LX)*  ...
cos(qC+qB)))*qBDDt + 0.5*LA*(LC*mC*cos(qC-qA)-2*mQ*(LY*sin(qC-qA)-(LC-LX)*cos(qC-qA)))*qADDt - g*mQ*(LY*sin(qC)-(LC-LX)*cos(qC))  ...
- 0.5*LC*mC*(LB*sin(qC+qB)*qBDt^2-LA*sin(qC-qA)*qADt^2) - mQ*(LB*LY*cos(qC+qB)*qBDt^2+LB*(LC-LX)*sin(qC+qB)*qBDt^2-LA*LY*cos(qC-qA)*qADt^2-  ...
LA*(LC-LX)*sin(qC-qA)*qADt^2);
xQ = (LC-LX)*cos(qC) + LA*cos(qA) - LY*sin(qC) - LB*cos(qB);
yQ = LY*cos(qC) + (LC-LX)*sin(qC) + LA*sin(qA) + LB*sin(qB);

Output = zeros( 1, 10 );
Output(1) = xQ;
Output(2) = yQ;

Output(3) = t;
Output(4) = qA;
Output(5) = qB;
Output(6) = qC;

Output(7) = t;
Output(8) = TA;
Output(9) = TB;
Output(10) = TC;
end


%===========================================================================
function OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile )
%===========================================================================
persistent FileIdentifier hasHeaderInformationBeenWritten;

if( isempty(Output) ),
   if( ~isempty(FileIdentifier) ),
      for( i = 1 : 3 ),  fclose( FileIdentifier(i) );  end
      clear FileIdentifier;
      fprintf( 1, '\n Output is in the files Squat.i  (i=1,2,3)\n\n' );
   end
   clear hasHeaderInformationBeenWritten;
   return;
end

if( isempty(hasHeaderInformationBeenWritten) ),
   if( shouldPrintToScreen ),
      fprintf( 1,                '%%      xQ             yQ\n' );
      fprintf( 1,                '%%      (m)            (m)\n\n' );
   end
   if( shouldPrintToFile && isempty(FileIdentifier) ),
      FileIdentifier = zeros( 1, 3 );
      FileIdentifier(1) = fopen('Squat.1', 'wt');   if( FileIdentifier(1) == -1 ), error('Error: unable to open file Squat.1'); end
      fprintf(FileIdentifier(1), '%% FILE: Squat.1\n%%\n' );
      fprintf(FileIdentifier(1), '%%      xQ             yQ\n' );
      fprintf(FileIdentifier(1), '%%      (m)            (m)\n\n' );
      FileIdentifier(2) = fopen('Squat.2', 'wt');   if( FileIdentifier(2) == -1 ), error('Error: unable to open file Squat.2'); end
      fprintf(FileIdentifier(2), '%% FILE: Squat.2\n%%\n' );
      fprintf(FileIdentifier(2), '%%       t             qA             qB             qc\n' );
      fprintf(FileIdentifier(2), '%%      (s)           (rad)          (rad)          (rad)\n\n' );
      FileIdentifier(3) = fopen('Squat.3', 'wt');   if( FileIdentifier(3) == -1 ), error('Error: unable to open file Squat.3'); end
      fprintf(FileIdentifier(3), '%% FILE: Squat.3\n%%\n' );
      fprintf(FileIdentifier(3), '%%       t             TA             TB             TC\n' );
      fprintf(FileIdentifier(3), '%%      (s)           (N*m)          (N*m)          (N*m)\n\n' );
   end
   hasHeaderInformationBeenWritten = 1;
end

if( shouldPrintToScreen ), WriteNumericalData( 1,                 Output(1:2) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(1), Output(1:2) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(2), Output(3:6) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(3), Output(7:10) );  end
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
data = load( 'Squat.1' ); 
plot( data(:,1),data(:,2),'-b', 'LineWidth',3 );
legend( 'yQ (m)' );
xlabel('xQ (m)');   ylabel('yQ (m)');   % title('Some plot title');
clear data;

figure;
data = load( 'Squat.2' ); 
plot( data(:,1),data(:,2),'-b', data(:,1),data(:,3),'-.g', data(:,1),data(:,4),'--r', 'LineWidth',3 );
legend( 'qA (rad)', 'qB (rad)', 'qc (rad)' );
xlabel('t (s)');   % ylabel('Some y-axis label');   title('Some plot title');
clear data;

figure;
data = load( 'Squat.3' ); 
plot( data(:,1),data(:,2),'-b', data(:,1),data(:,3),'-.g', data(:,1),data(:,4),'--r', 'LineWidth',3 );
legend( 'TA (N*m)', 'TB (N*m)', 'TC (N*m)' );
xlabel('t (s)');   % ylabel('Some y-axis label');   title('Some plot title');
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


%==============================
end    % End of function Squat
%==============================
