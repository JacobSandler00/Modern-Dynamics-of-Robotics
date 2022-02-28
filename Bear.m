function [t,VAR,Output] = Bear
%===========================================================================
% File: Bear.m created Nov 02 2021 by MotionGenesis 6.1.
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
Fx=0; qADDt=0; qCDDt=0; xDDt=0; TBC=0;


%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
g                               =  9.8;                    % m/s^2               Constant
I                               =  3.4;                    % kg*m^2              Constant
J                               =  3.2;                    % kg*m^2              Constant
K                               =  2.8;                    % kg*m^2              Constant
Lc                              = -35;                     % cm                  Constant
m                               =  2;                      % kg                  Constant
r                               =  30;                     % cm                  Constant
xDesired                        =  10;                     % m                   Constant

qA                              =  10;                     % deg                 Initial Value
qC                              =  0;                      % deg                 Initial Value
x                               =  0;                      % m                   Initial Value
qADt                            =  0;                      % rad/sec             Initial Value
qCDt                            =  0;                      % rad/sec             Initial Value
xDt                             =  0;                      % m/s                 Initial Value

tInitial                        =  0.0;                    % second              Initial Time
tFinal                          =  12;                     % sec                 Final Time
tStep                           =  0.05;                   % sec                 Integration Step
printIntScreen                  =  1;                      % 0 or +integer       0 is NO screen output
printIntFile                    =  1;                      % 0 or +integer       0 is NO file   output
absError                        =  1.0E-07;                %                     Absolute Error
relError                        =  1.0E-08;                %                     Relative Error
%-------------------------------+--------------------------+-------------------+-----------------

% Unit conversions.  UnitSystem: kg, meter, second.
DEGtoRAD = pi / 180.0;
RADtoDEG = 180.0 / pi;
Lc = Lc * 0.01;                                            %  Converted from cm 
r = r * 0.01;                                              %  Converted from cm 
qA = qA * DEGtoRAD;                                        %  Converted from deg 
qC = qC * DEGtoRAD;                                        %  Converted from deg 

VAR = SetMatrixFromNamedQuantities;
[t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, tStep, absError, relError, VAR, printIntScreen, printIntFile );
OutputToScreenOrFile( [], 0, 0 );   % Close output files.
if( printIntFile ~= 0 ),  PlotOutputFiles;  end


%===========================================================================
function sys = mdlDerivatives( t, VAR, uSimulink )
%===========================================================================
SetNamedQuantitiesFromMatrix( VAR );
% Quantities previously specified in MotionGenesis.
TBC = 0.3*xDesired - 0.3*x - 0.6*xDt;

COEF = zeros( 4, 4 );
COEF(1,1) = -K*sin(qC)^2 - (I+m*Lc^2)*cos(qC)^2 - m*r*(r+2*Lc*cos(qC));
COEF(2,4) = r;
COEF(3,2) = J + m*Lc^2;
COEF(3,3) = m*Lc*cos(qC);
COEF(4,2) = m*Lc*cos(qC);
COEF(4,3) = m;
COEF(4,4) = -1;
RHS = zeros( 1, 4 );
RHS(1) = -2*m*Lc^2*sin(qC)*cos(qC)*qADt*qCDt - 2*(I-K)*sin(qC)*cos(qC)*qADt*qCDt - m*(g*sin(qA)*(r+Lc*cos(qC))+2*Lc*r*sin(qC)*qADt*  ...
qCDt);
RHS(2) = TBC;
RHS(3) = m*g*Lc*sin(qC)*cos(qA) - TBC - m*Lc*r*sin(qC)*qADt^2 - m*Lc^2*sin(qC)*cos(qC)*qADt^2 - (I-K)*sin(qC)*cos(qC)*qADt^2;
RHS(4) = m*Lc*sin(qC)*qCDt^2;
SolutionToAlgebraicEquations = COEF \ transpose(RHS);

% Update variables after uncoupling equations
qADDt = SolutionToAlgebraicEquations(1);
qCDDt = SolutionToAlgebraicEquations(2);
xDDt = SolutionToAlgebraicEquations(3);
Fx = SolutionToAlgebraicEquations(4);

sys = transpose( SetMatrixOfDerivativesPriorToIntegrationStep );
end



%===========================================================================
function VAR = SetMatrixFromNamedQuantities
%===========================================================================
VAR = zeros( 1, 6 );
VAR(1) = qA;
VAR(2) = qC;
VAR(3) = x;
VAR(4) = qADt;
VAR(5) = qCDt;
VAR(6) = xDt;
end


%===========================================================================
function SetNamedQuantitiesFromMatrix( VAR )
%===========================================================================
qA = VAR(1);
qC = VAR(2);
x = VAR(3);
qADt = VAR(4);
qCDt = VAR(5);
xDt = VAR(6);
end


%===========================================================================
function VARp = SetMatrixOfDerivativesPriorToIntegrationStep
%===========================================================================
VARp = zeros( 1, 6 );
VARp(1) = qADt;
VARp(2) = qCDt;
VARp(3) = xDt;
VARp(4) = qADDt;
VARp(5) = qCDDt;
VARp(6) = xDDt;
end



%===========================================================================
function Output = mdlOutputs( t, VAR, uSimulink )
%===========================================================================
Output = zeros( 1, 10 );
Output(1) = t;
Output(2) = x;
Output(3) = qA*RADtoDEG;                              % Converted to deg
Output(4) = qC*RADtoDEG;                              % Converted to deg

Output(5) = t;
Output(6) = x;

Output(7) = t;
Output(8) = qA*RADtoDEG;

Output(9) = t;
Output(10) = qC*RADtoDEG;
end


%===========================================================================
function OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile )
%===========================================================================
persistent FileIdentifier hasHeaderInformationBeenWritten;

if( isempty(Output) ),
   if( ~isempty(FileIdentifier) ),
      for( i = 1 : 4 ),  fclose( FileIdentifier(i) );  end
      clear FileIdentifier;
      fprintf( 1, '\n Output is in the files Bear.i  (i=1,2,3,4)\n\n' );
   end
   clear hasHeaderInformationBeenWritten;
   return;
end

if( isempty(hasHeaderInformationBeenWritten) ),
   if( shouldPrintToScreen ),
      fprintf( 1,                '%%       t              x             qA             qC\n' );
      fprintf( 1,                '%%     (sec)           (m)           (deg)          (deg)\n\n' );
   end
   if( shouldPrintToFile && isempty(FileIdentifier) ),
      FileIdentifier = zeros( 1, 4 );
      FileIdentifier(1) = fopen('Bear.1', 'wt');   if( FileIdentifier(1) == -1 ), error('Error: unable to open file Bear.1'); end
      fprintf(FileIdentifier(1), '%% FILE: Bear.1\n%%\n' );
      fprintf(FileIdentifier(1), '%%       t              x             qA             qC\n' );
      fprintf(FileIdentifier(1), '%%     (sec)           (m)           (deg)          (deg)\n\n' );
      FileIdentifier(2) = fopen('Bear.2', 'wt');   if( FileIdentifier(2) == -1 ), error('Error: unable to open file Bear.2'); end
      fprintf(FileIdentifier(2), '%% FILE: Bear.2\n%%\n' );
      fprintf(FileIdentifier(2), '%%       t              x\n' );
      fprintf(FileIdentifier(2), '%%     (sec)           (m)\n\n' );
      FileIdentifier(3) = fopen('Bear.3', 'wt');   if( FileIdentifier(3) == -1 ), error('Error: unable to open file Bear.3'); end
      fprintf(FileIdentifier(3), '%% FILE: Bear.3\n%%\n' );
      fprintf(FileIdentifier(3), '%%       t             qa\n' );
      fprintf(FileIdentifier(3), '%%     (sec)          (deg)\n\n' );
      FileIdentifier(4) = fopen('Bear.4', 'wt');   if( FileIdentifier(4) == -1 ), error('Error: unable to open file Bear.4'); end
      fprintf(FileIdentifier(4), '%% FILE: Bear.4\n%%\n' );
      fprintf(FileIdentifier(4), '%%       t             qc\n' );
      fprintf(FileIdentifier(4), '%%     (sec)          (deg)\n\n' );
   end
   hasHeaderInformationBeenWritten = 1;
end

if( shouldPrintToScreen ), WriteNumericalData( 1,                 Output(1:4) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(1), Output(1:4) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(2), Output(5:6) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(3), Output(7:8) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(4), Output(9:10) );  end
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
data = load( 'Bear.2' ); 
plot( data(:,1),data(:,2),'-b', 'LineWidth',3 );
legend( 'x (m)' );
xlabel('t (sec)');   ylabel('x (m)');   % title('Some plot title');
clear data;

figure;
data = load( 'Bear.3' ); 
plot( data(:,1),data(:,2),'-b', 'LineWidth',3 );
legend( 'qa (deg)' );
xlabel('t (sec)');   ylabel('qa (deg)');   % title('Some plot title');
clear data;

figure;
data = load( 'Bear.4' ); 
plot( data(:,1),data(:,2),'-b', 'LineWidth',3 );
legend( 'qc (deg)' );
xlabel('t (sec)');   ylabel('qc (deg)');   % title('Some plot title');
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


%=============================
end    % End of function Bear
%=============================
