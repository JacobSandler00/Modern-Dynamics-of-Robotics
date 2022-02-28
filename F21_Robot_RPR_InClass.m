function [t,VAR,Output] = F21_Robot_RPR_InClass
%===========================================================================
% File: F21_Robot_RPR_InClass.m created Dec 01 2021 by MotionGenesis 6.1.
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
FB=0; TA=0; TC=0; dbDt=0; qaDt=0; qcDt=0; dbDDt=0; qaDDt=0; qcDDt=0; BCx=0; BCz=0; IAYY=0; ICXX=0; ICYY=0; Qx=0; Qy=0; Qz=0;


%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
g                               =  9.81;                   % m/s^2               Constant
LA                              =  1.0;                    % m                   Constant
LC                              =  0.8;                    % m                   Constant
mA                              =  1.0;                    % kg                  Constant
mB                              =  0.5;                    % kg                  Constant
mC                              =  0.5;                    % kg                  Constant
omega                           =  0.7853981633974483;     % rad/sec             Constant
r                               =  0.025;                  % m                   Constant
radiusCircle                    =  0.3;                    % rad                 Constant

db                              =  1.514188317724914;      % m                   Initial Value
qa                              = -5.495323031258774E-16;  % deg                 Initial Value
qc                              =  10.80692287486033;      % deg                 Initial Value

tInitial                        =  0.0;                    % seconds             Initial Time
tFinal                          =  8.0;                    % seconds             Final Time
tStep                           =  0.01;                   % seconds             Integration Step
printIntScreen                  =  0;                      % 0 or +integer       0 is NO screen output
printIntFile                    =  1;                      % 0 or +integer       0 is NO file   output
absError                        =  1.0E-05;                %                     Absolute Error
relError                        =  1.0E-08;                %                     Relative Error
%-------------------------------+--------------------------+-------------------+-----------------

% Unit conversions
DEGtoRAD = pi / 180.0;
qa = qa * DEGtoRAD;
qc = qc * DEGtoRAD;

VAR = SetMatrixFromNamedQuantities;
[t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, tStep, absError, relError, VAR, printIntScreen, printIntFile );
OutputToScreenOrFile( [], 0, 0 );   % Close output files.
if( printIntFile ~= 0 ),  PlotOutputFiles;  end


%===========================================================================
function sys = mdlDerivatives( t, VAR, uSimulink )
%===========================================================================
SetNamedQuantitiesFromMatrix( VAR );
% Quantities previously specified in MotionGenesis.
qaDt = -omega*radiusCircle*cos(omega*t+qa)/(db+LC*cos(qc));
dbDt = -omega*radiusCircle*sin(omega*t+qa);
qcDt = 0;
qaDDt = -(2*dbDt*qaDt-radiusCircle*omega^2*sin(omega*t+qa)-LC*sin(qc)*qaDt*qcDt)/(db+LC*cos(qc));
dbDDt = (db+LC*cos(qc))*qaDt^2 - radiusCircle*omega^2*cos(omega*t+qa);
qcDDt = 0;

sys = transpose( SetMatrixOfDerivativesPriorToIntegrationStep );
end



%===========================================================================
function VAR = SetMatrixFromNamedQuantities
%===========================================================================
VAR = zeros( 1, 3 );
VAR(1) = db;
VAR(2) = qa;
VAR(3) = qc;
end


%===========================================================================
function SetNamedQuantitiesFromMatrix( VAR )
%===========================================================================
db = VAR(1);
qa = VAR(2);
qc = VAR(3);
end


%===========================================================================
function VARp = SetMatrixOfDerivativesPriorToIntegrationStep
%===========================================================================
VARp = zeros( 1, 3 );
VARp(1) = dbDt;
VARp(2) = qaDt;
VARp(3) = qcDt;
end



%===========================================================================
function Output = mdlOutputs( t, VAR, uSimulink )
%===========================================================================
IAYY = 0.5*mA*r^2;
ICXX = 0.5*mC*r^2;
ICYY = 0.08333333333333333*mC*(LC^2+3*r^2);
TA = IAYY*qaDDt + mC*db*(2*dbDt*qaDt+db*qaDDt) + 0.5*mB*db*(2*dbDt*qaDt+db*qaDDt) + ICXX*sin(qc)*(2*cos(qc)*qaDt*qcDt+sin(qc)*  ...
qaDDt) + 0.25*cos(qc)*(4*ICYY*cos(qc)*qaDDt+4*LC*mC*(dbDt*qaDt+db*qaDDt)-8*ICYY*sin(qc)*qaDt*qcDt-mC*LC^2*(2*sin(qc)*qaDt*qcDt-cos(  ...
qc)*qaDDt)) - LC*mC*db*sin(qc)*qaDt*qcDt;
FB = -0.5*LC*mC*cos(qc)*qcDt^2 - 0.5*LC*mC*sin(qc)*qcDDt - 0.5*mB*(db*qaDt^2-dbDDt) - 0.5*mC*(2*db*qaDt^2+LC*cos(qc)*qaDt^2-2*  ...
dbDDt);
TC = 0.5*g*LC*mC*cos(qc) + ICYY*qcDDt + 0.25*mC*LC^2*qcDDt + 0.25*LC*mC*sin(qc)*(2*db*qaDt^2+LC*cos(qc)*qaDt^2-2*dbDDt)  ...
- (ICXX-ICYY)*sin(qc)*cos(qc)*qaDt^2;
BCx = db*cos(qa);
BCz = -db*sin(qa);
Qx = cos(qa)*(db+LC*cos(qc));
Qy = LA + LC*sin(qc);
Qz = -sin(qa)*(db+LC*cos(qc));

Output = zeros( 1, 28 );
Output(1) = t;
Output(2) = 0.0;
Output(3) = LA;
Output(4) = 0.0;
Output(5) = BCx;
Output(6) = LA;
Output(7) = BCz;
Output(8) = Qx;
Output(9) = Qy;
Output(10) = Qz;

Output(11) = t;
Output(12) = qa;
Output(13) = db;
Output(14) = qc;

Output(15) = t;
Output(16) = qa;
Output(17) = db;
Output(18) = qc;

Output(19) = t;
Output(20) = Qx;
Output(21) = Qy;
Output(22) = Qz;

Output(23) = Qz;
Output(24) = Qx;

Output(25) = t;
Output(26) = TA;
Output(27) = FB;
Output(28) = TC;
end


%===========================================================================
function OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile )
%===========================================================================
persistent FileIdentifier hasHeaderInformationBeenWritten;

if( isempty(Output) ),
   if( ~isempty(FileIdentifier) ),
      for( i = 1 : 6 ),  fclose( FileIdentifier(i) );  end
      clear FileIdentifier;
      fprintf( 1, '\n Output is in the files F21_Robot_RPR_InClass.i  (i=1, ..., 6)\n\n' );
   end
   clear hasHeaderInformationBeenWritten;
   return;
end

if( isempty(hasHeaderInformationBeenWritten) ),
   if( shouldPrintToScreen ),
      fprintf( 1,                '%%       t             BAx            BAy            BAz            BCx            BCy            BCz            Qx             Qy             Qz\n' );
      fprintf( 1,                '%%     (sec)           (m)            (m)            (m)            (m)            (m)            (m)            (m)            (m)            (m)\n\n' );
   end
   if( shouldPrintToFile && isempty(FileIdentifier) ),
      FileIdentifier = zeros( 1, 6 );
      FileIdentifier(1) = fopen('F21_Robot_RPR_InClass.1', 'wt');   if( FileIdentifier(1) == -1 ), error('Error: unable to open file F21_Robot_RPR_InClass.1'); end
      fprintf(FileIdentifier(1), '%% FILE: F21_Robot_RPR_InClass.1\n%%\n' );
      fprintf(FileIdentifier(1), '%%       t             BAx            BAy            BAz            BCx            BCy            BCz            Qx             Qy             Qz\n' );
      fprintf(FileIdentifier(1), '%%     (sec)           (m)            (m)            (m)            (m)            (m)            (m)            (m)            (m)            (m)\n\n' );
      FileIdentifier(2) = fopen('F21_Robot_RPR_InClass.2', 'wt');   if( FileIdentifier(2) == -1 ), error('Error: unable to open file F21_Robot_RPR_InClass.2'); end
      fprintf(FileIdentifier(2), '%% FILE: F21_Robot_RPR_InClass.2\n%%\n' );
      fprintf(FileIdentifier(2), '%%       t             qA             dB             qC\n' );
      fprintf(FileIdentifier(2), '%%     (sec)          (rad)           (m)           (rad)\n\n' );
      FileIdentifier(3) = fopen('F21_Robot_RPR_InClass.3', 'wt');   if( FileIdentifier(3) == -1 ), error('Error: unable to open file F21_Robot_RPR_InClass.3'); end
      fprintf(FileIdentifier(3), '%% FILE: F21_Robot_RPR_InClass.3\n%%\n' );
      fprintf(FileIdentifier(3), '%%       t             qA             dB             qC\n' );
      fprintf(FileIdentifier(3), '%%     (sec)          (rad)           (m)           (rad)\n\n' );
      FileIdentifier(4) = fopen('F21_Robot_RPR_InClass.4', 'wt');   if( FileIdentifier(4) == -1 ), error('Error: unable to open file F21_Robot_RPR_InClass.4'); end
      fprintf(FileIdentifier(4), '%% FILE: F21_Robot_RPR_InClass.4\n%%\n' );
      fprintf(FileIdentifier(4), '%%       t             Qx             Qy             Qz\n' );
      fprintf(FileIdentifier(4), '%%     (sec)           (m)            (m)            (m)\n\n' );
      FileIdentifier(5) = fopen('F21_Robot_RPR_InClass.5', 'wt');   if( FileIdentifier(5) == -1 ), error('Error: unable to open file F21_Robot_RPR_InClass.5'); end
      fprintf(FileIdentifier(5), '%% FILE: F21_Robot_RPR_InClass.5\n%%\n' );
      fprintf(FileIdentifier(5), '%%      Qz             Qx\n' );
      fprintf(FileIdentifier(5), '%%      (m)            (m)\n\n' );
      FileIdentifier(6) = fopen('F21_Robot_RPR_InClass.6', 'wt');   if( FileIdentifier(6) == -1 ), error('Error: unable to open file F21_Robot_RPR_InClass.6'); end
      fprintf(FileIdentifier(6), '%% FILE: F21_Robot_RPR_InClass.6\n%%\n' );
      fprintf(FileIdentifier(6), '%%       t             TA             FB             TC\n' );
      fprintf(FileIdentifier(6), '%%     (sec)          (N*m)           (N)           (N*m)\n\n' );
   end
   hasHeaderInformationBeenWritten = 1;
end

if( shouldPrintToScreen ), WriteNumericalData( 1,                 Output(1:10) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(1), Output(1:10) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(2), Output(11:14) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(3), Output(15:18) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(4), Output(19:22) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(5), Output(23:24) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(6), Output(25:28) );  end
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
data = load( 'F21_Robot_RPR_InClass.3' ); 
plot( data(:,1),data(:,2),'-b', data(:,1),data(:,3),'-.g', data(:,1),data(:,4),'--r', 'LineWidth',3 );
legend( 'qA (rad)', 'dB (m)', 'qC (rad)' );
xlabel('t (sec)');   % ylabel('Some y-axis label');   title('Some plot title');
clear data;

figure;
data = load( 'F21_Robot_RPR_InClass.4' ); 
plot( data(:,1),data(:,2),'-b', data(:,1),data(:,3),'-.g', data(:,1),data(:,4),'--r', 'LineWidth',3 );
legend( 'Qx (m)', 'Qy (m)', 'Qz (m)' );
xlabel('t (sec)');   % ylabel('Some y-axis label');   title('Some plot title');
clear data;

figure;
data = load( 'F21_Robot_RPR_InClass.5' ); 
plot( data(:,1),data(:,2),'-b', 'LineWidth',3 );
legend( 'Qx (m)' );
xlabel('Qz (m)');   ylabel('Qx (m)');   % title('Some plot title');
clear data;

figure;
data = load( 'F21_Robot_RPR_InClass.6' ); 
plot( data(:,1),data(:,2),'-b', data(:,1),data(:,3),'-.g', data(:,1),data(:,4),'--r', 'LineWidth',3 );
legend( 'TA (N*m)', 'FB (N)', 'TC (N*m)' );
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


%==============================================
end    % End of function F21_Robot_RPR_InClass
%==============================================
