function [t,VAR,Output] = Robotics_Quiz2
%===========================================================================
% File: Robotics_Quiz2.m created Oct 10 2021 by MotionGenesis 6.1.
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
qa=0; qb=0; x=0; qbDt=0; qcDt=0; qdDt=0; y=0;


%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
LA                              =  0.85;                   % m                   Constant
LB                              =  0.92;                   % m                   Constant
LC                              =  0.65;                   % m                   Constant

qc                              =  3;                      % UNITS               Initial Value
qd                              =  2;                      % UNITS               Initial Value

tInitial                        =  0.0;                    % second              Initial Time
tFinal                          =  6.0;                    % s                   Final Time
tStep                           =  0.01;                   % s                   Integration Step
printIntScreen                  =  0;                      % 0 or +integer       0 is NO screen output
printIntFile                    =  1;                      % 0 or +integer       0 is NO file   output
absError                        =  1.0E-05;                %                     Absolute Error
relError                        =  1.0E-08;                %                     Relative Error
%-------------------------------+--------------------------+-------------------+-----------------

% Unit conversions
RADtoDEG = 180.0 / pi;

% Evaluate constants
qa = 1.570796326794897;
qb = 0;
x = 0.5*LA*cos(qa) + LB*cos(qa+qb);
qbDt = 0;


VAR = SetMatrixFromNamedQuantities;
[t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, tStep, absError, relError, VAR, printIntScreen, printIntFile );
OutputToScreenOrFile( [], 0, 0 );   % Close output files.
if( printIntFile ~= 0 ),  PlotOutputFiles;  end


%===========================================================================
function sys = mdlDerivatives( t, VAR, uSimulink )
%===========================================================================
SetNamedQuantitiesFromMatrix( VAR );
qcDt = LB*(sin(qa+qb)*(LA*cos(qa)+2*LC*cos(qa-qc))-cos(qa+qb)*(LA*sin(qa)+2*LC*sin(qa-qc)))*qbDt/(LC*(sin(qa-qc)*(LA*cos(qa)+2*LB*  ...
cos(qa+qb))-cos(qa-qc)*(LA*sin(qa)+2*LB*sin(qa+qb))));
qdDt = LA*LB*sin(qb)*qbDt/(LC*(sin(qa-qc)*(LA*cos(qa)+2*LB*cos(qa+qb))-cos(qa-qc)*(LA*sin(qa)+2*LB*sin(qa+qb))));

% Quantities previously specified in MotionGenesis.
y = 0.5*LA*sin(qa) + LB*sin(qa+qb);

sys = transpose( SetMatrixOfDerivativesPriorToIntegrationStep );
end



%===========================================================================
function VAR = SetMatrixFromNamedQuantities
%===========================================================================
VAR = zeros( 1, 2 );
VAR(1) = qc;
VAR(2) = qd;
end


%===========================================================================
function SetNamedQuantitiesFromMatrix( VAR )
%===========================================================================
qc = VAR(1);
qd = VAR(2);
end


%===========================================================================
function VARp = SetMatrixOfDerivativesPriorToIntegrationStep
%===========================================================================
VARp = zeros( 1, 2 );
VARp(1) = qcDt;
VARp(2) = qdDt;
end



%===========================================================================
function Output = mdlOutputs( t, VAR, uSimulink )
%===========================================================================
Output = zeros( 1, 10 );
Output(1) = x;
Output(2) = y;

Output(3) = t;
Output(4) = qa*RADtoDEG;

Output(5) = t;
Output(6) = qb*RADtoDEG;

Output(7) = t;
Output(8) = qc*RADtoDEG;

Output(9) = t;
Output(10) = qd*RADtoDEG;
end


%===========================================================================
function OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile )
%===========================================================================
persistent FileIdentifier hasHeaderInformationBeenWritten;

if( isempty(Output) ),
   if( ~isempty(FileIdentifier) ),
      for( i = 1 : 5 ),  fclose( FileIdentifier(i) );  end
      clear FileIdentifier;
      fprintf( 1, '\n Output is in the files Robotics_Quiz2.i  (i=1, ..., 5)\n\n' );
   end
   clear hasHeaderInformationBeenWritten;
   return;
end

if( isempty(hasHeaderInformationBeenWritten) ),
   if( shouldPrintToScreen ),
      fprintf( 1,                '%%       x              y\n' );
      fprintf( 1,                '%%      (m)            (m)\n\n' );
   end
   if( shouldPrintToFile && isempty(FileIdentifier) ),
      FileIdentifier = zeros( 1, 5 );
      FileIdentifier(1) = fopen('Robotics_Quiz2.1', 'wt');   if( FileIdentifier(1) == -1 ), error('Error: unable to open file Robotics_Quiz2.1'); end
      fprintf(FileIdentifier(1), '%% FILE: Robotics_Quiz2.1\n%%\n' );
      fprintf(FileIdentifier(1), '%%       x              y\n' );
      fprintf(FileIdentifier(1), '%%      (m)            (m)\n\n' );
      FileIdentifier(2) = fopen('Robotics_Quiz2.2', 'wt');   if( FileIdentifier(2) == -1 ), error('Error: unable to open file Robotics_Quiz2.2'); end
      fprintf(FileIdentifier(2), '%% FILE: Robotics_Quiz2.2\n%%\n' );
      fprintf(FileIdentifier(2), '%%       t             qa\n' );
      fprintf(FileIdentifier(2), '%%      (s)         (degrees)\n\n' );
      FileIdentifier(3) = fopen('Robotics_Quiz2.3', 'wt');   if( FileIdentifier(3) == -1 ), error('Error: unable to open file Robotics_Quiz2.3'); end
      fprintf(FileIdentifier(3), '%% FILE: Robotics_Quiz2.3\n%%\n' );
      fprintf(FileIdentifier(3), '%%       t             qb\n' );
      fprintf(FileIdentifier(3), '%%      (s)         (degrees)\n\n' );
      FileIdentifier(4) = fopen('Robotics_Quiz2.4', 'wt');   if( FileIdentifier(4) == -1 ), error('Error: unable to open file Robotics_Quiz2.4'); end
      fprintf(FileIdentifier(4), '%% FILE: Robotics_Quiz2.4\n%%\n' );
      fprintf(FileIdentifier(4), '%%       t             qc\n' );
      fprintf(FileIdentifier(4), '%%      (s)         (degrees)\n\n' );
      FileIdentifier(5) = fopen('Robotics_Quiz2.5', 'wt');   if( FileIdentifier(5) == -1 ), error('Error: unable to open file Robotics_Quiz2.5'); end
      fprintf(FileIdentifier(5), '%% FILE: Robotics_Quiz2.5\n%%\n' );
      fprintf(FileIdentifier(5), '%%       t             qd\n' );
      fprintf(FileIdentifier(5), '%%      (s)         (degrees)\n\n' );
   end
   hasHeaderInformationBeenWritten = 1;
end

if( shouldPrintToScreen ), WriteNumericalData( 1,                 Output(1:2) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(1), Output(1:2) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(2), Output(3:4) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(3), Output(5:6) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(4), Output(7:8) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(5), Output(9:10) );  end
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
data = load( 'Robotics_Quiz2.1' ); 
plot( data(:,1),data(:,2),'-b', 'LineWidth',3 );
legend( 'y (m)' );
xlabel('x (m)');   ylabel('y (m)');   % title('Some plot title');
clear data;

figure;
data = load( 'Robotics_Quiz2.2' ); 
plot( data(:,1),data(:,2),'-b', 'LineWidth',3 );
legend( 'qa (degrees)' );
xlabel('t (s)');   ylabel('qa (degrees)');   % title('Some plot title');
clear data;

figure;
data = load( 'Robotics_Quiz2.3' ); 
plot( data(:,1),data(:,2),'-b', 'LineWidth',3 );
legend( 'qb (degrees)' );
xlabel('t (s)');   ylabel('qb (degrees)');   % title('Some plot title');
clear data;

figure;
data = load( 'Robotics_Quiz2.4' ); 
plot( data(:,1),data(:,2),'-b', 'LineWidth',3 );
legend( 'qc (degrees)' );
xlabel('t (s)');   ylabel('qc (degrees)');   % title('Some plot title');
clear data;

figure;
data = load( 'Robotics_Quiz2.5' ); 
plot( data(:,1),data(:,2),'-b', 'LineWidth',3 );
legend( 'qd (degrees)' );
xlabel('t (s)');   ylabel('qd (degrees)');   % title('Some plot title');
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


%=======================================
end    % End of function Robotics_Quiz2
%=======================================
