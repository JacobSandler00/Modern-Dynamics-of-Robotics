function [t,VAR,Output] = quizinclass
%===========================================================================
% File: quizinclass.m created Oct 19 2021 by MotionGenesis 6.1.
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
y=0; qaDt=0; qbDt=0; qcDt=0; qdDt=0; xDt=0; yDt=0; armangle=0;


%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
LA                              =  0.85;                   % m                   Constant
LB                              =  0.92;                   % m                   Constant
LC                              =  0.65;                   % UNITS               Constant
R                               =  0.3;                    % m                   Constant

qa                              =  39.1956341347808;       % deg                 Initial Value
qb                              =  9.914781131427836E-07;  % deg                 Initial Value
qc                              =  96.99140663080878;      % deg                 Initial Value
qd                              =  57.79577347220548;      % deg                 Initial Value
x                               =  1.042365099185502;      % m                   Initial Value

tInitial                        =  0.0;                    % second              Initial Time
tFinal                          =  6;                      % sec                 Final Time
tStep                           =  .01;                    % sec                 Integration Step
printIntScreen                  =  0;                      % 0 or +integer       0 is NO screen output
printIntFile                    =  1;                      % 0 or +integer       0 is NO file   output
absError                        =  1.0E-05;                %                     Absolute Error
relError                        =  1.0E-08;                %                     Relative Error
%-------------------------------+--------------------------+-------------------+-----------------

% Unit conversions
DEGtoRAD = pi / 180.0;
RADtoDEG = 180.0 / pi;
qa = qa * DEGtoRAD;
qb = qb * DEGtoRAD;
qc = qc * DEGtoRAD;
qd = qd * DEGtoRAD;

% Evaluate constants
qbDt = 0;


VAR = SetMatrixFromNamedQuantities;
[t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, tStep, absError, relError, VAR, printIntScreen, printIntFile );
OutputToScreenOrFile( [], 0, 0 );   % Close output files.
if( printIntFile ~= 0 ),  PlotOutputFiles;  end


%===========================================================================
function sys = mdlDerivatives( t, VAR, uSimulink )
%===========================================================================
SetNamedQuantitiesFromMatrix( VAR );
yDt = -0.6283185307179586*sin(pi*t);

xDt = -(LA*sin(qa)+2*LB*sin(qa+qb))*yDt/(LA*cos(qa)+2*LB*cos(qa+qb));
qaDt = 2*yDt/(LA*cos(qa)+2*LB*cos(qa+qb));
qcDt = (1+2*LC*cos(qa-qc)/(LA*cos(qa)+2*LB*cos(qa+qb)))*yDt/(LC*cos(qa-qc));
qdDt = ((R-LC*sin(qa-qc))/(LC*cos(qa-qc))+(LA*sin(qa)+2*LB*sin(qa+qb))/(LA*cos(qa)+2*LB*cos(qa+qb)))*yDt/R;

sys = transpose( SetMatrixOfDerivativesPriorToIntegrationStep );
end



%===========================================================================
function VAR = SetMatrixFromNamedQuantities
%===========================================================================
VAR = zeros( 1, 5 );
VAR(1) = qa;
VAR(2) = qb;
VAR(3) = qc;
VAR(4) = qd;
VAR(5) = x;
end


%===========================================================================
function SetNamedQuantitiesFromMatrix( VAR )
%===========================================================================
qa = VAR(1);
qb = VAR(2);
qc = VAR(3);
qd = VAR(4);
x = VAR(5);
end


%===========================================================================
function VARp = SetMatrixOfDerivativesPriorToIntegrationStep
%===========================================================================
VARp = zeros( 1, 5 );
VARp(1) = qaDt;
VARp(2) = qbDt;
VARp(3) = qcDt;
VARp(4) = qdDt;
VARp(5) = xDt;
end



%===========================================================================
function Output = mdlOutputs( t, VAR, uSimulink )
%===========================================================================
y = 0.65 + 0.2*cos(pi*t);
armangle = acos(cos(qa-qc));

Output = zeros( 1, 10 );
Output(1) = t;
Output(2) = qa*RADtoDEG;
Output(3) = qb*RADtoDEG;
Output(4) = qc*RADtoDEG;
Output(5) = qd*RADtoDEG;

Output(6) = x;
Output(7) = y;

Output(8) = t;
Output(9) = qa*RADtoDEG;
Output(10) = armangle*RADtoDEG;
end


%===========================================================================
function OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile )
%===========================================================================
persistent FileIdentifier hasHeaderInformationBeenWritten;

if( isempty(Output) ),
   if( ~isempty(FileIdentifier) ),
      for( i = 1 : 3 ),  fclose( FileIdentifier(i) );  end
      clear FileIdentifier;
      fprintf( 1, '\n Output is in the files quizinclass.i  (i=1,2,3)\n\n' );
   end
   clear hasHeaderInformationBeenWritten;
   return;
end

if( isempty(hasHeaderInformationBeenWritten) ),
   if( shouldPrintToScreen ),
      fprintf( 1,                '%%       t             qa             qb             qc             qd\n' );
      fprintf( 1,                '%%   (second)         (deg)          (deg)          (deg)          (deg)\n\n' );
   end
   if( shouldPrintToFile && isempty(FileIdentifier) ),
      FileIdentifier = zeros( 1, 3 );
      FileIdentifier(1) = fopen('quizinclass.1', 'wt');   if( FileIdentifier(1) == -1 ), error('Error: unable to open file quizinclass.1'); end
      fprintf(FileIdentifier(1), '%% FILE: quizinclass.1\n%%\n' );
      fprintf(FileIdentifier(1), '%%       t             qa             qb             qc             qd\n' );
      fprintf(FileIdentifier(1), '%%   (second)         (deg)          (deg)          (deg)          (deg)\n\n' );
      FileIdentifier(2) = fopen('quizinclass.2', 'wt');   if( FileIdentifier(2) == -1 ), error('Error: unable to open file quizinclass.2'); end
      fprintf(FileIdentifier(2), '%% FILE: quizinclass.2\n%%\n' );
      fprintf(FileIdentifier(2), '%%       x              y\n' );
      fprintf(FileIdentifier(2), '%%      (m)            (m)\n\n' );
      FileIdentifier(3) = fopen('quizinclass.3', 'wt');   if( FileIdentifier(3) == -1 ), error('Error: unable to open file quizinclass.3'); end
      fprintf(FileIdentifier(3), '%% FILE: quizinclass.3\n%%\n' );
      fprintf(FileIdentifier(3), '%%       t             qa          armangle\n' );
      fprintf(FileIdentifier(3), '%%     (sec)          (deg)          (deg)\n\n' );
   end
   hasHeaderInformationBeenWritten = 1;
end

if( shouldPrintToScreen ), WriteNumericalData( 1,                 Output(1:5) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(1), Output(1:5) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(2), Output(6:7) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(3), Output(8:10) );  end
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
data = load( 'quizinclass.1' ); 
plot( data(:,1),data(:,2),'-b', data(:,1),data(:,3),'-.g', data(:,1),data(:,4),'--r', data(:,1),data(:,5),'-m', 'LineWidth',3 );
legend( 'qa (deg)', 'qb (deg)', 'qc (deg)', 'qd (deg)' );
xlabel('t (second)');   % ylabel('Some y-axis label');   title('Some plot title');
clear data;

figure;
data = load( 'quizinclass.2' ); 
plot( data(:,1),data(:,2),'-b', 'LineWidth',3 );
legend( 'y (m)' );
xlabel('x (m)');   ylabel('y (m)');   % title('Some plot title');
clear data;

figure;
data = load( 'quizinclass.3' ); 
plot( data(:,1),data(:,2),'-b', data(:,1),data(:,3),'-.g', 'LineWidth',3 );
legend( 'qa (deg)', 'armangle (deg)' );
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


%====================================
end    % End of function quizinclass
%====================================
