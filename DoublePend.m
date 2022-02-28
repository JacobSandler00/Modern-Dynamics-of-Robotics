function [t,VAR,Output] = DoublePend
%===========================================================================
% File: DoublePend.m created Nov 02 2021 by MotionGenesis 6.1.
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
qADDt=0; qBDDt=0; KE=0; MechanicalEnergy=0; PE=0;


%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
g                               =  9.81;                   % m/s^2               Constant
IAx                             =  50;                     % g*cm^2              Constant
IBx                             =  2500;                   % g*cm^2              Constant
IBy                             =  500;                    % g*cm^2              Constant
IBz                             =  2000;                   % g*cm^2              Constant
LA                              =  7.5;                    % cm                  Constant
LB                              =  20;                     % cm                  Constant
mA                              =  10;                     % grams               Constant
mB                              =  100;                    % grams               Constant

qA                              =  90;                     % deg                 Initial Value
qB                              =  1.0;                    % deg                 Initial Value
qADt                            =  0.0;                    % rad/sec             Initial Value
qBDt                            =  0.0;                    % rad/sec             Initial Value

tInitial                        =  0.0;                    % second              Initial Time
tFinal                          =  10;                     % sec                 Final Time
tStep                           =  0.02;                   % sec                 Integration Step
printIntScreen                  =  1;                      % 0 or +integer       0 is NO screen output
printIntFile                    =  1;                      % 0 or +integer       0 is NO file   output
absError                        =  1.0E-07;                %                     Absolute Error
relError                        =  1.0E-07;                %                     Relative Error
%-------------------------------+--------------------------+-------------------+-----------------

% Unit conversions.  UnitSystem: kg, meter, second.
DEGtoRAD = pi / 180.0;
RADtoDEG = 180.0 / pi;
IAx = IAx * 1.0E-07;                                       %  Converted from g*cm^2 
IBx = IBx * 1.0E-07;                                       %  Converted from g*cm^2 
IBy = IBy * 1.0E-07;                                       %  Converted from g*cm^2 
IBz = IBz * 1.0E-07;                                       %  Converted from g*cm^2 
LA = LA * 0.01;                                            %  Converted from cm 
LB = LB * 0.01;                                            %  Converted from cm 
mA = mA * 0.001;                                           %  Converted from grams 
mB = mB * 0.001;                                           %  Converted from grams 
qA = qA * DEGtoRAD;                                        %  Converted from deg 
qB = qB * DEGtoRAD;                                        %  Converted from deg 

VAR = SetMatrixFromNamedQuantities;
[t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, tStep, absError, relError, VAR, printIntScreen, printIntFile );
OutputToScreenOrFile( [], 0, 0 );   % Close output files.
if( printIntFile ~= 0 ),  PlotOutputFiles;  end


%===========================================================================
function sys = mdlDerivatives( t, VAR, uSimulink )
%===========================================================================
SetNamedQuantitiesFromMatrix( VAR );
qADDt = -(g*(mA*LA+mB*LB)*sin(qA)-2*(IBx-IBy)*sin(qB)*cos(qB)*qADt*qBDt)/(IAx+mA*LA^2+mB*LB^2+IBx*cos(qB)^2+IBy*sin(qB)^2);
qBDDt = -(IBx-IBy)*sin(qB)*cos(qB)*qADt^2/IBz;

sys = transpose( SetMatrixOfDerivativesPriorToIntegrationStep );
end



%===========================================================================
function VAR = SetMatrixFromNamedQuantities
%===========================================================================
VAR = zeros( 1, 4 );
VAR(1) = qA;
VAR(2) = qB;
VAR(3) = qADt;
VAR(4) = qBDt;
end


%===========================================================================
function SetNamedQuantitiesFromMatrix( VAR )
%===========================================================================
qA = VAR(1);
qB = VAR(2);
qADt = VAR(3);
qBDt = VAR(4);
end


%===========================================================================
function VARp = SetMatrixOfDerivativesPriorToIntegrationStep
%===========================================================================
VARp = zeros( 1, 4 );
VARp(1) = qADt;
VARp(2) = qBDt;
VARp(3) = qADDt;
VARp(4) = qBDDt;
end



%===========================================================================
function Output = mdlOutputs( t, VAR, uSimulink )
%===========================================================================
KE = 0.5*IAx*qADt^2 + 0.5*IBx*qADt^2 + 0.5*IBz*qBDt^2 + 0.5*mA*LA^2*qADt^2 + 0.5*mB*LB^2*qADt^2 - 0.5*(IBx-IBy)*sin(qB)^2*qADt^2;
PE = -g*(mA*LA+mB*LB)*cos(qA);
MechanicalEnergy = PE + KE;

Output = zeros( 1, 6 );
Output(1) = t;
Output(2) = qA*RADtoDEG;                              % Converted to deg
Output(3) = qB*RADtoDEG;                              % Converted to deg
Output(4) = MechanicalEnergy;

Output(5) = t;
Output(6) = qB*RADtoDEG;
end


%===========================================================================
function OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile )
%===========================================================================
persistent FileIdentifier hasHeaderInformationBeenWritten;

if( isempty(Output) ),
   if( ~isempty(FileIdentifier) ),
      for( i = 1 : 2 ),  fclose( FileIdentifier(i) );  end
      clear FileIdentifier;
      fprintf( 1, '\n Output is in the files DoublePend.i  (i=1,2)\n\n' );
   end
   clear hasHeaderInformationBeenWritten;
   return;
end

if( isempty(hasHeaderInformationBeenWritten) ),
   if( shouldPrintToScreen ),
      fprintf( 1,                '%%       t             qA             qB       MechanicalEnergy\n' );
      fprintf( 1,                '%%     (sec)          (deg)          (deg)         (Joules)\n\n' );
   end
   if( shouldPrintToFile && isempty(FileIdentifier) ),
      FileIdentifier = zeros( 1, 2 );
      FileIdentifier(1) = fopen('DoublePend.1', 'wt');   if( FileIdentifier(1) == -1 ), error('Error: unable to open file DoublePend.1'); end
      fprintf(FileIdentifier(1), '%% FILE: DoublePend.1\n%%\n' );
      fprintf(FileIdentifier(1), '%%       t             qA             qB       MechanicalEnergy\n' );
      fprintf(FileIdentifier(1), '%%     (sec)          (deg)          (deg)         (Joules)\n\n' );
      FileIdentifier(2) = fopen('DoublePend.2', 'wt');   if( FileIdentifier(2) == -1 ), error('Error: unable to open file DoublePend.2'); end
      fprintf(FileIdentifier(2), '%% FILE: DoublePend.2\n%%\n' );
      fprintf(FileIdentifier(2), '%%       t             qb\n' );
      fprintf(FileIdentifier(2), '%%     (sec)          (deg)\n\n' );
   end
   hasHeaderInformationBeenWritten = 1;
end

if( shouldPrintToScreen ), WriteNumericalData( 1,                 Output(1:4) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(1), Output(1:4) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(2), Output(5:6) );  end
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
data = load( 'DoublePend.2' ); 
plot( data(:,1),data(:,2),'-b', 'LineWidth',3 );
legend( 'qb (deg)' );
xlabel('t (sec)');   ylabel('qb (deg)');   % title('Some plot title');
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


%===================================
end    % End of function DoublePend
%===================================
