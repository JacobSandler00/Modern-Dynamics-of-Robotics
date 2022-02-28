function [t,VAR,Output] = FInalComputation
%===========================================================================
% File: FInalComputation.m created Dec 08 2021 by MotionGenesis 6.1.
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
qa1=0; qa2=0; qb1=0; qb2=0; qc1=0; qc2=0; qd1=0; qd2=0; qa=0; qb=0; qc=0; qd=0; qe=0; TAN=0; TBA=0; TCB=0; TDC=0; TED=0; dumbDt=0;
qaDt=0; qbDt=0; qcDt=0; qdDt=0; qeDt=0; qaDDt=0; qbDDt=0; qcDDt=0; qdDDt=0; qeDDt=0;


%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
g                               =  9.81;                   % m/s^2               Constant
Ia                              =  0.05567153725950489;    % UNITS               Constant
Ib                              =  0.2740537379284354;     % UNITS               Constant
Ic                              =  3.743559541263018;      % UNITS               Constant
Id                              =  0.01225765153256133;    % UNITS               Constant
Ie                              =  0.002279841968795617;   % UNITS               Constant
LA                              =  432.3;                  % mm                  Constant
LB                              =  368.5;                  % mm                  Constant
LC                              =  529.3;                  % mm                  Constant
LD                              =  275.1;                  % mm                  Constant
LE                              =  264.3;                  % mm                  Constant
ma                              =  2.97739;                % UNITS               Constant
mb                              =  9.14882;                % UNITS               Constant
mBox                            =  5;                      % kg                  Constant
mc                              =  26.35083;               % UNITS               Constant
md                              =  1.57845;                % UNITS               Constant
me                              =  0.85422;                % UNITS               Constant

dumb                            =  0;                      % UNITS               Initial Value

tInitial                        =  0.0;                    % second              Initial Time
tFinal                          =  3;                      % s                   Final Time
tStep                           =  .01;                    % s                   Integration Step
printIntScreen                  =  1;                      % 0 or +integer       0 is NO screen output
printIntFile                    =  1;                      % 0 or +integer       0 is NO file   output
absError                        =  1.0E-05;                %                     Absolute Error
relError                        =  1.0E-08;                %                     Relative Error
%-------------------------------+--------------------------+-------------------+-----------------

% Unit conversions.  UnitSystem: kg, meter, second.
RADtoDEG = 180.0 / pi;
LA = LA * 0.001;                                           %  Converted from mm 
LB = LB * 0.001;                                           %  Converted from mm 
LC = LC * 0.001;                                           %  Converted from mm 
LD = LD * 0.001;                                           %  Converted from mm 
LE = LE * 0.001;                                           %  Converted from mm 

% Evaluate constants
qa1 = 1.396263401595464;
qa2 = 0.2617993877991494;
qb1 = 1.396263401595464;
qb2 = 0.3084996989990423;
qc1 = 0.8726646259971648;
qc2 = 0.174532925199433;
qd1 = 4.235771860988153;
qd2 = -1.6483472150552;
dumbDt = 0.1234;


VAR = SetMatrixFromNamedQuantities;
[t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, tStep, absError, relError, VAR, printIntScreen, printIntFile );
OutputToScreenOrFile( [], 0, 0 );   % Close output files.
if( printIntFile ~= 0 ),  PlotOutputFiles;  end


%===========================================================================
function sys = mdlDerivatives( t, VAR, uSimulink )
%===========================================================================
SetNamedQuantitiesFromMatrix( VAR );
sys = transpose( SetMatrixOfDerivativesPriorToIntegrationStep );
end



%===========================================================================
function VAR = SetMatrixFromNamedQuantities
%===========================================================================
VAR = zeros( 1, 1 );
VAR(1) = dumb;
end


%===========================================================================
function SetNamedQuantitiesFromMatrix( VAR )
%===========================================================================
dumb = VAR(1);
end


%===========================================================================
function VARp = SetMatrixOfDerivativesPriorToIntegrationStep
%===========================================================================
VARp = zeros( 1, 1 );
VARp(1) = dumbDt;
end



%===========================================================================
function Output = mdlOutputs( t, VAR, uSimulink )
%===========================================================================
qa = qa1 + 0.1591549430918954*(qa1-qa2)*sin(2.094395102393195*t) - 0.3333333333333333*(qa1-qa2)*t;
qaDt = 0.3333333333333333*qa2 + 0.3333333333333333*(qa1-qa2)*cos(2.094395102393195*t) - 0.3333333333333333*qa1;
qaDDt = -0.6981317007977317*(qa1-qa2)*sin(2.094395102393195*t);
qb = qb1 + 0.1591549430918954*(qb1-qb2)*sin(2.094395102393195*t) - 0.3333333333333333*(qb1-qb2)*t;
qbDt = 0.3333333333333333*qb2 + 0.3333333333333333*(qb1-qb2)*cos(2.094395102393195*t) - 0.3333333333333333*qb1;
qbDDt = -0.6981317007977317*(qb1-qb2)*sin(2.094395102393195*t);
qc = qc1 + 0.1591549430918954*(qc1-qc2)*sin(2.094395102393195*t) - 0.3333333333333333*(qc1-qc2)*t;
qcDt = 0.3333333333333333*qc2 + 0.3333333333333333*(qc1-qc2)*cos(2.094395102393195*t) - 0.3333333333333333*qc1;
qcDDt = -0.6981317007977317*(qc1-qc2)*sin(2.094395102393195*t);
qd = qd1 + 0.1591549430918954*(qd1-qd2)*sin(2.094395102393195*t) - 0.3333333333333333*(qd1-qd2)*t;
qdDt = 0.3333333333333333*qd2 + 0.3333333333333333*(qd1-qd2)*cos(2.094395102393195*t) - 0.3333333333333333*qd1;
qdDDt = -0.6981317007977317*(qd1-qd2)*sin(2.094395102393195*t);
qe = 4.235771860988153 + 0.5512983908219449*sin(2.094395102393195*t) - 1.154636649694731*t;
qeDt = -1.154636649694731 + 1.154636649694731*cos(2.094395102393195*t);
qeDDt = -2.418265344164333*sin(2.094395102393195*t);
TAN = g*LA*mb*sin(qa) + g*LA*mBox*sin(qa) + g*LA*mc*sin(qa) + g*LA*md*sin(qa) + g*LA*me*sin(qa) + g*LC*mBox*sin(qc) + g*LC*md*sin(qc)  ...
+ g*LC*me*sin(qc) + g*LD*mBox*cos(qd) + g*LD*me*cos(qd) + g*LE*mBox*cos(qe) + 0.5*g*LA*ma*sin(qa) + 0.5*g*LC*mc*sin(qc)  ...
+ 0.5*g*LD*md*cos(qd) + 0.5*g*LE*me*cos(qe) + 0.5*LB*(LA*mb*sin(qa+qb)*qaDt^2+mc*(LC*sin(qb+qc)*qcDt^2+2*LA*sin(qa+qb)*qaDt^2)+md*(  ...
2*LA*sin(qa+qb)*qaDt^2+2*LC*sin(qb+qc)*qcDt^2+LD*cos(qb-qd)*qdDt^2)+2*mBox*(LA*sin(qa+qb)*qaDt^2+LC*sin(qb+qc)*qcDt^2+LD*cos(qb-qd)*qdDt^2+  ...
LE*cos(qb-qe)*qeDt^2)+me*(2*LA*sin(qa+qb)*qaDt^2+2*LC*sin(qb+qc)*qcDt^2+LE*cos(qb-qe)*qeDt^2+2*LD*cos(qb-qd)*qdDt^2)) + 0.5*LA*LE*(  ...
me+2*mBox)*sin(qa+qe)*qeDDt + 0.5*LC*LE*(me+2*mBox)*sin(qc+qe)*qeDDt + 0.25*(4*Ie+me*LE^2+4*mBox*LE^2)*qeDDt + 0.5*LD*LE*(me+2*mBox)*  ...
cos(qd-qe)*qdDDt + 0.5*LD*LE*(me+2*mBox)*cos(qd-qe)*qeDDt + 0.5*LA*LD*(md+2*mBox+2*me)*sin(qa+qd)*qdDDt + 0.5*LC*LD*(md+2*mBox+2*  ...
me)*sin(qc+qd)*qdDDt + 0.5*LB*LC*(mc+2*mBox+2*md+2*me)*cos(qb+qc)*qbDDt + 0.25*(4*Id+md*LD^2+4*mBox*LD^2+4*me*LD^2)*qdDDt  ...
+ 0.5*LA*LB*(mb+2*mBox+2*mc+2*md+2*me)*cos(qa+qb)*qbDDt + 0.25*(4*Ib+mb*LB^2+4*mBox*LB^2+4*mc*LB^2+4*md*LB^2+4*me*LB^2)*qbDDt  ...
- g*LB*mBox*sin(qb) - g*LB*mc*sin(qb) - g*LB*md*sin(qb) - g*LB*me*sin(qb) - 0.5*g*LB*mb*sin(qb) - 0.5*LE*(me+2*mBox)*(LA*cos(qa+qe)*qaDt^2+  ...
LC*cos(qc+qe)*qcDt^2+LB*cos(qb-qe)*qbDt^2+LD*sin(qd-qe)*qdDt^2) - 0.5*LD*(md*(LA*cos(qa+qd)*qaDt^2+LC*cos(qc+qd)*qcDt^2+LB*cos(qb-  ...
qd)*qbDt^2)+2*mBox*(LA*cos(qa+qd)*qaDt^2+LC*cos(qc+qd)*qcDt^2+LB*cos(qb-qd)*qbDt^2-LE*sin(qd-qe)*qeDt^2)+me*(2*LA*cos(qa+qd)*qaDt^2+  ...
2*LC*cos(qc+qd)*qcDt^2+2*LB*cos(qb-qd)*qbDt^2-LE*sin(qd-qe)*qeDt^2)) - 0.5*LC*(mc*(LB*sin(qb+qc)*qbDt^2-LA*sin(qa-qc)*qaDt^2)+md*(  ...
2*LB*sin(qb+qc)*qbDt^2-LD*cos(qc+qd)*qdDt^2-2*LA*sin(qa-qc)*qaDt^2)+2*mBox*(LB*sin(qb+qc)*qbDt^2-LD*cos(qc+qd)*qdDt^2-LE*cos(qc+qe)*qeDt^2-  ...
LA*sin(qa-qc)*qaDt^2)+me*(2*LB*sin(qb+qc)*qbDt^2-2*LD*cos(qc+qd)*qdDt^2-LE*cos(qc+qe)*qeDt^2-2*LA*sin(qa-qc)*qaDt^2)) - 0.5*LA*(LB*  ...
mb*sin(qa+qb)*qbDt^2+mc*(2*LB*sin(qa+qb)*qbDt^2+LC*sin(qa-qc)*qcDt^2)-md*(LD*cos(qa+qd)*qdDt^2-2*LB*sin(qa+qb)*qbDt^2-2*LC*sin(qa-  ...
qc)*qcDt^2)-2*mBox*(LD*cos(qa+qd)*qdDt^2+LE*cos(qa+qe)*qeDt^2-LB*sin(qa+qb)*qbDt^2-LC*sin(qa-qc)*qcDt^2)-me*(LE*cos(qa+qe)*qeDt^2+  ...
2*LD*cos(qa+qd)*qdDt^2-2*LB*sin(qa+qb)*qbDt^2-2*LC*sin(qa-qc)*qcDt^2)) - 0.5*LA*LE*(me+2*mBox)*sin(qa+qe)*qaDDt - 0.5*LC*LE*(me+2*  ...
mBox)*sin(qc+qe)*qcDDt - 0.5*LB*LE*(me+2*mBox)*sin(qb-qe)*qbDDt - 0.5*LB*LE*(me+2*mBox)*sin(qb-qe)*qeDDt - 0.5*LA*LD*(md+2*mBox+2*  ...
me)*sin(qa+qd)*qaDDt - 0.5*LC*LD*(md+2*mBox+2*me)*sin(qc+qd)*qcDDt - 0.5*LB*LD*(md+2*mBox+2*me)*sin(qb-qd)*qbDDt - 0.5*LB*LD*(md+  ...
2*mBox+2*me)*sin(qb-qd)*qdDDt - 0.5*LB*LC*(mc+2*mBox+2*md+2*me)*cos(qb+qc)*qcDDt - 0.5*LA*LC*(mc+2*mBox+2*md+2*me)*cos(qa-qc)*qaDDt  ...
- 0.5*LA*LC*(mc+2*mBox+2*md+2*me)*cos(qa-qc)*qcDDt - 0.5*LA*LB*(mb+2*mBox+2*mc+2*md+2*me)*cos(qa+qb)*qaDDt - 0.25*(4*Ic+mc*LC^2+4*  ...
mBox*LC^2+4*md*LC^2+4*me*LC^2)*qcDDt - 0.25*(4*Ia+ma*LA^2+4*mb*LA^2+4*mBox*LA^2+4*mc*LA^2+4*md*LA^2+4*me*LA^2)*qaDDt;
TBA = g*LC*mBox*sin(qc) + g*LC*md*sin(qc) + g*LC*me*sin(qc) + g*LD*mBox*cos(qd) + g*LD*me*cos(qd) + g*LE*mBox*cos(qe) + 0.5*g*LC*  ...
mc*sin(qc) + 0.5*g*LD*md*cos(qd) + 0.5*g*LE*me*cos(qe) + 0.5*LB*(LA*mb*sin(qa+qb)*qaDt^2+mc*(LC*sin(qb+qc)*qcDt^2+2*LA*sin(qa+qb)*qaDt^2)+  ...
md*(2*LA*sin(qa+qb)*qaDt^2+2*LC*sin(qb+qc)*qcDt^2+LD*cos(qb-qd)*qdDt^2)+2*mBox*(LA*sin(qa+qb)*qaDt^2+LC*sin(qb+qc)*qcDt^2+LD*cos(  ...
qb-qd)*qdDt^2+LE*cos(qb-qe)*qeDt^2)+me*(2*LA*sin(qa+qb)*qaDt^2+2*LC*sin(qb+qc)*qcDt^2+LE*cos(qb-qe)*qeDt^2+2*LD*cos(qb-qd)*qdDt^2))  ...
+ 0.5*LC*LE*(me+2*mBox)*sin(qc+qe)*qeDDt + 0.25*(4*Ie+me*LE^2+4*mBox*LE^2)*qeDDt + 0.5*LD*LE*(me+2*mBox)*cos(qd-qe)*qdDDt  ...
+ 0.5*LD*LE*(me+2*mBox)*cos(qd-qe)*qeDDt + 0.5*LC*LD*(md+2*mBox+2*me)*sin(qc+qd)*qdDDt + 0.5*LB*LC*(mc+2*mBox+2*md+2*me)*cos(qb+qc)*qbDDt  ...
+ 0.25*(4*Id+md*LD^2+4*mBox*LD^2+4*me*LD^2)*qdDDt + 0.25*(4*Ib+mb*LB^2+4*mBox*LB^2+4*mc*LB^2+4*md*LB^2+4*me*LB^2)*qbDDt  ...
- g*LB*mBox*sin(qb) - g*LB*mc*sin(qb) - g*LB*md*sin(qb) - g*LB*me*sin(qb) - 0.5*g*LB*mb*sin(qb) - 0.5*LE*(me+2*mBox)*(LA*cos(qa+qe)*qaDt^2+  ...
LC*cos(qc+qe)*qcDt^2+LB*cos(qb-qe)*qbDt^2+LD*sin(qd-qe)*qdDt^2) - 0.5*LD*(md*(LA*cos(qa+qd)*qaDt^2+LC*cos(qc+qd)*qcDt^2+LB*cos(qb-  ...
qd)*qbDt^2)+2*mBox*(LA*cos(qa+qd)*qaDt^2+LC*cos(qc+qd)*qcDt^2+LB*cos(qb-qd)*qbDt^2-LE*sin(qd-qe)*qeDt^2)+me*(2*LA*cos(qa+qd)*qaDt^2+  ...
2*LC*cos(qc+qd)*qcDt^2+2*LB*cos(qb-qd)*qbDt^2-LE*sin(qd-qe)*qeDt^2)) - 0.5*LC*(mc*(LB*sin(qb+qc)*qbDt^2-LA*sin(qa-qc)*qaDt^2)+md*(  ...
2*LB*sin(qb+qc)*qbDt^2-LD*cos(qc+qd)*qdDt^2-2*LA*sin(qa-qc)*qaDt^2)+2*mBox*(LB*sin(qb+qc)*qbDt^2-LD*cos(qc+qd)*qdDt^2-LE*cos(qc+qe)*qeDt^2-  ...
LA*sin(qa-qc)*qaDt^2)+me*(2*LB*sin(qb+qc)*qbDt^2-2*LD*cos(qc+qd)*qdDt^2-LE*cos(qc+qe)*qeDt^2-2*LA*sin(qa-qc)*qaDt^2)) - 0.5*LA*LE*(  ...
me+2*mBox)*sin(qa+qe)*qaDDt - 0.5*LC*LE*(me+2*mBox)*sin(qc+qe)*qcDDt - 0.5*LB*LE*(me+2*mBox)*sin(qb-qe)*qbDDt - 0.5*LB*LE*(me+2*mBox)*  ...
sin(qb-qe)*qeDDt - 0.5*LA*LD*(md+2*mBox+2*me)*sin(qa+qd)*qaDDt - 0.5*LC*LD*(md+2*mBox+2*me)*sin(qc+qd)*qcDDt - 0.5*LB*LD*(md+2*mBox+  ...
2*me)*sin(qb-qd)*qbDDt - 0.5*LB*LD*(md+2*mBox+2*me)*sin(qb-qd)*qdDDt - 0.5*LB*LC*(mc+2*mBox+2*md+2*me)*cos(qb+qc)*qcDDt  ...
- 0.5*LA*LC*(mc+2*mBox+2*md+2*me)*cos(qa-qc)*qaDDt - 0.5*LA*LB*(mb+2*mBox+2*mc+2*md+2*me)*cos(qa+qb)*qaDDt - 0.25*(4*Ic+mc*LC^2+4*  ...
mBox*LC^2+4*md*LC^2+4*me*LC^2)*qcDDt;
TCB = g*LC*mBox*sin(qc) + g*LC*md*sin(qc) + g*LC*me*sin(qc) + g*LD*mBox*cos(qd) + g*LD*me*cos(qd) + g*LE*mBox*cos(qe) + 0.5*g*LC*  ...
mc*sin(qc) + 0.5*g*LD*md*cos(qd) + 0.5*g*LE*me*cos(qe) + 0.5*LC*LE*(me+2*mBox)*sin(qc+qe)*qeDDt + 0.25*(4*Ie+me*LE^2+4*mBox*LE^2)*qeDDt  ...
+ 0.5*LD*LE*(me+2*mBox)*cos(qd-qe)*qdDDt + 0.5*LD*LE*(me+2*mBox)*cos(qd-qe)*qeDDt + 0.5*LC*LD*(md+2*mBox+2*me)*sin(qc+qd)*qdDDt  ...
+ 0.5*LB*LC*(mc+2*mBox+2*md+2*me)*cos(qb+qc)*qbDDt + 0.25*(4*Id+md*LD^2+4*mBox*LD^2+4*me*LD^2)*qdDDt - 0.5*LE*(me+2*mBox)*(LA*cos(  ...
qa+qe)*qaDt^2+LC*cos(qc+qe)*qcDt^2+LB*cos(qb-qe)*qbDt^2+LD*sin(qd-qe)*qdDt^2) - 0.5*LD*(md*(LA*cos(qa+qd)*qaDt^2+LC*cos(qc+qd)*qcDt^2+  ...
LB*cos(qb-qd)*qbDt^2)+2*mBox*(LA*cos(qa+qd)*qaDt^2+LC*cos(qc+qd)*qcDt^2+LB*cos(qb-qd)*qbDt^2-LE*sin(qd-qe)*qeDt^2)+me*(2*LA*cos(qa+  ...
qd)*qaDt^2+2*LC*cos(qc+qd)*qcDt^2+2*LB*cos(qb-qd)*qbDt^2-LE*sin(qd-qe)*qeDt^2)) - 0.5*LC*(mc*(LB*sin(qb+qc)*qbDt^2-LA*sin(qa-qc)*qaDt^2)+  ...
md*(2*LB*sin(qb+qc)*qbDt^2-LD*cos(qc+qd)*qdDt^2-2*LA*sin(qa-qc)*qaDt^2)+2*mBox*(LB*sin(qb+qc)*qbDt^2-LD*cos(qc+qd)*qdDt^2-LE*cos(  ...
qc+qe)*qeDt^2-LA*sin(qa-qc)*qaDt^2)+me*(2*LB*sin(qb+qc)*qbDt^2-2*LD*cos(qc+qd)*qdDt^2-LE*cos(qc+qe)*qeDt^2-2*LA*sin(qa-qc)*qaDt^2))  ...
- 0.5*LA*LE*(me+2*mBox)*sin(qa+qe)*qaDDt - 0.5*LC*LE*(me+2*mBox)*sin(qc+qe)*qcDDt - 0.5*LB*LE*(me+2*mBox)*sin(qb-qe)*qbDDt  ...
- 0.5*LA*LD*(md+2*mBox+2*me)*sin(qa+qd)*qaDDt - 0.5*LC*LD*(md+2*mBox+2*me)*sin(qc+qd)*qcDDt - 0.5*LB*LD*(md+2*mBox+2*me)*sin(qb-qd)*qbDDt  ...
- 0.5*LA*LC*(mc+2*mBox+2*md+2*me)*cos(qa-qc)*qaDDt - 0.25*(4*Ic+mc*LC^2+4*mBox*LC^2+4*md*LC^2+4*me*LC^2)*qcDDt;
TDC = g*LD*mBox*cos(qd) + g*LD*me*cos(qd) + g*LE*mBox*cos(qe) + 0.5*g*LD*md*cos(qd) + 0.5*g*LE*me*cos(qe) + 0.25*(4*Ie+me*LE^2+4*  ...
mBox*LE^2)*qeDDt + 0.5*LD*LE*(me+2*mBox)*cos(qd-qe)*qdDDt + 0.5*LD*LE*(me+2*mBox)*cos(qd-qe)*qeDDt + 0.25*(4*Id+md*LD^2+4*mBox*LD^2+  ...
4*me*LD^2)*qdDDt - 0.5*LE*(me+2*mBox)*(LA*cos(qa+qe)*qaDt^2+LC*cos(qc+qe)*qcDt^2+LB*cos(qb-qe)*qbDt^2+LD*sin(qd-qe)*qdDt^2)  ...
- 0.5*LD*(md*(LA*cos(qa+qd)*qaDt^2+LC*cos(qc+qd)*qcDt^2+LB*cos(qb-qd)*qbDt^2)+2*mBox*(LA*cos(qa+qd)*qaDt^2+LC*cos(qc+qd)*qcDt^2+LB*  ...
cos(qb-qd)*qbDt^2-LE*sin(qd-qe)*qeDt^2)+me*(2*LA*cos(qa+qd)*qaDt^2+2*LC*cos(qc+qd)*qcDt^2+2*LB*cos(qb-qd)*qbDt^2-LE*sin(qd-qe)*qeDt^2))  ...
- 0.5*LA*LE*(me+2*mBox)*sin(qa+qe)*qaDDt - 0.5*LC*LE*(me+2*mBox)*sin(qc+qe)*qcDDt - 0.5*LB*LE*(me+2*mBox)*sin(qb-qe)*qbDDt  ...
- 0.5*LA*LD*(md+2*mBox+2*me)*sin(qa+qd)*qaDDt - 0.5*LC*LD*(md+2*mBox+2*me)*sin(qc+qd)*qcDDt - 0.5*LB*LD*(md+2*mBox+2*me)*sin(qb-qd)*qbDDt;
TED = g*LE*mBox*cos(qe) + 0.5*g*LE*me*cos(qe) + 0.25*(4*Ie+me*LE^2+4*mBox*LE^2)*qeDDt + 0.5*LD*LE*(me+2*mBox)*cos(qd-qe)*qdDDt  ...
- 0.5*LE*(me+2*mBox)*(LA*cos(qa+qe)*qaDt^2+LC*cos(qc+qe)*qcDt^2+LB*cos(qb-qe)*qbDt^2+LD*sin(qd-qe)*qdDt^2) - 0.5*LA*LE*(me+2*mBox)*  ...
sin(qa+qe)*qaDDt - 0.5*LC*LE*(me+2*mBox)*sin(qc+qe)*qcDDt - 0.5*LB*LE*(me+2*mBox)*sin(qb-qe)*qbDDt;

Output = zeros( 1, 12 );
Output(1) = t;
Output(2) = TAN;
Output(3) = TBA;
Output(4) = TCB;
Output(5) = TDC;
Output(6) = TED;

Output(7) = t;
Output(8) = qa*RADtoDEG;                              % Converted to deg
Output(9) = qb*RADtoDEG;                              % Converted to deg
Output(10) = qc*RADtoDEG;                             % Converted to deg
Output(11) = qd*RADtoDEG;                             % Converted to deg
Output(12) = qe*RADtoDEG;                             % Converted to deg
end


%===========================================================================
function OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile )
%===========================================================================
persistent FileIdentifier hasHeaderInformationBeenWritten;

if( isempty(Output) ),
   if( ~isempty(FileIdentifier) ),
      for( i = 1 : 2 ),  fclose( FileIdentifier(i) );  end
      clear FileIdentifier;
      fprintf( 1, '\n Output is in the files FInalComputation.i  (i=1,2)\n\n' );
   end
   clear hasHeaderInformationBeenWritten;
   return;
end

if( isempty(hasHeaderInformationBeenWritten) ),
   if( shouldPrintToScreen ),
      fprintf( 1,                '%%       t             tan            tba            tcb            tdc            ted\n' );
      fprintf( 1,                '%%     (sec)         (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)\n\n' );
   end
   if( shouldPrintToFile && isempty(FileIdentifier) ),
      FileIdentifier = zeros( 1, 2 );
      FileIdentifier(1) = fopen('FInalComputation.1', 'wt');   if( FileIdentifier(1) == -1 ), error('Error: unable to open file FInalComputation.1'); end
      fprintf(FileIdentifier(1), '%% FILE: FInalComputation.1\n%%\n' );
      fprintf(FileIdentifier(1), '%%       t             tan            tba            tcb            tdc            ted\n' );
      fprintf(FileIdentifier(1), '%%     (sec)         (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)\n\n' );
      FileIdentifier(2) = fopen('FInalComputation.2', 'wt');   if( FileIdentifier(2) == -1 ), error('Error: unable to open file FInalComputation.2'); end
      fprintf(FileIdentifier(2), '%% FILE: FInalComputation.2\n%%\n' );
      fprintf(FileIdentifier(2), '%%       t             qa             qb             qc             qd             qe\n' );
      fprintf(FileIdentifier(2), '%%     (sec)          (deg)          (deg)          (deg)          (deg)          (deg)\n\n' );
   end
   hasHeaderInformationBeenWritten = 1;
end

if( shouldPrintToScreen ), WriteNumericalData( 1,                 Output(1:6) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(1), Output(1:6) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(2), Output(7:12) );  end
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
data = load( 'FInalComputation.1' ); 
plot( data(:,1),data(:,2),'-b', data(:,1),data(:,3),'-.g', data(:,1),data(:,4),'--r', data(:,1),data(:,5),'-m', data(:,1),data(:,6),'-.k', 'LineWidth',3 );
legend( 'tan', 'tba', 'tcb', 'tdc', 'ted' );
xlabel('t (sec)');   % ylabel('Some y-axis label');   title('Some plot title');
clear data;

figure;
data = load( 'FInalComputation.2' ); 
plot( data(:,1),data(:,2),'-b', data(:,1),data(:,3),'-.g', data(:,1),data(:,4),'--r', data(:,1),data(:,5),'-m', data(:,1),data(:,6),'-.k', 'LineWidth',3 );
legend( 'qa (deg)', 'qb (deg)', 'qc (deg)', 'qd (deg)', 'qe (deg)' );
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


%=========================================
end    % End of function FInalComputation
%=========================================
