function [t,VAR,Output] = FinalComp
%===========================================================================
% File: FinalComp.m created Dec 07 2021 by MotionGenesis 6.1.
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
qa1=0; qa2=0; qb1=0; qb2=0; qc1=0; qc2=0; qd1=0; qd2=0; qa=0; qb=0; qc=0; qd=0; TAN=0; TBA=0; TCB=0; TDC=0; TED=0; qaDt=0; qbDt=0;
qcDt=0; qdDt=0; qeDt=0; qaDDt=0; qbDDt=0; qcDDt=0; qdDDt=0; Ia=0; Ib=0; Ic=0; Id=0; Ie=0; LACM=0; LAL=0; LBCM=0; LBL=0; LCCM=0; LCL=0;
LDCM=0; LDL=0; LECM=0; LEL=0; ma=0; mb=0; mc=0; md=0; me=0; rA=0; rB=0; rC=0; rD=0; rE=0;


%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
g                               =  9.81;                   % m/s^2               Constant
LA                              =  432.3;                  % mm                  Constant
LB                              =  368.5;                  % mm                  Constant
LC                              =  529.3;                  % mm                  Constant
LD                              =  275.1;                  % mm                  Constant
LE                              =  264.3;                  % mm                  Constant
LPACM                           =  .4416;                  % UNITS               Constant
LPBCM                           =  .3612;                  % UNITS               Constant
LPCCM                           =  .4151;                  % UNITS               Constant
LPDCM                           =  .5754;                  % UNITS               Constant
LPECM                           =  .4559;                  % UNITS               Constant
mbody                           =  61.9;                   % kg                  Constant
mBox                            =  5;                      % kg                  Constant
mpa                             =  .0481;                  % UNITS               Constant
mpb                             =  .1478;                  % UNITS               Constant
mpc                             =  .4257;                  % UNITS               Constant
mpd                             =  .0255;                  % UNITS               Constant
mpe                             =  .0138;                  % UNITS               Constant
rpa                             =  .093;                   % UNITS               Constant
rpb                             =  .162;                   % UNITS               Constant
rpc                             =  .171;                   % UNITS               Constant
rpd                             =  .148;                   % UNITS               Constant
rpe                             =  .094;                   % UNITS               Constant

qe                              = -62.69185061489564;      % deg                 Initial Value

tInitial                        =  0.0;                    % second              Initial Time
tFinal                          =  3;                      % s                   Final Time
tStep                           =  .01;                    % s                   Integration Step
printIntScreen                  =  1;                      % 0 or +integer       0 is NO screen output
printIntFile                    =  1;                      % 0 or +integer       0 is NO file   output
absError                        =  1.0E-05;                %                     Absolute Error
relError                        =  1.0E-08;                %                     Relative Error
%-------------------------------+--------------------------+-------------------+-----------------

% Unit conversions.  UnitSystem: kg, meter, second.
DEGtoRAD = pi / 180.0;
RADtoDEG = 180.0 / pi;
LA = LA * 0.001;                                           %  Converted from mm 
LB = LB * 0.001;                                           %  Converted from mm 
LC = LC * 0.001;                                           %  Converted from mm 
LD = LD * 0.001;                                           %  Converted from mm 
LE = LE * 0.001;                                           %  Converted from mm 
qe = qe * DEGtoRAD;                                        %  Converted from deg 

% Evaluate constants
qa1 = 1.396263401595464;
qa2 = 0.2617993877991495;
qb1 = 1.396263401595464;
qb2 = 0.3084996989990425;
qc1 = 0.8726646259971648;
qc2 = 0.174532925199433;
qd1 = -1.09417920739836;
qd2 = 0.6657969425451742;
LACM = LA*LPACM;
LAL = LA - LACM;
ma = mbody*mpa;
me = mbody*mpe;
LBCM = LB*LPBCM;
LBL = LB - LBCM;
LECM = LE*LPECM;
LEL = LE - LECM;
LCCM = LC*LPCCM;
mc = mbody*mpc;
LDCM = LD*LPDCM;
md = mbody*mpd;
rA = LA*rpa;
rB = LB*rpb;
rC = LC*rpc;
rD = LD*rpd;
rE = LE*rpe;
LCL = LC - LCCM;
LDL = LD - LDCM;
mb = mbody*mpb;
Ia = mbody*LA^2*ma*rA^2;
Ib = mbody*LB^2*mb*rB^2;
Ic = mbody*LC^2*mc*rC^2;
Id = mbody*LD^2*md*rD^2;
Ie = mbody*LE^2*me*rE^2;


VAR = SetMatrixFromNamedQuantities;
[t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, tStep, absError, relError, VAR, printIntScreen, printIntFile );
OutputToScreenOrFile( [], 0, 0 );   % Close output files.
if( printIntFile ~= 0 ),  PlotOutputFiles;  end


%===========================================================================
function sys = mdlDerivatives( t, VAR, uSimulink )
%===========================================================================
SetNamedQuantitiesFromMatrix( VAR );
qa = qa1 + 0.1591549430918954*(qa1-qa2)*sin(2.094395102393195*t) - 0.3333333333333333*(qa1-qa2)*t;
qaDt = 0.3333333333333333*qa2 + 0.3333333333333333*(qa1-qa2)*cos(2.094395102393195*t) - 0.3333333333333333*qa1;
qb = qb1 + 0.1591549430918954*(qb1-qb2)*sin(2.094395102393195*t) - 0.3333333333333333*(qb1-qb2)*t;
qbDt = 0.3333333333333333*qb2 + 0.3333333333333333*(qb1-qb2)*cos(2.094395102393195*t) - 0.3333333333333333*qb1;
qc = qc1 + 0.1591549430918954*(qc1-qc2)*sin(2.094395102393195*t) - 0.3333333333333333*(qc1-qc2)*t;
qcDt = 0.3333333333333333*qc2 + 0.3333333333333333*(qc1-qc2)*cos(2.094395102393195*t) - 0.3333333333333333*qc1;
qd = qd1 + 0.1591549430918954*(qd1-qd2)*sin(2.094395102393195*t) - 0.3333333333333333*(qd1-qd2)*t;
qdDt = 0.3333333333333333*qd2 + 0.3333333333333333*(qd1-qd2)*cos(2.094395102393195*t) - 0.3333333333333333*qd1;

qeDt = 0.4*(1+2.5*(LACM+LAL)*sin(qa)*qaDt+2.5*(LBCM+LBL)*sin(qb)*qbDt+2.5*(LCCM+LCL)*sin(qc)*qcDt-cos(2.094395102393195*t)-2.5*(LDCM+  ...
LDL)*cos(qd)*qdDt)/((LECM+LEL)*cos(qe));

sys = transpose( SetMatrixOfDerivativesPriorToIntegrationStep );
end



%===========================================================================
function VAR = SetMatrixFromNamedQuantities
%===========================================================================
VAR = zeros( 1, 1 );
VAR(1) = qe;
end


%===========================================================================
function SetNamedQuantitiesFromMatrix( VAR )
%===========================================================================
qe = VAR(1);
end


%===========================================================================
function VARp = SetMatrixOfDerivativesPriorToIntegrationStep
%===========================================================================
VARp = zeros( 1, 1 );
VARp(1) = qeDt;
end



%===========================================================================
function Output = mdlOutputs( t, VAR, uSimulink )
%===========================================================================
qaDDt = -0.6981317007977317*(qa1-qa2)*sin(2.094395102393195*t);
qbDDt = -0.6981317007977317*(qb1-qb2)*sin(2.094395102393195*t);
qcDDt = -0.6981317007977317*(qc1-qc2)*sin(2.094395102393195*t);
qdDDt = -0.6981317007977317*(qd1-qd2)*sin(2.094395102393195*t);
TAN = g*LACM*ma*sin(qa) + g*mb*(LACM+LAL)*sin(qa) + g*mc*(LACM+LAL)*sin(qa) + g*md*(LACM+LAL)*sin(qa) + 0.8377580409572779*(LACM+  ...
LAL)*(1.193662073189215*md*(LDCM*cos(qa+qd)*qdDt^2-(LBCM+LBL)*sin(qa+qb)*qbDt^2-(LCCM+LCL)*sin(qa-qc)*qcDt^2)-1.193662073189215*LBCM*  ...
mb*sin(qa+qb)*qbDt^2-1.193662073189215*mc*(LCCM*sin(qa-qc)*qcDt^2+(LBCM+LBL)*sin(qa+qb)*qbDt^2)-Ie*sin(qa)*(sin(2.094395102393195*  ...
t)+1.193662073189215*(LACM+LAL)*cos(qa)*qaDt^2+1.193662073189215*(LBCM+LBL)*cos(qb)*qbDt^2+1.193662073189215*(LCCM+LCL)*cos(qc)*qcDt^2+  ...
1.193662073189215*(LDCM+LDL)*sin(qd)*qdDt^2+1.193662073189215*(LECM+LEL)*sin(qe)*qeDt^2)/((LECM+LEL)^2*cos(qe)^2)-me*(1.193662073189216*(  ...
LBCM+LBL)*sin(qa+qb)*qbDt^2+1.193662073189216*(LCCM+LCL)*sin(qa-qc)*qcDt^2-1.193662073189216*(LDCM+LDL)*cos(qa+qd)*qdDt^2-1.193662073189216*  ...
LECM*cos(qa+qe)*qeDt^2-LECM*(1*sin(qa+qe)-LECM*sin(qa)/((LECM+LEL)*cos(qe)))*(sin(2.094395102393195*t)+1.193662073189215*(LACM+LAL)*  ...
cos(qa)*qaDt^2+1.193662073189215*(LBCM+LBL)*cos(qb)*qbDt^2+1.193662073189215*(LCCM+LCL)*cos(qc)*qcDt^2+1.193662073189215*(LDCM+LDL)*  ...
sin(qd)*qdDt^2+1.193662073189215*(LECM+LEL)*sin(qe)*qeDt^2)/((LECM+LEL)*cos(qe)))-mBox*(1.193662073189216*(LBCM+LBL)*sin(qa+qb)*qbDt^2+  ...
1.193662073189216*(LCCM+LCL)*sin(qa-qc)*qcDt^2-1.193662073189216*(LDCM+LDL)*cos(qa+qd)*qdDt^2-1.193662073189216*(LECM+LEL)*cos(qa+  ...
qe)*qeDt^2-(1*sin(qa+qe)*(sin(2.094395102393195*t)+1.193662073189215*(LACM+LAL)*cos(qa)*qaDt^2+1.193662073189215*(LBCM+LBL)*cos(qb)*qbDt^2+  ...
1.193662073189215*(LCCM+LCL)*cos(qc)*qcDt^2+1.193662073189215*(LDCM+LDL)*sin(qd)*qdDt^2+1.193662073189215*(LECM+LEL)*sin(qe)*qeDt^2)+  ...
sin(qa)*(1.193662073189215*(LACM+LAL)*cos(qa+qe)*qaDt^2+1.193662073189215*(LCCM+LCL)*cos(qc+qe)*qcDt^2+1.193662073189215*(LBCM+LBL)*  ...
cos(qb-qe)*qbDt^2+1.193662073189215*(LDCM+LDL)*sin(qd-qe)*qdDt^2-(sin(2.094395102393195*t)+1.193662073189215*(LACM+LAL)*cos(qa)*qaDt^2+  ...
1.193662073189215*(LBCM+LBL)*cos(qb)*qbDt^2+1.193662073189215*(LCCM+LCL)*cos(qc)*qcDt^2+1.193662073189215*(LDCM+LDL)*sin(qd)*qdDt^2+  ...
1.193662073189215*(LECM+LEL)*sin(qe)*qeDt^2)/cos(qe)))/cos(qe))) + 1*(LACM+LAL)*(1*LDCM*md*sin(qa+qd)+Ie*(LDCM+LDL)*sin(qa)*cos(qd)/((  ...
LECM+LEL)^2*cos(qe)^2)+mBox*(LDCM+LDL)*(1*sin(qa+qd)-(1*cos(qd)*sin(qa+qe)+sin(qa)*(cos(qd-qe)-cos(qd)/cos(qe)))/cos(qe))+me*(LDCM+  ...
LDL)*(1*sin(qa+qd)-LECM*cos(qd)*(1*sin(qa+qe)-LECM*sin(qa)/((LECM+LEL)*cos(qe)))/((LECM+LEL)*cos(qe))))*qdDDt + 0.8377580409572779*(  ...
LACM+LAL)*sin(qa)*(1.193662073189215*g*LECM*me*cos(qe)+(Ie+me*LECM^2)*(sin(2.094395102393195*t)+1.193662073189215*(LACM+LAL)*cos(  ...
qa)*qaDt^2+1.193662073189215*(LBCM+LBL)*cos(qb)*qbDt^2+1.193662073189215*(LCCM+LCL)*cos(qc)*qcDt^2+1.193662073189215*(LDCM+LDL)*sin(  ...
qd)*qdDt^2+1.193662073189215*(LECM+LEL)*sin(qe)*qeDt^2+1.193662073189215*(LACM+LAL)*sin(qa)*qaDDt+1.193662073189215*(LBCM+LBL)*sin(  ...
qb)*qbDDt+1.193662073189215*(LCCM+LCL)*sin(qc)*qcDDt-1.193662073189215*(LDCM+LDL)*cos(qd)*qdDDt)/((LECM+LEL)*cos(qe)))/((LECM+LEL)*  ...
cos(qe)) - g*me*(LACM+LAL)*(-1+LECM/(LECM+LEL))*sin(qa) - (LACM+LAL)*(1*LCCM*mc*cos(qa-qc)+1*md*(LCCM+LCL)*cos(qa-qc)+Ie*(LCCM+LCL)*  ...
sin(qa)*sin(qc)/((LECM+LEL)^2*cos(qe)^2)+mBox*(LCCM+LCL)*(1*cos(qa-qc)-(1*sin(qc)*sin(qa+qe)+sin(qa)*(sin(qc+qe)-sin(qc)/cos(qe)))/  ...
cos(qe))+me*(LCCM+LCL)*(1*cos(qa-qc)-LECM*sin(qc)*(1*sin(qa+qe)-LECM*sin(qa)/((LECM+LEL)*cos(qe)))/((LECM+LEL)*cos(qe))))*qcDDt  ...
- (1*Ia+1*ma*LACM^2+1*mb*(LACM+LAL)^2+1*mc*(LACM+LAL)^2+1*md*(LACM+LAL)^2+Ie*(LACM+LAL)^2*sin(qa)^2/((LECM+LEL)^2*cos(qe)^2)+mBox*(  ...
LACM+LAL)*(1*LACM+1*LAL+(LACM+LAL)*sin(qa)^2/cos(qe)^2-2*(LACM+LAL)*sin(qa)*sin(qa+qe)/cos(qe))-me*(LACM+LAL)*(LECM*(LACM+LAL)*sin(  ...
qa)*sin(qa+qe)/((LECM+LEL)*cos(qe))-LACM-LAL-LECM^2*(LACM+LAL)*sin(qa)^2/((LECM+LEL)^2*cos(qe)^2)))*qaDDt - (LACM+LAL)*(Ie*(LBCM+  ...
LBL)*sin(qa)*sin(qb)/((LECM+LEL)^2*cos(qe)^2)-LBCM*mb*cos(qa+qb)-mc*(LBCM+LBL)*cos(qa+qb)-md*(LBCM+LBL)*cos(qa+qb)-mBox*(LBCM+LBL)*(  ...
1*cos(qa+qb)+(1*sin(qb)*sin(qa+qe)+sin(qa)*(sin(qb-qe)-sin(qb)/cos(qe)))/cos(qe))-me*(LBCM+LBL)*(1*cos(qa+qb)+LECM*sin(qb)*(1*sin(  ...
qa+qe)-LECM*sin(qa)/((LECM+LEL)*cos(qe)))/((LECM+LEL)*cos(qe))))*qbDDt;
TBA = 1*g*me*(LBCM+LBL)*(-1+LECM/(LECM+LEL))*sin(qb) + LBCM*mb*(LACM+LAL)*sin(qa+qb)*qaDt^2 + mc*(LBCM+LBL)*(LCCM*sin(qb+qc)*qcDt^2+(  ...
LACM+LAL)*sin(qa+qb)*qaDt^2) + md*(LBCM+LBL)*(LDCM*cos(qb-qd)*qdDt^2+(LACM+LAL)*sin(qa+qb)*qaDt^2+(LCCM+LCL)*sin(qb+qc)*qcDt^2)  ...
+ 0.8377580409572779*Ie*(LBCM+LBL)*sin(qb)*(sin(2.094395102393195*t)+1.193662073189215*(LACM+LAL)*cos(qa)*qaDt^2+1.193662073189215*(  ...
LBCM+LBL)*cos(qb)*qbDt^2+1.193662073189215*(LCCM+LCL)*cos(qc)*qcDt^2+1.193662073189215*(LDCM+LDL)*sin(qd)*qdDt^2+1.193662073189215*(  ...
LECM+LEL)*sin(qe)*qeDt^2)/((LECM+LEL)^2*cos(qe)^2) + 0.8377580409572779*me*(LBCM+LBL)*(1.193662073189216*(LACM+LAL)*sin(qa+qb)*qaDt^2+  ...
1.193662073189216*(LCCM+LCL)*sin(qb+qc)*qcDt^2+1.193662073189216*(LDCM+LDL)*cos(qb-qd)*qdDt^2+1.193662073189216*LECM*cos(qb-qe)*qeDt^2-  ...
LECM*(1*sin(qb-qe)-LECM*sin(qb)/((LECM+LEL)*cos(qe)))*(sin(2.094395102393195*t)+1.193662073189215*(LACM+LAL)*cos(qa)*qaDt^2+1.193662073189215*(  ...
LBCM+LBL)*cos(qb)*qbDt^2+1.193662073189215*(LCCM+LCL)*cos(qc)*qcDt^2+1.193662073189215*(LDCM+LDL)*sin(qd)*qdDt^2+1.193662073189215*(  ...
LECM+LEL)*sin(qe)*qeDt^2)/((LECM+LEL)*cos(qe))) + 0.8377580409572779*mBox*(LBCM+LBL)*(1.193662073189216*(LACM+LAL)*sin(qa+qb)*qaDt^2+  ...
1.193662073189216*(LCCM+LCL)*sin(qb+qc)*qcDt^2+1.193662073189216*(LDCM+LDL)*cos(qb-qd)*qdDt^2+1.193662073189216*(LECM+LEL)*cos(qb-  ...
qe)*qeDt^2-(1*sin(qb-qe)*(sin(2.094395102393195*t)+1.193662073189215*(LACM+LAL)*cos(qa)*qaDt^2+1.193662073189215*(LBCM+LBL)*cos(qb)*qbDt^2+  ...
1.193662073189215*(LCCM+LCL)*cos(qc)*qcDt^2+1.193662073189215*(LDCM+LDL)*sin(qd)*qdDt^2+1.193662073189215*(LECM+LEL)*sin(qe)*qeDt^2)+  ...
sin(qb)*(1.193662073189215*(LACM+LAL)*cos(qa+qe)*qaDt^2+1.193662073189215*(LCCM+LCL)*cos(qc+qe)*qcDt^2+1.193662073189215*(LBCM+LBL)*  ...
cos(qb-qe)*qbDt^2+1.193662073189215*(LDCM+LDL)*sin(qd-qe)*qdDt^2-(sin(2.094395102393195*t)+1.193662073189215*(LACM+LAL)*cos(qa)*qaDt^2+  ...
1.193662073189215*(LBCM+LBL)*cos(qb)*qbDt^2+1.193662073189215*(LCCM+LCL)*cos(qc)*qcDt^2+1.193662073189215*(LDCM+LDL)*sin(qd)*qdDt^2+  ...
1.193662073189215*(LECM+LEL)*sin(qe)*qeDt^2)/cos(qe)))/cos(qe)) + 1*(LBCM+LBL)*(Ie*(LCCM+LCL)*sin(qb)*sin(qc)/((LECM+LEL)^2*cos(qe)^2)-  ...
LCCM*mc*cos(qb+qc)-md*(LCCM+LCL)*cos(qb+qc)-mBox*(LCCM+LCL)*(1*cos(qb+qc)+(1*sin(qc)*sin(qb-qe)+sin(qb)*(sin(qc+qe)-sin(qc)/cos(qe)))/  ...
cos(qe))-me*(LCCM+LCL)*(1*cos(qb+qc)+LECM*sin(qc)*(1*sin(qb-qe)-LECM*sin(qb)/((LECM+LEL)*cos(qe)))/((LECM+LEL)*cos(qe))))*qcDDt  ...
+ 1*(1*Ib+1*mb*LBCM^2+1*mc*(LBCM+LBL)^2+1*md*(LBCM+LBL)^2+Ie*(LBCM+LBL)^2*sin(qb)^2/((LECM+LEL)^2*cos(qe)^2)+mBox*(LBCM+LBL)*(1*LBCM+  ...
1*LBL+(LBCM+LBL)*sin(qb)^2/cos(qe)^2-2*(LBCM+LBL)*sin(qb)*sin(qb-qe)/cos(qe))+me*(LBCM+LBL)*(1*LBCM+1*LBL+LECM^2*(LBCM+LBL)*sin(qb)^2/((  ...
LECM+LEL)^2*cos(qe)^2)-LECM*(LBCM+LBL)*sin(qb)*sin(qb-qe)/((LECM+LEL)*cos(qe))))*qbDDt + 1*(LACM+LAL)*(Ie*(LBCM+LBL)*sin(qa)*sin(  ...
qb)/((LECM+LEL)^2*cos(qe)^2)-LBCM*mb*cos(qa+qb)-mc*(LBCM+LBL)*cos(qa+qb)-md*(LBCM+LBL)*cos(qa+qb)-mBox*(LBCM+LBL)*(1*cos(qa+qb)+(  ...
1*sin(qa)*sin(qb-qe)+sin(qb)*(sin(qa+qe)-sin(qa)/cos(qe)))/cos(qe))-me*(LBCM+LBL)*(1*cos(qa+qb)+LECM*sin(qa)*(1*sin(qb-qe)-LECM*sin(  ...
qb)/((LECM+LEL)*cos(qe)))/((LECM+LEL)*cos(qe))))*qaDDt - g*LBCM*mb*sin(qb) - g*mc*(LBCM+LBL)*sin(qb) - g*md*(LBCM+LBL)*sin(qb)  ...
- g*mBox*(1.110223024625157E-16*LBCM+1.110223024625157E-16*LBL)*sin(qb) - 0.8377580409572779*(LBCM+LBL)*sin(qb)*(1.193662073189215*  ...
g*LECM*me*cos(qe)+(Ie+me*LECM^2)*(sin(2.094395102393195*t)+1.193662073189215*(LACM+LAL)*cos(qa)*qaDt^2+1.193662073189215*(LBCM+LBL)*  ...
cos(qb)*qbDt^2+1.193662073189215*(LCCM+LCL)*cos(qc)*qcDt^2+1.193662073189215*(LDCM+LDL)*sin(qd)*qdDt^2+1.193662073189215*(LECM+LEL)*  ...
sin(qe)*qeDt^2+1.193662073189215*(LACM+LAL)*sin(qa)*qaDDt+1.193662073189215*(LBCM+LBL)*sin(qb)*qbDDt+1.193662073189215*(LCCM+LCL)*  ...
sin(qc)*qcDDt-1.193662073189215*(LDCM+LDL)*cos(qd)*qdDDt)/((LECM+LEL)*cos(qe)))/((LECM+LEL)*cos(qe)) - (LBCM+LBL)*(1*LDCM*md*sin(  ...
qb-qd)+Ie*(LDCM+LDL)*sin(qb)*cos(qd)/((LECM+LEL)^2*cos(qe)^2)+mBox*(LDCM+LDL)*(1*sin(qb-qd)-(1*cos(qd)*sin(qb-qe)+sin(qb)*(cos(qd-  ...
qe)-cos(qd)/cos(qe)))/cos(qe))+me*(LDCM+LDL)*(1*sin(qb-qd)-LECM*cos(qd)*(1*sin(qb-qe)-LECM*sin(qb)/((LECM+LEL)*cos(qe)))/((LECM+LEL)*  ...
cos(qe))))*qdDDt;
TCB = g*LCCM*mc*sin(qc) + g*md*(LCCM+LCL)*sin(qc) + 1*(LCCM+LCL)*(1*LDCM*md*sin(qc+qd)+Ie*(LDCM+LDL)*sin(qc)*cos(qd)/((LECM+LEL)^2*  ...
cos(qe)^2)+mBox*(LDCM+LDL)*(1*sin(qc+qd)-(1*cos(qd)*sin(qc+qe)+sin(qc)*(cos(qd-qe)-cos(qd)/cos(qe)))/cos(qe))+me*(LDCM+LDL)*(1*sin(  ...
qc+qd)-LECM*cos(qd)*(1*sin(qc+qe)-LECM*sin(qc)/((LECM+LEL)*cos(qe)))/((LECM+LEL)*cos(qe))))*qdDDt + 0.8377580409572779*(LCCM+LCL)*  ...
sin(qc)*(1.193662073189215*g*LECM*me*cos(qe)+(Ie+me*LECM^2)*(sin(2.094395102393195*t)+1.193662073189215*(LACM+LAL)*cos(qa)*qaDt^2+  ...
1.193662073189215*(LBCM+LBL)*cos(qb)*qbDt^2+1.193662073189215*(LCCM+LCL)*cos(qc)*qcDt^2+1.193662073189215*(LDCM+LDL)*sin(qd)*qdDt^2+  ...
1.193662073189215*(LECM+LEL)*sin(qe)*qeDt^2+1.193662073189215*(LACM+LAL)*sin(qa)*qaDDt+1.193662073189215*(LBCM+LBL)*sin(qb)*  ...
qbDDt+1.193662073189215*(LCCM+LCL)*sin(qc)*qcDDt-1.193662073189215*(LDCM+LDL)*cos(qd)*qdDDt)/((LECM+LEL)*cos(qe)))/((LECM+LEL)*cos(  ...
qe)) - g*me*(LCCM+LCL)*(-1+LECM/(LECM+LEL))*sin(qc) - LCCM*mc*((LBCM+LBL)*sin(qb+qc)*qbDt^2-(LACM+LAL)*sin(qa-qc)*qaDt^2)  ...
- md*(LCCM+LCL)*((LBCM+LBL)*sin(qb+qc)*qbDt^2-LDCM*cos(qc+qd)*qdDt^2-(LACM+LAL)*sin(qa-qc)*qaDt^2) - 0.8377580409572779*Ie*(LCCM+  ...
LCL)*sin(qc)*(sin(2.094395102393195*t)+1.193662073189215*(LACM+LAL)*cos(qa)*qaDt^2+1.193662073189215*(LBCM+LBL)*cos(qb)*qbDt^2+1.193662073189215*(  ...
LCCM+LCL)*cos(qc)*qcDt^2+1.193662073189215*(LDCM+LDL)*sin(qd)*qdDt^2+1.193662073189215*(LECM+LEL)*sin(qe)*qeDt^2)/((LECM+LEL)^2*cos(  ...
qe)^2) - 0.8377580409572779*me*(LCCM+LCL)*(1.193662073189216*(LBCM+LBL)*sin(qb+qc)*qbDt^2-1.193662073189216*(LDCM+LDL)*cos(qc+qd)*qdDt^2-  ...
1.193662073189216*(LACM+LAL)*sin(qa-qc)*qaDt^2-1.193662073189216*LECM*cos(qc+qe)*qeDt^2-LECM*(1*sin(qc+qe)-LECM*sin(qc)/((LECM+LEL)*  ...
cos(qe)))*(sin(2.094395102393195*t)+1.193662073189215*(LACM+LAL)*cos(qa)*qaDt^2+1.193662073189215*(LBCM+LBL)*cos(qb)*qbDt^2+1.193662073189215*(  ...
LCCM+LCL)*cos(qc)*qcDt^2+1.193662073189215*(LDCM+LDL)*sin(qd)*qdDt^2+1.193662073189215*(LECM+LEL)*sin(qe)*qeDt^2)/((LECM+LEL)*cos(  ...
qe))) - 0.8377580409572779*mBox*(LCCM+LCL)*(1.193662073189216*(LBCM+LBL)*sin(qb+qc)*qbDt^2-1.193662073189216*(LDCM+LDL)*cos(qc+qd)*qdDt^2-  ...
1.193662073189216*(LACM+LAL)*sin(qa-qc)*qaDt^2-1.193662073189216*(LECM+LEL)*cos(qc+qe)*qeDt^2-(1*sin(qc+qe)*(sin(2.094395102393195*  ...
t)+1.193662073189215*(LACM+LAL)*cos(qa)*qaDt^2+1.193662073189215*(LBCM+LBL)*cos(qb)*qbDt^2+1.193662073189215*(LCCM+LCL)*cos(qc)*qcDt^2+  ...
1.193662073189215*(LDCM+LDL)*sin(qd)*qdDt^2+1.193662073189215*(LECM+LEL)*sin(qe)*qeDt^2)+sin(qc)*(1.193662073189215*(LACM+LAL)*cos(  ...
qa+qe)*qaDt^2+1.193662073189215*(LCCM+LCL)*cos(qc+qe)*qcDt^2+1.193662073189215*(LBCM+LBL)*cos(qb-qe)*qbDt^2+1.193662073189215*(LDCM+  ...
LDL)*sin(qd-qe)*qdDt^2-(sin(2.094395102393195*t)+1.193662073189215*(LACM+LAL)*cos(qa)*qaDt^2+1.193662073189215*(LBCM+LBL)*cos(qb)*qbDt^2+  ...
1.193662073189215*(LCCM+LCL)*cos(qc)*qcDt^2+1.193662073189215*(LDCM+LDL)*sin(qd)*qdDt^2+1.193662073189215*(LECM+LEL)*sin(qe)*qeDt^2)/  ...
cos(qe)))/cos(qe)) - (1*Ic+1*mc*LCCM^2+1*md*(LCCM+LCL)^2+Ie*(LCCM+LCL)^2*sin(qc)^2/((LECM+LEL)^2*cos(qe)^2)+mBox*(LCCM+LCL)*(1*LCCM+  ...
1*LCL+(LCCM+LCL)*sin(qc)^2/cos(qe)^2-2*(LCCM+LCL)*sin(qc)*sin(qc+qe)/cos(qe))-me*(LCCM+LCL)*(LECM*(LCCM+LCL)*sin(qc)*sin(qc+qe)/((  ...
LECM+LEL)*cos(qe))-LCCM-LCL-LECM^2*(LCCM+LCL)*sin(qc)^2/((LECM+LEL)^2*cos(qe)^2)))*qcDDt - (LBCM+LBL)*(Ie*(LCCM+LCL)*sin(qb)*sin(  ...
qc)/((LECM+LEL)^2*cos(qe)^2)-LCCM*mc*cos(qb+qc)-md*(LCCM+LCL)*cos(qb+qc)-mBox*(LCCM+LCL)*(1*cos(qb+qc)+(1*sin(qb)*sin(qc+qe)+sin(  ...
qc)*(sin(qb-qe)-sin(qb)/cos(qe)))/cos(qe))-me*(LCCM+LCL)*(1*cos(qb+qc)+LECM*sin(qb)*(1*sin(qc+qe)-LECM*sin(qc)/((LECM+LEL)*cos(qe)))/((  ...
LECM+LEL)*cos(qe))))*qbDDt - (LACM+LAL)*(1*LCCM*mc*cos(qa-qc)+1*md*(LCCM+LCL)*cos(qa-qc)+Ie*(LCCM+LCL)*sin(qa)*sin(qc)/((LECM+LEL)^2*  ...
cos(qe)^2)+mBox*(LCCM+LCL)*(1*cos(qa-qc)-(1*sin(qa)*sin(qc+qe)+sin(qc)*(sin(qa+qe)-sin(qa)/cos(qe)))/cos(qe))+me*(LCCM+LCL)*(1*cos(  ...
qa-qc)-LECM*sin(qa)*(1*sin(qc+qe)-LECM*sin(qc)/((LECM+LEL)*cos(qe)))/((LECM+LEL)*cos(qe))))*qaDDt;
TDC = g*LDCM*md*cos(qd) + 1*(1*Id+1*md*LDCM^2+Ie*(LDCM+LDL)^2*cos(qd)^2/((LECM+LEL)^2*cos(qe)^2)+mBox*(LDCM+LDL)*(1*LDCM+1*LDL+(LDCM+  ...
LDL)*cos(qd)^2/cos(qe)^2-2*(LDCM+LDL)*cos(qd)*cos(qd-qe)/cos(qe))+me*(LDCM+LDL)*(1*LDCM+1*LDL+LECM^2*(LDCM+LDL)*cos(qd)^2/((LECM+  ...
LEL)^2*cos(qe)^2)-LECM*(LDCM+LDL)*cos(qd)*cos(qd-qe)/((LECM+LEL)*cos(qe))))*qdDDt + 0.8377580409572779*(LDCM+LDL)*cos(qd)*(1.193662073189215*  ...
g*LECM*me*cos(qe)+(Ie+me*LECM^2)*(sin(2.094395102393195*t)+1.193662073189215*(LACM+LAL)*cos(qa)*qaDt^2+1.193662073189215*(LBCM+LBL)*  ...
cos(qb)*qbDt^2+1.193662073189215*(LCCM+LCL)*cos(qc)*qcDt^2+1.193662073189215*(LDCM+LDL)*sin(qd)*qdDt^2+1.193662073189215*(LECM+LEL)*  ...
sin(qe)*qeDt^2+1.193662073189215*(LACM+LAL)*sin(qa)*qaDDt+1.193662073189215*(LBCM+LBL)*sin(qb)*qbDDt+1.193662073189215*(LCCM+LCL)*  ...
sin(qc)*qcDDt-1.193662073189215*(LDCM+LDL)*cos(qd)*qdDDt)/((LECM+LEL)*cos(qe)))/((LECM+LEL)*cos(qe)) - g*me*(LDCM+LDL)*(-1+LECM/(  ...
LECM+LEL))*cos(qd) - LDCM*md*((LACM+LAL)*cos(qa+qd)*qaDt^2+(LCCM+LCL)*cos(qc+qd)*qcDt^2+(LBCM+LBL)*cos(qb-qd)*qbDt^2) - 0.8377580409572779*  ...
Ie*(LDCM+LDL)*cos(qd)*(sin(2.094395102393195*t)+1.193662073189215*(LACM+LAL)*cos(qa)*qaDt^2+1.193662073189215*(LBCM+LBL)*cos(qb)*qbDt^2+  ...
1.193662073189215*(LCCM+LCL)*cos(qc)*qcDt^2+1.193662073189215*(LDCM+LDL)*sin(qd)*qdDt^2+1.193662073189215*(LECM+LEL)*sin(qe)*qeDt^2)/((  ...
LECM+LEL)^2*cos(qe)^2) - 0.8377580409572779*me*(LDCM+LDL)*(1.193662073189216*(LACM+LAL)*cos(qa+qd)*qaDt^2+1.193662073189216*(LCCM+  ...
LCL)*cos(qc+qd)*qcDt^2+1.193662073189216*(LBCM+LBL)*cos(qb-qd)*qbDt^2-1.193662073189216*LECM*sin(qd-qe)*qeDt^2-LECM*(1*cos(qd-qe)-  ...
LECM*cos(qd)/((LECM+LEL)*cos(qe)))*(sin(2.094395102393195*t)+1.193662073189215*(LACM+LAL)*cos(qa)*qaDt^2+1.193662073189215*(LBCM+  ...
LBL)*cos(qb)*qbDt^2+1.193662073189215*(LCCM+LCL)*cos(qc)*qcDt^2+1.193662073189215*(LDCM+LDL)*sin(qd)*qdDt^2+1.193662073189215*(LECM+  ...
LEL)*sin(qe)*qeDt^2)/((LECM+LEL)*cos(qe))) - 0.8377580409572779*mBox*(LDCM+LDL)*(1.193662073189216*(LACM+LAL)*cos(qa+qd)*qaDt^2+  ...
1.193662073189216*(LCCM+LCL)*cos(qc+qd)*qcDt^2+1.193662073189216*(LBCM+LBL)*cos(qb-qd)*qbDt^2-1.193662073189216*(LECM+LEL)*sin(qd-  ...
qe)*qeDt^2-(1*cos(qd-qe)*(sin(2.094395102393195*t)+1.193662073189215*(LACM+LAL)*cos(qa)*qaDt^2+1.193662073189215*(LBCM+LBL)*cos(qb)*qbDt^2+  ...
1.193662073189215*(LCCM+LCL)*cos(qc)*qcDt^2+1.193662073189215*(LDCM+LDL)*sin(qd)*qdDt^2+1.193662073189215*(LECM+LEL)*sin(qe)*qeDt^2)+  ...
cos(qd)*(1.193662073189215*(LACM+LAL)*cos(qa+qe)*qaDt^2+1.193662073189215*(LCCM+LCL)*cos(qc+qe)*qcDt^2+1.193662073189215*(LBCM+LBL)*  ...
cos(qb-qe)*qbDt^2+1.193662073189215*(LDCM+LDL)*sin(qd-qe)*qdDt^2-(sin(2.094395102393195*t)+1.193662073189215*(LACM+LAL)*cos(qa)*qaDt^2+  ...
1.193662073189215*(LBCM+LBL)*cos(qb)*qbDt^2+1.193662073189215*(LCCM+LCL)*cos(qc)*qcDt^2+1.193662073189215*(LDCM+LDL)*sin(qd)*qdDt^2+  ...
1.193662073189215*(LECM+LEL)*sin(qe)*qeDt^2)/cos(qe)))/cos(qe)) - (LACM+LAL)*(1*LDCM*md*sin(qa+qd)+Ie*(LDCM+LDL)*sin(qa)*cos(qd)/((  ...
LECM+LEL)^2*cos(qe)^2)+mBox*(LDCM+LDL)*(1*sin(qa+qd)-(1*sin(qa)*cos(qd-qe)+cos(qd)*(sin(qa+qe)-sin(qa)/cos(qe)))/cos(qe))+me*(LDCM+  ...
LDL)*(1*sin(qa+qd)-LECM*sin(qa)*(1*cos(qd-qe)-LECM*cos(qd)/((LECM+LEL)*cos(qe)))/((LECM+LEL)*cos(qe))))*qaDDt - (LCCM+LCL)*(1*LDCM*  ...
md*sin(qc+qd)+Ie*(LDCM+LDL)*sin(qc)*cos(qd)/((LECM+LEL)^2*cos(qe)^2)+mBox*(LDCM+LDL)*(1*sin(qc+qd)-(1*sin(qc)*cos(qd-qe)+cos(qd)*(  ...
sin(qc+qe)-sin(qc)/cos(qe)))/cos(qe))+me*(LDCM+LDL)*(1*sin(qc+qd)-LECM*sin(qc)*(1*cos(qd-qe)-LECM*cos(qd)/((LECM+LEL)*cos(qe)))/((  ...
LECM+LEL)*cos(qe))))*qcDDt - (LBCM+LBL)*(1*LDCM*md*sin(qb-qd)+Ie*(LDCM+LDL)*sin(qb)*cos(qd)/((LECM+LEL)^2*cos(qe)^2)+mBox*(LDCM+LDL)*(  ...
1*sin(qb-qd)-(1*sin(qb)*cos(qd-qe)+cos(qd)*(sin(qb-qe)-sin(qb)/cos(qe)))/cos(qe))+me*(LDCM+LDL)*(1*sin(qb-qd)-LECM*sin(qb)*(1*cos(  ...
qd-qe)-LECM*cos(qd)/((LECM+LEL)*cos(qe)))/((LECM+LEL)*cos(qe))))*qbDDt;
TED = g*LECM*me*cos(qe) + 0.837758040957278*(Ie+me*LECM^2)*(sin(2.094395102393195*t)+1.193662073189215*(LACM+LAL)*cos(qa)*qaDt^2+  ...
1.193662073189215*(LBCM+LBL)*cos(qb)*qbDt^2+1.193662073189215*(LCCM+LCL)*cos(qc)*qcDt^2+1.193662073189215*(LDCM+LDL)*sin(qd)*qdDt^2+  ...
1.193662073189215*(LECM+LEL)*sin(qe)*qeDt^2+1.193662073189215*(LACM+LAL)*sin(qa)*qaDDt+1.193662073189215*(LBCM+LBL)*sin(qb)*  ...
qbDDt+1.193662073189215*(LCCM+LCL)*sin(qc)*qcDDt-1.193662073189215*(LDCM+LDL)*cos(qd)*qdDDt)/((LECM+LEL)*cos(qe)) - LECM*me*(LACM+  ...
LAL)*(cos(qa+qe)*qaDt^2+sin(qa+qe)*qaDDt) - LECM*me*(LCCM+LCL)*(cos(qc+qe)*qcDt^2+sin(qc+qe)*qcDDt) - LECM*me*(LBCM+LBL)*(cos(qb-  ...
qe)*qbDt^2+sin(qb-qe)*qbDDt) - LECM*me*(LDCM+LDL)*(sin(qd-qe)*qdDt^2-cos(qd-qe)*qdDDt);

Output = zeros( 1, 12 );
Output(1) = t;
Output(2) = qa*RADtoDEG;                              % Converted to deg
Output(3) = qb*RADtoDEG;                              % Converted to deg
Output(4) = qc*RADtoDEG;                              % Converted to deg
Output(5) = qd*RADtoDEG;                              % Converted to deg
Output(6) = qe*RADtoDEG;                              % Converted to deg

Output(7) = t;
Output(8) = TAN;
Output(9) = TBA;
Output(10) = TCB;
Output(11) = TDC;
Output(12) = TED;
end


%===========================================================================
function OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile )
%===========================================================================
persistent FileIdentifier hasHeaderInformationBeenWritten;

if( isempty(Output) ),
   if( ~isempty(FileIdentifier) ),
      for( i = 1 : 2 ),  fclose( FileIdentifier(i) );  end
      clear FileIdentifier;
      fprintf( 1, '\n Output is in the files FinalComp.i  (i=1,2)\n\n' );
   end
   clear hasHeaderInformationBeenWritten;
   return;
end

if( isempty(hasHeaderInformationBeenWritten) ),
   if( shouldPrintToScreen ),
      fprintf( 1,                '%%       t             qa             qb             qc             qd             qe\n' );
      fprintf( 1,                '%%     (sec)          (deg)          (deg)          (deg)          (deg)          (deg)\n\n' );
   end
   if( shouldPrintToFile && isempty(FileIdentifier) ),
      FileIdentifier = zeros( 1, 2 );
      FileIdentifier(1) = fopen('FinalComp.1', 'wt');   if( FileIdentifier(1) == -1 ), error('Error: unable to open file FinalComp.1'); end
      fprintf(FileIdentifier(1), '%% FILE: FinalComp.1\n%%\n' );
      fprintf(FileIdentifier(1), '%%       t             qa             qb             qc             qd             qe\n' );
      fprintf(FileIdentifier(1), '%%     (sec)          (deg)          (deg)          (deg)          (deg)          (deg)\n\n' );
      FileIdentifier(2) = fopen('FinalComp.2', 'wt');   if( FileIdentifier(2) == -1 ), error('Error: unable to open file FinalComp.2'); end
      fprintf(FileIdentifier(2), '%% FILE: FinalComp.2\n%%\n' );
      fprintf(FileIdentifier(2), '%%       t             TAN            TBA            TCB            TDC            TED\n' );
      fprintf(FileIdentifier(2), '%%     (sec)         (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)\n\n' );
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
data = load( 'FinalComp.1' ); 
plot( data(:,1),data(:,2),'-b', data(:,1),data(:,3),'-.g', data(:,1),data(:,4),'--r', data(:,1),data(:,5),'-m', data(:,1),data(:,6),'-.k', 'LineWidth',3 );
legend( 'qa (deg)', 'qb (deg)', 'qc (deg)', 'qd (deg)', 'qe (deg)' );
xlabel('t (sec)');   % ylabel('Some y-axis label');   title('Some plot title');
clear data;

figure;
data = load( 'FinalComp.2' ); 
plot( data(:,1),data(:,2),'-b', data(:,1),data(:,3),'-.g', data(:,1),data(:,4),'--r', data(:,1),data(:,5),'-m', data(:,1),data(:,6),'-.k', 'LineWidth',3 );
legend( 'TAN', 'TBA', 'TCB', 'TDC', 'TED' );
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
end    % End of function FinalComp
%==================================
