
FILENAME DataFile '/home/markgewhite0/JoB2/Model Outputs v1.xlsx';
FILENAME RegFile '/home/markgewhite0/JoB2/Registration Outputs v1.xlsx';
FILENAME TFile '/home/markgewhite0/JoB2/TStat Outputs v1.xlsx';

FILENAME ParmFile '/home/markgewhite0/JoB2/Parameter Estimates.xlsx';

/* loading data */

proc import datafile=DataFile replace	DBMS=XLSX Out=PerfData;
	GetNames = Yes;
	Sheet = "Perf";
run;

proc import datafile=RegFile replace	DBMS=XLSX Out=RegData;
	GetNames = Yes;
	Sheet = "Registrations";
run;

proc import datafile=TFile replace	DBMS=XLSX Out=TStatData;
	GetNames = Yes;
	Sheet = "tStat";
run;


/* data processing */

data PerfData;
	set PerfData;
	Loglikelihood = 2000 - 2*loglikelihood;
run;

proc sort data=PerfData;
	by Outcome;
run;

data RegData;
	set RegData;
	where RSq>0 and Row not in (13, 72);
	LogPhaVar = log( PhaVar );
	LogAmpVar = log( AmpVar );
run;



data PerfData1;
	set PerfData;
	Partition = 'Train';
	RSq = TrainRSq;
	RMSE = TrainRMSE;
	Accuracy = TrainAccuracy;
run;

data PerfData2;
	set PerfData;
	Partition = 'Test';
	RSq = TestRSq;
	RMSE = TestRMSE;
	Accuracy = TestAccuracy;
run;

data OutcomeData;
	set PerfData1 PerfData2;
run;

proc sort data=OutcomeData;
	by Outcome Partition;
run;
	

/* Loglikelihood models */

proc mixed data=PerfData plots=studentpanel;
    by Outcome;
	class Norm (ref='PAD') LMReg (ref='0000') CTReg (ref='N') Predictor (ref='PCAU');
	model Loglikelihood = Predictor | Norm | LMReg | CTReg @2 / solution;
run;

ods pdf file = '/home/markgewhite0/JoB2/Loglikelihood Models.pdf';

proc glimmix data=PerfData plots=studentpanel ;
    by Outcome;
	class Norm (ref='PAD') LMReg (ref='0000') CTReg (ref='N') Predictor (ref='PCAU');
	model Loglikelihood = Predictor | Norm | LMReg | CTReg @2 / solution dist=gamma link=identity;
	ods output ParameterEstimates = JHParamEst;
run;

ods pdf close;



/* Registration models */

proc mixed data=RegData plots=studentpanel;
	class Norm (ref='PAD') LMReg (ref='0000') CTReg (ref='N');
	model AmpVar = Norm | LMReg | CTReg @2 / solution;
run;

proc mixed data=RegData plots=studentpanel;
	class Norm (ref='PAD') LMReg (ref='0000') CTReg (ref='N');
	model PhaVar = Norm | LMReg | CTReg @2 / solution;
run;

proc mixed data=RegData plots=studentpanel;
	class LM1-LM4 (ref='0') CTReg (ref='N');
	model RSq = LM1 | LM2 | LM3 | LM4 | CTReg / noint solution;
run;


/* T-statistic models */

proc mixed data=TStatData plots=studentpanel;
	class Norm (ref='PAD') LMReg (ref='0000') CTReg (ref='N') Predictor (ref='PCAU') Outcome (ref='JHtov');
	model T_U1 = Predictor | Norm | LMReg | CTReg | Outcome @2 / solution;
run;

proc mixed data=TStatData plots=studentpanel;
	class Norm (ref='PAD') LMReg (ref='0000') CTReg (ref='N') Predictor (ref='PCAU') Outcome (ref='JHtov');
	model T_U2 = Predictor | Norm | LMReg | CTReg | Outcome @2 / solution;
run;

proc mixed data=TStatData plots=studentpanel;
	class Norm (ref='PAD') LMReg (ref='0000') CTReg (ref='N') Predictor (ref='PCAU') Outcome (ref='JHtov');
	model T_U3 = Predictor | Norm | LMReg | CTReg | Outcome @2 / solution;
run;

proc mixed data=TStatData plots=studentpanel;
	class Norm (ref='PAD') LMReg (ref='0000') CTReg (ref='N') Predictor (ref='PCAU') Outcome (ref='JHtov');
	model T_UW1 = Predictor | Norm | LMReg | CTReg | Outcome @2 / solution;
run;

proc mixed data=TStatData plots=studentpanel;
	class Norm (ref='PAD') LMReg (ref='0000') CTReg (ref='N') Predictor (ref='PCAU') Outcome (ref='JHtov');
	model T_UW2 = Predictor | Norm | LMReg | CTReg | Outcome @2 / solution;
run;


/* Outcome variable models */

ods pdf file = '/home/markgewhite0/JoB2/Outcome Variable Models.pdf';

proc glimmix data=OutcomeData plots=studentpanel;
	by Outcome;
	where Outcome in ('JHtov' 'PP');
	class Norm (ref='PAD') LMReg (ref='0000') CTReg (ref='N') Predictor (ref='PCAU') Partition (ref='Train');
	model RMSE = Norm LMReg CTReg Predictor Partition 
				Norm*Partition LMReg*Partition CTReg*Partition Predictor*Partition
				/ solution dist=gamma link=identity;
run;

proc glimmix data=OutcomeData plots=studentpanel ;
	where Outcome='jumpType';
	class Norm (ref='PAD') LMReg (ref='0000') CTReg (ref='N') Predictor (ref='PCAU') Partition (ref='Train');
	model Accuracy = Norm LMReg CTReg Predictor Partition 
				Norm*Partition LMReg*Partition CTReg*Partition Predictor*Partition
				/ solution dist=gamma link=identity;
run;


ods pdf close;


