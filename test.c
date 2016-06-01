#include"../shared/include_all.h"
#include "../shared/branches.h"
#include "../shared/tagcorrection.h"
#include "../../Dss/BF/DssWeightsParam.h"
#include "../shared/BrCorrection.h"


void test(){
enum{signal_decay,charm_lep, charm_s_lep, charm_nonres_lep, charm_ss_lep, double_charm, charm_other, FakeLep, ulnu, rare,  otherB, badB, continuum, data, nContributions};

TFile *f = new TFile("~/rootFiles/smallFiles/charged_s0.root","read");
TTree *t = (TTree*)f->Get("tree");

t->SetBranchAddress("run", &run);
t->SetBranchAddress("btag", &btag);
t->SetBranchAddress("lep1", &lep1);
t->SetBranchAddress("gx", &gx);
t->SetBranchAddress("event", &event);
t->SetBranchAddress("truth", &truth);
t->SetBranchAddress("gmiss", &gmiss);
t->SetBranchAddress("sBDecay", &sBDecay);
t->SetBranchAddress("sDDecay", &sDDecay);
int nEvents = 0;
int nEventsBG = 0;

TH1F **hist = new TH1F*[2];
hist[0] = new TH1F("h1","h1",8,1,9);
hist[1] = new TH1F("h2","h2",8,1,9);

CorrectBf CorrBf = CorrectBf();

for(int i = 0; i< t->GetEntries(); i++)
{
	t->GetEntry(i);
	int iLep = lep1.fl_lep%10;
	
	
	
	int cont = -1;

	if(event.lclass==0)
		cont = badB;

	else if((lep1.fl_lep != 10 && lep1.fl_lep!=21)  ){
			cont = FakeLep;
	}
	else if(abs(event.lclass)==1 )
		cont = signal_decay;
	else if(abs(event.lclass)==2)
	{
		if(event.dclass==1)
			cont = charm_lep;
		else if (event.dclass ==3)
			cont = charm_s_lep;
		else if (event.dclass ==2 && event.dclass ==4)
			cont = charm_nonres_lep;
		else
			cont = charm_ss_lep;
	
	}	
	else if(abs(event.lclass)==4)
		cont = double_charm;
	else if(abs(event.lclass)>=8)
		cont = otherB;
	else if(abs(event.lclass)>4)
		cont = charm_other;
	else if( lep1.fl_mother != 3  )
		cont = rare;


	

	if( btag.m_bc<5.27) continue;
	if( btag.pcode_b*lep1.q>0) continue;
	//if( iLep == 0 && acos(lep1.th)/M_PI*180.<60 ) continue;
	if( log(btag.NB)<-4) continue;
	if( gmiss.m2<-5. ) continue;
	if( gmiss.m2>20. ) continue;
	//if( iLep == 0 && (lep1.ps<0.2||lep1.ps>2.5) ) continue;
	//if( iLep == 1 && (lep1.ps<0.4||lep1.ps>2.5) ) continue;
	
	if(event.cos_thrAm > 0.8) continue;
	//if(abs(event.lclass) <= 2) continue;
	hist[0]->Fill(event.nch,tagcorr(btag.b_mode, (double)btag.NB));
   	 hist[1]->Fill(event.nch);
	
	if(cont == signal_decay) nEvents++;
	else  nEventsBG++;
}

cout<<nEvents<<", "<<nEventsBG<<endl;
//hist[1]->Divide(hist[0]);
hist[0]->Divide(hist[1]);
hist[0]->Draw();
}
