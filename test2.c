#include"../shared/include_all.h"
#include "../shared/branches.h"
#include "../shared/tagcorrection.h"
#include "../shared/BrCorrection.h"
#include "../shared/pid/LepFake.h"
#include "../shared/pid/LepEff.h"
#include "../shared/SemiLepWeights2D.h"
#include "../shared/SemiLepDssWeights2D.h"
#include <map>

using namespace std;
CorrectBf CorrBf = CorrectBf();

void test2(){

SemiLepWeights2D slWeight = SemiLepWeights2D();
SemiLepDssWeights2D slDssWeight = SemiLepDssWeights2D();
LepFake MuonFake1 = LepFake(1,"muon fake1");
LepFake ElecFake1 = LepFake(1,"elec fake1");

LepEff MuonEff1 = LepEff(1,"muon eff1",0.9,1);
LepEff ElecEff1 = LepEff(1,"elec eff1",0.9,0);


CorrectBf CorrBf = CorrectBf();

int dclass2lclass[12] = {0,2,0,1,0,3,4,5,6,0,0,0};



enum{signal_decay,charm_lep, charm_s_lep, charm_nonres_lep, charm_ss_lep, double_charm, charm_other, FakeLep, ulnu, rare,  otherB, badB, continuum, data, nContributions};

TFile *f = new TFile("~/rootFiles/smallFiles/DssMC.root","read");
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
hist[0] = new TH1F("h1","h1",20,0.2,2.3);
hist[1] = new TH1F("h2","h2",20,0.2,2.3);

std::map<std::string,TH1F*> mymap;
cout<<"start"<<endl;
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
	//if( gmiss.m2<0. ) continue;
	if( gmiss.m2>20. ) continue;
	//if( iLep == 0 && (lep1.ps<0.2||lep1.ps>2.5) ) continue;
	//if( iLep == 1 && (lep1.ps<0.4||lep1.ps>2.5) ) continue;
    if(lep1.ps<0.5) continue;
	
	if(event.cos_thrAm > 0.8) continue;
	
	nEvents++;
	
	auto DDecay = NoCharge(sDDecay);
	string BDecay = NoCharge(sBDecay);
	
	
	
		double MC_Correction = tagcorr(btag.b_mode, (double)btag.NB);
		//	if(it_f->Type == 'o' || it_f->Type == 'c')
	//			MC_Correction *= genMCCorr(run.exp);
	
	
	


		if(abs(event.lclass) == 2)
		{
			if(event.dclass<5)
			{
		       		MC_Correction*= slWeight.weight(dclass2lclass[event.dclass],truth.w, truth.costh,truth.q2,truth.p_l);
			}
			else
			{
				MC_Correction*= slDssWeight.weight(event.dclass-5,truth.w, truth.costh);

			}
		}
	//  MC_Correction *= BrCorrection(event.lclass, event.dclass, sBDecay, sDDecay);
		if(abs(event.lclass) == 2 )
		{
			MC_Correction*= CorrBf.weightB(event.lclass, event.dclass);
	//cout<<event.dclass<<" "<<CorrBf.weightB(event.lclass, event.dclass)<<endl;
	}
	else
	{
		MC_Correction*=CorrBf.weightB(sBDecay);

	}
	MC_Correction*= CorrBf.weightD(sDDecay);
	//cout<<CorrBf.weightD(sDDecay)<<endl;
	if(abs(event.lclass) == 2 && event.dclass >4){
		//if(it_f->Type == 's')
			MC_Correction*= CorrBf.FixDss(sDDecay);
	
		}

	//lep1
	//electrons
	if(lep1.fl_lep==10)//true e
		{MC_Correction*= ElecEff1.weight(lep1.fl_lep%10,run.exp,acos(lep1.th)/M_PI*180., lep1.p);}// cout<<Corrections[4]<<endl;}
	else if(lep1.fl_lep%10 == 0 )//fake e
		MC_Correction*= ElecFake1.weight(lep1.fl_lep%10,lep1.fl_lep/10-1,acos(lep1.th)/M_PI*180., lep1.p);

	//muons
	else if(lep1.fl_lep==21)//true mu
		MC_Correction*= MuonEff1.weight(lep1.fl_lep%10,run.exp,acos(lep1.th)/M_PI*180., lep1.p);
	else if(lep1.fl_lep%10 == 1 )//fake mu
		MC_Correction*= MuonFake1.weight(lep1.fl_lep%10,lep1.fl_lep/10-1,acos(lep1.th)/M_PI*180., lep1.p);
	else
		cout<<"Error: No Lepton\n";




	

    //BDecay[1] = 'x';
    int last = BDecay.length()-1;
  /*  BDecay[last ] = 'x';
    BDecay[last -3] = 'x';
    if(BDecay[last -7] !=3)
        BDecay[last -7] = 'x';
	*/
	//if(abs(event.lclass) !=2) continue;
	//if(event.dclass !=7) continue;
	
	vector<string> *DDec = CorrectBf::SplitDString(DDecay);
	for(int i = 0; i<DDec->size(); i++){
	
		//cout<<(*DDec)[i]<<endl;
		/*if((*DDec)[i].find("411")!= string::npos && ((*DDec)[i].find("11_12")!= string::npos || (*DDec)[i].find("13_14")!= string::npos)){
			event.dclass += 100;
			break;
		} else if((*DDec)[i].find("421")!= string::npos && ((*DDec)[i].find("11_12")!= string::npos || (*DDec)[i].find("13_14")!= string::npos)){
		event.dclass += 200;
			break;
		}*/
		if(/*(*DDec)[i].find("_411_221") == string::npos &&(*DDec)[i].find("_421_221") == string::npos &&(*DDec)[i].find("_10411_221") == string::npos &&(*DDec)[i].find("_10421_221") == string::npos
                */ // (*DDec)[i].find("100423") != 0 &&(*DDec)[i].find("100421") != 0 && (*DDec)[i].find("100413") != 0 &&(*DDec)[i].find("100411") != 0
                  (*DDec)[i].find("10423") != 0 && (*DDec)[i].find("10413") != 0) continue;
        //cout<<(*DDec)[i].find("221")<<" "<<BDecay<<endl;
		if(mymap.count((*DDec)[i]))
		{
			mymap[(*DDec)[i]]->Fill(gmiss.m2,MC_Correction);
		}else{
			mymap[(*DDec)[i]] = new TH1F("mm2","mm2",30,-5,15);
			mymap[(*DDec)[i]]->Fill(gmiss.m2,MC_Correction);
		}
	}
	delete DDec;
	/*
	if(mymap.count(BDecay))
	{
		mymap[BDecay]->Fill(gmiss.m2);
	}else{
		mymap[BDecay] = new TH1F("mm2","mm2",40,-5,15);
		mymap[BDecay]->Fill(gmiss.m2);
	}
	*/
}
cout<<nEvents<<endl;

TCanvas *cv = new TCanvas("cv","cv",600,600);
cv->Print("mm2.pdf[");
int count = 0;
for(auto it = mymap.begin(); it!=mymap.end(); ++it)
{
	//if(it->second->Integral()<10) continue;
    cout<<it->first<<" "<<it->second->Integral()<<endl;
	//if(it->second->GetMean()>1) continue;
	count++;
	if(count > 1000 )break;
	it->second->Draw();
	it->second->GetXaxis()->SetTitle(it->first.c_str());
	cv->Print("mm2.pdf");
}
cv->Print("mm2.pdf]");

cv->Clear();
cv->cd(0);

count = 0;
unsigned colors[] = {kBlack, kRed, kBlue, kGreen+1, kOrange};
TLegend leg(0.6,0.7,0.9,0.9);
for(auto it = mymap.begin(); it!=mymap.end(); ++it)
{
	if(it->second->Integral()<10) continue;
	//if(it->second->GetMean()>1) continue;
	
	if(count > 5 )break;
	it->second->SetLineColor(colors[count]);
	it->second->SetLineWidth(2);
	it->second->Scale(1./it->second->Integral());
	it->second->Draw(count>0?"same":"");
	if(count == 0) it->second->GetYaxis()->SetRangeUser(0, 1.2*it->second->GetMaximum());
	leg.AddEntry(it->second,it->first.c_str(),"l");
	count++;
	//it->second->GetXaxis()->SetTitle(it->first.c_str());
	
}
leg.Draw();
cv->Print("mm2_overlay.pdf");
}
