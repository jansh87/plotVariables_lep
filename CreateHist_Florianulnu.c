
#define _BrLog 1
#define nbins_const 80
#define _pilep 0
#define _drawData 1
#define Scale_MC 0
#define _newTree 0


#include"../shared/include_all.h"
#include "../shared/CVar.h"
#include "../shared/CFile.h"
#include "../shared/CCont.h"
#include "../shared/tagcorrection.h"
//#include "../shared/q2_correction.h"

#include "../shared/branches.h"
#include "../shared/MyParticle.h"
#include "../shared/BrCorrection.h"
#include "../shared/SemiLepWeights2D.h"
#include "../shared/SemiLepDssWeights2D.h"
#include "../shared/pid/LepFake.h"
#include "../shared/pid/LepEff.h"
#include "../shared/pid/KPiEff.h"
#include <algorithm>
enum{signal_decay,charm_lep, charm_s_lep, charm_nonres_lep, charm_ss_lep, double_charm, charm_other, FakeLep, ulnu, rare,  otherB, badB, continuum, data, nContributions};


int sign( int val )
{ if(val>=0) return 1; else return -1;}


float sqr( float x)
{return x*x;}
double sqr(double x)
{return x*x; }


std::vector<CVar> vVar;
float mysqrt(float x)
{
	if(x>=0) return sqrt(x);
	else return -sqrt(-x);
}

TH1F* newHist(char* sName, CVar v1)
{
	TH1F *out;

			if(v1.xBins == NULL )
				out =  new TH1F(sName,sName,v1.nBins,v1.xMin,v1.xMax);
			else
				out =  new TH1F(sName,sName,v1.nBins,v1.xBins);
			out->GetXaxis()->SetTitle(v1.xTitle.c_str());
			out->GetYaxis()->SetTitle("Events");
			return out;
}

std::vector<CCont> vCont;
std::vector<CFile> vFile;
int nVariables;
int nExtraVariables;



struct t_momInfo {float e; float m2; float pt; float p3; float p; float costh;} t_gamma,t_track,t_btag,t_PmisswTag,t_PmissWoGamma, t_Vis;
float fMy4s;
float ll_theta, fnkss, fnkos, fnk, fq, BDT1,BDT2, pl1_bdt, fnch,fnn;
float flclass,fl_lep1, fl_lep2, fl_mother1, fl_mother2, fl_mother12;
float sqrt_gxm2, sqrt_xlepm2;
float test;
float vis_miss;
float cosgmissthc, coslep1thc,coslep2thc;
float fNFS, fb_mode, fexp;
float fpi0;
float frunno;
float vtx2_vtx1;
float fmm;
float fEgamma, fPgamma;
float fMD0Cand;
float fMDpCand;
float fMOmegaCand;
float fMKK,fMKPi,fMPiPi=-1;
float gmissm2Trans;
float fThLepMiss,fThLepX, fThXMiss;
float px1, px2;
float SumPsig;
int nVarReconMode;
int nVarPs1Elec;
int nVarPs2Elec;
int iVarMx,iVarMb,igPmiss;
int ihistKK;
int useDilep = 0;
float xps1el[] = {.3,	0.5,	0.7, 	0.9,	1.1, 	1.3,	1.5,	1.7,	 	2.5};
float xps1mu[] = {	0.5,	0.7, 	0.9,	1.1, 	1.3,	1.5,	1.7,	 	2.5};
float Corrections[8];
void DefVariables(){

	//variables
	CVar tmp;
	vVar.clear();

	/*
	tmp = CVar("Mbc","Mbc/GeV",&btag.m_bc, 60,5.24, 5.29);
	vVar.push_back(tmp);
	tmp = CVar("Mbc1","Mbc/GeV",&btag.m_bc, 60,5.24, 5.29);
	vVar.push_back(tmp);
	tmp = CVar("Mbc2","Mbc/GeV",&btag.m_bc, 60,5.24, 5.29);
	vVar.push_back(tmp);
	tmp = CVar("Mbc3","Mbc/GeV",&btag.m_bc, 60,5.24, 5.29);
	vVar.push_back(tmp);
	tmp = CVar("Mbc4","Mbc/GeV",&btag.m_bc, 60,5.24, 5.29);
	vVar.push_back(tmp);
	tmp = CVar("Mbc5","Mbc/GeV",&btag.m_bc, 60,5.24, 5.29);
	vVar.push_back(tmp);
*/

/*
	tmp = CVar("Mx0","MX/GeV",&gx.m2, 60,0, 10);
	vVar.push_back(tmp);
	tmp = CVar("Mx1","MX/GeV",&gx.m2, 60,0, 10);
	vVar.push_back(tmp);
	tmp = CVar("Mx2","MX/GeV",&gx.m2, 60,0, 10);
	vVar.push_back(tmp);
	tmp = CVar("Mx3","MX/GeV",&gx.m2, 60,0, 10);
	vVar.push_back(tmp);
	tmp = CVar("Mx4","MX/GeV",&gx.m2, 60,0, 10);
	vVar.push_back(tmp);
	tmp = CVar("Mx5","MX/GeV",&gx.m2, 60,0, 10);
	vVar.push_back(tmp);
	*//*
	tmp = CVar("Mx6","MX/GeV",&gx.m2, 60,0, 10);
	vVar.push_back(tmp);
	tmp = CVar("Mx7","MX/GeV",&gx.m2, 60,0, 10);
	vVar.push_back(tmp);
	tmp = CVar("Mx8","MX/GeV",&gx.m2, 60,0, 10);
	vVar.push_back(tmp);
	tmp = CVar("Mx9","MX/GeV",&gx.m2, 60,0, 10);
	vVar.push_back(tmp);
	tmp = CVar("Mx10","MX/GeV",&gx.m2, 60,0, 10);
	vVar.push_back(tmp);
	tmp = CVar("Mx11","MX/GeV",&gx.m2, 60,0, 10);
	vVar.push_back(tmp);
	*/

	//tmp = CVar("MM2","MM2/GeV2",&gmiss.m2, 60,-5, 20);
	//vVar.push_back(tmp);
	//tmp = CVar("MX","MX/GeV",&gx.m2, 60,0, 6);
	//vVar.push_back(tmp);



//	tmp = CVar("P track0","P track/GeV",&miss.m2, 60,-2, 2);
//	vVar.push_back(tmp);
//	tmp = CVar("P track1","P track/GeV",&miss.m2, 60,-2, 2);
//	vVar.push_back(tmp);
//	tmp = CVar("P track2","P track/GeV",&miss.m2, 60,-2, 2);
//	vVar.push_back(tmp);
//	tmp = CVar("P track3","P track/GeV",&miss.m2, 60,-2, 2);
//	vVar.push_back(tmp);
//	tmp = CVar("P track4","P track/GeV",&miss.m2, 60,-2, 2);
//	vVar.push_back(tmp);
//	tmp = CVar("P track5","P track/GeV",&miss.m2, 60,-2, 2);
//	vVar.push_back(tmp);
//
//	tmp = CVar("P track6","P track/GeV",&miss.m2, 60,-2, 2);
//	vVar.push_back(tmp);
//	tmp = CVar("P track7","P track/GeV",&miss.m2, 60,-2, 2);
//	vVar.push_back(tmp);
//	tmp = CVar("P track8","P track/GeV",&miss.m2, 60,-2, 2);
//	vVar.push_back(tmp);
//	tmp = CVar("P track9","P track/GeV",&miss.m2, 60,-2, 2);
//	vVar.push_back(tmp);
//	tmp = CVar("P track10","P track/GeV",&miss.m2, 60,-2, 2);
//	vVar.push_back(tmp);
//	tmp = CVar("P track11","P track/GeV",&miss.m2, 60,-2, 2);
//	vVar.push_back(tmp);

	tmp = CVar("P track0","P track/GeV",&miss.m2, 60,0, 2);
	vVar.push_back(tmp);
	tmp = CVar("P track1","P track/GeV",&miss.m2, 60,0, 2);
	vVar.push_back(tmp);
	tmp = CVar("P track2","P track/GeV",&miss.m2, 60,0, 2);
	vVar.push_back(tmp);
	tmp = CVar("P track3","P track/GeV",&miss.m2, 60,0, 2);
	vVar.push_back(tmp);
	tmp = CVar("P track4","P track/GeV",&miss.m2, 60,0, 2);
	vVar.push_back(tmp);
	tmp = CVar("P track5","P track/GeV",&miss.m2, 60,0, 2);
	vVar.push_back(tmp);

	tmp = CVar("P track6","P track/GeV",&miss.m2, 60,0, 2);
	vVar.push_back(tmp);
	tmp = CVar("P track7","P track/GeV",&miss.m2, 60,0, 2);
	vVar.push_back(tmp);
	tmp = CVar("P track8","P track/GeV",&miss.m2, 60,0, 2);
	vVar.push_back(tmp);
	tmp = CVar("P track9","P track/GeV",&miss.m2, 60,0, 2);
	vVar.push_back(tmp);
	tmp = CVar("P track10","P track/GeV",&miss.m2, 60,0, 2);
	vVar.push_back(tmp);
	tmp = CVar("P track11","P track/GeV",&miss.m2, 60,0, 2);
	vVar.push_back(tmp);
/*
	tmp = CVar("MM20","MM2/GeV2",&miss.m2, 60,-15, 15);
	vVar.push_back(tmp);
	tmp = CVar("MM21","MM2/GeV2",&miss.m2, 60,-15, 15);
	vVar.push_back(tmp);
	tmp = CVar("MM22","MM2/GeV2",&miss.m2, 60,-15, 15);
	vVar.push_back(tmp);
	tmp = CVar("MM23","MM2/GeV2",&miss.m2, 60,-15, 15);
	vVar.push_back(tmp);
	tmp = CVar("MM24","MM2/GeV2",&miss.m2, 60,-15, 15);
	vVar.push_back(tmp);
	tmp = CVar("MM25","MM2/GeV2",&miss.m2, 60,-15, 15);
	vVar.push_back(tmp);

	tmp = CVar("MM26","MM2/GeV2",&miss.m2, 60,-15, 15);
	vVar.push_back(tmp);
	tmp = CVar("MM27","MM2/GeV2",&miss.m2, 60,-15, 15);
	vVar.push_back(tmp);
	tmp = CVar("MM28","MM2/GeV2",&miss.m2, 60,-15, 15);
	vVar.push_back(tmp);
	tmp = CVar("MM29","MM2/GeV2",&miss.m2, 60,-15, 15);
	vVar.push_back(tmp);
	tmp = CVar("MM210","MM2/GeV2",&miss.m2, 60,-15, 15);
	vVar.push_back(tmp);
	tmp = CVar("MM211","MM2/GeV2",&miss.m2, 60,-15, 15);
	vVar.push_back(tmp);
	*/
	tmp = CVar("MM20","MM2/GeV2",&miss.m2, 60,0, 4);
	vVar.push_back(tmp);
	tmp = CVar("MM21","MM2/GeV2",&miss.m2, 60,0, 4);
	vVar.push_back(tmp);
	tmp = CVar("MM22","MM2/GeV2",&miss.m2, 60,0, 4);
	vVar.push_back(tmp);
	tmp = CVar("MM23","MM2/GeV2",&miss.m2, 60,0, 4);
	vVar.push_back(tmp);
	tmp = CVar("MM24","MM2/GeV2",&miss.m2, 60,0, 4);
	vVar.push_back(tmp);
	tmp = CVar("MM25","MM2/GeV2",&miss.m2, 60,0, 4);
	vVar.push_back(tmp);

	tmp = CVar("MM26","MM2/GeV2",&miss.m2, 60,0, 4);
	vVar.push_back(tmp);
	tmp = CVar("MM27","MM2/GeV2",&miss.m2, 60,0, 4);
	vVar.push_back(tmp);
	tmp = CVar("MM28","MM2/GeV2",&miss.m2, 60,0, 4);
	vVar.push_back(tmp);
	tmp = CVar("MM29","MM2/GeV2",&miss.m2, 60,0, 4);
	vVar.push_back(tmp);
	tmp = CVar("MM210","MM2/GeV2",&miss.m2, 60,0, 4);
	vVar.push_back(tmp);
	tmp = CVar("MM211","MM2/GeV2",&miss.m2, 60,0, 4);
	vVar.push_back(tmp);

	/*
	tmp = CVar("Ptag20","P_{tag}/GeV",&btag.p, 60,1.7, 3);
	vVar.push_back(tmp);
	tmp = CVar("Ptag21","P_{tag}/GeV",&btag.p, 60,1.7, 3);
	vVar.push_back(tmp);
	tmp = CVar("Ptag22","P_{tag}/GeV",&btag.p, 60,1.7, 3);
	vVar.push_back(tmp);
	tmp = CVar("Ptag23","P_{tag}/GeV",&btag.p, 60,1.7, 3);
	vVar.push_back(tmp);
	tmp = CVar("Ptag24","P_{tag}/GeV",&btag.p, 60,1.7, 3);
	vVar.push_back(tmp);
	tmp = CVar("Ptag25","P_{tag}/GeV",&btag.p, 60,1.7, 3);
	vVar.push_back(tmp);
	tmp = CVar("Ptag26","P_{tag}/GeV",&btag.p, 60,1.7, 3);
	vVar.push_back(tmp);
	tmp = CVar("Ptag27","P_{tag}/GeV",&btag.p, 60,1.7, 3);
	vVar.push_back(tmp);
	tmp = CVar("Ptag28","P_{tag}/GeV",&btag.p, 60,1.7, 3);
	vVar.push_back(tmp);
	tmp = CVar("Ptag29","P_{tag}/GeV",&btag.p, 60,1.7, 3);
	vVar.push_back(tmp);
	tmp = CVar("Ptag210","P_{tag}/GeV",&btag.p, 60,1.7, 3);
	vVar.push_back(tmp);
	tmp = CVar("Ptag211","P_{tag}/GeV",&btag.p, 60,1.7, 3);
	vVar.push_back(tmp);

	*/

	tmp = CVar("M_{bc} / GeV",&btag.m_bc, 5.26, 5.29);
	vVar.push_back(tmp);

	tmp = CVar("contNB",&btag.contNB, -5., 0.);
	vVar.push_back(tmp);

	tmp = CVar("NB",&btag.NB, -6., 0.);
	vVar.push_back(tmp);

	tmp = CVar("deltaE","#Delta E / GeV",&btag.delta_e, -0.18, 0.15);
	vVar.push_back(tmp);

	tmp = CVar("Ebeam","E_{CMS} / GeV",&run.etot, 20,11.49, 11.5);
	vVar.push_back(tmp);
	tmp = CVar("Eher","E_{HER} / GeV",&run.eher, 20,7.99,8. );
	vVar.push_back(tmp);
	tmp = CVar("Eler","E_{LER} / GeV",&run.eler, 20,3.495, 3.5);
	vVar.push_back(tmp);
	tmp = CVar("Mtot","M_{Y4S} / GeV",&fMy4s, 10.57, 10.58);
	vVar.push_back(tmp);



	tmp = CVar("runno","run number", &frunno, 0, 3000);
	vVar.push_back(tmp);


	tmp = CVar("btagcosthr","cos(#theta_{thr})_{tag}",&btag.cos_thr, -1, 1);
	//vVar.push_back(tmp);

	tmp = CVar("tag thrust",&btag.thr,0.6, 1);
	//vVar.push_back(tmp);

	tmp = CVar("ptag","p_{tag} / GeV",&t_btag.p,0, 0.5);
	vVar.push_back(tmp);

	tmp = CVar("tagpz","p_{z} / GeV",&t_btag.p3, -1,1);
	vVar.push_back(tmp);

	tmp = CVar("tagptrans","p_{tag trans} / GeV",&t_btag.pt,0.,0.5 );
	vVar.push_back(tmp);

	tmp = CVar("tagpe","e_{tag} / GeV",&t_btag.e,5.2, 5.4);
	vVar.push_back(tmp);

	tmp = CVar("tagmb","M_{tag} / GeV",&t_btag.m2,5.24, 5.4);
	vVar.push_back(tmp);

	tmp = CVar("tagth","th_{tag} ",&t_btag.costh,-1, 1);
	vVar.push_back(tmp);

	tmp = CVar("tag NFS",&fNFS,12,1.5, 13.5);
	vVar.push_back(tmp);

	tmp = CVar("B+tag recon mode",&fb_mode,17,20.5, 37.5);
	vVar.push_back(tmp);
	nVarReconMode = vVar.size() -1;

	tmp = CVar("B0tag recon mode",&fb_mode,15,0.5, 15.5);
	vVar.push_back(tmp);

	tmp = CVar("ExpNum",&fexp,65-6,6.5, 65.5);
	vVar.push_back(tmp);

	tmp = CVar("qtot","|q_{tot}|",&fq,6,-0.5,5.5);
	vVar.push_back(tmp);



	tmp = CVar("event.cos_thrA",&event.cos_thrA,0.,1);
	vVar.push_back(tmp);

	tmp = CVar("event.cos_thrAm",&event.cos_thrAm,0.,1);
	vVar.push_back(tmp);

	tmp = CVar("event.cos_thrB",&event.cos_thrB,0.,1);
	vVar.push_back(tmp);

	tmp = CVar("thrA2",&event.cos_thrA2,-1,1);
	vVar.push_back(tmp);

	tmp = CVar("thrC",&event.cos_thrC,-1,1);
	vVar.push_back(tmp);

	tmp = CVar("event.thrX",&event.thrX,0.5,1);
	vVar.push_back(tmp);

	tmp = CVar("event.thrSigX",&event.thrSigX,0.5,1);
	vVar.push_back(tmp);

	tmp = CVar("event.thrSigXm",&event.thrSigXm,0.5,1);
	vVar.push_back(tmp);

	tmp = CVar("thrSig2X",&event.thrSig2X,0.5,1);
	vVar.push_back(tmp);

	//tmp = CVar("ps1","p*_{#tau} / GeV",&lep1.ps, 0, 2.5);
	tmp = CVar("ps1","p*_{#tau} / GeV",&lep1.ps, 7,xps1mu);
//	vVar.push_back(tmp);

	tmp = CVar("p1","p_{#tau} / GeV",&lep1.p, 0.3, 2.5);
//	vVar.push_back(tmp);


	tmp = CVar("ps1e","electron p*_{#tau} / GeV",&lep1.ps, 0.3, 2.5);
	vVar.push_back(tmp);
	nVarPs1Elec = vVar.size();

	tmp = CVar("ps1mu","muon p*_{#tau} / GeV",&lep1.ps, 0.5, 2.5);
	vVar.push_back(tmp);

	/*
	tmp = CVar("ps1e","electron p*_{#tau} / GeV",&lep1.ps, 0, 2.5);
	vVar.push_back(tmp);
	nVarPs2Elec = vVar.size();

	tmp = CVar("ps1mu","muon p*_{#tau} / GeV",&lep1.ps, 0, 2.5);
	vVar.push_back(tmp);
*/
	tmp = CVar("labcostheta","lab cos#theta_{#tau}",&lep1.th, -1, 1);
	vVar.push_back(tmp);

	tmp = CVar("mm2","M_{miss}^{2} / GeV^{2}",&miss.m2, -5, 20);
	vVar.push_back(tmp);

	tmp = CVar("emiss","E*_{miss} / GeV",&miss.es, 0, 2.5);
	vVar.push_back(tmp);

	tmp = CVar("pmiss","p*_{miss} / GeV",&miss.ps, 0, 2.5);
	vVar.push_back(tmp);

	tmp = CVar("mm2g","M_{miss}^{2} / GeV^{2}",&gmiss.m2, 60,-5, 15);
	vVar.push_back(tmp);

	tmp = CVar("mm2gnPeak","M_{miss}^{2} / GeV^{2} no Peak",&gmiss.m2, 60,-5, 15);
	vVar.push_back(tmp);

	tmp = CVar("mm2gtrans","M_{miss Trans}^{2} / GeV^{2}",&gmissm2Trans,60, -5, 15);
	vVar.push_back(tmp);

	tmp = CVar("emissg","E*_{miss} / GeV",&gmiss.es, -2, 4.5);
	vVar.push_back(tmp);

	tmp = CVar("pmissg","p*_{miss} / GeV",&gmiss.ps, 0, 4.5);
	vVar.push_back(tmp);
	igPmiss = vVar.size();


	tmp = CVar("cosgmiss.ths","cos(#theta*_{miss})",&gmiss.ths, -1, 1);
	vVar.push_back(tmp);

	tmp = CVar("cosgmiss.th","cos(#theta_{miss})",&gmiss.th, -1, 1);
	vVar.push_back(tmp);

	tmp = CVar("mg","M / GeV",&fmm, -5, 5);
	vVar.push_back(tmp);


	tmp = CVar("eecl","Eecl / GeV",&event.Eecl, 0,2.5);
	//vVar.push_back(tmp);

	tmp = CVar("nkss","N_{kaon, ss}",&fnkss,5,-0.5,4.5);
	vVar.push_back(tmp);

	tmp = CVar("nkos","N_{kaon, os}",&fnkos,5,-0.5,4.5);
	vVar.push_back(tmp);

	tmp = CVar("nktot","N_{kaon}",&fnk,5,-0.5,4.5);
	vVar.push_back(tmp);
	tmp = CVar("sumP","#SIGMA |p|",&SumPsig,0.,5);
	vVar.push_back(tmp);

	tmp = CVar("M2gamma","M2_{gamma} / GeV^{2}",&t_gamma.m2, 0,2.);
	vVar.push_back(tmp);
	tmp = CVar("Egamma_","E_{gamma} / GeV",&t_gamma.e, 0,2.5);
	vVar.push_back(tmp);
	tmp = CVar("Ptgamma","Pt_{gamma} / GeV",&t_gamma.pt, 0,1.5);
	vVar.push_back(tmp);
	tmp = CVar("Pzgamma","Pz_{gamma} / GeV",&t_gamma.p3, -1.,1.);
	vVar.push_back(tmp);
	tmp = CVar("Pgamma","P_{gamma} / GeV",&t_gamma.p, 0,2.);
	vVar.push_back(tmp);
	tmp = CVar("thgamma","th_{gamma} ",&t_gamma.costh, -1,1.);
	vVar.push_back(tmp);

	tmp = CVar("M2track","M2_{track} / GeV",&t_track.m2, 0.,5);
	vVar.push_back(tmp);
	tmp = CVar(" Etrack","E_{track} / GeV",&t_track.e, 0,5);
	vVar.push_back(tmp);
	tmp = CVar("Pttrack","Pt_{track} / GeV",&t_track.pt, 0,3.5);
	vVar.push_back(tmp);
	tmp = CVar("Pztrack","Pz_{track} / GeV",&t_track.p3, -1.5,2.5);
	vVar.push_back(tmp);
	tmp = CVar(" Ptrack","P_{track} / GeV",&t_track.p, 0,3.5);
	vVar.push_back(tmp);
	tmp = CVar("thtrack","th_{track} ",&t_track.costh, -1,1.);
	vVar.push_back(tmp);
	tmp = CVar("px1","P_{X}^{1} ",&px1, 0,10);
	vVar.push_back(tmp);
	tmp = CVar("px2","P_{X}^{2} ",&px2, 0,10);
	vVar.push_back(tmp);

	tmp = CVar("missM","M2_{miss} / GeV",&t_PmisswTag.m2, 80,-5,15);
	vVar.push_back(tmp);
	tmp = CVar("missE","E_{miss} / GeV",&t_PmisswTag.e, -1,5);
	vVar.push_back(tmp);
	tmp = CVar("misspt","Pt_{miss} / GeV",&t_PmisswTag.pt, 0,1.5);
	vVar.push_back(tmp);
	tmp = CVar("misspz","Pz_{miss} / GeV",&t_PmisswTag.p3, -2,2);
	vVar.push_back(tmp);
	tmp = CVar("missp","P_{miss} / GeV",&t_PmisswTag.p, 0,3);
	vVar.push_back(tmp);
	tmp = CVar("thmiss","th_{miss} ",&t_PmisswTag.costh, -1,1.);
	vVar.push_back(tmp);
	
	tmp = CVar("missM","M2_{miss} / GeV",&t_PmissWoGamma.m2, 80,-5,15);
	vVar.push_back(tmp);
	tmp = CVar("missE","E_{miss} / GeV",&t_PmissWoGamma.e, -1,5);
	vVar.push_back(tmp);
	tmp = CVar("misspt","Pt_{miss} / GeV",&t_PmissWoGamma.pt, 0,1.5);
	vVar.push_back(tmp);
	tmp = CVar("misspz","Pz_{miss} / GeV",&t_PmissWoGamma.p3, -2,2);
	vVar.push_back(tmp);
	tmp = CVar("missp","P_{miss} / GeV",&t_PmissWoGamma.p, 0,3);
	vVar.push_back(tmp);
	tmp = CVar("thmiss","th_{miss} ",&t_PmissWoGamma.costh, -1,1.);
	vVar.push_back(tmp);

	tmp = CVar("visM","M2_{vis} / GeV",&t_Vis.m2, 60,0,4);
	vVar.push_back(tmp);
	tmp = CVar("visE","E_{vis} / GeV",&t_Vis.e, 0,7);
	vVar.push_back(tmp);
	tmp = CVar("vispt","Pt_{vis} / GeV",&t_Vis.pt, 0,5);
	vVar.push_back(tmp);
	tmp = CVar("vispz","Pz_{vis} / GeV",&t_Vis.p3, -1,4);
	vVar.push_back(tmp);
	tmp = CVar("visp","P_{vis} / GeV",&t_Vis.p, 0,4);
	vVar.push_back(tmp);
	tmp = CVar("thvis","th_{vis} ",&t_Vis.costh, -2,2.);
	vVar.push_back(tmp);

	tmp = CVar("thLepMiss","#Theta_{l,Miss} ",&fThLepMiss, -1,1.);
	vVar.push_back(tmp);
	tmp = CVar("thLepX","#Theta_{l,X} ",&fThLepX, -1,1.);
	vVar.push_back(tmp);
	tmp = CVar("thXMiss","#Theta_{X,Miss} ",&fThXMiss, -1,1.);
	vVar.push_back(tmp);

	tmp = CVar("gx.ps","p*_{X} / GeV",&gx.ps,0,3);
	vVar.push_back(tmp);

	tmp = CVar("gx.es","E*_{X} / GeV",&gx.es,0,6);
	vVar.push_back(tmp);

	tmp = CVar("gx.m2","M_{X} / GeV",&gx.m2,0,6);
	vVar.push_back(tmp);

	tmp = CVar("gx.m2.zoom","M_{X} / GeV",&gx.m2,40,1.8,2.1);
	vVar.push_back(tmp);

	tmp = CVar("gx.ths","cos(#theta_{X})",&gx.ths,-1,1);
	vVar.push_back(tmp);

	tmp = CVar("x.ps","p*_{X tracks} / GeV",&x.ps,0,3);
	vVar.push_back(tmp);

	tmp = CVar("x.es","E*_{X tracks} / GeV",&x.es,0,6);
	vVar.push_back(tmp);

	tmp = CVar("x.m2","M_{X tracks} / GeV",&x.m2,80,0,4);
	vVar.push_back(tmp);
	iVarMx = vVar.size();

	tmp = CVar("x.ths","cos(#theta_{X tracks})",&x.ths,-1,1);
	vVar.push_back(tmp);

	tmp = CVar("Egamma","E*_{$gamma} / GeV",&fEgamma,-1,6);
	vVar.push_back(tmp);

	tmp = CVar("Pgamma","p*_{$gamma} / GeV",&fPgamma,0,6);
	vVar.push_back(tmp);

	tmp = CVar("MD0","M_{D0} / GeV",&fMD0Cand,1.83,1.9);
	vVar.push_back(tmp);

	tmp = CVar("MPhi","M_{#Phi} / GeV",&fMDpCand,0.99,1.05);
//	tmp = CVar("MDp","M_{D+} / GeV",&fMDpCand,1.83,1.9);
	vVar.push_back(tmp);
	ihistKK = vVar.size();
	tmp = CVar("MKK","M_{KK} / GeV",&fMKK,100,0.95,1.9);
	vVar.push_back(tmp);
	tmp = CVar("MKPi","M_{K#pi} / GeV",&fMKPi,200,0.6,2.);
	vVar.push_back(tmp);
	tmp = CVar("MPiPi","M_{#pi#pi} / GeV",&fMPiPi,50,0.5,1.);
	vVar.push_back(tmp);
    tmp = CVar("MOmega","M_{#Omega} / GeV",&fMOmegaCand,0.49,0.51);
	vVar.push_back(tmp);

	tmp = CVar("xlep.th_xl",&xlep.th_xl,0.,180);
	//vVar.push_back(tmp);

	tmp = CVar("xlep.ths_xl",&xlep.ths_xl,-1.,1);
//	vVar.push_back(tmp);

	tmp = CVar("xlep.th_xl",&xlep.th_xl,-1.,1);
//	vVar.push_back(tmp);



	tmp = CVar("lclass",&flclass,7,8,15);
	//vVar.push_back(tmp);

	tmp = CVar("fl_lep1",&fl_lep1,35,9.5,44.5);
	vVar.push_back(tmp);

	tmp = CVar("fl_mother1",&fl_mother1,8,0.5,8.5);
//	tmp.fl_logy = true;
	vVar.push_back(tmp);

	tmp = CVar("fl_lep2",&fl_lep2,35,9.5,44.5);
	//vVar.push_back(tmp);

	tmp = CVar("fl_mother2",&fl_mother2,8,0.5,8.5);
	//vVar.push_back(tmp);

	tmp = CVar("fl_mother12",&fl_mother12,7,0,7);
	//vVar.push_back(tmp);

	tmp = CVar("truth q2",&truth.q2,0,15);
	//vVar.push_back(tmp);

	tmp = CVar("BDT1",&BDT1,-0.2,0.6);
	//vVar.push_back(tmp);

	tmp = CVar("BDT2",&BDT2,-0.7,0.);
	//vVar.push_back(tmp);

	tmp = CVar("q2",&event.q2,-5,25);
	vVar.push_back(tmp);
	tmp = CVar("ntrack","N_{track}",&fnch,10,0.5,10.5);
	vVar.push_back(tmp);
	tmp = CVar("ngamma","N_{gamma}",&fnn,16,-0.5,15.5);
	vVar.push_back(tmp);

	tmp = CVar("tagcorr","Tag Correction",&Corrections[0],60,0.,2);
	vVar.push_back(tmp);

	tmp = CVar("BrCorr","BR Correction",&Corrections[1],60,0.,2);
	vVar.push_back(tmp);

	tmp = CVar("sl","SL Correction",&Corrections[2],60,0.,2);
	vVar.push_back(tmp);

	tmp = CVar("Dss","Dss Correction",&Corrections[3],60,0.,2);
	vVar.push_back(tmp);

	tmp = CVar("eEff","e Eff",&Corrections[4],60,0.,2);
	vVar.push_back(tmp);

	tmp = CVar("eFake","e Fake",&Corrections[5],60,0.,2);
	vVar.push_back(tmp);


	tmp = CVar("muEff","mu Eff",&Corrections[6],60,0.,2);
	vVar.push_back(tmp);

	tmp = CVar("muFake","mu Fake",&Corrections[7],60,0.,2);
	vVar.push_back(tmp);

	nVariables = vVar.size();
	tmp = CVar("trk_sep",&event.q2,7,-0.5,6.5);
	vVar.push_back(tmp);
	tmp = CVar("trk_p_e",&event.q2,0,3);
	vVar.push_back(tmp);
	tmp = CVar("trk_th_e",&event.q2,-1,1);
	vVar.push_back(tmp);
	tmp = CVar("trk_dr_e",&event.q2,-0.2,0.2);
	vVar.push_back(tmp);
	tmp = CVar("trk_dz_e",&event.q2,-2,2);
	vVar.push_back(tmp);
	tmp = CVar("trk_p_mu",&event.q2,0,3);
	vVar.push_back(tmp);
	tmp = CVar("trk_th_mu",&event.q2,-1,1);
	vVar.push_back(tmp);
	tmp = CVar("trk_dr_mu",&event.q2,-0.2,0.2);
	vVar.push_back(tmp);
	tmp = CVar("trk_dz_mu",&event.q2,-2,2);
	vVar.push_back(tmp);
	tmp = CVar("trk_p_pi",&event.q2,0,3);
	vVar.push_back(tmp);
	tmp = CVar("trk_th_pi",&event.q2,-1,1);
	vVar.push_back(tmp);
	tmp = CVar("trk_dr_pi",&event.q2,-0.2,0.2);
	vVar.push_back(tmp);
	tmp = CVar("trk_dz_pi",&event.q2,-2,2);
	vVar.push_back(tmp);
	tmp = CVar("trk_N_pi",&event.q2,10,-0.5,9.5);
	vVar.push_back(tmp);
	tmp = CVar("trk_p_k",&event.q2,0,3);
	vVar.push_back(tmp);
	tmp = CVar("trk_th_k",&event.q2,-1,1);
	vVar.push_back(tmp);
	tmp = CVar("trk E -q ",&event.q2,-3,3);
	vVar.push_back(tmp);
	tmp = CVar("trk E q",&event.q2,-3,3);
	vVar.push_back(tmp);
	tmp = CVar("trk p -q",&event.q2,-3,3);
	vVar.push_back(tmp);
	tmp = CVar("trk p q",&event.q2,-3,3);
	vVar.push_back(tmp);
	tmp = CVar("trk th",&event.q2,-1,1);
	vVar.push_back(tmp);
	tmp = CVar("trk dr",&event.q2,-2,2);
	vVar.push_back(tmp);
	tmp = CVar("trk dz",&event.q2,-3,3);
	vVar.push_back(tmp);
	tmp = CVar("trk DeltaTh ss",&event.q2,-1,1);
	vVar.push_back(tmp);
	tmp = CVar("trk DeltaTh os",&event.q2,-2,2);
	vVar.push_back(tmp);
	tmp = CVar("trk DeltaTh ss",&event.q2,-1,1);
	vVar.push_back(tmp);
	tmp = CVar("trk DeltaTh os",&event.q2,0,2);
	vVar.push_back(tmp);
	tmp = CVar("trk DeltaTh ss",&event.q2,0,2);
	vVar.push_back(tmp);
	tmp = CVar("trk DeltaTh os",&event.q2,-1,1);
	vVar.push_back(tmp);
	tmp = CVar("trk DeltaTh ss",&event.q2,-1,1);
	vVar.push_back(tmp);
	tmp = CVar("trk DeltaTh os",&event.q2,-1,1);
	vVar.push_back(tmp);
	tmp = CVar("dz Kaon from D0",&event.q2,-0.1,0.1);
	vVar.push_back(tmp);
	tmp = CVar("dz Kaon from Phi",&event.q2,-0.1,0.1);
	vVar.push_back(tmp);
    tmp = CVar("dr Kaon from D0",&event.q2,-0.1,0.1);
	vVar.push_back(tmp);
	tmp = CVar("dr Kaon from Phi",&event.q2,-0.1,0.1);
	vVar.push_back(tmp);
    tmp = CVar("th Kaon from D0",&event.q2,-1,1);
	vVar.push_back(tmp);
	tmp = CVar("th Kaon from Phi",&event.q2,-1,1);
	vVar.push_back(tmp);
	tmp = CVar("P Kaon from D0",&event.q2,0.05,2);
	vVar.push_back(tmp);
	tmp = CVar("P Kaon from Phi",&event.q2,0.05,2);
	vVar.push_back(tmp);

	nExtraVariables = vVar.size()-nVariables;
	nVariables = vVar.size();

	//contributions
	CCont tmp_cont;

	tmp_cont = CCont("signal", kRed);
	vCont.push_back(tmp_cont);

	tmp_cont = CCont("Dlnu", kGreen );
	vCont.push_back(tmp_cont);

	tmp_cont = CCont("Dslnu", kRed-2 );
	vCont.push_back(tmp_cont);

	tmp_cont = CCont("Dpilnu", kBlue );
	vCont.push_back(tmp_cont);


	tmp_cont = CCont("Dsslnu", kGreen +2);
	//tmp_cont.Weight = 0.7;
	vCont.push_back(tmp_cont);

	tmp_cont = CCont("XcXc", kBlue+1);
	vCont.push_back(tmp_cont);

	tmp_cont = CCont("XcXu", kCyan+1);
	vCont.push_back(tmp_cont);

	tmp_cont = CCont("FakeLeptons", kBlue -1);
	vCont.push_back(tmp_cont);

	tmp_cont = CCont("ulnu", kOrange +2);
	vCont.push_back(tmp_cont);

	tmp_cont = CCont("rare", kGray+2);
	//tmp_cont.Weight = 0.8;
	vCont.push_back(tmp_cont);

	tmp_cont = CCont("otherB", kOrange -2);
	vCont.push_back(tmp_cont);

	tmp_cont = CCont("badB", kOrange-6);
	vCont.push_back(tmp_cont);

	tmp_cont = CCont("Continuum", kMagenta -3);
	vCont.push_back(tmp_cont);

	tmp_cont = CCont("Data", kBlack);
	vCont.push_back(tmp_cont);

	CFile tmp_file;
/*	if(_pilep)
		tmp_file = CFile("/scratch/bctaunu/onres.root","/scratch/bctaunu/onres_pilep_bdt2.root",'o');
	else
		tmp_file = CFile("/scratch/bctaunu/onres_s0.root",'o');
	vFile.push_back(tmp_file );


	if(_pilep)
	tmp_file = CFile("/scratch/bctaunu/continuum.root","/scratch/bctaunu/continuum_pilep_bdt2.root",'c');
	else
	tmp_file = CFile("/scratch/bctaunu/continuum_s0.root",'c');
	vFile.push_back(tmp_file );

	if(_pilep)
	tmp_file = CFile("/scratch/bctaunu/data.root","/scratch/bctaunu/data_pilep_bdt2.root",'d');
	else
	tmp_file = CFile("/scratch/bctaunu/data.root",'d');
	vFile.push_back(tmp_file );
*/
	float nStreams = 2. ;
	float nContiStreams = 2. ;
	//string path = string("/scratch/bctaunu/kakuno/robin/new");
//	string path = string("/scratch/bctaunu/kakuno/bugfree2/");
	string path = string("/scratch/bctaunu/kakuno/new_0815/");


	tmp_file = CFile(path,"DssMC.root","lep/tree",'s');
	tmp_file.weight = 1./2.8;
	//vFile.push_back(tmp_file );

	tmp_file = CFile(path,"ulnu.root","lep/tree",'u');
	tmp_file.weight = 1./20;
	//vFile.push_back(tmp_file );
	
	tmp_file = CFile(path,"rare.root","lep/tree",'r');
	tmp_file.weight = 1./50;
	//vFile.push_back(tmp_file );
	tmp_file =CFile(path, "charged_s0.root", "lep/tree",'o');
	tmp_file.weight = 1./nStreams;
	vFile.push_back(tmp_file );
	tmp_file =CFile(path, "mixed_s0.root", "lep/tree",'o');
	tmp_file.weight = 1./nStreams;
	vFile.push_back(tmp_file );
	tmp_file =CFile(path, "continuum_s0.root", "lep/tree",'c');
	tmp_file.weight = 1./nContiStreams;
	//vFile.push_back(tmp_file );
/*
	   tmp_file =CFile(path, "charged_s1.root", "lep/tree",'o');
	tmp_file.weight = 1./nStreams;
	vFile.push_back(tmp_file );
	tmp_file =CFile(path, "mixed_s1.root", "lep/tree",'o');
	tmp_file.weight = 1./nStreams;
	vFile.push_back(tmp_file );
	tmp_file =CFile(path, "continuum_s1.root", "lep/tree",'c');
	tmp_file.weight = 1./nContiStreams;
	vFile.push_back(tmp_file );

	tmp_file =CFile(path, "data.root", "lep/tree",'d');
	vFile.push_back(tmp_file );
*/
/*	tmp_file =CFile(path, "charged_s0.root", "dilep/tree",'o');
	tmp_file.weight = 1./nStreams;
	vFile.push_back(tmp_file );
	tmp_file =CFile(path, "mixed_s0.root", "dilep/tree",'o');
	tmp_file.weight = 1./nStreams;
	vFile.push_back(tmp_file );
	tmp_file =CFile(path, "continuum_s0.root", "dilep/tree",'c');
	tmp_file.weight = 1./nContiStreams;
	vFile.push_back(tmp_file );
	tmp_file =CFile(path, "data.root", "dilep/tree",'d');
	vFile.push_back(tmp_file );
*/

}

TH1F ***hists;
THStack **StackMC;
TLegend *lgd;
TH1F ***hists_uw;//unweighted

/*
void CatEvents(vector<MyParticle> &v)
{
	static TAxis** axis = null;
	int nMaxTracks = 6;
	if(!TAxis)
	{
		axis = new TAxis[nMaxTracks];
		for(int i = 0; i<nMaxTracks; i++)
			axis[i] = new TAxis("",3,0,2);
	}
	for(auto it = v.begin(); it!=v.end(); ++it)
	{
		
	}
	
}*/

void DefHist()
{
	char sHistName[200];
	char sStackName[200];
    hists 	= new TH1F**[nContributions];
	hists_uw= new TH1F**[nContributions];
	lgd	= new TLegend(0,0.5,.5,1);
	for(int i = 0; i<nContributions; i++)
	{
		hists[i] = new TH1F*[nVariables];
		hists_uw[i] = new TH1F*[nVariables];
		for( int j =0; j<nVariables;j++)
		{
			sprintf(sHistName,"hist_%s_%s",vVar[j].Name.c_str(),vCont[i].Name.c_str());
			hists[i][j] = newHist(sHistName,vVar[j]);
			hists[i][j] -> SetLineColor(vCont[i].Color);
			hists[i][j] -> SetFillColor(vCont[i].Color);
			sprintf(sHistName,"hist_uw_%s_%s",vVar[j].Name.c_str(),vCont[i].Name.c_str());
			hists_uw[i][j] = newHist(sHistName,vVar[j]);

			//hists[i][j] -> SetFillStyle(3354);
			hists[i][j] -> SetLineWidth( 2 );

		}
		lgd->AddEntry(hists[i][0],vCont[i].Title.c_str(),"lf");
		lgd->SetFillColor( kWhite );


	}


	StackMC = new THStack*[nVariables];
	for(int j =0; j<nVariables; ++j)
	{
		sprintf(sStackName, "StackMC_%s",vVar[j].Name.c_str());
		StackMC[j] = new THStack(sStackName,sStackName);
		sprintf(sHistName,"hist_div_%s",vVar[j].Name.c_str());
		//hists_div[j] = new TH1F(sHistName,sHistName,nBins[j],fXmin[j],fXmax[j]);

		StackMC[j] -> Add(hists[continuum][j]);
		StackMC[j] -> Add(hists[badB][j]);
		StackMC[j] -> Add(hists[otherB][j]);
		StackMC[j] -> Add(hists[rare][j]);
		StackMC[j] -> Add(hists[ulnu][j]);
		StackMC[j] -> Add(hists[FakeLep][j]);
		StackMC[j] -> Add(hists[charm_other][j]);
		StackMC[j] -> Add(hists[double_charm][j]);
		StackMC[j] -> Add(hists[charm_ss_lep][j]);
		StackMC[j] -> Add(hists[charm_nonres_lep][j]);
		StackMC[j] -> Add(hists[charm_s_lep][j]);
		StackMC[j] -> Add(hists[charm_lep][j]);
		//StackMC[j] -> Add(hists[signal_cf][j]);
		StackMC[j] -> Add(hists[signal_decay][j]);

		hists[data][j] -> SetMarkerStyle(20);
	}


}

void Divide(TH1F* h1, TH1F* h2)
{
	for(int i= 1; i<=h1->GetNbinsX();i++)
	{
		float a = h1->GetBinContent(i);
		float b = h2->GetBinContent(i);
		if(b>0 && a>0)
			a = a/b;
		else a=1;
		h1->SetBinContent(i,a);
	}

}

void SmearHist(TH1F* h, float width)
{

	int xbins = h->GetNbinsX();
		//Make Matrix
	float** m = new float*[xbins];
	for(int i =0; i<xbins; i++)
	{
		m[i] = new float[xbins];
		float Norm = 0;
		for(int j = 0; j<xbins;j++)
		{
			m[i][j] = TMath::Gaus(j,i,width);
			Norm+=m[i][j];
		}
		for(int j = 0; j<xbins;j++)
		{
			m[i][j]/=Norm;
		}
	}

	//apply matrix
	float *val = new float[xbins];
	for(int i =0; i<xbins; i++)
	{
		val[i] = 0;
		for(int j = 0; j<xbins;j++)
		{
			val[i]+=m[i][j]*h->GetBinContent(j+1);
		}
	}
	for(int i =0; i<xbins; i++)
		h->SetBinContent(i+1,val[i]);

	for(int i =0; i<xbins; i++)
		delete[] m[i];
	delete[] m;
	delete[] val;

}
void SimpleFit(TH1F* hdata, TH1F* hc1, TH1F* hc2, double &x, double &dx, double &yield, TH1F *hc1_uw = NULL, TH1F *hc2_uw = NULL)
{

	TH1F *w1 = NULL, *w2 = NULL;
	TObjArray *mc = new TObjArray(2);        // MC histograms are put in this array
	if(hc1_uw!= NULL && hc2_uw !=NULL)
	{
		w1 = (TH1F*)hc1->Clone();
		Divide(w1,hc1_uw);
		w2 = (TH1F*)hc2->Clone();
		Divide(w2,hc2_uw);
		mc->Add(hc1_uw);
		mc->Add(hc2_uw);

	}
	else
	{
		mc->Add(hc1);
		mc->Add(hc2);
	}
	TFractionFitter* fit = new TFractionFitter(hdata, mc); // initialise
	if(hc1_uw!= NULL && hc2_uw !=NULL)
	{
		fit->SetWeight(0,w1);
		fit->SetWeight(1,w2);

	}

	//fit->Constrain(1,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
	Int_t status = fit->Fit();               // perform the fit
	//fit->ErrorAnalysis(0);
	fit->GetResult(0,x,dx);
	yield = fit->GetPlot()->Integral();
	delete fit;
	delete mc;
	delete w1;
	delete w2;

}
void PerfFit(TH1F* hdata, TH1F** hmc, int nHist, double *x, double *dx, double &yield, TH1F **hmc_weight = NULL)
{
	TObjArray *mc = new TObjArray(nHist);        // MC histograms are put in this array
	for(int i =0; i<nHist; ++i)
		mc->Add(hmc[i]);
	TFractionFitter* fit = new TFractionFitter(hdata, mc); // initialise
	if(hmc_weight != NULL )
	{
		for(int i =0; i<nHist; ++i)
		{
			fit->SetWeight(i,hmc_weight[i]);
			//fit->Constrain(i,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
		}

	}

	Int_t status = fit->Fit();               // perform the fit
	//fit->ErrorAnalysis(0);
	for(int i =0; i<nHist; ++i)
		fit->GetResult(i,x[i],dx[i]);
	yield = fit->GetPlot()->Integral();
	delete fit;
	delete mc;
}

void SimpleFit(TH1F* hdata, TH1F** hmc, int nHist, double *x, double *dx, double &yield, TH1F **hmc_uw = NULL)
{

	TH1F **w = new TH1F*[nHist] ;
	if(hmc_uw!= NULL)
	{
		for(int i =0; i<nHist; ++i)
		{
			w[i] = (TH1F*)hmc[i]->Clone();
			Divide(w[i],hmc_uw[i]);

		}
		PerfFit(hdata,hmc_uw,nHist,x,dx,yield,w);

	}
	else
	{
		PerfFit(hdata,hmc,nHist,x,dx,yield);
	}
	for(int i =0; i<nHist; ++i)
		delete w[i];
	delete[] w;


}
void PrintMaxima(TH1F* hist, int nMax)
{
	cout<<"Print Max: "<<endl;
	int maxBin;
	for(int i = 0; i<nMax; i++)
	{
		if(i>hist->GetNbinsX())continue;
		maxBin = hist->GetMaximumBin();
		if(hist->GetBinContent(maxBin)<=0) break;
		cout.width(50);
		cout
			<<left
			<<hist->GetXaxis()->GetBinLabel(maxBin)
			<<right
			<<hist->GetBinContent(maxBin)
			<<endl;
		hist->SetBinContent(maxBin,0);
	}
	cout<<endl<<endl;
}
char sCut[500];
char sCut_tmp[500];

int dclass2lclass[12] = {0,2,0,1,0,3,4,5,6,0,0,0};



float mvis;
int GetCont(std::vector<CFile>::iterator it_f)
{
			int cont = -1;
			if(it_f -> Type == 'c')
				cont = continuum;
			//else if( it_f -> Type == 'u' || it_f->Type == 'r' )
			else if( it_f -> Type == 'u' )
				cont = signal_decay;
			else if(event.lclass==0 || btag.fl_tag<0)
				cont = badB;
            else cont = charm_lep;
            /*
			else if(abs(event.lclass)==1 )
				cont = signal_decay;
			else if((lep1.fl_lep != 10 && lep1.fl_lep!=21)  )
				cont = FakeLep;
			else if( lep1.fl_mother != 3  )
				cont = rare;
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
                */
		//	else cont = charm_lep;
			if(cont == -1 ) cout<<"err"<<endl;
		/*
			if(abs(event.lclass) ==1 || abs(event.lclass) == 2)
			{
				cont = (truth.w -1)*10;
			}
			else cont = continuum;
*/

//			if(btag.pcode_b<0)
//				cont = event.nk%10;
//			else
//				cont = event.nk/10;
		/*cont = fnkss + 3*fnkos;
			if(fnk >=5 ||fnkos>1 || fnkss >2) cont = 5;
			cont = 5-cont;*/
			/*if(fnkss == 1) cont = 0;
			else cont = 1;
			if(fnk == 0) cont = 2;
		*///	cont = 4-btag.fl_tag;
//			if(btag.fl_tag<2) cont =1;
//			else cont = 0;
			/*if(abs(event.lclass) == 2)
			{
				if(event.dclass ==1 || event.dclass == 3) cont = 0;
			//	else if(event.dclass ==3) cont = 1;
				else cont = 1;
			}
			else cont = 2;
		*/
			/*if(cont<double_charm || cont == ulnu)
				cont = 0;
			else cont = 1;
			*/

			/*
			if(lep1.fl_mother == 3) cont = 0;
			else cont = 1;
			if(lep1.fl_lep != 10 && lep1.fl_lep !=21) cont = 1;

			*/

//            if(gmiss.e<0) cont = 1;
//            else cont = 2;

			if(it_f -> Type == 'd')
				cont = data;
 //			else if(gmiss.m2<-0.01) cont = 1;
//			else cont = 3;
			return cont;
}

float RecoDp2KPiPi(std::vector<MyParticle> &v)
{
    float fD = 999.;
    float mD0 = 1.86961;
    for(auto it = v.begin(); it!= v.end(); ++it)
    {
        if(it->pid != 2) continue;
        for(auto kt = v.begin(); kt!= v.end(); ++kt)
    	{
        	if(kt->pid != 3) continue;
        	for(auto jt = v.begin(); jt!= v.end(); ++jt)
        	{
        	    if(jt->pid != 3) continue;
       		    float m = (it->P + jt->P + kt->P).M();
       		    if(fabs(fD-mD0) > fabs(m-mD0))
          	    fD = m;
          	}
        }
    }
    return fD;
}
float RecoDp2PiPiPi(std::vector<MyParticle> &v)
{
    float fD = 999.;
    float mD0 = 1.86961;
    for(auto it = v.begin(); it!= v.end(); ++it)
    {
        if(it->pid != 3) continue;
        for(auto kt = v.begin(); kt!= v.end(); ++kt)
    	{
        	if(kt->pid != 3) continue;
        	for(auto jt = v.begin(); jt!= v.end(); ++jt)
        	{
        	    if(jt->pid != 3) continue;
       		    float m = (it->P + jt->P + kt->P).M();
       		    if(fabs(fD-mD0) > fabs(m-mD0))
          	    fD = m;
          	}
        }
    }
    return fD;
}
float RecoD2KPi(std::vector<MyParticle> &v)
{
    float fD = 999.;
    float mD0 = 1.86486;
    for(auto it = v.begin(); it!= v.end(); ++it)
    {
        if(it->pid != 2) continue;
        for(auto jt = v.begin(); jt!= v.end(); ++jt)
        {
            if(jt->pid != 3) continue;
            float m = (it->P + jt->P).M();
            if(fabs(fD-mD0) > fabs(m-mD0))
                fD = m;
        }
    }
    return fD;
}

float RecoD2KPiPiPi(std::vector<MyParticle> &v)
{
    float fD = 999.;
    float mD0 = 1.86486;
    for(auto it = v.begin(); it!= v.end(); ++it)
    {
        if(it->pid != 2) continue;
        for(auto kt = v.begin(); kt!= v.end(); ++kt)
    	{
        	if(kt->pid != 3) continue;
        	for(auto jt = v.begin(); jt!= v.end(); ++jt)
        	{
			if(jt->pid != 3) continue;
        		for(auto lt = v.begin(); lt!= v.end(); ++lt)
        		{
        	    		if(lt->pid != 3) continue;
       		    		float m = (it->P + jt->P + kt->P + lt->P).M();
       		    		if(fabs(fD-mD0) > fabs(m-mD0))
          	    		fD = m;
			}
          	}
        }
    }
    return fD;
}

float RecoD2KK(std::vector<MyParticle> &v, MyParticle *p)
{
    float fD = 999.;
    float mD0 = 1.86486;
    for(auto it = v.begin(); it!= v.end(); ++it)
    {
        if(it->pid != 2) continue;
        for(auto jt = v.begin(); jt!= v.end(); ++jt)
        {
            if(jt->pid != 2) continue;
            float m = (it->P + jt->P).M();
            if(fabs(fD-mD0) > fabs(m-mD0))
            {
                fD = m;
                p[0] = *it;
                p[1] = *jt;
            }
        }
    }
    return fD;
}
float RecoD24Pi(std::vector<MyParticle> &v)
{
    float fD = 999.;
    float mD0 = 1.86486;
    for(auto it = v.begin(); it!= v.end(); ++it)
    {
        if(it->pid != 3) continue;
        for(auto jt = it+1; jt!= v.end(); ++jt)
        {
            if(jt->pid != 3) continue;
            for(auto kt = jt+1; kt!= v.end(); ++kt)
            {
                if(kt->pid != 3) continue;
                for(auto lt = kt+1; lt!= v.end(); ++lt)
                {
                    if(lt->pid != 3) continue;
                    float m = (it->P + jt->P + kt->P + lt->P).M();
                    if(fabs(fD-mD0) > fabs(m-mD0))
                        fD = m;
                }
            }
        }
    }
    return fD;
}
float RecoPhi2KK(std::vector<MyParticle> &v, MyParticle *p)
{
    float fCand = 999.;
    float mPDG = 1.019455;
    for(auto it = v.begin(); it!= v.end(); ++it)
    {
        if(it->pid != 2) continue;
        for(auto jt = v.begin(); jt!= v.end(); ++jt)
        {
            if(jt->pid != 2) continue;
            float m = (it->P + jt->P).M();
            if(fabs(fCand-mPDG) > fabs(m-mPDG))
            {
                fCand = m;
                p[0] = *it;
                p[1] = *jt;
            }
        }
    }
    return fCand;
}
float RecoOmega2PiPi(std::vector<MyParticle> &v)
{
    float fCand = 999.;
    float mPDG = 0.78265;
    for(auto it = v.begin(); it!= v.end(); ++it)
    {
        if(it->pid != 3) continue;
        for(auto jt = v.begin(); jt!= v.end(); ++jt)
        {
            if(jt->pid != 3) continue;
            float m = (it->P + jt->P).M();
            if(fabs(fCand-mPDG) > fabs(m-mPDG))
                fCand = m;
        }
    }
    return fCand;
}
float RecoK0S2PiPi(std::vector<MyParticle> &v)
{
    float fCand = 999.;
    float mPDG = 0.4976;
    for(auto it = v.begin(); it!= v.end(); ++it)
    {
        if(it->pid != 3) continue;
        for(auto jt = v.begin(); jt!= v.end(); ++jt)
        {
            if(jt->pid != 3) continue;
            float m = (it->P + jt->P).M();
            if(fabs(fCand-mPDG) > fabs(m-mPDG))
                fCand = m;
        }
    }
    return fCand;
}

float Corr2Track[] = {
1	,
1	,
1	,
0.976176	,
1.04735	,
0.976953	,
0.994644	,
0.989452	,
0.973591	,
0.961919	,
1.05482	,
1.00791	,
1.00377	,
1.02367	,
1.00038	,
0.975737	,
1.00572	,
1.02795	,
1.02633	,
1.05808	,
1.01465	,
0.988255	,
0.998573	,
1.01384	,
1.03947	,
0.978806	,
0.971848	,
1.02078	,
0.979847	,
0.998557	,
0.994802	,
0.952634	,
0.967487	,
0.992261	,
1.03152	,
0.926549	,
0.944721	,
0.952123	,
1.00385	,
0.932115	,
1.01823	,
0.965955	,
0.875781	,
1.0698	,
0.935974	,
0.927772	,
0.945634	,
0.954734	,
0.897369	,
0.991722	,
1.01488	,
1.09588	,
1.10338	,
1.05723	,
0.98945	,
1.01682	,
0.882137	,
1.16991	,
0.962177	,
1.04746	};

float Corr3Track[] = {
1	,
1	,
1	,
0.932988	,
0.971332	,
0.967803	,
0.998411	,
1.00286	,
0.985913	,
0.998678	,
1.01715	,
0.999114	,
1.03027	,
1.0385	,
1.03235	,
1.04138	,
0.985212	,
1.01646	,
1.03468	,
1.03675	,
0.981736	,
1.01894	,
1.02457	,
0.94936	,
0.972879	,
1.03418	,
0.95896	,
0.97411	,
1.01172	,
0.918873	,
0.853356	,
0.998491	,
0.953021	,
0.958134	,
1.11864	,
1.05095	,
0.925757	,
0.902148	,
0.921314	,
1.04296	,
0.805661	,
0.952761	,
1.2673	,
0.722042	,
0.915587	,
0.82451	,
0.776314	,
1.77708	,
1.32065	,
2.14385	,
0.601328	,
1.05825	,
1.3322	,
1.16132	,
0.764365	,
2.42273	,
1.03079	,
0.434396	,
1	,
2.12052	

};

float Corr4Track[] = {
1	,
1	,
1	,
0.944771	,
0.987531	,
0.983396	,
0.976368	,
0.999179	,
0.995388	,
1.03351	,
1.05507	,
1.00514	,
1.02704	,
1.03461	,
1.06843	,
1.04964	,
1.06351	,
1.00861	,
0.976955	,
1.07457	,
1.06586	,
0.899726	,
0.88821	,
0.970023	,
0.755301	,
0.803871	,
0.742683	,
1.17025	,
0.530409	,
0.861632	,
0.446962	,
0.43016	,
0.687488	,
0.627267	,
0.0859432	,
1.18629	,
1.13134	,
0.321104	,
1.02146	,
1.93415	,
1	,
1	,
1	,
2.32344	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1
};

float Corr5Track[] = {
1	,
1	,
1	,
0.911828	,
0.941147	,
0.946753	,
0.998633	,
1.03932	,
1.04848	,
1.03774	,
1.10476	,
1.05934	,
1.07284	,
1.1272	,
1.03002	,
1.06315	,
0.956717	,
0.870378	,
1.06345	,
1.00548	,
0.844438	,
0.532839	,
1.32511	,
0.928205	,
0.712776	,
1.19205	,
0.751014	,
1	,
1.69325	,
1	,
1.31352	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	,
1	
};

/*
float nTrackCorr[40][10] = {
{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
{ 1.72477, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
{ 1.29083, 9.79779, 1, 1, 1, 1, 1, 1, 1, 1 },
{ 1.10466, 0.903693, 1.08578, 1, 1, 1, 1, 1, 1, 1 },
{ 1.26463, 1.01405, 0.771969, 1, 1, 1, 1, 1, 1, 1 },
{ 1.25032, 0.935346, 0.560523, 2.01477, 1, 1, 1, 1, 1, 1 },
{ 1.49034, 1.13473, 1.04862, 0.670222, 1, 1, 1, 1, 1, 1 },
{ 1.7236, 1.1053, 0.811924, 0.987567, 1.07225, 1, 1, 1, 1, 1 },
{ 1.74061, 1.12476, 0.969379, 1.06006, 0.520922, 1, 1, 1, 1, 1 },
{ 1.64563, 1.11501, 1.01731, 0.949831, 0.862003, 0.354656, 1, 1, 1, 1 },
{ 2.01821, 1.10407, 1.05088, 0.996699, 0.826702, 0.169928, 1, 1, 1, 1 },
{ 2.10913, 1.0669, 1.02359, 1.04122, 0.978318, 0.475114, 1, 1, 1, 1 },
{ 2.14428, 1.12397, 1.0072, 1.06948, 0.985333, 0.64721, 1.03535, 1, 1, 1 },
{ 2.42264, 1.07792, 1.07104, 1.11663, 1.03412, 0.86356, 1.34582, 1, 1, 1 },
{ 2.7831, 1.0975, 1.0158, 1.14179, 1.06617, 1.24752, 0.964107, 4.35723, 1, 1 },
{ 4.2222, 1.03903, 1.0298, 1.13311, 1.13906, 0.991812, 0.705759, 1, 1, 1 },
{ 7.39916, 1.0853, 1.00883, 1.11337, 1.14024, 0.998995, 0.987065, 0.243378, 1, 1 },
{ 8.04158, 0.957391, 1.10669, 1.11715, 1.09331, 0.998498, 0.924893, 0.863001, 1, 1 },
{ 6.75711, 1.11049, 0.997528, 1.09511, 1.11129, 1.10966, 0.953847, 1.17968, 0.512371, 1 },
{ 14.9811, 0.85001, 0.991722, 1.10455, 1.18285, 1.04591, 1.09581, 0.870413, 0.57122, 1 },
{ 53.3768, 1.04568, 0.945121, 1.13359, 1.13309, 1.12285, 1.13199, 0.881789, 0.75468, 1 },
{ 1, 1.10303, 0.901466, 1.10877, 1.15027, 1.08535, 1.08986, 0.884481, 0.580455, 0.624027 },
{ 1, 1.10525, 0.898987, 1.06594, 1.20718, 1.08729, 1.06416, 0.91033, 0.824801, 0.936538 },
{ 1, 1.28548, 0.667661, 0.982544, 1.11084, 1.15372, 1.00884, 1.03453, 1.13968, 1.31645 },
{ 1, 1.33259, 0.758789, 0.990572, 1.03775, 1.12121, 1.06173, 1.48547, 1.18727, 0.873535 },
{ 1, 1.20582, 0.712522, 1.01351, 1.01762, 1.08228, 1.2149, 1.07187, 1.495, 1.62321 },
{ 1, 2.01959, 0.648214, 1.00365, 1.06421, 1.08246, 1.19983, 0.999012, 1.48173, 1.10496 },
{ 1, 1.37273, 0.348468, 0.830307, 0.967738, 1.11327, 1.12795, 1.37337, 1.59982, 1.41933 },
{ 1, 0.769993, 0.431925, 0.961029, 1.08311, 1.08084, 1.04536, 1.19801, 1.03083, 1.3397 },
{ 1, -0.331335, 0.361081, 1.02223, 0.960148, 1.08532, 1.14796, 1.16087, 1.4358, 0.476358 },
{ 1, 1.46088, -0.341141, 0.971384, 0.960308, 1.02626, 1.15422, 1.12056, 1.18145, 2.16642 },
{ 1, 0.330474, -0.219279, 0.72436, 0.893102, 1.03265, 1.1989, 1.36756, 1.25148, 1.64353 },
{ 1, 8.41737, 0.763261, 0.862649, 0.763999, 0.940366, 1.15551, 1.26082, 1.04662, 1.43863 },
{ 1, 3.34894, -0.71023, 0.781938, 0.845961, 0.86681, 0.893497, 1.16374, 1.15513, 1.82119 },
{ 1, 10.374, 1.1431, 0.947535, 0.96169, 0.727219, 0.891126, 1.10614, 2.11751, 1.25404 },
{ 1, 1, 0.707446, 0.844532, 0.75976, 0.908288, 0.925949, 1.15277, 1.47738, 1.24462 },
{ 1, -0.829344, 0.722382, 1.21832, 0.870082, 0.908464, 0.944844, 1.27002, 1.50302, 0.727807 },
{ 1, 3.72042, 2.47058, 0.991239, 0.905043, 0.888918, 1.09385, 1.4183, 2.02362, 2.08913 },
{ 1, 2.25982, 2.63686, 0.953181, 1.16829, 0.887548, 0.853564, 1.31672, 2.25636, 3.20746 },
{ 1, 1, 0.547715, 1.01639, 0.708948, 0.705313, 1.11718, 0.926838, 25.1704, 1.40895 }
};
*/
float nTrackCorr[40][30] = {
{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
{ 1, 2.23017, 1, 1, 1, 1, 1, 1.69306, 0.975091, 1.15069, 3.97529, 0.650689, 2.76227, 1.98488, 1.00154, 12.3129, 1.51535, 1, 0.46781, 5.7726, 1, 1, 2.68493, 1, 1, 1, 1, 1, 1, 1 },
{ 2.66905, 1, 0.93734, 1.57298, 0.458013, 1.27237, 2.78462, 0.438955, 1.24781, 1.24064, 0.930701, 1.24892, 2.01898, 1.73804, 1.47506, 1.14959, 1.39162, 1.27095, 2.29593, 1.11681, 1.05503, 0.86805, 0.689581, 0.785288, 8.36663, -1.62101, 0.944571, 1.12069, 1, 1 },
{ 1.99871, 0.0803983, 1.82899, 1.67475, 0.180479, 1.06515, 0.645474, 1.57004, 1.09973, 0.918369, 0.997164, 1.06179, 1.21731, 0.711921, 1.16836, 1.13635, 1.06544, 0.991982, 1.35473, 0.882152, 0.540863, 2.49356, 2.00226, 0.948838, 2.28405, 0.155179, 1.84119, 1.10915, 1, 1 },
{ 1.04156, 1.02718, 1.40583, 0.697813, 0.481316, 0.573089, 0.288591, 0.739858, 1.0188, 1.54772, 1.43032, 1.15419, 1.40738, 0.999839, 0.874045, 0.961913, 1.45834, 1.209, 1.83785, 1.43404, 1.64866, 0.855156, 0.721156, 1.26039, 0.888006, 1.93, 0.517283, 0.879001, 1, 1 },
{ 0.889731, 0.607019, 0.639195, 0.712941, 0.619153, 0.638272, 1.28463, 0.902472, 0.950162, 1.09325, 1.12862, 0.964537, 0.990211, 0.962309, 0.760681, 1.06416, 0.8636, 0.952407, 1.18674, 1.07298, 1.09773, 1.09944, 1.43493, 1.34338, 0.901385, 2.2744, 3.24304, 21.7546, 1, 1 },
{ 0.895483, 1.21694, 1.58825, 0.451607, 2.31901, 0.751027, 0.995733, 1.04511, 1.19839, 1.12815, 0.78901, 1.01618, 1.21431, 1.18144, 1.15757, 1.50077, 1.51267, 1.73521, 1.07113, 0.925381, 1.35874, 0.918259, 1.20101, 1.20745, 1.14326, 1.44312, 1.19774, 2.11212, 3.70476, 1 },
{ 1.12326, 0.668525, 0.878282, 0.834241, 1.10339, 0.9234, 0.995359, 0.869291, 0.867829, 1.30467, 1.00068, 0.875687, 1.25248, 0.921387, 1.02791, 1.15773, 1.12742, 0.963685, 1.1952, 1.59841, 1.30843, 1.56449, 1.27503, 0.777339, 1.38465, 0.948422, 1.6623, 1.46562, 1, 1 },
{ 1.55946, 0.894488, 0.719161, 1.16171, 0.820408, 0.962854, 1.09794, 0.851475, 1.08238, 1.02246, 0.971132, 0.970156, 1.22652, 1.22187, 1.11826, 1.02959, 1.24049, 1.27405, 1.15645, 0.904423, 1.32959, 1.0224, 1.42468, 1.19347, 1.68254, 0.934618, 1.99347, 1.08089, 1.43858, 1 },
{ 1.04951, 1.4636, 1.06154, 0.899051, 1.04932, 0.631111, 0.980467, 0.995527, 1.05695, 0.980391, 0.997697, 1.12322, 1.10451, 1.07044, 1.06699, 1.04037, 1.16519, 1.18667, 1.05412, 0.936167, 1.23671, 1.14809, 1.06852, 1.22281, 1.39947, 1.93296, 0.676461, 1.56061, 0.372632, 1 },
{ 1.33489, 0.795252, 0.972412, 0.94614, 0.737085, 0.725158, 1.00114, 1.09189, 0.831902, 1.21372, 1.29653, 1.03065, 1.19613, 1.02276, 0.93549, 0.95277, 1.02549, 1.06918, 1.41434, 1.024, 1.06631, 1.18517, 1.36758, 1.73765, 1.43606, 1.73225, 1.79506, 2.57247, 3.46892, 1 },
{ 1.10486, 0.978087, 0.892557, 1.19725, 0.94921, 1.17001, 0.929187, 1.23708, 1.0113, 1.07632, 1.09243, 1.00391, 0.792069, 1.02451, 1.03255, 1.117, 0.945157, 0.896579, 1.13751, 1.30914, 1.29165, 1.23647, 1.48589, 1.37539, 1.40677, 1.21986, 1.42516, 2.52961, 0.842845, 1 },
{ 1.0563, 0.883468, 0.93902, 0.984908, 1.00172, 0.935056, 0.934102, 0.907719, 0.864298, 1.0667, 1.01985, 1.02386, 1.1486, 1.07017, 1.05136, 1.15205, 0.962109, 1.34078, 1.06696, 1.34241, 1.07704, 1.16988, 1.16989, 1.12806, 1.10487, 0.848545, 0.995165, 1.31321, 0.587767, 1 },
{ 1.40159, 1.16059, 0.981103, 0.853363, 0.823852, 1.00747, 0.932628, 0.866965, 1.01275, 1.11716, 1.01643, 1.2519, 1.04179, 1.12883, 1.01788, 1.05554, 1.04445, 1.03354, 0.978113, 1.46924, 1.26874, 1.3246, 1.67875, 1.49001, 1.92282, 1.25925, 1.00322, 1.3323, 3.67021, 1 },
{ 1.18864, 1.11558, 0.995875, 0.904852, 0.773244, 0.913567, 0.904438, 1.02789, 1.02894, 1.00999, 0.94673, 1.17939, 1.15761, 1.0336, 1.04706, 1.09501, 1.00049, 1.10362, 1.07703, 1.35333, 1.53816, 1.03735, 1.63853, 1.52839, 1.31778, 1.85942, 1.96064, 1.6554, 5.64975, 1.06937 },
{ 1.18509, 1.39, 1.00102, 1.0373, 1.07231, 0.949113, 0.941183, 0.922163, 1.14979, 1.10183, 1.03953, 1.0619, 1.15074, 0.984025, 1.17604, 0.968733, 0.964389, 1.10274, 1.21868, 1.25452, 1.17255, 1.20597, 1.1085, 1.69185, 1.09196, 1.51471, 0.873209, 1.06374, 2.20896, 1 },
{ 1.0028, 1.20964, 0.972831, 0.80376, 0.962239, 1.1042, 0.882368, 0.756932, 1.14763, 1.00627, 1.0787, 1.10818, 0.976354, 0.981381, 1.1278, 1.06109, 1.01677, 1.07764, 1.20787, 1.1378, 1.24053, 1.56515, 1.35585, 1.41061, 1.62997, 0.883439, 1.30649, 1.58568, 1.91209, 1 },
{ 1.13096, 0.805966, 1.04795, 0.879954, 0.912062, 0.745231, 0.949462, 1.06183, 0.940619, 1.02456, 1.18896, 1.01759, 1.00746, 1.05736, 1.02931, 1.15586, 1.23385, 1.16338, 0.971354, 1.35476, 1.391, 1.24478, 1.36458, 1.37355, 1.12163, 1.12976, 1.38501, 1.06052, 0.27025, 1 },
{ 1.10193, 0.849379, 0.60802, 0.852398, 0.956145, 0.897476, 0.88525, 1.00474, 0.999108, 1.00203, 1.02052, 1.04811, 1.09618, 1.08749, 1.10486, 1.05315, 1.07784, 1.06946, 1.02599, 1.25281, 1.31745, 1.34183, 1.33061, 1.63331, 1.26817, 1.92456, 1.73085, 0.943338, 0.453117, 1 },
{ 0.786976, 0.796478, 0.717648, 0.991826, 0.821735, 0.950999, 1.00855, 0.884657, 1.01625, 1.03533, 1.0442, 1.05644, 1.14176, 1.02554, 1.0721, 1.03563, 1.0135, 1.01979, 1.19291, 1.26035, 1.43829, 1.39371, 1.33569, 1.38897, 1.79792, 1.25384, 1.47493, 1.58926, 1.22689, 1 },
{ 1.01265, 0.98698, 0.649373, 1.08557, 1.2146, 1.06611, 1.01971, 0.973621, 0.835304, 0.909653, 1.12116, 1.07194, 1.11374, 1.08588, 0.969273, 1.05339, 1.09444, 1.18047, 1.19151, 1.26369, 1.25726, 1.28096, 1.22206, 1.35743, 1.38603, 1.69748, 1.6033, 1.57362, 1.14764, 1 },
{ 0.993384, 0.910583, 1.07059, 1.09986, 0.880983, 0.856126, 0.964023, 0.93046, 0.94488, 1.14666, 1.17482, 0.965768, 0.940729, 1.02123, 1.02475, 1.20135, 0.986373, 1.25149, 0.979497, 1.44737, 1.15867, 1.24923, 1.15353, 1.43877, 1.18232, 1.1786, 1.24581, 1.34011, 0.385812, 1 },
{ 1.01648, 0.882821, 0.920423, 0.903265, 1.09529, 0.970349, 0.915489, 0.972953, 0.949632, 0.836879, 0.986931, 1.22905, 0.963988, 1.09156, 1.28141, 1.11032, 1.04545, 1.08678, 1.13207, 1.11565, 1.10627, 1.22213, 1.56404, 1.24451, 1.13122, 0.927497, 1.69173, 0.526548, 1.04507, 1 },
{ 0.657198, 1.06157, 0.87655, 1.00251, 1.04508, 0.768253, 0.769836, 0.776018, 1.07229, 0.883676, 1.08778, 0.9849, 1.0876, 1.05401, 0.783315, 0.977122, 0.910492, 1.12293, 1.04923, 1.25154, 1.32066, 1.14438, 1.16597, 1.05217, 1.36869, 1.25, 1.40417, 1.06336, -0.0507276, 1 },
{ 0.934537, 0.918278, 0.923641, 1.03266, 0.837168, 0.848878, 0.741994, 0.686517, 1.02489, 0.729337, 1.10628, 1.11698, 0.986989, 1.09672, 0.933008, 0.829321, 1.18216, 0.975896, 1.15174, 1.33009, 1.558, 1.21765, 1.16996, 1.22959, 1.38682, 1.85026, 1.58242, 1.71844, 0.697296, 1 },
{ 0.858734, 0.97751, 0.853295, 1.09774, 0.803811, 0.911339, 0.965327, 0.860426, 1.059, 0.939341, 1.14688, 0.849513, 0.995323, 1.00031, 1.05644, 1.03784, 0.889334, 1.07129, 1.0061, 1.26827, 1.25896, 1.17924, 1.21949, 1.31228, 1.3857, 1.50544, 1.94359, 1.19194, 0.484175, 1 },
{ 0.765657, 0.804183, 0.781777, 1.0195, 0.825488, 1.011, 0.867032, 0.72796, 0.870855, 1.10794, 1.14139, 1.13268, 1.15378, 0.957316, 1.05582, 0.948713, 1.03662, 1.10816, 0.793349, 1.14924, 1.09931, 1.3997, 1.61443, 1.55269, 1.33661, 1.34875, 1.33903, 0.962627, 0.0284541, 1 },
{ 0.907382, 0.931784, 1.13279, 0.691581, 0.916433, 0.618067, 0.788561, 0.707094, 0.834406, 0.97098, 1.0611, 1.2015, 0.860019, 0.928175, 0.909447, 0.80239, 0.97608, 0.951857, 0.971387, 0.888884, 1.34385, 0.912143, 1.32212, 1.45296, 1.92123, 0.983976, 0.939178, 1.73744, 1.75826, 1 },
{ 0.78153, 1.43538, 1.24228, 0.750424, 0.939412, 0.680748, 1.01648, 0.835815, 1.09839, 1.03668, 1.10123, 1.08138, 1.1155, 0.873493, 0.787814, 0.968528, 0.886903, 1.05036, 0.971161, 1.33896, 1.12316, 1.2123, 1.11016, 1.51182, 1.2559, 1.55083, 1.23959, 0.976591, 1.805, 1 },
{ 1.09069, 0.701507, 0.74124, 0.815021, 0.935916, 0.847031, 1.08517, 0.832034, 0.856937, 0.879777, 1.30117, 1.09589, 1.06318, 1.14986, 0.772011, 0.967189, 0.891932, 1.08888, 0.992017, 1.21122, 0.990114, 1.22789, 1.39881, 1.61552, 0.727558, 1.05979, 0.405245, 1.29183, 1.18609, -0.358535 },
{ 0.875387, 0.754856, 0.821365, 0.395063, 0.947101, 1.13079, 0.717402, 0.884681, 1.00426, 0.957824, 0.707513, 1.06531, 1.15407, 0.948619, 0.745968, 0.951327, 0.959118, 1.01944, 1.25783, 1.03904, 1.32914, 0.873359, 1.22042, 0.865437, 1.70116, 1.70147, 1.80971, 3.15834, 3.22856, 1 },
{ 2.77584, 0.701679, 0.910546, 0.83129, 0.654015, 0.738041, 0.847115, 0.852816, 0.779262, 0.837936, 0.980448, 0.823362, 0.88058, 0.729422, 0.939892, 1.04922, 1.04416, 0.729113, 1.15327, 1.23526, 1.00741, 1.27616, 0.782594, 1.61955, 1.50586, 1.18822, 1.5261, 1.81522, -0.273015, 1 },
{ 0.691531, 1.11288, 0.615842, 0.574422, 0.858997, 0.942352, 0.670554, 0.810021, 0.641311, 0.901952, 1.04334, 0.991078, 0.903071, 0.818815, 0.969196, 0.837101, 0.968136, 0.978101, 1.14892, 1.02856, 1.10002, 1.58612, 1.06567, 1.30024, 1.5052, 0.682549, 1.58697, 2.63386, 1.30547, 1 },
{ 1.24067, 0.867285, 0.821762, 0.51288, 0.767689, 0.592559, 0.531702, 0.444864, 0.984586, 0.982384, 1.38009, 0.776989, 0.628948, 0.575735, 0.928039, 0.728469, 0.963052, 0.681601, 0.825048, 1.07125, 1.20579, 1.44582, 1.10105, 1.43485, 2.5686, 1.40602, 1.51556, 1.3149, 0.871623, 1 },
{ 2.03763, 1.17927, 0.518766, 0.580949, 0.413369, 1.47607, 0.777546, 0.783802, 0.659492, 1.26274, 1.19845, 0.877409, 0.712508, 0.994451, 0.839, 0.822803, 1.19968, 0.825563, 1.0741, 0.900913, 1.12728, 1.19617, 2.24896, 1.45542, 1.28361, 0.373693, 0.701897, 0.785298, 1.54321, 1 },
{ 0.531422, 0.68418, 1.15416, 0.43741, 0.956526, 0.794399, 1.29894, 0.910657, 0.851456, 0.996965, 1.11662, 0.882106, 1.12877, 0.824031, 1.01866, 1.26515, 0.750627, 0.575754, 1.02355, 0.617995, 1.51638, 0.92236, -0.0523698, 0.48621, 2.23445, 0.49622, 2.15834, 1, 1, 1 },
{ 0.934124, 1.20761, 1.44011, 0.647147, 0.51338, 0.693244, 0.523223, 1.11143, 1.00966, 0.595848, 0.889168, 0.693047, 1.02195, 1.03551, 1.18437, 0.933441, 1.4622, 1.08945, 1.55386, 0.303153, 1.19353, 4.86698, 0.919036, 0.606542, 0.886229, 22.2344, 1.01314, 1, 1, 1 },
{ 0.19024, 1.41257, 1.52232, 1.02229, 0.684931, 0.774686, 1.15751, 1.28045, 1.19835, 1.24585, 0.880512, 1.02214, 1.02753, 0.84065, 0.78469, 0.813695, 0.822055, 2.01435, 1.81722, 1.7649, 0.765381, 1.69023, 0.968618, 1.51585, 2.30977, -0.176007, 4.71334, 1, 1, 1 },
{ 6, 0.728038, 2.19011, 0.969365, 1.4465, 1.7343, 1.58503, 2.6738, 1.08906, 0.989379, 0.933308, 1.11793, 0.708995, 1.06829, 0.923042, 1.47591, 1.02025, 0.898011, 0.465877, 1.2576, 1.01761, 1.09238, 0.821909, 0.546691, 1.33033, 1.05111, 1.70984, 1, 1.78571, 1 },
{ 8.10616, 0.728915, 1.11498, 1.89806, 4.1635, 0.182029, 1.6917, 0.972322, 0.796127, 0.835079, 1.04974, 1.38912, 0.431852, 0.443401, 1.11179, 1.38886, 0.560667, 0.519029, 1.66363, 2.35747, 3.60853, 4.1104, 1.95871, 1, 3.13119, 1, 1, 1, 1, 1 }};

float nTrackKaonCorr[10][5] = {
{ 1, 1, 1, 1, 1 },
{ 1.14204, 1, 1, 1, 1 },
{ 1.0242, 1.01284, 1, 1, 1 },
{ 0.964421, 0.980556, 1.05239, 1, 1 },
{ 0.947115, 0.963309, 0.986594, 1.03021, 1 },
{ 0.926197, 1.00214, 1.03134, 0.916187, 1.03245 },
{ 0.968412, 1.03123, 1.03675, 0.992196, 0.91806 },
{ 0.968853, 1.05726, 1.08222, 0.986642, 1.22805 },
{ 1.06092, 1.08837, 1.05504, 1.20801, 0.852954 },
{ 1.06717, 1.1358, 1.15011, 0.982436, 1.96327 }
};

float mytagcorr[30][30] = {
{ 1, 1, 1, 0.957851, 0.762565, 0.308832, 2.82979, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
{ 1, 1, 1, 1.04128, 1.32864, 0.460588, 0.395338, 0.294994, 1, 3.10522, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
{ 1, 1, 1, 3.89942, 0.302, 0.693749, 0.888405, 0.975708, 4.76059, 0.667499, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
{ 1, 1, 1, 0.820558, 0.990975, 1.16532, 1.07905, 0.751849, 1.10287, 0.938791, 0.627963, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
{ 1, 1, 1, 1.21713, 0.714611, 0.95437, 1.26466, 0.846641, 0.929532, 0.91447, 0.687965, 1.99946, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
{ 1, 1, 1, 1.94447, 0.761155, 0.640169, 0.97011, 0.856743, 0.867143, 0.955434, 0.98349, 1.07779, 1.35173, 0.298429, 1.57321, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
{ 1, 1, 1, 2.0261, 0.98186, 1.2227, 0.878563, 0.937386, 0.988487, 0.937718, 0.924071, 0.895853, 1.06326, 1.07499, 0.34473, 1.60531, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
{ 1, 1, 1, 1, 0.682935, 0.716778, 1.11226, 0.718046, 1.09568, 0.964922, 0.975243, 0.848627, 0.896181, 0.813142, 1.1638, 1.40832, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
{ 1, 1, 1, 1, 0.854559, 1.01545, 1.03956, 0.919203, 1.07059, 0.962376, 1.00874, 1.06683, 1.03572, 0.987928, 0.906489, 0.969922, 1.19436, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
{ 1, 1, 1, 1, 0.732747, 1.06845, 0.975406, 0.92788, 0.868253, 0.899763, 0.910892, 1.02517, 0.908673, 0.914922, 1.00906, 0.947738, 1.00682, 0.705516, 0.124844, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
{ 1, 1, 1, 1, 1, 0.685801, 0.911891, 0.903289, 0.979109, 1.05455, 1.03241, 1.00236, 0.983108, 0.926264, 0.943233, 0.945488, 0.949854, 0.970668, 0.808121, 4.00132, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
{ 1, 1, 1, 1, 1, 0.707279, 0.815855, 0.853189, 0.925017, 0.972558, 1.00536, 1.00471, 1.03817, 0.927479, 0.961747, 0.925044, 0.970784, 0.880202, 0.954319, 0.892344, 0.380883, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
{ 1, 1, 1, 1, 1, 2.062, 2.57314, 0.863909, 0.988228, 0.929596, 0.926323, 1.02616, 1.03168, 1.0268, 1.00793, 0.98355, 0.98829, 1.0101, 0.847943, 0.875608, 0.897768, 2.71475, 1, 1, 1, 1, 1, 1, 1, 1 },
{ 1, 1, 1, 1, 1, 1, 0.24657, 0.490541, 0.744362, 1.0339, 0.969763, 0.943875, 1.01059, 1.07484, 1.07148, 0.971739, 0.969007, 1.02416, 1.03554, 0.990036, 0.855166, 0.822697, 1.32949, 1, 1, 1, 1, 1, 1, 1 },
{ 1, 1, 1, 1, 1, 1, 0.544248, 0.524124, 0.987559, 0.992573, 1.12012, 0.9982, 0.963562, 1.0148, 1.07455, 1.04141, 1.00234, 0.982616, 1.04362, 0.913125, 0.935396, 0.798635, 0.986336, 1.50938, 1, 1, 1, 1, 1, 1 },
{ 1, 1, 1, 1, 1, 1, 1, 0.954721, 1.66658, 0.976613, 0.860534, 1.09697, 0.974189, 0.98888, 0.972604, 1.06635, 1.08213, 0.97413, 1.0148, 0.982081, 0.975279, 0.914971, 1.00607, 0.682009, 1.25128, 1, 1, 1, 1, 1 },
{ 1, 1, 1, 1, 1, 1, 1, 0.588484, 1.57916, 1.16832, 1.2891, 1.33079, 1.2097, 1.09922, 1.05896, 0.979221, 1.0363, 1.09918, 1.00424, 0.910763, 0.990076, 1.0363, 0.834168, 0.989907, 0.666844, 1.56962, 1, 1, 1, 1 },
{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.971606, 0.861461, 1.12812, 1.01724, 1.13819, 1.07286, 0.979057, 1.01496, 1.07276, 1.02426, 0.975056, 0.824724, 0.884919, 0.771595, 0.775165, 0.762535, 1, 1, 1, 1 },
{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.724111, 1.2673, 1.29987, 1.16436, 1.12448, 1.14339, 1.05009, 0.976588, 0.998951, 1.06582, 1.01071, 0.991485, 0.960445, 0.927855, 1.2094, 0.57015, 0.125019, 1, 1, 1 },
{ 1, 1, 1, 1, 1, 1, 1, 4.86907, 1, 1, 1, 1, 0.873486, 0.989943, 1.01609, 1.1346, 1.11449, 1.0602, 0.998078, 0.986498, 1.04615, 1.03742, 0.976413, 1.02349, 0.877824, 0.874217, 0.500611, 1, 1, 1 },
{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.944452, 0.787161, 1.10582, 1.04748, 1.0327, 1.10717, 1.08182, 0.988003, 0.969119, 1.01929, 1.06024, 0.943775, 1.01739, 0.880665, 0.991439, 1, 1, 1 },
{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 6.18012, 1.27244, 1, 1.62396, 1.32071, 1.01583, 1.21617, 1.08081, 1.08147, 1.04377, 0.942071, 1.00959, 0.99575, 0.940494, 0.87069, 0.568619, 0.754432, 1, 1 },
{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.542, 1.61521, 0.885397, 1.10562, 1.07483, 1.14119, 1.02504, 0.936339, 0.946152, 1.01715, 0.934864, 1.00256, 1, 1, 1 },
{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.265412, 1, 1, 1, 1.60231, 2.16649, 0.768723, 1.08281, 1.00487, 1.07933, 0.988354, 0.941322, 0.951246, 0.961957, 0.756765, 4.61486, 1, 1 },
{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.542818, 0.856728, 1.365, 1.25365, 0.953549, 1.06271, 1.09395, 0.919543, 0.99057, 0.99725, 0.592598, 1, 1 },
{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.687257, 1, 0.616769, 0.0797017, 1.59058, 0.98666, 1.04845, 0.996098, 1.04918, 0.886095, 0.831867, 1.17522, 1.28021, 1 },
{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.105557, 0.770693, 1.24728, 0.895165, 0.948774, 1.15405, 0.900326, 0.946786, 1.18703, 1 },
{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.70714, 0.820833, 0.781167, 0.963826, 0.910623, 0.606199, 1.16651, 1.64882, 1 },
{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.28517, 1.64651, 0.979808, 0.84367, 0.73945, 0.661985, 1 },
{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.41254, 1, 4.2208, 0.754966, 0.3952, 5.29143, 1 }};


float nTrack[] = {
  1.16408,
 1.00155,
 0.960201,
 0.980238,
 0.979162,
 0.997814,
 1.03451,
 1.11787,
 1.30359,
 1.25997
};


float tagCorr1[] = {
0.90839,
 1.0753,
 1.27674,
 1.00508,
 0.943151,
 1.13422,
 1,
 1,
 1,
 1,
 1,
 1,
 0.936078,
1,
 0.96263,
1,
 0.910089
};

float tagCorr0[] = {
1, 
 1.07067,
 1.0752,
 1,
 1.00911,
 0.830767,
 1,
 1.19175,
 1.04741,
 1.35003,
 1,
 1.0366,
 1,
 1,
 1.43072
};
void FillHist()
{

cout<<"Fill Hists..."<<endl;

	TTree *t = NULL;
	TTree *t_friend = NULL;
	TTree *t_new = NULL;
	TFile *f_new = NULL;
	std::vector<float> *vtrk_p = new std::vector<float>;
	std::vector<float> *vtrk_p1 = new std::vector<float>;
	std::vector<float> *vtrk_p2 = new std::vector<float>;
	std::vector<float> *vtrk_p3 = new std::vector<float>;
	std::vector<float> *vtrk_e = new std::vector<float>;
	std::vector<float> *vtrk_th = new std::vector<float>;
	std::vector<float> *vtrk_dr = new std::vector<float>;
	std::vector<float> *vtrk_dz = new std::vector<float>;
	std::vector<float> *vtrk_kpi_id = new std::vector<float>;
	std::vector<int> *vtrk_sep = new std::vector<int>;
	std::vector<int> *vtrk_q = new std::vector<int>;
	std::vector<float> *cal_e = new std::vector<float>;
	std::vector<float> *cal_th = new std::vector<float>;
	std::vector<double> *efc_e = new std::vector<double>;
	std::vector<double> *efc_th = new std::vector<double>;
	std::vector<float> *gam_p1 = new std::vector<float>;
	std::vector<float> *gam_p2 = new std::vector<float>;
	std::vector<float> *gam_p3 = new std::vector<float>;
	std::vector<float> *vtrk_gen_p1 = new std::vector<float>;
	std::vector<float> *vtrk_gen_p2 = new std::vector<float>;
	std::vector<float> *vtrk_gen_p3 = new std::vector<float>;
	std::vector<float> *vtrk_gen_e = new std::vector<float>;
	std::vector<int> *vtrk_gen_id = new std::vector<int>;
	std::vector<int> *vtrk_gen_mo_id = new std::vector<int>;
	std::vector<int> *vtrk_gen_ist = new std::vector<int>;


	TH1F *BDecays = new TH1F("bdecays","bdecays",10000,0,10000);
	TH1F *DDecays = new TH1F("ddecays","ddecays",10000,0,10000);

	char sName[200];
	TH1F **BDecaysBins = new TH1F*[12];
	TH1F **DDecaysBins = new TH1F*[12];
	for(int i = 0;i<12;i++)
	{

		sprintf(sName,"BdecaysBins_%i",i);
		BDecaysBins[i] = new TH1F(sName,sName,10000,0,10000);
		sprintf(sName,"DdecaysBins_%i",i);
		DDecaysBins[i] = new TH1F(sName,sName,10000,0,10000);
	}
	TH1F *hWeight1 = new TH1F("hWeight","hWeight",16,0,15);
	TH1F *hWeight2 = new TH1F("hWeight2","hWeight2",16,0,15);

	float nEvents[nContributions];
	float nEventsSplit[2][2][nContributions];
	float nEventsSigReg[nContributions];
	float nEventsPerFile[nContributions];
	for(int i=0;i<nContributions;++i) nEvents[i] = 0;
	for(int i=0;i<nContributions;++i) {
		nEvents[i] = 0;
		nEventsSigReg[i] = 0;
		for(int j=0;j<2;j++){
			for(int k=0; k<2;k++){
				nEventsSplit[j][k][i] = 0;
			}
		}
	}
	
	cout<<"test"<<endl;
	float nEventsBCuts[2][4][nContributions];
	for(int i=0;i<2;++i) {
		for(int j=0;j<4;++j) {
			for(int k=0;k<nContributions;++k) {
				nEventsBCuts[i][j][k] = 0;
			}
		}
	}

	int double_events = 0;
	float MC_Correction;


	SemiLepWeights2D slWeight = SemiLepWeights2D();
	SemiLepDssWeights2D slDssWeight = SemiLepDssWeights2D();
	LepFake MuonFake1 = LepFake(1,"muon fake1");
	LepFake ElecFake1 = LepFake(1,"elec fake1");
	LepFake MuonFake2 = LepFake(1,"muon fake2");
	LepFake ElecFake2 = LepFake(1,"elec fake2");

	LepEff MuonEff1 = LepEff(1,"muon eff1",0.9,1);
	LepEff ElecEff1 = LepEff(1,"elec eff1",0.5,0);


	KPiEff KPiEffFake = KPiEff(0.1,1);
	KPiEff PiKEffFake = KPiEff(0.9,1);

	CorrectBf CorrBf = CorrectBf();

	int countD1=0;

	int fCount =0;

	TH2F** h2d = new TH2F*[5];
	struct {int bin; float min; float max;} xa,ya;
	xa.bin = 30; xa.min = 0.; xa.max = 3;
	xa.bin = 30; xa.min = 0; xa.max = 0.5;
	//xa.bin = 20; xa.min = -1.; xa.max = 1;
//	xa.bin = 50; xa.min = 0; xa.max = 2.5;
//	xa.bin = 10; xa.min = -.5; xa.max = 9.5;
//	xa.bin = 30; xa.min = 5.4; xa.max = 6.;
	ya.bin = 20; ya.min = 0.; ya.max = 7;
	ya.bin = 20; ya.min = 5.27; ya.max = 5.29;
	//ya.bin = 5; ya.min = -.0; ya.max = 4.5;
	ya.bin = 50; ya.min = 1; ya.max = 2.5;
	ya.bin = 50; ya.min =-0.15; ya.max =.15;
	ya.bin = 50; ya.min =0; ya.max =5;
	ya.bin = 10; ya.min = 0.5; ya.max = 10.5;
	ya.bin = 30; ya.min = -1; ya.max = 10;
	//ya.bin = 20; ya.min = -1.; ya.max = 1;
	h2d[0]= new TH2F("h2d0","h2d0",xa.bin,xa.min,xa.max,ya.bin,ya.min,ya.max);
	h2d[1]= new TH2F("h2d1","h2d1",xa.bin,xa.min,xa.max,ya.bin,ya.min,ya.max);
	h2d[2]= new TH2F("h2d2","h2d2",xa.bin,xa.min,xa.max,ya.bin,ya.min,ya.max);
	h2d[3]= new TH2F("h2d3","h2d3",xa.bin,xa.min,xa.max,ya.bin,ya.min,ya.max);
	h2d[4]= new TH2F("h2d4","h2d4",xa.bin,xa.min,xa.max,ya.bin,ya.min,ya.max);

	for(vector<CFile>::iterator it_f = vFile.begin(); it_f != vFile.end(); ++it_f)
	{

		for(int i = 0; i<nContributions; i++)
			nEventsPerFile[i] = 0;
		/*TFile f( it_f -> FileName.c_str() ,"read");
		if(fCount>=4 && useDilep)
			t = (TTree*)f.Get("dilep/tree");
		else
			t = (TTree*)f.Get("lep/tree");
		if(!t) cout<<"no tree"<<endl;
		*/fCount++;
        cout<<it_f->Path<<" "<<it_f->Name<<" "<<it_f->Tree<<" "<<it_f->FileName<<" "<<endl;

        TFile f( it_f -> FileName.c_str() ,"read");
        t = (TTree*)f.Get(it_f->Tree.c_str());
        cout<<t<<endl;

		TTreePerfStats *ps =
			   new TTreePerfStats(it_f->Name.c_str(),t);


        if(_newTree)
        {
            f_new = new TFile(it_f->Name.c_str(),"RECREATE","0");
            t_new = new TTree("tree","tree");
            t_new->Branch("run", &run, "etot/F:eher:eler:eth:exp/i:num/i");
            t_new->Branch("evtnr", &evtnr, "evtnr/i");

            t_new->Branch("btag",   &btag,   	"pcode_b/f:b_mode:sub1_mod:sub2_mod:sub3_mod:sub4_mod:MCinfo:NB:contNB:m_bc:delta_e:cos_thr:thr:p:p1:p2:p3:e:NFS/I:fl_tag/I");
            t_new->Branch("lep1", &lep1,		"elid/f:muid:atckpi:atcpk:atcppi:e/f:p/f:th:dr:dz:es:ps:ths:ec:pc:thc:id_truth/I:q/I:mother_id/I:fl_lep:fl_mother");
            t_new->Branch("gmiss", &gmiss,	"e:p:th:m2:es:ps:ths:ec:pc:thc");
            t_new->Branch("gx", &gx,		"e:p:th:m2:es:ps:ths:ec:pc:thc");
            t_new->Branch("event", &event,	"Eecl:q2:cos_thrA:cos_thrA2:cos_thrAm:cos_thrB:cos_thrC:thrX:thrSigX:thrSigX2:thrSigXm:nch/I:nn/I:npi0/I:q/I:nk:lclass:dclass");
            t_new->Branch("truth", &truth,	"q2/F:p_l:w:costh:m:th:thc:ths:p:pc:ps:miss_m2:miss_e:miss_p:miss_th:miss_es:miss_ps:miss_ths:nkl/I:nch/I");
            t_new->Branch("sBDecay",		&sBDecay, "sBDecay/C");
            t_new->Branch("sDDecay",		&sDDecay, "sDDecay/C");
            t_new->Branch("thXMiss", &fThXMiss);

        }

		t->SetBranchAddress("run", &run);
		t->SetBranchAddress("evtnr", &evtnr);
		t->SetBranchAddress("btag", &btag);
		t->SetBranchAddress("lep1", &lep1);
//		t->SetBranchAddress("lep2", &lep2);
//		t->SetBranchAddress("dilep", &dilep);
		t->SetBranchAddress("x", &x);
		t->SetBranchAddress("gx", &gx);
		t->SetBranchAddress("xlep", &xlep);
		t->SetBranchAddress("event", &event);
		t->SetBranchAddress("gmiss", &gmiss);
		t->SetBranchAddress("miss", &miss);
	//	t->SetBranchAddress("BDT1", &BDT1);
	//	t->SetBranchAddress("BDT2", &BDT2);
		t->SetBranchAddress("truth", &truth);

		t->SetBranchAddress("trk_p", &vtrk_p);
		t->SetBranchAddress("trk_e", &vtrk_e);
		t->SetBranchAddress("trk_p1", &vtrk_p1);
		t->SetBranchAddress("trk_p2", &vtrk_p2);
		t->SetBranchAddress("trk_p3", &vtrk_p3);
		t->SetBranchAddress("trk_th", &vtrk_th);
		t->SetBranchAddress("trk_sep", &vtrk_sep);
		t->SetBranchAddress("trk_kpi_id", &vtrk_kpi_id);
		t->SetBranchAddress("trk_dr", &vtrk_dr);
		t->SetBranchAddress("trk_dz", &vtrk_dz);
		t->SetBranchAddress("cal_e", &cal_e);
		t->SetBranchAddress("cal_th", &cal_th);
		t->SetBranchAddress("cal_p1", &gam_p1);
		t->SetBranchAddress("cal_p2", &gam_p2);
		t->SetBranchAddress("cal_p3", &gam_p3);
		t->SetBranchAddress("trk_gen_e", &vtrk_gen_e);
		t->SetBranchAddress("trk_gen_p1", &vtrk_gen_p1);
		t->SetBranchAddress("trk_gen_p2", &vtrk_gen_p2);
		t->SetBranchAddress("trk_gen_p3", &vtrk_gen_p3);
		t->SetBranchAddress("trk_gen_is", &vtrk_gen_ist);
		t->SetBranchAddress("trk_gen_id", &vtrk_gen_id);
		t->SetBranchAddress("trk_gen_mo_is", &vtrk_gen_mo_id);
		/*
		t->SetBranchAddress("EFC_energy", &efc_e);
		t->SetBranchAddress("EFC_theta", &efc_th);
		t->SetBranchAddress("Ecl_energy", &cal_e);
		t->SetBranchAddress("Ecl_theta", &cal_th);
*/
		t->SetBranchAddress("sBDecay", &sBDecay);
		t->SetBranchAddress("sDDecay", &sDDecay);

		t->SetBranchStatus("*",0);
		t->SetBranchStatus("run",1);
		t->SetBranchStatus("evtnr",1);
		t->SetBranchStatus("btag",1);
		t->SetBranchStatus("lep1",1);
		t->SetBranchStatus("gx",1);
		t->SetBranchStatus("x",1);
		t->SetBranchStatus("event",1);
		t->SetBranchStatus("gmiss",1);
		t->SetBranchStatus("miss",1);
		t->SetBranchStatus("truth",1);
		t->SetBranchStatus("trk_p",1);
		t->SetBranchStatus("trk_p1",1);
		t->SetBranchStatus("trk_p2",1);
		t->SetBranchStatus("trk_p3",1);
		t->SetBranchStatus("trk_e",1);
		t->SetBranchStatus("trk_th",1);
		t->SetBranchStatus("trk_sep",1);
		t->SetBranchStatus("trk_kpi_id",1);
		t->SetBranchStatus("trk_dr",1);
		t->SetBranchStatus("trk_dz",1);
		t->SetBranchStatus("cal_e",1);
		t->SetBranchStatus("cal_th",1);
		t->SetBranchStatus("cal_p1",1);
		t->SetBranchStatus("cal_p2",1);
		t->SetBranchStatus("cal_p3",1);
		t->SetBranchStatus("sBDecay",1);
		t->SetBranchStatus("sDDecay",1);
		t->SetBranchStatus("trk_gen_p1",1);
		t->SetBranchStatus("trk_gen_p2",1);
		t->SetBranchStatus("trk_gen_p3",1);
		t->SetBranchStatus("trk_gen_e",1);
		t->SetBranchStatus("trk_gen_id",1);
		t->SetBranchStatus("trk_gen_mo_is",1);
		t->SetBranchStatus("trk_gen_is",1);

		t->SetCacheSize(500000000);
		t->AddBranchToCache("*");

		int evtnr_old =-1,num_old=-1;

		cout<<"Read "<<it_f->FileName<<" with "<<t->GetEntries()<<" Events and weight "<<it_f->weight<<endl;
		for(int i =0; i<t->GetEntries(); ++i)
		{
			t->GetEntry(i);
		//	if(i>100000) break;

			if(i%10000 == 0) cout<<"\r"<<i*100/t->GetEntries()<<flush;
		//	if(run.exp>37) continue;
			if( btag.m_bc<5.27) continue;
			if( log(btag.NB)<-4 ) continue;
			//if(lep1.fl_lep%10 == 0 && lep1.elid<0.9) continue;
			//if( log(btag.contNB)<-3 ) continue;
//			if(!((sqrt(gx.m2)>1.85 && sqrt(gx.m2)<1.89) || (sqrt(gx.m2)>1.98 && sqrt(gx.m2)<2.03)))continue;
			//if( log(btag.NB)<-4) continue;
			//if( run.etot>11.495) continue;
			//if(btag.p<0.3) continue;
			//if( abs(event.q)==0 ) continue;
			if(lep1.ps <1) continue;
			//if( gmiss.m2>0. ) continue;
			//if(event.cos_thrA2 > 0.9) continue;
			

				//if(lep1.elid<0) continue;
			//if(miss.e < 0) continue;

			//if(event.nch + btag.NFS < 10 || event.nch + btag.NFS > 13) continue;
			vtrk_q->clear();
			for(vector<int>::iterator itsep = vtrk_sep->begin(); itsep!=vtrk_sep->end(); ++itsep)
			{
				//cout<<*itsep<<" ";
				vtrk_q->push_back(sign(*itsep));
				if(it_f->Type != 'd')
					*itsep = abs(*itsep) - 11;
				else *itsep = abs(*itsep) -1;
				//cout<<*itsep<<endl;
			}



			int iLep = lep1.fl_lep%10;
			//if(iLep!=0) continue;

			//if(abs(btag.pcode_b) == 511) continue;



			fb_mode = 2*(fabs(btag.pcode_b) - 511) + fabs(getnum(btag.b_mode));
			/*if(
					fb_mode== 1 ||
					fb_mode== 4 ||
					fb_mode== 7 ||
					fb_mode== 11 ||
					fb_mode== 13 ||
					fb_mode== 14 ||
					fb_mode== 27 ||
					fb_mode== 28 ||
					fb_mode== 29 ||
					fb_mode== 30 ||
					fb_mode== 31 ||
					fb_mode== 32 ||
					fb_mode== 34 ||
					fb_mode== 36
					) continue;*/
	//excl. signal region
//		if(lep1.ps<1.2 ) continue;

		//	if(sqrt(xlep.m2)<2.) continue;
			//if(lep1.ps<1.2 && sqrt(xlep.m2)<2.) continue;
			//signal region
			//if(lep1.ps>1.1) continue;
			//if(sqrt(xlep.m2)<2.3) continue;
			//if(sqrt(xlep.m2)>6) continue;

			string BDecay(sBDecay);
			string DDecay(sDDecay);
			string bad_part("10042");
			//if(BDecay.find(bad_part)!=string::npos) continue;





			
			
			


			//remove D**lnu from generic
			if( it_f->Type == 'o'){ if(  abs(event.lclass) == 2 && event.dclass > 4) continue; }
			// add only true D**lnu events
			else if( it_f->Type == 's'){ if(!(abs(event.lclass) == 2 && event.dclass > 4))continue;}
			if(abs(event.lclass == 2) && event.dclass == 5) countD1++;





			int nkp = event.nk%10;
			int nkm = event.nk/10;
			if(lep1.q > 0)
			{
				fnkss = (float)nkp;
				fnkos = (float)nkm;
			}else
			{
				fnkos = (float)nkp;
				fnkss = (float)nkm;
			}
			fnk = fnkos + fnkss;
			

			int cont = GetCont(it_f);
			int cont_old =cont;

			//if (cont != charm_lep  && cont != charm_s_lep && cont != charm_ss_lep) continue;

			MC_Correction = 1;
            for(int i =0; i<8; i++)
				 Corrections[i] = -1;
			if(cont != data )
			{



				if(it_f->Type != 'c')
					Corrections[0] = tagcorr(btag.b_mode, (double)btag.NB);
				else Corrections[0]=0.7;


				if(abs(event.lclass) == 2)
				{
					if(event.dclass<5)
					{
				       		Corrections[2]= slWeight.weight(dclass2lclass[event.dclass],truth.w, truth.costh,truth.q2,truth.p_l);
					}
					else
					{
						Corrections[3] = slDssWeight.weight(event.dclass-5,truth.w, truth.costh);

					}
				}
			//  MC_Correction *= BrCorrection(event.lclass, event.dclass, sBDecay, sDDecay);
				if(abs(event.lclass) == 2 && event.dclass<9)
				{
				       	Corrections[1] = CorrBf.weightB(event.lclass, event.dclass);
			//cout<<event.dclass<<" "<<CorrBf.weightB(event.lclass, event.dclass)<<endl;
			}
			else
			{
				Corrections[1] = CorrBf.weightB(sBDecay);

			}
			Corrections[1] *= CorrBf.weightD(sDDecay);


			//lep1
			//electrons
			if(lep1.fl_lep==10)//true e
				{Corrections[4]= ElecEff1.weight(lep1.fl_lep%10,run.exp,acos(lep1.th)/M_PI*180., lep1.p);}// cout<<Corrections[4]<<endl;}
			else if(lep1.fl_lep%10 == 0 )//fake e
				Corrections[5]= ElecFake1.weight(lep1.fl_lep%10,lep1.fl_lep/10-1,acos(lep1.th)/M_PI*180., lep1.p);

			//muons
			else if(lep1.fl_lep==21)//true mu
				Corrections[6]= MuonEff1.weight(lep1.fl_lep%10,run.exp,acos(lep1.th)/M_PI*180., lep1.p);
			else if(lep1.fl_lep%10 == 1 )//fake mu
				Corrections[7]= MuonFake1.weight(lep1.fl_lep%10,lep1.fl_lep/10-1,acos(lep1.th)/M_PI*180., lep1.p);
			else
				cout<<"Error: No Lepton\n";




            //kid
			std::vector<float>::iterator itp1 = vtrk_p1->begin();
			std::vector<float>::iterator itp2 = vtrk_p2->begin();
			std::vector<float>::iterator itp3 = vtrk_p3->begin();
			std::vector<float>::iterator ite = vtrk_e->begin();
			std::vector<int>::iterator itsep = vtrk_sep->begin();
			TLorentzVector lvTmp;
            		float fKid=1;
			for(; itp1!=vtrk_p1->end();++itp1, ++itp2, ++itp3, ++ite, ++itsep)
			{
				lvTmp =  TLorentzVector(*itp1, *itp2, *itp3, *ite);

                if(*itsep%10 >= 5) continue;
			    if(*itsep==22)
				fKid *= KPiEffFake.weight(0,lvTmp.Theta()/M_PI*180., lvTmp.Rho(),run.exp);
			    else if(*itsep%10 == 2)
				fKid *= KPiEffFake.weight(1,lvTmp.Theta()/M_PI*180., lvTmp.Rho(),run.exp);
			    else if(*itsep == 33 )
				fKid *= PiKEffFake.weight(2,lvTmp.Theta()/M_PI*180., lvTmp.Rho(),run.exp);
			    else if(*itsep%10 == 3)
				fKid *= PiKEffFake.weight(3,lvTmp.Theta()/M_PI*180., lvTmp.Rho(),run.exp);
            }


//            MC_Correction*=Corrections[0];
				for(int i = 0; i<8; i++)
				{
					if(Corrections[i]>0)
						MC_Correction*=Corrections[i];
				}
       //        		 MC_Correction*=fKid;



			} //end of Mc corection
MC_Correction = 1;
			float theta = 0.022;
//			TLorentzVector lvBeam( run.eher*sin( theta ), 0.0, run.eher*cos(theta ) - run.eler, run.eher + run.eler );
			TLorentzVector lvBeam( run.eher*sin( theta ), 0.0, run.eher*cos(theta ) - run.eler, run.eher + run.eler );
			fMy4s = lvBeam.M();
			TLorentzVector lvBreco(btag.p1,btag.p2,btag.p3,btag.e);
//			lvBreco.Boost(-lvBeam.BoostVector());
//            lvBreco.SetRho(sqrt(-sqr(5.2796)+sqr(fMy4s/2.)));
//            lvBreco.SetE(fMy4s/2.);
//            lvBreco.Boost(lvBeam.BoostVector());

//			btag.delta_e-=lvBreco.E();
//			if(cont!=data) lvBreco.SetE( lvBreco.E()-0.0009);
//			btag.delta_e +=lvBreco.E();
//			lvBreco.SetRho(lvBreco)

			TRandom3 gRand;
			TLorentzVector lvPmisswTag = lvBeam - lvBreco;
			TLorentzVector lvTracks;
			std::vector<float>::iterator itp1 = vtrk_p1->begin();
			std::vector<float>::iterator itp2 = vtrk_p2->begin();
			std::vector<float>::iterator itp3 = vtrk_p3->begin();
			std::vector<float>::iterator ite = vtrk_e->begin();
			std::vector<float>::iterator itdr = vtrk_dr->begin();
			std::vector<float>::iterator itdz = vtrk_dz->begin();
			std::vector<int>::iterator itsep = vtrk_sep->begin();
			std::vector<float>::iterator itkpid = vtrk_kpi_id->begin();
			std::vector<int>::iterator itq = vtrk_q->begin();
            		std::vector<MyParticle> vParticles;
			vParticles.clear();
			std::vector<float>::iterator itp1_gen = vtrk_gen_p1->begin();
			std::vector<float>::iterator itp2_gen = vtrk_gen_p2->begin();
			std::vector<float>::iterator itp3_gen = vtrk_gen_p3->begin();
			std::vector<float>::iterator ite_gen = vtrk_gen_e->begin();
			std::vector<int>::iterator itid_gen = vtrk_gen_id->begin();
			std::vector<int>::iterator itis_gen = vtrk_gen_ist->begin();
			std::vector<int>::iterator it_moid_gen = vtrk_gen_mo_id->begin();
			int nTracks = 0;
			int nTracksPoor = 0;
			TLorentzVector lvTmp;
			int nSlow = 0;
			int nFake=0;
			int nFakeK=0;
			event.nch = 0;
			int SignalLep = 0;
			float fKid = 1;
			enum {t_SigLep = 1, t_Ghost};
			fnk = fnkos = fnkss = 0;
			TLorentzVector lvSignal;
			SumPsig =0;
			int nPi = 0;
			//cout<<event.nch<<endl;
			for(; itp1!=vtrk_p1->end();++itp1, ++itp2, ++itp3, ++ite, ++itsep,++itdr, ++itdz, ++itq,++itp1_gen, ++itp2_gen, ++itp3_gen, ++ite_gen, ++itid_gen, ++itis_gen, ++it_moid_gen,++itkpid)
			{

				//if(cont!=data && gRand.Uniform(1.)<0.37e-2) continue;

				MyParticle mp_tmp;
				mp_tmp.pid = *itsep%10;
				mp_tmp.pid_true = int(*itsep/10);
				mp_tmp.charge = *itq;
				mp_tmp.dr = *itdr;
				mp_tmp.dz = *itdz;
				mp_tmp.P = TLorentzVector(*itp1,*itp2,*itp3,*ite);
				
//mp_tmp.P = TLorentzVector(*itp1,*itp2, 0 , sqrt(mp_tmp.P.M() + sqr(*itp1) + sqr(*itp2) );
//                if(cont!=data && mp_tmp.P.Rho()<0.5)
//                {
//                    float m2 = mp_tmp.P.M2();
//                    mp_tmp.P.SetRho(mp_tmp.P.Rho()+0.01);
////                    mp_tmp.P.SetZ(mp_tmp.P.Z()+0.00/sqr(mp_tmp.P.Rho()));
//                    mp_tmp.P.SetE(sqrt(sqr(mp_tmp.P.Rho()) + m2));
//                }
				
				if(mp_tmp.pid >=5) continue;
				//if(mp_tmp.P.Perp()<0.1) continue;
				//if(abs(mp_tmp.dr)>.5) continue;
				//if(abs(mp_tmp.dz)>1.5) continue;
				if(mp_tmp.pid>1 && mp_tmp.pid<5)
				{
				if(*itkpid>0.9) {
					mp_tmp.pid = 2;
					mp_tmp.P = TLorentzVector(*itp1,*itp2,*itp3,sqrt(sqr(mp_tmp.P.Rho()) + sqr(.493)));
				}
				else if(*itkpid<0.9) {
					mp_tmp.pid = 3;
					mp_tmp.P = TLorentzVector(*itp1,*itp2,*itp3,sqrt(sqr(mp_tmp.P.Rho()) + sqr(.139)));
					nPi++;
					}
				else continue;
				}
				//if( (mp_tmp.P.CosTheta()) > 0.956 || mp_tmp.P.CosTheta()<-0.866 ) continue;


					TLorentzVector tmp = mp_tmp.P;
					//tmp.Boost(-lvBeam.BoostVector());
					SumPsig += tmp.Rho();
				mp_tmp.gen_id = *itid_gen;
				mp_tmp.gen_ist = *itis_gen;
				mp_tmp.gen_mo_id = *it_moid_gen;
				mp_tmp.genP = TLorentzVector(*itp1_gen, *itp2_gen,0, *ite_gen);
				
				
				
				if(mp_tmp.pid == 3 && mp_tmp.pid_true != 3 && mp_tmp.P.Rho()<.3 ) nFake++;
				//if(mp_tmp.pid != 2 && mp_tmp.pid_true == 2 && mp_tmp.P.Rho()<1 ) MC_Correction*=1.1;

				if(mp_tmp.pid == 0 || mp_tmp.pid ==1)
				{
					mp_tmp.flag  = 1;
					SignalLep++;
					lvSignal = mp_tmp.P;
				}
				else
				{
					
					mp_tmp.flag = 0;
					//if(*itkpid<0.95 && *itkpid > 0.05) continue;
					//if(cont!=data && *itkpid) mp_tmp.P = TLorentzVector(*itp1,*itp2,*itp3,sqrt(pow(mp_tmp.P.Rho(),2) + pow(0.49,2)));
				}
				if(mp_tmp.pid == 2)
				{
                    fnk++;
                    mp_tmp.charge*lep1.q>0?fnkss++:fnkos++;
				}

				event.nch++;
				vParticles.push_back(mp_tmp);

				lvTmp =  mp_tmp.P;
				lvTracks += lvTmp;
				nTracks++;
				    if(*itsep==22)
					fKid *= KPiEffFake.weight(0,mp_tmp.P.Theta()/M_PI*180., mp_tmp.P.Rho(),run.exp);
				    else if(*itsep%10 == 2)
					fKid *= KPiEffFake.weight(1,mp_tmp.P.Theta()/M_PI*180., mp_tmp.P.Rho(),run.exp);
				    else if(*itsep == 33 )
					fKid *= PiKEffFake.weight(2,mp_tmp.P.Theta()/M_PI*180., mp_tmp.P.Rho(),run.exp);
				    else if(*itsep%10 == 3)
					fKid *= PiKEffFake.weight(3,mp_tmp.P.Theta()/M_PI*180., mp_tmp.P.Rho(),run.exp);

			}
			if(!SignalLep) continue;
			if(fnk>0) continue;
			//if(nPi<2) continue;
			//if(event.nch != 3) continue;
			
			//if(event.nch<11) MC_Correction *= nTrack[event.nch-1];
			

//			if(event.nch < 10 && fnk<5 && cont != data)
//			{
//                MC_Correction*=nTrackKaonCorr[event.nch][(int)fnk];
//			}

			//cont = (cont == data?data:nFake);
			//if(event.nch<10) MC_Correction*=nTrackCorr[event.nch-1];
			/*if(cont != data)
			{
				if(SumPsig<5 && cont != charm_lep && cont != charm_s_lep)
				{
					MC_Correction*=nTrackCorr[h2d[0]->GetXaxis()->FindBin(SumPsig)-1][h2d[0]->GetYaxis()->FindBin(btag.m_bc)-1];
					//if(nTrackCorr[h2d[0]->GetXaxis()->FindBin(SumPsig)-1][h2d[0]->GetXaxis()->FindBin(SumPsig)-1] > 10)
					//cout<<h2d[0]->GetXaxis()->FindBin(SumPsig)-1<<" "<<event.nch<< " "<<nTrackCorr[h2d[0]->GetXaxis()->FindBin(SumPsig)-1][h2d[0]->GetXaxis()->FindBin(SumPsig)-1]<<endl;
				}	
			}
*/
			for(auto it = vParticles.begin(); it!=vParticles.end();++it) {
				for(auto jt = it+1; jt!=vParticles.end();++jt) {
					if(it->pid_true == jt->pid_true && it->P == jt->P)
					{
						if(it->flag&2==false) it->flag+=2;
						if(jt->flag&2==false) jt->flag+=2;
						cout<<"equal "<< it->flag<<endl;
					}
				}
			}

			TLorentzVector lvGamma;
			itp1 = gam_p1->begin();
			itp2 = gam_p2->begin();
			itp3 = gam_p3->begin();
			int nGamma = 0;
			for(; itp1!=gam_p1->end();++itp1, ++itp2, ++itp3)
			{
				//if(cont!=data && gRand.Uniform(1.)<0.02) continue;				
				TLorentzVector tmp  = TLorentzVector(*itp1, *itp2, *itp3,
						sqrt(pow(*itp1,2)+ pow(*itp2,2)+ pow(*itp3,2)) );
		                if(tmp.E()<0.15) continue;
				lvGamma += tmp;
				nGamma++;
			}
            event.nn = nGamma;
			
			lvPmisswTag -= lvTracks;
			TLorentzVector lvPmissWoGamma = lvPmisswTag;
			lvPmisswTag -= lvGamma;
			/*if(cont!=data){
				lvPmisswTag.SetE(lvPmisswTag.E()+0.01);
				lvPmissWoGamma.SetE(lvPmissWoGamma.E()+0.01);
			}
*/



		//SumPsig -= lvBreco.Rho();
t_track.costh = lvTracks.CosTheta();

		TLorentzVector lvVis = lvTracks + lvGamma-lvSignal;
		TLorentzVector lvVis2 = lvTracks + lvGamma;
        event.q2 = (lvBeam-lvBreco-lvVis).M2();
        sort(vParticles.begin(),vParticles.end(),[](const MyParticle &a, const MyParticle &b){return a.P.Rho() < b.P.Rho();});
        /*if(event.nch>2)
        {	
        	if(vParticles[0].flag ==0)
        		lvPmisswTag+=vParticles[0].P;
        	else
        		lvPmisswTag+=vParticles[1].P;
        }*/
        
/*
		lvTracks.Boost(-lvBreco.BoostVector());
		lvGamma.Boost(-lvBreco.BoostVector());
		lvVis.Boost(-lvBreco.BoostVector());
		lvPmisswTag.Boost(-lvBreco.BoostVector());
		lvPmissWoGamma.Boost(-lvBreco.BoostVector());
*/		
		lvTracks.Boost(-lvBeam.BoostVector());
		lvGamma.Boost(-lvBeam.BoostVector());
		lvVis.Boost(-lvBeam.BoostVector());
		lvPmisswTag.Boost(-lvBeam.BoostVector());
		lvPmissWoGamma.Boost(-lvBeam.BoostVector());
		
//		lvBreco.Boost(-lvBreco.BoostVector());
        gmiss.e = lvPmisswTag.E();
        cont = GetCont(it_f);



        fThLepMiss = cos(lvPmisswTag.Angle(lvSignal.Vect()));
        fThLepX = cos(lvVis.Angle(lvSignal.Vect()));
        fThXMiss = cos(lvPmisswTag.Angle(lvVis.Vect()));
        //if(fThXMiss>0) continue;

		t_PmisswTag.e = lvPmisswTag.E();
		t_PmisswTag.p3 = lvPmisswTag.Pz();
		t_PmisswTag.pt = lvPmisswTag.Perp();	
		t_PmisswTag.m2 = lvPmisswTag.M2();
		//if(cont != data) t_PmisswTag.m2 -= 0.04;
		//t_PmisswTag.m2 = lvBeam.M()/2. - ((lvTracks + lvGamma).E() + (lvTracks + lvGamma).Rho());
		//t_PmisswTag.m2 = sqr(5.278) + (lvTracks + lvGamma).M2() - 2*(lvTracks + lvGamma).Vect()*lvBreco.Vect() - 2*5.278*(lvTracks + lvGamma).E();
		//t_PmisswTag.m2 = SumPsig;
		//cout<<t_PmisswTag.m2<<endl;
		t_PmisswTag.p = lvPmisswTag.Rho();
		//if(t_PmisswTag.pt<0.35) continue;
		t_PmisswTag.costh = lvPmisswTag.CosTheta();
		gmiss.m2 = lvPmisswTag.M2();
			gmiss.ps = lvPmisswTag.Rho();
			gmiss.es = lvPmisswTag.E();
			//if(gmiss.es<1.75) continue;
		//if(gmiss.m2<-4 || gmiss.m2>4) continue;

		t_PmissWoGamma.e = lvPmissWoGamma.E();
		t_PmissWoGamma.p3 = lvPmissWoGamma.Pz();
		t_PmissWoGamma.pt = lvPmissWoGamma.Perp();
		t_PmissWoGamma.m2 = lvPmissWoGamma.M2();
		//if(cont != data) t_PmissWoGamma.m2 -= 0.04;
		//t_PmissWoGamma.m2 = lvBeam.M()/2. - ((lvTracks).E() + (lvTracks).Rho());
		//t_PmissWoGamma.m2 = sqr(5.278) + (lvTracks).M2() - 2*(lvTracks).Vect()*lvBreco.Vect() - 2*5.278*(lvTracks).E();
		t_PmissWoGamma.p = lvPmissWoGamma.Rho();
		t_PmissWoGamma.costh = lvPmissWoGamma.CosTheta();
		
			t_Vis.e = lvVis.E();
			t_Vis.p3 = lvVis.Pz();
			t_Vis.pt = lvVis.Perp();
			t_Vis.m2 = mysqrt(lvVis.M2());
			t_Vis.p = lvVis.Rho();
			t_Vis.costh = lvVis.CosTheta();
			mvis = t_track.e;
			gx.m2 = lvVis.M2();
			gx.ps = lvVis.Rho();
			gx.es = lvVis.E();
			//cont = GetCont(it_f);
            //if(t_Vis.p>2) continue;
           // if(t_Vis.m2>2.8) continue;
           lvBreco.Boost(-lvBeam.BoostVector());	
			t_btag.e = lvBreco.E();
			t_btag.p3 = lvBreco.Pz();
			t_btag.pt = lvBreco.Perp();
			t_btag.m2 = mysqrt(lvBreco.M2());
			t_btag.p = lvBreco.Rho();
			
			lvVis2.Boost(-lvBeam.BoostVector());
			t_btag.costh = cos(lvBreco.Vect().Angle(lvVis2.Vect()));
			if(nTracks > 0){
				t_track.e = lvTracks.E();
				t_track.p3 = lvTracks.Pz();
				t_track.pt = lvTracks.Perp();
				t_track.m2 = mysqrt(lvTracks.M2());
				t_track.p = lvTracks.Rho();
				//t_track.costh = lvTracks.CosTheta();
			}else{
				t_track.e = t_track.p3 = t_track.pt = t_track.m2 = t_track.p = t_track.costh=-3;
				//continue;
			}
			if(nGamma > 0){
				t_gamma.e = lvGamma.E();
				t_gamma.p3 = lvGamma.Pz();
				t_gamma.pt = lvGamma.Perp();
				t_gamma.m2 = mysqrt(lvGamma.M2());
				t_gamma.p = lvGamma.Rho();
				t_gamma.costh = lvGamma.CosTheta();
			}else {
				t_gamma.e = t_gamma.p3 = t_gamma.pt = t_gamma.m2 = t_gamma.p = t_gamma.costh =-3;
				//continue;
			}
			//cont = GetCont(it_f);



			btag.contNB = log(btag.contNB);
			btag.NB = log(btag.NB);
			fpi0 = (float) event.npi0;
			fnch = (float) event.nch;
			fnn = (float) event.nn;
			fq = (float) abs(event.q);
			fNFS = (float) btag.NFS;
			fb_mode = 2*(fabs(btag.pcode_b) - 511) + fabs(getnum(btag.b_mode));
			/*if(cont!=data){
				if(fb_mode>20 && fb_mode<38)
					MC_Correction*=tagCorr1[(int)fb_mode-21];
				else if(fb_mode>0 && fb_mode <16)
					MC_Correction *= tagCorr0[(int)fb_mode-1];
			}*/

			gmissm2Trans = gmiss.m2*sin(acos(gmiss.th));

//            {
//            int ix = h2d[0]->GetXaxis()->FindBin(t_btag.e)-1;
//            int iy = h2d[0]->GetYaxis()->FindBin(t_btag.p)-1;
//            if(ix<30 && iy<30)
//                MC_Correction*=mytagcorr[ix][iy];
//            }



			fl_lep1 = (float)lep1.fl_lep;
			fl_lep2 = (float)lep2.fl_lep;
			fl_mother1 = (float)lep1.fl_mother;
			fl_mother2 = (float)lep2.fl_mother;
			fl_mother12 = fl_mother1+fl_mother2;
			btag.p1 = sqrt(btag.p1*btag.p1 + btag.p2*btag.p2);

			fexp = (int)run.exp;
			frunno = (float)run.num;

			fmm = gmiss.es- gmiss.ps;

			fEgamma = gx.es - x.es;
			fPgamma = 0;



			if(lep1.fl_lep%10 == 0) lep1.muid = 2;
			else if(lep1.fl_lep%10 == 1) lep1.elid = 2;
			if(lep2.fl_lep%10 == 0) lep2.muid = 2;
			else if(lep2.fl_lep%10 == 1) lep2.elid = 2;

			lep1.atckpi = log(lep1.atckpi);
			lep2.atckpi = log(lep2.atckpi);



			gx.m2 = mysqrt(gx.m2);
			x.m2 = mysqrt(x.m2);

            MyParticle p[2];
        /*    fMD0Cand = RecoD2KPi(vParticles);
	    float fMD0Cand2 = RecoD2KPiPiPi(vParticles);
           // if(fabs(fMD0Cand-1.86)>0.1) continue;
            float fMDpCand = RecoDp2KPiPi(vParticles);
	    float fMDpCand2 = RecoDp2PiPiPi(vParticles);
	    int iMassWindow = 0;
	    int nTrk = 0;
	    if(fabs(fMDpCand-1.869)<0.1) {iMassWindow++; nTrk =3; }
	   // if(fabs(fMDpCand2-1.869)<0.1) {iMassWindow++; nTrk =5; }
	    if(fabs(fMD0Cand-1.86)<0.1) {iMassWindow++; nTrk =4; }
	    if(fabs(fMD0Cand2-1.86)<0.1) {iMassWindow++; nTrk =4; }

            if(iMassWindow!=1) continue;
	  */  //if(nTrk != event.nch) continue;

            if(p[0].P.Rho()>1.e-3)
            {
                hists[cont][nVariables - 8]->Fill(p[0].dz,MC_Correction);
                hists[cont][nVariables - 8]->Fill(p[1].dz,MC_Correction);
                hists_uw[cont][nVariables - 8]->Fill(p[0].dz,1);
                hists_uw[cont][nVariables - 8]->Fill(p[1].dz,1);
                hists[cont][nVariables - 6]->Fill(p[0].dr,MC_Correction);
                hists[cont][nVariables - 6]->Fill(p[1].dr,MC_Correction);
                hists_uw[cont][nVariables - 6]->Fill(p[0].dr,1);
                hists_uw[cont][nVariables - 6]->Fill(p[1].dr,1);
                hists[cont][nVariables - 4]->Fill(p[0].P.CosTheta(),MC_Correction);
                hists[cont][nVariables - 4]->Fill(p[1].P.CosTheta(),MC_Correction);
                hists_uw[cont][nVariables - 4]->Fill(p[0].P.CosTheta(),1);
                hists_uw[cont][nVariables - 4]->Fill(p[1].P.CosTheta(),1);
                hists[cont][nVariables - 2]->Fill(p[0].P.Rho(),MC_Correction);
                hists[cont][nVariables - 2]->Fill(p[1].P.Rho(),MC_Correction);
                hists_uw[cont][nVariables - 2]->Fill(p[0].P.Rho(),1);
                hists_uw[cont][nVariables - 2]->Fill(p[1].P.Rho(),1);
            }
            fMDpCand = RecoPhi2KK(vParticles,p);
            if(p[0].P.Rho()>1e-3)
            {
                hists[cont][nVariables - 7]->Fill(p[0].dz,MC_Correction);
                hists[cont][nVariables - 7]->Fill(p[1].dz,MC_Correction);
                hists_uw[cont][nVariables - 7]->Fill(p[0].dz,1);
                hists_uw[cont][nVariables - 7]->Fill(p[1].dz,1);
                hists[cont][nVariables - 5]->Fill(p[0].dr,MC_Correction);
                hists[cont][nVariables - 5]->Fill(p[1].dr,MC_Correction);
                hists_uw[cont][nVariables - 5]->Fill(p[0].dr,1);
                hists_uw[cont][nVariables - 5]->Fill(p[1].dr,1);
                hists[cont][nVariables - 3]->Fill(p[0].P.CosTheta(),MC_Correction);
                hists[cont][nVariables - 3]->Fill(p[1].P.CosTheta(),MC_Correction);
                hists_uw[cont][nVariables - 3]->Fill(p[0].P.CosTheta(),1);
                hists_uw[cont][nVariables - 3]->Fill(p[1].P.CosTheta(),1);
                hists[cont][nVariables - 1]->Fill(p[0].P.Rho(),MC_Correction);
                hists[cont][nVariables - 1]->Fill(p[1].P.Rho(),MC_Correction);
                hists_uw[cont][nVariables - 1]->Fill(p[0].P.Rho(),1);
                hists_uw[cont][nVariables - 1]->Fill(p[1].P.Rho(),1);
            }
            fMOmegaCand = RecoK0S2PiPi(vParticles);


            for(auto it = vParticles.begin(); it!=vParticles.end(); ++it)
            {
                if(it->pid != 2) continue;
                for(auto jt = it+1; jt!=vParticles.end(); ++jt)
                {
                    if(jt->pid != 2) continue;
                    if(it->charge*jt->charge > 0) continue;
                    hists[cont][ihistKK]->Fill((it->P + jt->P).M(),MC_Correction*it_f->weight);
                    hists_uw[cont][ihistKK]->Fill((it->P + jt->P).M());
                }

            }
            for(auto it = vParticles.begin(); it!=vParticles.end(); ++it)
            {
                if(it->pid != 2) continue;
                for(auto jt = vParticles.begin(); jt!=vParticles.end(); ++jt)
                {
                    if(jt->pid != 3) continue;
                    if(it->charge*jt->charge > 0) continue;
                    hists[cont][ihistKK+1]->Fill((it->P + jt->P).M(),MC_Correction*it_f->weight);
                    hists_uw[cont][ihistKK+1]->Fill((it->P + jt->P).M());
                }

            }
            for(auto it = vParticles.begin(); it!=vParticles.end(); ++it)
            {
                if(it->pid != 3) continue;
                for(auto jt = it+1; jt!=vParticles.end(); ++jt)
                {
                    if(jt->pid != 3) continue;
                    if(it->charge*jt->charge > 0) continue;
                    float M = (it->P + jt->P).M();
                    if(M<0.5 || M > 1.) continue;
                    hists[cont][ihistKK+2]->Fill(M,1);//MC_Correction*it_f->weight);
                    hists_uw[cont][ihistKK+2]->Fill((it->P + jt->P).M());
                }

            }

			if(cont == data) MC_Correction=1;

			std::vector<float>::iterator it_p = vtrk_p->begin();
			std::vector<float>::iterator it_e = vtrk_e->begin();
			std::vector<float>::iterator it_th = vtrk_th->begin();
			float ptot = 0;
		    for(auto iPar= vParticles.begin(); iPar!=vParticles.end();++iPar)
		    {
			if(iPar->pid>5) continue;
			//if(iPar->P.Perp()>0.5) continue;
			ptot+= iPar->P.Perp();
		    }
		    //if(ptot>2.5) continue;
		   // for(auto it= vParticles.begin(); it!=vParticles.end();++it)
		    {
			//if(it->pid>5) continue;
			//if(iPar->P.Perp()>0.5) continue;
			//for(auto jPar = iPar; jPar!=vParticles.end();++jPar)
			{
				//if(iPar!=jPar) continue;
				//(jPar->pid>5) continue;
				//if(jPar->P.Perp()>0.5) continue;
				/*float th = iPar->charge*jPar->charge*cos(iPar->P.Angle(jPar->P.Vect()));
				float p1 = iPar->charge*jPar->charge*iPar->P.Perp();
				float p2 = jPar->P.Perp();
				float p = iPar->charge*jPar->charge*fabs(iPar->P.Perp()-jPar->P.Perp());
				float Dr2 = sqr(iPar->dr) + sqr(jPar->dr) - 2.*fabs(iPar->dr*jPar->dr)*cos(iPar->P.Vect().DeltaPhi(jPar->P.Vect()));
			//	cout<<jPar->P.Perp()<<" "<<sqrt(Dr2)<<endl;
				float dy  = sqrt(Dr2 + sqr(iPar->dz - jPar->dz));
				*/
				//TLorentzVector pi = it->P;
				//pi.Boost(-lvBeam.BoostVector());
				if(cont == data)
					h2d[1]->Fill(t_btag.p,t_PmissWoGamma.m2);
				else{
					h2d[0]->Fill(t_btag.p,t_PmissWoGamma.m2,MC_Correction*it_f->weight);
					if(cont == charm_lep || cont == charm_s_lep )
						h2d[3]->Fill(t_btag.p,t_PmissWoGamma.m2,MC_Correction*it_f->weight);
					else 
						h2d[4]->Fill(t_btag.p,t_PmissWoGamma.m2,MC_Correction*it_f->weight);
						
				}
			}
		    }
			//if(nSlow>1) continue;


			
			/*if(cont!=data)
			{
				if(event.nch ==2) MC_Correction *= Corr2Track[hists[0][0]->FindBin(vParticles[0].P.Rho())-1];
				else if(event.nch ==3) MC_Correction *= Corr3Track[hists[0][0]->FindBin(vParticles[0].P.Rho())-1];
				else if(event.nch ==4) MC_Correction *= Corr4Track[hists[0][0]->FindBin(vParticles[0].P.Rho())-1];
				else if(event.nch ==5) MC_Correction *= Corr5Track[hists[0][0]->FindBin(vParticles[0].P.Rho())-1];
			}*/
			for(int j = 0; j<nVariables; ++j)
			{

				if((j==nVarPs1Elec-1) && lep1.fl_lep%10 != 0) continue;
				if((j==nVarPs1Elec) && lep1.fl_lep%10 != 1) continue;
				if(j<(nVariables -nExtraVariables))
				{
					if(j<12)
					{
						if(event.nch == 2 && j == 0)
						{
							int iTrack = 0;
							for( auto it = vParticles.begin(); it!=vParticles.end(); ++it)
							{
								if(it->flag ==1 ) continue;
								if(it_f->Type != 'd') cont = iTrack;
								else cont = data;
//								hists[cont][j]->Fill(lep1.q*it->charge*it->P.Rho(),MC_Correction*it_f->weight );
//								hists_uw[cont][j]->Fill(lep1.q*it->charge*it->P.Rho(),1. );
                                				hists[cont][j+iTrack]->Fill(it->P.Rho(),MC_Correction*it_f->weight );
								hists_uw[cont][j+iTrack]->Fill(it->P.Rho(),1. );
								iTrack++;
							}
						}
						else if(event.nch == 3 && j == 1)
						{
							int iTrack = 0;
							for( auto it = vParticles.begin(); it!=vParticles.end(); ++it)
							{
								if(it->flag ==1 ) continue;
								if(it_f->Type != 'd') cont = iTrack;
								else cont = data;
//								hists[cont][j]->Fill(lep1.q*it->charge*it->P.Rho(),MC_Correction*it_f->weight );
//								hists_uw[cont][j]->Fill(lep1.q*it->charge*it->P.Rho(),1. );
                                				hists[cont][j+iTrack]->Fill(it->P.Rho(),MC_Correction*it_f->weight );
								hists_uw[cont][j+iTrack]->Fill(it->P.Rho(),1. );
								iTrack++;
							}
						}
						else if(event.nch == 4 && j == 3)
						{
							int iTrack = 0;
							for( auto it = vParticles.begin(); it!=vParticles.end(); ++it)
							{
								if(it->flag ==1 ) continue;
								if(it_f->Type != 'd') cont = iTrack;
								else cont = data;
//								hists[cont][j]->Fill(lep1.q*it->charge*it->P.Rho(),MC_Correction*it_f->weight );
//								hists_uw[cont][j]->Fill(lep1.q*it->charge*it->P.Rho(),1. );
                                				hists[cont][j+iTrack]->Fill(it->P.Rho(),MC_Correction*it_f->weight );
								hists_uw[cont][j+iTrack]->Fill(it->P.Rho(),1. );
								iTrack++;
							}
						}
						else if(event.nch == 5 && j == 6)
						{
							int iTrack = 0;
							for( auto it = vParticles.begin(); it!=vParticles.end(); ++it)
							{
								if(it->flag ==1 ) continue;
								if(it_f->Type != 'd') cont = iTrack;
								else cont = data;
//								hists[cont][j]->Fill(lep1.q*it->charge*it->P.Rho(),MC_Correction*it_f->weight );
//								hists_uw[cont][j]->Fill(lep1.q*it->charge*it->P.Rho(),1. );
                                				hists[cont][j+iTrack]->Fill(it->P.Rho(),MC_Correction*it_f->weight );
								hists_uw[cont][j+iTrack]->Fill(it->P.Rho(),1. );
								iTrack++;
							}
						}





					}else if(j<24){
						cont = GetCont(it_f);
						if(event.nch == j-12+1)
						{
						hists[cont][j]->Fill(SumPsig,MC_Correction*it_f->weight );
						hists_uw[cont][j]->Fill(SumPsig,1. );
						}
					}else{

						cont = GetCont(it_f);
					hists[cont][j]->Fill(vVar[j].val(),MC_Correction*it_f->weight );
					hists_uw[cont][j]->Fill(vVar[j].val(),1. );
					}
				}
				else if(j == nVariables -nExtraVariables)
				{
					for(std::vector<int>::iterator it = vtrk_sep->begin(); it!=vtrk_sep->end(); ++it)
					{
						hists[cont][j]->Fill(*it%10,MC_Correction*it_f->weight );
						hists_uw[cont][j]->Fill(*it%10,1. );
					}
				}
				else if(j == nVariables -  30)
				{
					for( auto it = vParticles.begin(); it!=vParticles.end(); ++it)
					{
						if(it->pid != 0 ) continue;
						if(it_f->Type != 'd') cont = it->pid_true;
						else cont = data;

						hists[cont][j]->Fill(it->P.Rho(),MC_Correction*it_f->weight );
						hists_uw[cont][j]->Fill(it->P.Rho(),1. );
					}
				}
				else if(j == nVariables -  29)
				{
					std::vector<int>::iterator it_sep = vtrk_sep->begin();
					for(std::vector<float>::iterator it = vtrk_th->begin(); it!=vtrk_th->end(); ++it,++it_sep)
					{
						if(*it_sep%10 !=0) continue;
						if(cont!=data) cont = *it_sep/10;
						hists[cont][j]->Fill(cos(*it/180*3.1415),MC_Correction*it_f->weight );
						hists_uw[cont][j]->Fill(cos(*it/180*3.1415),1. );
					}
				}
				else if(j == nVariables -  28)
				{
					std::vector<int>::iterator it_sep = vtrk_sep->begin();
					for(std::vector<float>::iterator it = vtrk_dr->begin(); it!=vtrk_dr->end(); ++it,++it_sep)
					{
						if(*it_sep%10 !=0) continue;
						if(cont!=data) cont = *it_sep/10;
						hists[cont][j]->Fill(*it,MC_Correction*it_f->weight );
						hists_uw[cont][j]->Fill(*it,1. );
					}
				}
				else if(j == nVariables -  27)
				{
					std::vector<int>::iterator it_sep = vtrk_sep->begin();
					for(std::vector<float>::iterator it = vtrk_dz->begin(); it!=vtrk_dz->end(); ++it,++it_sep)
					{
						if(*it_sep%10 !=0) continue;
						if(cont!=data) cont = *it_sep/10;
						hists[cont][j]->Fill(*it,MC_Correction*it_f->weight );
						hists_uw[cont][j]->Fill(*it,1. );
					}
				}
				else if(j == nVariables -  26)
				{
					for( auto it = vParticles.begin(); it!=vParticles.end(); ++it)
					{
						if(it->pid != 1 ) continue;
						if(it_f->Type != 'd') cont = it->pid_true;
						else cont = data;

						hists[cont][j]->Fill(it->P.Rho(),MC_Correction*it_f->weight );
						hists_uw[cont][j]->Fill(it->P.Rho(),1. );
					}
				}
				else if(j == nVariables -  25)
				{
					std::vector<int>::iterator it_sep = vtrk_sep->begin();
					for(std::vector<float>::iterator it = vtrk_th->begin(); it!=vtrk_th->end(); ++it,++it_sep)
					{
						if(*it_sep%10 !=1) continue;
						if(cont!=data) cont = *it_sep/10;
						hists[cont][j]->Fill(cos(*it/180*3.1415),MC_Correction*it_f->weight );
						hists_uw[cont][j]->Fill(cos(*it/180*3.1415),1. );
					}
				}
				else if(j == nVariables -  24)
				{
					std::vector<int>::iterator it_sep = vtrk_sep->begin();
					for(std::vector<float>::iterator it = vtrk_dr->begin(); it!=vtrk_dr->end(); ++it,++it_sep)
					{
						if(*it_sep%10 !=1) continue;
						if(cont!=data) cont = *it_sep/10;
						hists[cont][j]->Fill(*it,MC_Correction*it_f->weight );
						hists_uw[cont][j]->Fill(*it,1. );
					}
				}
				else if(j == nVariables -  23)
				{
					std::vector<int>::iterator it_sep = vtrk_sep->begin();
					for(std::vector<float>::iterator it = vtrk_dz->begin(); it!=vtrk_dz->end(); ++it,++it_sep)
					{
						if(*it_sep%10 !=1) continue;
						if(cont!=data) cont = *it_sep/10;
						hists[cont][j]->Fill(*it,MC_Correction*it_f->weight );
						hists_uw[cont][j]->Fill(*it,1. );
					}
				}
				else if(j == nVariables -  22)
				{
					int n = 0;
					std::vector<int>::iterator it_sep = vtrk_sep->begin();
					for(std::vector<float>::iterator it = vtrk_p->begin(); it!=vtrk_p->end(); ++it,++it_sep)
					{
						if(cont!=data) cont = *it_sep/10;
						if(*it_sep%10 != 2 ) continue;
						hists[cont][j]->Fill(*it,MC_Correction*it_f->weight );
						hists_uw[cont][j]->Fill(*it,1. );
						n++;
					}
					cont = GetCont(it_f);
					hists[cont][j+4]->Fill(n,MC_Correction*it_f->weight );
					hists_uw[cont][j+4]->Fill(n,1. );
				}
				else if(j == nVariables -  21)
				{
					std::vector<int>::iterator it_sep = vtrk_sep->begin();
					for(std::vector<float>::iterator it = vtrk_th->begin(); it!=vtrk_th->end(); ++it,++it_sep)
					{
						if(*it_sep%10 !=2 ) continue;
						if(cont!=data) cont = *it_sep/10;
						hists[cont][j]->Fill(cos(*it/180*3.1415),MC_Correction*it_f->weight );
						hists_uw[cont][j]->Fill(cos(*it/180*3.1415),1. );
					}
				}
				else if(j == nVariables -  20)
				{
					std::vector<int>::iterator it_sep = vtrk_sep->begin();
					for(std::vector<float>::iterator it = vtrk_dr->begin(); it!=vtrk_dr->end(); ++it,++it_sep)
					{
						if(*it_sep%10 !=2) continue;
						if(cont!=data) cont = *it_sep/10;
						hists[cont][j]->Fill(*it,MC_Correction*it_f->weight );
						hists_uw[cont][j]->Fill(*it,1. );
					}
				}
				else if(j == nVariables -  19)
				{
					std::vector<int>::iterator it_sep = vtrk_sep->begin();
					for(std::vector<float>::iterator it = vtrk_dz->begin(); it!=vtrk_dz->end(); ++it,++it_sep)
					{
						if(*it_sep%10 !=2) continue;
						if(cont!=data) cont = *it_sep/10;
						hists[cont][j]->Fill(*it,MC_Correction*it_f->weight );
						hists_uw[cont][j]->Fill(*it,1. );
					}
				}
				else if(j == nVariables -  18)
				{
					int n=0;

					cont = GetCont(it_f);

					for(std::vector<int>::iterator it = vtrk_sep->begin(); it!=vtrk_sep->end(); ++it)
					{
						if(*it%10 !=2) continue;
						n++;
					}
			//		hists[cont][j]->Fill(n,MC_Correction*it_f->weight );
			//		hists_uw[cont][j]->Fill(n,1. );
				}
				else if(j == nVariables -  17)
				{
					std::vector<int>::iterator it_sep = vtrk_sep->begin();
					for(std::vector<float>::iterator it = vtrk_p->begin(); it!=vtrk_p->end(); ++it,++it_sep)
					{
						if(cont!=data) cont = *it_sep/10;
						if(*it_sep%10 != 3) continue;
						hists[cont][j]->Fill(*it,MC_Correction*it_f->weight );
						hists_uw[cont][j]->Fill(*it,1. );
					}
				}
				else if(j == nVariables -  16)
				{
					std::vector<int>::iterator it_sep = vtrk_sep->begin();
					for(std::vector<float>::iterator it = vtrk_th->begin(); it!=vtrk_th->end(); ++it,++it_sep)
					{
						if(*it_sep%10 !=3) continue;
						if(cont!=data) cont = *it_sep/10;
						hists[cont][j]->Fill(cos(*it/180*3.1415),MC_Correction*it_f->weight );
						hists_uw[cont][j]->Fill(cos(*it/180*3.1415),1. );
					}
				}
				else if(j == nVariables -  15)
				{
					for( auto it = vParticles.begin(); it!=vParticles.end(); ++it)
					{
						for( auto jt = it+1; jt!=vParticles.end(); ++jt)
						{
							//if(it->flag ==1 ) continue;
							if(it_f->Type != 'd') cont = it->pid_true;
							else cont = data;
							TLorentzVector pi = it->P, pj = jt->P;
							pi.Boost(-lvBeam.BoostVector()); pj.Boost(-lvBeam.BoostVector());
							hists[cont][j]->Fill(pi.Vect()*(pj.Vect()),MC_Correction*it_f->weight );
							hists_uw[cont][j]->Fill(pi.Vect()*(pj.Vect()),1. );
						}
					}
				}
				else if(j == nVariables -  14)
				{
					for( auto it = vParticles.begin(); it!=vParticles.end(); ++it)
					{
						for( auto jt = it+1; jt!=vParticles.end(); ++jt)
						{
							//if(it->flag ==1 ) continue;
							if(it_f->Type != 'd') cont = it->pid_true;
							else cont = data;
							TLorentzVector pi = it->P, pj = jt->P;
							pi.Boost(-lvBeam.BoostVector()); pj.Boost(-lvBeam.BoostVector());
							hists[cont][j]->Fill(cos(pi.Angle(pj.Vect())),MC_Correction*it_f->weight );
							hists_uw[cont][j]->Fill(cos(pi.Angle(pj.Vect())),1. );
						}
					}
				}
				else if(j == nVariables -  13)
				{
					//if(lep1.q>0){
					for( auto it = vParticles.begin(); it!=vParticles.end(); ++it)
					{
						//if(it->flag ==1 ) continue;
						if(it_f->Type != 'd') cont = it->pid_true;
						else cont = data;
						TLorentzVector pi = it->P;
							pi.Boost(-lvBeam.BoostVector()); 
						hists[cont][j]->Fill(/*it->charge*/pi.Rho(),MC_Correction*it_f->weight );
						hists_uw[cont][j]->Fill(/*it->charge*/pi.Rho(),1. );
					}
					//}
				}
				else if(j == nVariables -  12)
				{
					
					for( auto it = vParticles.begin(); it!=vParticles.end(); ++it)
					{
						for( auto jt = it+1; jt!=vParticles.end(); ++jt)
						{
							//if(it->flag ==1 ) continue;
							if(it_f->Type != 'd') cont = it->pid_true;
							else cont = data;
							TLorentzVector pi = it->P, pj = jt->P;
							pi.Boost(-lvBeam.BoostVector()); pj.Boost(-lvBeam.BoostVector());
							hists[cont][j]->Fill(pi.Rho()* pj.Rho(),MC_Correction*it_f->weight );
							hists_uw[cont][j]->Fill(pi.Rho()* pj.Rho(),1. );
						}
					}
					
				}
				else if(j == nVariables -  11)
				{
					for( auto it = vParticles.begin(); it!=vParticles.end(); ++it)
					{
					//	if(it->charge < 0 ) continue;
					//	if(it->P.Rho()>1) continue;
						if(cont!=data) cont = fabs( it->gen_ist-1 );
						if( cont > 50 ) cont = 12;
						hists[cont][j]->Fill(it->P.CosTheta(),MC_Correction*it_f->weight );
						hists_uw[cont][j]->Fill(it->P.CosTheta(),1. );
					}
				}
				else if(j == nVariables -  10)
				{
					for( auto it = vParticles.begin(); it!=vParticles.end(); ++it)
					{
					//	if(it->charge > 0 ) continue;
						if(cont!=data) cont = fabs( it->gen_ist-1 );
						if( cont > 50 ) cont = 12;
						hists[cont][j]->Fill(it->dr,MC_Correction*it_f->weight );
						hists_uw[cont][j]->Fill(it->dr,1. );
					}
				}
				else if(j == nVariables -  9)
				{
					for( auto it = vParticles.begin(); it!=vParticles.end(); ++it)
					{
					//	if(it->charge > 0 ) continue;
						if(cont!=data) cont = fabs( it->gen_ist-1 );
						if( cont > 50 ) cont = 12;
						hists[cont][j]->Fill(it->dz,MC_Correction*it_f->weight );
						hists_uw[cont][j]->Fill(it->dz,1. );
					}
				}
				else if(j == nVariables -  8)
				{
					for( auto it = vParticles.begin(); it!=vParticles.end(); ++it)
					{
						if(it->P.Perp()>0.200) continue;
						for( auto jt = it+1;jt!=vParticles.end(); ++jt)
						{
							if(jt->P.Rho()>0.200) continue;
							if(it->charge*jt->charge < 0)continue;
							if(cont!=data)
							{
								if(it->pid_true == jt->pid_true && (it->P - jt->P).Rho()<0.1)
									cont = 0;
								else cont = 1;
							}

							hists[cont][j]->Fill(cos(it->P.Angle(jt->P.Vect())),MC_Correction*it_f->weight );
							hists_uw[cont][j]->Fill(cos(it->P.Angle(jt->P.Vect())),1. );
						}
					}
				}
				else if(j == nVariables -  7)
				{
					for( auto it = vParticles.begin(); it!=vParticles.end(); ++it)
					{
						if(it->P.Rho()>0.200) continue;
						for( auto jt = it+1;jt!=vParticles.end(); ++jt)
						{
							if(jt->P.Rho()>0.200) continue;
							if(it->charge*jt->charge > 0)continue;
							if(cont!=data)
							{
								if(it->pid_true == jt->pid_true && (it->P - jt->P).Rho()<0.1)
									cont = 0;
								else cont = 1;
							}

							hists[cont][j]->Fill(cos(it->P.Angle(jt->P.Vect())),MC_Correction*it_f->weight );
							hists_uw[cont][j]->Fill(cos(it->P.Angle(jt->P.Vect())),1. );
						}
					}
				}
				else if(j == nVariables -  6)
				{
					for( auto it = vParticles.begin(); it!=vParticles.end(); ++it)
					{
						if(it->P.Rho()>0.300 || it->P.Perp()<0.200) continue;
						for( auto jt = it+1;jt!=vParticles.end(); ++jt)
						{
							if(it->P.Rho()>0.300 || it->P.Perp()<0.200) continue;
							if(it->charge*jt->charge < 0)continue;
							if(cont!=data)
							{
								if(it->pid_true == jt->pid_true && (it->P - jt->P).Rho()<0.1)
									cont = 0;
								else cont = 1;
							}

							hists[cont][j]->Fill(cos(it->P.Angle(jt->P.Vect())),MC_Correction*it_f->weight );
							hists_uw[cont][j]->Fill(cos(it->P.Angle(jt->P.Vect())),1. );
						}
					}
				}
				else if(j == nVariables -  5)
				{
					for( auto it = vParticles.begin(); it!=vParticles.end(); ++it)
					{
						if(it->P.Rho()>0.300 || it->P.Perp()<0.200) continue;
						for( auto jt = it+1;jt!=vParticles.end(); ++jt)
						{
							if(it->P.Rho()>0.300 || it->P.Perp()<0.200) continue;
							if(it->charge*jt->charge > 0)continue;
							if(cont!=data)
							{
								if(it->pid_true == jt->pid_true && (it->P - jt->P).Rho()<0.1)
									cont = 0;
								else cont = 1;
							}

							hists[cont][j]->Fill(cos(it->P.Angle(jt->P.Vect())),MC_Correction*it_f->weight );
							hists_uw[cont][j]->Fill(cos(it->P.Angle(jt->P.Vect())),1. );
						}
					}
				}
				else if(j == nVariables -  4)
				{
					for( auto it = vParticles.begin(); it!=vParticles.end(); ++it)
					{
						if(it->P.Rho()>0.400 || it->P.Perp()<0.300) continue;
						for( auto jt = it+1;jt!=vParticles.end(); ++jt)
						{
							if(it->P.Rho()>0.400 || it->P.Perp()<0.300) continue;
							if(it->charge*jt->charge < 0)continue;
							if(cont!=data)
							{
								if(it->pid_true == jt->pid_true && (it->P - jt->P).Rho()<0.1)
									cont = 0;
								else cont = 1;
							}

							hists[cont][j]->Fill(cos(it->P.Angle(jt->P.Vect())),MC_Correction*it_f->weight );
							hists_uw[cont][j]->Fill(cos(it->P.Angle(jt->P.Vect())),1. );
						}
					}
				}
				else if(j == nVariables -  3)
				{
					for( auto it = vParticles.begin(); it!=vParticles.end(); ++it)
					{
						if(it->P.Rho()>0.400 || it->P.Perp()<0.300) continue;
						for( auto jt = it+1;jt!=vParticles.end(); ++jt)
						{
							if(it->P.Rho()>0.400 || it->P.Perp()<0.300) continue;
							if(it->charge*jt->charge > 0)continue;
							if(cont!=data)
							{
								if(it->pid_true == jt->pid_true && (it->P - jt->P).Rho()<0.1)
									cont = 0;
								else cont = 1;
							}

							hists[cont][j]->Fill(cos(it->P.Angle(jt->P.Vect())),MC_Correction*it_f->weight );
							hists_uw[cont][j]->Fill(cos(it->P.Angle(jt->P.Vect())),1. );
						}
					}
				}
				else if(j == nVariables -  2)
				{
					for( auto it = vParticles.begin(); it!=vParticles.end(); ++it)
					{
						if( it->P.Perp()<0.400) continue;
						for( auto jt = it+1;jt!=vParticles.end(); ++jt)
						{
							if( it->P.Perp()<0.400) continue;
							if(it->charge*jt->charge < 0)continue;
							if(cont!=data)
							{
								if(it->pid_true == jt->pid_true && (it->P - jt->P).Rho()<0.1)
									cont = 0;
								else cont = 1;
							}

							hists[cont][j]->Fill(cos(it->P.Angle(jt->P.Vect())),MC_Correction*it_f->weight );
							hists_uw[cont][j]->Fill(cos(it->P.Angle(jt->P.Vect())),1. );
						}
					}
				}
				else if(j == nVariables -  1)
				{
					for( auto it = vParticles.begin(); it!=vParticles.end(); ++it)
					{
						if( it->P.Perp()<0.400) continue;
						for( auto jt = it+1;jt!=vParticles.end(); ++jt)
						{
							if( it->P.Perp()<0.400) continue;
							if(it->charge*jt->charge > 0)continue;
							if(cont!=data)
							{
								if(it->pid_true == jt->pid_true && (it->P - jt->P).Rho()<0.1)
									cont = 0;
								else cont = 1;
							}

							hists[cont][j]->Fill(cos(it->P.Angle(jt->P.Vect())),MC_Correction*it_f->weight );
							hists_uw[cont][j]->Fill(cos(it->P.Angle(jt->P.Vect())),1. );
						}
					}
				}
			}

			cont = GetCont(it_f);
			nEventsSplit[iLep][fabs(btag.pcode_b)==511?0:1][cont]+=MC_Correction*it_f->weight ;
			nEvents[cont]+=MC_Correction*it_f->weight ;
			nEventsPerFile[0]+=MC_Correction*it_f->weight;
			if(lep1.ps<1.2 && xlep.m2 <2.3)
				nEventsSigReg[cont]+=MC_Correction*it_f->weight;
			if(1)
			{
				if(BDecay[0] != '0')
				{
					BDecays->Fill(NoCharge(BDecay).c_str(),1);
					if(DDecay[0] != '0')
						DDecays->Fill(NoSubDecay(NoCharge(DDecay)).c_str(),1);
				}
			}

			if(_newTree) t_new->Fill();

		}

		cout<<"\rDone" <<endl;

		float nRho = 0;
		for(int i = 0; i< nContributions; i++)
			nRho += hists[i][ihistKK+2]->Integral();
		cout<<"//////////////////////////////////\n"<<nRho<<endl;

		ps->Print();
		char sNam[100];
		sprintf(sNam,"ReadPerf_%s.root",it_f->Name.c_str());
		ps->SaveAs();
		delete ps;
		t->DropBranchFromCache("*");
		t->SetCacheSize(0);
		t->Delete();
		f.Close();
		cout<<"Events Per File: \n";
		for(int i=0;i<nContributions;++i)
			cout<<nEventsPerFile[i]<<"\t";
		cout<<endl;
		if(_newTree)
		{

            f_new->Write();
            t_new->Delete();
            f_new->Close();
		}


	}


	for(int i =1; i<=h2d[0]->GetNbinsX();i++)
	{
		cout<<"{ ";
		for(int j =1; j<=h2d[0]->GetNbinsY();j++)
		{
			float dat = h2d[1]->GetBinContent(i,j);
			float mc = h2d[0]->GetBinContent(i,j);
			float mcc = h2d[3]->GetBinContent(i,j); //const mc
			float mcf = h2d[4]->GetBinContent(i,j); //scale mc

			float tmp=0;
			if( dat )
			{
				tmp = (dat-mc)/sqrt(dat);
			}
			h2d[2]->SetBinContent(i,j,tmp);
			tmp = 1;
			if(dat&&mcf)
				tmp = (dat-mcc)/mcf;

			cout<<tmp;
			if(j+1<=h2d[0]->GetNbinsY())
				cout<<", ";
		} 
		cout<<" },"<<endl;
	}


	TCanvas *c2 = new TCanvas ("v","v",1000,1000);
	c2 -> Divide(2,2);
	c2->cd(1);
	h2d[0]->Draw("colz");
	c2->cd(2);
	h2d[2]->Draw("colz");
	c2->cd(4);
	h2d[1]->Draw("colz");
	c2->cd(0);
	c2->Print("plot2D.pdf");

	cout<<endl<<"Events: "<<endl;
	float sumB = 0;
	for(int i=0;i<nContributions;++i)
	{
		cout<<nEvents[i]<<"\t";
		if(i!=signal_decay)
			sumB += nEvents[i];
	}
	cout<<endl
		<<"S/B: "<<nEvents[0]/sumB
		<<endl
		<<"S/sqrt(S+B): " <<nEvents[0]/sqrt(nEvents[0]+sumB)
		<<endl;

	cout<<"in signal region:"<<endl;
	cout<<endl<<"Events: "<<endl;
	sumB = 0;
	for(int i=0;i<nContributions;++i)
	{
		cout<<nEventsSigReg[i]<<"\t";
		if(i!=signal_decay)
			sumB += nEventsSigReg[i];
	}
	cout<<endl
		<<"S/B: "<<nEventsSigReg[0]/sumB
		<<endl
		<<"S/sqrt(S+B): " <<nEventsSigReg[0]/sqrt(nEventsSigReg[0]+sumB)
		<<endl;

	cout<<"Events:"<<endl;
	for(int il = 0; il<2; il++){
	cout<<"e mu = "<<il<<endl;
	for(int iB = 0; iB<2; iB++){
	cout<<"B B0"<<iB<<endl;
	for(int i=0;i<nContributions;++i){
		cout<<nEventsSplit[il][iB][i]<<"\t";
	}
	cout<<endl;
	}
	cout<<endl;
	}
	cout<<endl;

	cout<<"s,b in brecon"<<endl;
	float sig, bg;
	cout<<"B+"<<endl;
	for(int i = 1; i<=hists[0][nVarReconMode]->GetNbinsX();i++)
	{
		sig = hists[0][nVarReconMode]->GetBinContent(i);
		bg = 0;
		for(int j = 1; j<nContributions;j++)
			bg += hists[j][nVarReconMode]->GetBinContent(i);

		cout<<i<<"\t"<<sig<<"\t"<<bg<<"\t"<<sig/bg*100<<endl;
	}

	cout<<"B0"<<endl;
	nVarReconMode++;
	for(int i = 1; i<=hists[0][nVarReconMode]->GetNbinsX();i++)
	{
		sig = hists[0][nVarReconMode]->GetBinContent(i);
		bg = 0;
		for(int j = 1; j<nContributions;j++)
			bg += hists[j][nVarReconMode]->GetBinContent(i);

		cout<<i<<"\t"<<sig<<"\t"<<bg<<"\t"<<sig/bg*100<<endl;
	}


	PrintMaxima(BDecays,40);
	PrintMaxima(DDecays,40);

	for(int i=0;i<12;i++)
	{
		PrintMaxima(BDecaysBins[i],20);
		PrintMaxima(DDecaysBins[i],20);
	}



	cout<<double_events<<endl;

	TCanvas *cv1 = new TCanvas("c","c",600,600);
	hWeight2 -> Divide(hWeight1);
	hWeight2 -> Draw();
	cv1 -> Print("weight.pdf");
	cv1->Delete();

	cout<< "ND1: "<<countD1<<endl;




	cout<< "N Events::"<<endl;
	for(int i=0;i<2;++i) {
		for(int j=0;j<4;++j) {
			cout<<i<<" "<<j<<" ";
			for(int k=0;k<nContributions;++k) {
				cout<<nEventsBCuts[i][j][k]<<" ";
			}
			cout<<endl;
		}
	}






}



void DrawHist()
{
	char sCanvName[200];
	char sPadName[200];
	char sHistTitle[200];

	sprintf(sCanvName,"test");

	TCanvas *cv1 = new TCanvas(sCanvName,sCanvName,1050,1485);
	TVirtualPad *pad;

	sprintf(sHistTitle," ");
	cv1 -> Print("plot.pdf[");

	TH1F *signal_tmp[nVariables];
	TH1F *sumMC[nVariables];
	TH1F *MCerr[nVariables];
	TH1F *sumMC_uw[nVariables];
	TLine *line = new TLine();
	line->SetLineWidth(1);
	line->SetLineColor(12);
	line->SetLineStyle(0);
	float MC2DataScale = 1;




	//fit mbc

	/*

	const int nFitBins = 12;
	const int nFitComp = 3;
	TH1F ***hc = new TH1F**[nFitBins];
	TH1F ***hc_uw = new TH1F**[nFitBins];
	for(int i = 0; i<nFitBins; i++)
	{
		hc[i] = new TH1F*[nFitComp];
		hc_uw[i] = new TH1F*[nFitComp];
	}

	double y[nFitBins][nFitComp], dy[nFitBins][nFitComp] ,x[nFitBins], dx[nFitBins], yield[nFitBins], fitScaling[nFitBins][nFitComp], fitEvents[nFitBins][nFitComp];
	for(int i = 0; i<nFitBins; i++)
	{
		for(int j = 0;j<nFitComp;j++)
		{
			sprintf(sHistTitle,"hist_%i_%i_%s",i,j,vVar[i].Name.c_str());
			hc[i][j] = newHist(sHistTitle,vVar[i]);
			sprintf(sHistTitle,"hist_uw_%i_%i_%s",i,j,vVar[i].Name.c_str());
			hc_uw[i][j] = newHist(sHistTitle,vVar[i]);
			hc[i][j]->Add(hists[j][i]);
			hc_uw[i][j]->Add(hists_uw[j][i]);
		}
		x[i]=i;
		dx[i]=0.5;
		SimpleFit(hists[data][i],hc[i],nFitComp,y[i],dy[i],yield[i],hc_uw[i]);
	}

	cout<<"x dx well/poor d scale well d scale poor d"<<endl;
	for(int i = 0; i<nFitBins; i++)
	{


		for(int j = 0;j<nFitComp;j++)
		{
			fitEvents[i][j] = y[i][j]*yield[i];
			fitScaling[i][j] = y[i][j]*yield[i]/hists[j][i]->Integral();
			hists[j][i]->Scale(fitScaling[i][j]);

			cout<< i << " " << x[i]<<" "<< dx[i]<< " "<< yield[i] << " "<<j << " "<< y[i][j]<< " "<< dy[i][j] <<" "<<fitEvents[i][j] << " " << fitScaling[i][j]<< endl;

		}
	}*/
	/*for(int i = 0; i<nFitBins; i++)
	{
		for(int j = 0;j<nFitComp;j++)
		{
			hists[j][i]->Scale(fitScaling[i][j]);
		}
	}*/
	/*TGraphErrors *gr_mbc = new TGraphErrors(nFitBins, x,y,dx,dy);
	TGraphErrors *gr_mbc2 = new TGraphErrors(nFitBins, x,y2,dx,dy2);
	TGraphErrors *gr_mbcr = new TGraphErrors(nFitBins, x,yr,dx,dyr);
	cv1->Divide(2,2);
	cv1->cd(1);
	gr_mbc->SetMarkerStyle(20);
	gr_mbc->Draw("AP");
	cv1->cd(2);
	gr_mbc2->SetMarkerStyle(20);
	gr_mbc2->SetMarkerColor(kRed);
	gr_mbc2->SetFillColor(kRed);
	gr_mbc2->Draw("AP");
	cv1->cd(3);
	gr_mbcr->SetMarkerStyle(20);
	gr_mbcr->SetMarkerColor(kBlue);
	gr_mbcr->SetFillColor(kBlue);
	gr_mbcr->Draw("AP");
	cv1->cd(0);
	cv1->Print("graphMbc.pdf");
*/
	for(int j = 0; j<nVariables; ++j)
	{
		if( j%6 == 0)
		{
			cv1 -> Clear();
			cv1 -> Divide(2,3);
		}



		for(int i=0; i<nContributions; i++)
		{
			if(i == data) continue;
			//SmearHist(hists[i][j],0.9);
			//SmearHist(hists_uw[i][j],0.9);
		}

		sprintf(sHistTitle,"sumMC_%i",j);
		sumMC[j] = newHist(sHistTitle, vVar[j]);
		sprintf(sHistTitle,"MCerr_%i",j);
		MCerr[j] = newHist(sHistTitle, vVar[j]);
		sprintf(sHistTitle,"sumMC_uw_%i",j);
		sumMC_uw[j] = newHist(sHistTitle, vVar[j]);

		pad = cv1 -> cd( j%6 +1 );
		pad->cd();
		TPad *pad1;
        if(_drawData){
            pad1 = new TPad("p1","p1",0.,0.3,1,1);
            pad1->SetTopMargin(0.1);
            pad1->SetBottomMargin(0);
            pad1->Draw();
            pad1->cd();
        }

		if(vVar[j].fl_logy) {_drawData?pad1->SetLogy():pad->SetLogy();}

		StackMC[j]->Draw("HIST");
		StackMC[j]->GetXaxis()->SetTitle(vVar[j].xTitle.c_str());
		StackMC[j]->GetYaxis()->SetTitle("Events");
		signal_tmp[j] = (TH1F*)hists[signal_decay][j]->Clone();
		signal_tmp[j]->Draw("same HIST");
		signal_tmp[j]->SetLineStyle(2);
		signal_tmp[j]->SetFillStyle(0);

		for(int i = 0; i<nContributions;++i)
		{
			if(i==data) continue;
			//hists[i][j] -> Scale(fScale[i]);
			hists[i][j] -> Scale(vCont[i].Weight);
			sumMC[j]->Add(hists[i][j]);
			sumMC_uw[j]->Add(hists_uw[i][j]);
//			hists[i][j] -> Scale(1.0/hists[i][j]->Integral());
		}
		if(vVar[j].Name == string("missM"))
		{
			FILE *fp = fopen("mmiss.txt","w");

			for(int i = 1; i<= sumMC[j]->GetNbinsX(); i++)
                fprintf(fp,"%i\t%f\t%f\n",i,sumMC[j]->GetBinContent(i), hists[data][j]->GetBinContent(i));

            fclose(fp);
		}
		if(Scale_MC)
		{
			if(j==0) MC2DataScale = hists[data][j]->Integral()/sumMC[j]->Integral();
			cout<<MC2DataScale<<" "<<hists[data][j]->Integral()<<" "<<sumMC[j]->Integral()<<endl;
			sumMC[j] -> Scale(MC2DataScale);
			for(int i = 0; i<nContributions;++i)
			{
				if(i==data) continue;
				hists[i][j] -> Scale(MC2DataScale);
			}
		}
		if(signal_tmp[j]->Integral()>0)
			signal_tmp[j]->Scale(0.2*sumMC[j]->Integral()/signal_tmp[j]->Integral() );
		StackMC[j]->Modified();
		StackMC[j]->SetTitle(vVar[j].LeafName.c_str());
		if(j==0)
		StackMC[j]->SetTitle(sCut);
		StackMC[j]->GetXaxis()->SetTitle(vVar[j].xTitle.c_str());
		StackMC[j]->GetYaxis()->SetTitle("Events");


	if(_drawData)
		hists[data][j]->Draw("E1 X0 same");
		hists[data][j]->GetXaxis()->SetTitle(vVar[j].xTitle.c_str());
		hists[data][j]->GetYaxis()->SetTitle("Events");
		//hists[data][j]->GetYaxis()->SetTitleOffset(1.25);

		if(hists[data][j]->GetMaximum()+sqrt(hists[data][j]->GetMaximum()) > StackMC[j]->GetMaximum())
		{
			StackMC[j]->SetMaximum(hists[data][j]->GetMaximum()+sqrt(hists[data][j]->GetMaximum()));
		}
		StackMC[j]->SetMinimum(0.1);

		cv1->Update();

		pad->cd();
		TPad *pad2;
        if(_drawData)
        {
            pad2 = new TPad("p2","p2",0.,0.,1,0.3);
            pad2->SetTopMargin(0);
            pad2->SetBottomMargin(0.45);
            pad2->Draw();
            pad2->cd();


            sumMC[j]->Sumw2();
            //sumMC[j]->Divide(hists[data][j]);
            float val,err,dat;
	    cout<<vVar[j].xTitle<<endl;
            for(int i =1; i<=vVar[j].nBins;i++)
            {

                val = sumMC[j]->GetBinContent(i);
                dat = hists[data][j]->GetBinContent(i);
                err = sumMC_uw[j]->GetBinContent(i);
		if(err != 0)
			err = val/err * sqrt(err); // MC_error = w_i * sqrt(n_i)
                if((val != 0 && dat != 0) &&  _drawData)
                {
                    sumMC[j]->SetBinContent(i, (dat -val)/sqrt(dat/* + err*err*/) );
                    sumMC[j]->SetBinError(i, 1.);
                    //sumMC[j]->SetBinContent(i, dat/val );
                    //sumMC[j]->SetBinError(i, sqrt(dat)/val);
                }
                else
                {
                    sumMC[j]->SetBinContent(i, -50 );
                    sumMC[j]->SetBinError(i, 1 );
                }

                //if(val>0) MCerr[j]->SetBinError(i,dat*err/val/val);
                if(dat>0) MCerr[j]->SetBinError(i,err/sqrt(dat));
                else MCerr[j]->SetBinError(i,0);
                MCerr[j]->SetBinContent(i,0);

		cout<<i<<" "<<dat/val<<"\n";

            }
            //gPad->SetLogy();
            sumMC[j]->SetMarkerStyle(20);
            //sumMC[j]->SetAxisRange(0.667,1.5,"Y");
            sumMC[j]->SetAxisRange(-15.9,15.9,"Y");
            //sumMC[j]->SetAxisRange(-5.9,5.9,"Y");
            sumMC[j]->GetYaxis()->SetNdivisions(505);
            //sumMC[j]->GetYaxis()->SetTitle("#frac{Data}{MC}");
            sumMC[j]->GetYaxis()->SetTitle("#frac{Data - MC}{#Delta}");
            sumMC[j]->GetXaxis()->SetTitle(vVar[j].xTitle.c_str());
            sumMC[j]->GetXaxis()->SetNdivisions(505);
            sumMC[j]->GetYaxis()->SetTitleOffset(pad2->GetAbsHNDC()/pad->GetAbsHNDC());
            sumMC[j]->SetLabelSize(0.05/pad2->GetAbsHNDC()*pad->GetAbsHNDC(),"x");
            sumMC[j]->SetLabelSize(0.05/pad2->GetAbsHNDC()*pad->GetAbsHNDC(),"y");
            sumMC[j]->SetTitleSize(0.06/pad2->GetAbsHNDC()*pad->GetAbsHNDC(),"x");
            sumMC[j]->SetTitleSize(0.6*0.06/pad2->GetAbsHNDC()*pad->GetAbsHNDC(),"y");
            sumMC[j]->Draw("ep X0");
            MCerr[j]->SetFillStyle(3354);
            MCerr[j]->SetFillColor(12);
            MCerr[j]->Draw("e2 same");
            line->DrawLine(vVar[j].xMin,0.,vVar[j].xMax,0.);
            sumMC[j]->Draw("ep X0 same");


            //sumMC[j]->Draw("ep");

            StackMC[j]->GetYaxis()->SetTitleOffset(1.05*pad1->GetAbsHNDC()/pad->GetAbsHNDC());
            StackMC[j]->GetXaxis()->SetLabelSize(0.05/pad1->GetAbsHNDC()*pad->GetAbsHNDC());
            StackMC[j]->GetYaxis()->SetLabelSize(0.05/pad1->GetAbsHNDC()*pad->GetAbsHNDC());
            StackMC[j]->GetXaxis()->SetTitleSize(0.06/pad1->GetAbsHNDC()*pad->GetAbsHNDC());
            StackMC[j]->GetYaxis()->SetTitleSize(0.06/pad1->GetAbsHNDC()*pad->GetAbsHNDC());
            StackMC[j]->GetYaxis()->Draw();
        }
		pad->cd();


		sprintf(sPadName,"sp/%s.pdf",vVar[j].Name.c_str());
		//pad -> Print(sPadName);

		if((j + 1)%6 == 0 || j+1 ==nVariables)
		{
			cv1 -> cd(0);
			cv1 -> Print("plot.pdf");
		}

	}
	cv1 -> cd(0);
	cv1 -> Print("plot.pdf]");
	//delete cv1;

	/*
	cv1 -> Print("plot_ratio.pdf[");

	TH1F *htmp;

	for(int j = 0; j<nVariables; ++j)
	{
		if( j%6 == 0)
		{
			cv1 -> Clear();
			cv1 -> Divide(2,3);
		}

		pad = cv1 -> cd( j%6 +1 );

		htmp = (TH1F*)hists[0][j]->Clone();


		for(int i = 0; i<nContributions;++i)
		{
			if(i==data) continue;
			if(i==signal_decay) continue;
			htmp->Add(hists[i][j]);
		}

		htmp->Divide(hists[data][j]);

		htmp->Draw("E1 X0");
		htmp->GetXaxis()->SetTitle(vVar[j].xTitle.c_str());
		htmp->GetYaxis()->SetTitle("N_{MC}/N_{Data}");
		//hists[data][j]->GetYaxis()->SetTitleOffset(1.25);


		cv1->Update();


		if((j + 1)%6 == 0 || j+1 ==nVariables)
		{
			cv1 -> cd(0);
			cv1 -> Print("plot_ratio.pdf");
		}
		htmp->Delete();
	}
	cv1 -> cd(0);
	cv1 -> Print("plot_ratio.pdf]");
	*/
	cv1 -> Print("plot_NoStack.pdf[");

  gROOT->SetStyle("belleStyle");
  gROOT->ForceStyle();

	float Ymax;
	for(int j = 0; j<nVariables; ++j)
	{

		if( j%6 == 0)
		{
			cv1 -> Clear();
			cv1 -> Divide(2,3);
		}

		cv1 -> cd( j%6 +1 );

		Ymax =0;

			sumMC[j]->Reset("ICESM");
			sumMC_uw[j]->Reset("ICESM");
		for(int i = 0; i<nContributions;++i)
		{
			if(i==data) continue;
			hists[i][j] -> SetFillStyle(0);
			sumMC[j]->Add(hists[i][j]);
			sumMC[j]->SetFillStyle(0);
			sumMC[j]->SetLineColor(kBlue+2);
			sumMC_uw[j]->Add(hists_uw[i][j]);
		}
		StackMC[j]->Draw("HIST nostack");
		StackMC[j]->SetTitle(vVar[j].LeafName.c_str());
		StackMC[j]->GetXaxis()->SetTitle(vVar[j].xTitle.c_str());
		StackMC[j]->GetYaxis()->SetTitle("Events");
		StackMC[j]->GetYaxis()->SetTitleOffset(1.05);
		StackMC[j]->GetXaxis()->SetLabelSize(0.05);
		StackMC[j]->GetYaxis()->SetLabelSize(0.05);
		StackMC[j]->GetXaxis()->SetTitleSize(0.06);
		StackMC[j]->GetYaxis()->SetTitleSize(0.06);
		StackMC[j]->GetYaxis()->Draw();
		sumMC[j]->Draw("HIST same");
		if(_drawData)
			hists[data][j]->Draw("E1 X0 same");
		hists[data][j]->GetXaxis()->SetTitle(vVar[j].xTitle.c_str());
		hists[data][j]->GetYaxis()->SetTitle("Events");
		Ymax = hists[data][j]->GetMaximum();
		Ymax +=sqrt(Ymax);
		if(Ymax < sumMC[j]->GetMaximum()) Ymax = sumMC[j]->GetMaximum();
		StackMC[j]->SetMaximum(Ymax*1.05);
		cv1->Update();

		if((j + 1)%6 == 0 || j+1 ==nVariables)
		{
			cv1 -> cd(0);
			cv1 -> Print("plot_NoStack.pdf");
		}
	}
	cv1 -> cd(0);
	cv1 -> Print("plot_NoStack.pdf]");
	cout<<"\n\nMC 2 Data scale: "<<MC2DataScale<<endl;

}

void CreateHist()
{

	gStyle -> SetHistLineWidth( 1 );
	gStyle -> SetHistFillColor( kWhite );
//	gStyle -> SetLegendFillColor ( kWhite ); // only in root 5.30
	gStyle -> SetCanvasColor(kWhite);
	gStyle -> SetFrameFillColor(kWhite);
	gStyle -> SetTitleFillColor( kWhite );
	gStyle -> SetTitleBorderSize( 0 );
	gStyle -> SetOptStat( 0 );
	gStyle->SetHatchesSpacing(2);
	//gStyle -> SetTitleFontSize( 2 );
	//gStyle -> SetTitleAlign ( 0 );

 TStyle *belleStyle= new TStyle("belleStyle","Phill's  un-official plots style");

  // use helvetica-bold-r-normal, precision 2 (rotatable)
  Int_t belleFont = 62;
  // line thickness
  Double_t belleWidth = 2;

  // use plain black on white colors
  belleStyle->SetFrameBorderMode(0);
  belleStyle->SetCanvasBorderMode(0);
  belleStyle->SetPadBorderMode(0);
  belleStyle->SetPadColor(0);
  belleStyle->SetCanvasColor(0);
  belleStyle->SetStatColor(0);
  belleStyle->SetPalette(1);
  belleStyle->SetFillColor(0);

  // set the paper & margin sizes
  belleStyle->SetPaperSize(20,26);
  belleStyle->SetPadTopMargin(0.05);
  belleStyle->SetPadRightMargin(0.15); // increase for colz plots!!
  belleStyle->SetPadBottomMargin(0.16);
  belleStyle->SetPadLeftMargin(0.14);

  // use large fonts
  belleStyle->SetTextFont(belleFont);
  belleStyle->SetTextSize(0.08);
  belleStyle->SetLabelFont(belleFont,"x");
  belleStyle->SetLabelFont(belleFont,"y");
  belleStyle->SetLabelFont(belleFont,"z");
  belleStyle->SetLabelSize(0.05,"x");
  belleStyle->SetLabelSize(0.05,"y");
  belleStyle->SetLabelSize(0.05,"z");
  belleStyle->SetTitleFont(belleFont);
  belleStyle->SetTitleSize(0.06,"x");
  belleStyle->SetTitleSize(0.06,"y");
  belleStyle->SetTitleSize(0.06,"z");
	belleStyle->SetHatchesSpacing(2);

  // use bold lines and markers
  belleStyle->SetLineWidth(belleWidth);
  belleStyle->SetFrameLineWidth(belleWidth);
  belleStyle->SetHistLineWidth(belleWidth);
  belleStyle->SetFuncWidth(belleWidth);
  belleStyle->SetGridWidth(belleWidth);
  belleStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
 // belleStyle->SetMarkerStyle(8);
 // belleStyle->SetMarkerSize(1.0);

  // label offsets
  belleStyle->SetLabelOffset(0.015);

  // by default, do not display histogram decorations:
  belleStyle->SetOptStat(0);
  //belleStyle->SetOptStat(1110);  // show only nent, mean, rms
  belleStyle->SetOptTitle(0);
  belleStyle->SetOptFit(0);
//  belleStyle->SetOptFit(1011); // show probability, parameters and errors

  // look of the statistics box:
  belleStyle->SetStatBorderSize(1);
  belleStyle->SetStatFont(belleFont);
  belleStyle->SetStatFontSize(0.05);
  belleStyle->SetStatX(0.9);
  belleStyle->SetStatY(0.9);
  belleStyle->SetStatW(0.25);
  belleStyle->SetStatH(0.15);

  // put tick marks on top and RHS of plots
  belleStyle->SetPadTickX(1);
  belleStyle->SetPadTickY(1);

  // histogram divisions: only 5 in x to avoid label overlaps
  belleStyle->SetNdivisions(505,"x");
  belleStyle->SetNdivisions(510,"y");

  gROOT->SetStyle("belleStyle");
  gROOT->ForceStyle();

  TPaveText *belleName = new TPaveText(0.65,0.8,0.9,0.9,"BRNDC");
  belleName->SetFillColor(0);
  belleName->SetTextAlign(12);
  belleName->SetBorderSize(0);
  belleName->AddText("Belle");

  TText *belleLabel = new TText();
  belleLabel->SetTextFont(belleFont);
  belleLabel->SetTextColor(1);
  belleLabel->SetTextSize(0.04);
  belleLabel->SetTextAlign(12);

  TLatex *belleLatex = new TLatex();
  belleLatex->SetTextFont(belleFont);
  belleLatex->SetTextColor(1);
  belleLatex->SetTextSize(0.04);
  belleLatex->SetTextAlign(12);
  TGaxis::SetMaxDigits(3);

//Add Belle Preliminary, 710 fb-1 to plots
belleLatex->SetTextSize(0.07);
//belleLatex->DrawLatex((h_data[0]->GetXaxis()->GetXmax()-h_data[0]->GetXaxis()->GetXmin())*.1+h_data[0]->GetXaxis()->GetXmin(), 6./7.*h_data[0]->GetMaximum(), "#splitline{#splitline{Belle}{Preliminary}}{#scale[0.7]{710 fb^{-1}}}");


	DefVariables();
	DefHist();
	FillHist();
	DrawHist();

}
