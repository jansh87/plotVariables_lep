
#define _BrLog 1
#define nbins_const 80
#define _pilep 0
#define _drawData 0
#define Scale_MC 0
#define _newTree 0
#define smallFile


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
float flclass,fdclass,fl_lep1, fl_lep2, fl_mother1, fl_mother2, fl_mother12;
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
float fMKK,fMKPi,fMPiPi;
float gmissm2Trans;
float fThLepMiss,fThLepX, fThXMiss;
float px1, px2;
float SumPsig;
float fUmiss;
float EmissLepSys;
int nVarReconMode;
int hadevtmidx;
int nVarPs1Elec;
int nVarPs2Elec;
int iVarMx,iVarMb,igPmiss,iEgammaIdx;
int ihistKK;
float pmisspl;
float abspmisspl;
float fHadEvtE,fHadEvtP,fHadEvtM;
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
	
	tmp = CVar("MM20","MM2/GeV2",&gmiss.m2, 60,-5, 15);
	vVar.push_back(tmp);
	tmp = CVar("MM21","MM2/GeV2",&gmiss.m2, 60,-5, 15);
	vVar.push_back(tmp);
	tmp = CVar("MM22","MM2/GeV2",&gmiss.m2, 60,-5, 15);
	vVar.push_back(tmp);
	tmp = CVar("MM23","MM2/GeV2",&gmiss.m2, 60,-5, 15);
	vVar.push_back(tmp);
	tmp = CVar("MM24","MM2/GeV2",&gmiss.m2, 60,-5, 15);
	vVar.push_back(tmp);
	tmp = CVar("MM25","MM2/GeV2",&gmiss.m2, 60,-5, 15);
	vVar.push_back(tmp);

	tmp = CVar("MM26","MM2/GeV2",&gmiss.m2, 60,-5, 15);
	vVar.push_back(tmp);
	tmp = CVar("MM27","MM2/GeV2",&gmiss.m2, 60,-5, 15);
	vVar.push_back(tmp);
	tmp = CVar("MM28","MM2/GeV2",&gmiss.m2, 60,-5, 15);
	vVar.push_back(tmp);
	tmp = CVar("MM29","MM2/GeV2",&gmiss.m2, 60,-5, 15);
	vVar.push_back(tmp);
	tmp = CVar("MM210","MM2/GeV2",&gmiss.m2, 60,-5, 15);
	vVar.push_back(tmp);
	tmp = CVar("MM211","MM2/GeV2",&gmiss.m2, 60,-5, 15);
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

	tmp = CVar("mbc","M_{bc} / GeV",&btag.m_bc, 5.26, 5.29);
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


	tmp = CVar("ps1e","electron p*_{#tau} / GeV",&lep1.ps, 0.2, 2.5);
	vVar.push_back(tmp);
	nVarPs1Elec = vVar.size();

	tmp = CVar("ps1mu","muon p*_{#tau} / GeV",&lep1.ps, 0.5, 2.5);
	vVar.push_back(tmp);

	tmp = CVar("labcosthetael","electron cos#theta_{#tau}^{lab}",&lep1.th, -1, 1);
	vVar.push_back(tmp);
	tmp = CVar("labcosthetamu","muon cos#theta_{#tau}^{lab}",&lep1.th, -1, 1);
	vVar.push_back(tmp);

	tmp = CVar("pmisspl","p_{miss}+p_{l}",&pmisspl, -0, 2.5);
	vVar.push_back(tmp);
	tmp = CVar("abspmisspl","|p_{miss}|+|p_{l}|",&abspmisspl, 0, 5.5);
	vVar.push_back(tmp);
	
	tmp = CVar("Emisslepsys","E^{#tau}_{miss} / GeV",&EmissLepSys, 0, 5.5);
	vVar.push_back(tmp);
	
	tmp = CVar("mm2","M_{miss}^{2} / GeV^{2}",&miss.m2, -5, 20);
	vVar.push_back(tmp);

	tmp = CVar("emiss","E*_{miss} / GeV",&miss.es, 0, 2.5);
	vVar.push_back(tmp);

	tmp = CVar("pmiss","p*_{miss} / GeV",&miss.ps, 0, 2.5);
	vVar.push_back(tmp);

	tmp = CVar("mm2g","M_{miss}^{2} / GeV^{2}",&gmiss.m2, 71,-5, 15);
	vVar.push_back(tmp);

	tmp = CVar("mm2gnPeak","M_{miss}^{2} / GeV^{2}",&gmiss.m2, 30,-5, 15);
	vVar.push_back(tmp);

	tmp = CVar("mm2gtrans","E*_{miss} - P*_{miss}/ GeV",&fUmiss,60, -2, 4);
	vVar.push_back(tmp);
	
	tmp = CVar("mm2truth","truth M_{miss}^{2} / GeV^{2}",&truth.miss_m2, 0.0001, 15);
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
	tmp = CVar("HadEvtE","HadE/E",&fHadEvtE, 0, 1);
	vVar.push_back(tmp);
	
	tmp = CVar("HadEvtP","HadP/P",&fHadEvtP, 0, 1.3);
	vVar.push_back(tmp);
	tmp = CVar("HadEvtM","HadM/M",&fHadEvtM, 0, 1);
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
	iEgammaIdx = vVar.size()-1;
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

	tmp = CVar("gx.es","E*_{X} / GeV",&gx.es,5.1,5.5);
	vVar.push_back(tmp);
	hadevtmidx = vVar.size()-1;
	tmp = CVar("gx.m2","M_{X} / GeV",&gx.m2,0,6);
	vVar.push_back(tmp);

	tmp = CVar("gx.m2.zoom","M_{X} / GeV",&gx.m2,40,1.5,2.1);
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
	tmp = CVar("MPiPi","M_{#pi#pi} / GeV",&fMPiPi,200,0.25,1.25);
	vVar.push_back(tmp);
    tmp = CVar("MOmega","M_{#Omega} / GeV",&fMOmegaCand,0.49,0.51);
	vVar.push_back(tmp);

	tmp = CVar("xlep.th_xl",&xlep.th_xl,0.,180);
	//vVar.push_back(tmp);

	tmp = CVar("xlep.ths_xl",&xlep.ths_xl,-1.,1);
//	vVar.push_back(tmp);

	tmp = CVar("xlep.th_xl",&xlep.th_xl,-1.,1);
//	vVar.push_back(tmp);



	tmp = CVar("lclass",&flclass,12,0,12);
	vVar.push_back(tmp);
	tmp = CVar("dclass",&fdclass,12,0,12);
	vVar.push_back(tmp);
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
	//{signal_decay,charm_lep, charm_s_lep, charm_nonres_lep, charm_ss_lep, double_charm, charm_other, FakeLep, ulnu, rare,  otherB, badB, continuum, data, nContributions};
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
	
	float nStreams = 5. ;
	float nContiStreams = 5. ;
	//string path = string("/scratch/bctaunu/kakuno/robin/new");
	string path = string("/scratch/bctaunu/kakuno/new_0815/");
	//string path = string("~");
	
	#ifdef smallFile 
		path = string("~/rootFiles/smallFiles/");
	#endif

	tmp_file = CFile(path,"DssMC.root","lep/tree",'s');
	tmp_file.weight = 0.405;
	vFile.push_back(tmp_file );

	tmp_file = CFile(path,"ulnu.root","lep/tree",'u');
	tmp_file.weight = 1./20;
	//vFile.push_back(tmp_file );
	
	tmp_file = CFile(path,"rare.root","lep/tree",'r');
	tmp_file.weight = 1./50;
	//vFile.push_back(tmp_file );
	
	tmp_file =CFile(path, "charged_s2.root", "lep/tree",'o');
	tmp_file.weight = 1./nStreams;
	vFile.push_back(tmp_file );
	tmp_file =CFile(path, "mixed_s2.root", "lep/tree",'o');
	tmp_file.weight = 1./nStreams;
	vFile.push_back(tmp_file );
	tmp_file =CFile(path, "continuum_s2.root", "lep/tree",'c');
	tmp_file.weight = 1./nContiStreams;
	//vFile.push_back(tmp_file );
	
	tmp_file =CFile(path, "charged_s3.root", "lep/tree",'o');
	tmp_file.weight = 1./nStreams;
	vFile.push_back(tmp_file );
	tmp_file =CFile(path, "mixed_s3.root", "lep/tree",'o');
	tmp_file.weight = 1./nStreams;
	vFile.push_back(tmp_file );
	tmp_file =CFile(path, "continuum_s3.root", "lep/tree",'c');
	tmp_file.weight = 1./nContiStreams;
	//vFile.push_back(tmp_file );
	
	tmp_file =CFile(path, "charged_s4.root", "lep/tree",'o');
	tmp_file.weight = 1./nStreams;
	vFile.push_back(tmp_file );
	tmp_file =CFile(path, "mixed_s4.root", "lep/tree",'o');
	tmp_file.weight = 1./nStreams;
	vFile.push_back(tmp_file );
	tmp_file =CFile(path, "continuum_s4.root", "lep/tree",'c');
	tmp_file.weight = 1./nContiStreams;
	vFile.push_back(tmp_file );
	
	tmp_file =CFile(path, "charged_s1.root", "lep/tree",'o');
	tmp_file.weight = 1./nStreams;
	vFile.push_back(tmp_file );
	tmp_file =CFile(path, "mixed_s1.root", "lep/tree",'o');
	tmp_file.weight = 1./nStreams;
	vFile.push_back(tmp_file );
	tmp_file =CFile(path, "continuum_s1.root", "lep/tree",'c');
	tmp_file.weight = 1./nContiStreams;
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
	tmp_file = CFile(path,"dilepDssMC.root","tree",'s');
	tmp_file.weight = 0.4012;
	vFile.push_back(tmp_file );

	tmp_file = CFile(path,"dilepulnu.root","tree",'u');
	tmp_file.weight = 1./20;
	vFile.push_back(tmp_file );
	
	tmp_file = CFile(path,"dileprare.root","tree",'r');
	tmp_file.weight = 1./50;
	//vFile.push_back(tmp_file );
	tmp_file =CFile(path, "dilepcharged_s0.root", "tree",'o');
	tmp_file.weight = 1./nStreams;
	vFile.push_back(tmp_file );
	tmp_file =CFile(path, "dilepmixed_s0.root", "tree",'o');
	tmp_file.weight = 1./nStreams;
	vFile.push_back(tmp_file );
	
	tmp_file =CFile(path, "dilepcontinuum_s0.root", "tree",'c');
	tmp_file.weight = 1./nContiStreams;
	vFile.push_back(tmp_file );
	*/
	
	
	
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
*/
	tmp_file =CFile(path, "data.root", "lep/tree",'d');
	//vFile.push_back(tmp_file );

	tmp_file =CFile(path, "dilepdata.root", "tree",'d');
	//vFile.push_back(tmp_file );

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
TH1F **hOverlap;


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
		StackMC[j] -> Add(hists[signal_decay][j]);
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
		

		hists[data][j] -> SetMarkerStyle(20);
	}
	
	hOverlap= new TH1F*[2];
	
	for(int i = 0; i<2; i++)
	{
		sprintf(sHistName,"hist_%i",i);
		hOverlap[i] = newHist(sHistName,CVar("ps1e","p*_{#tau} / GeV",&lep1.ps, 100,0.1, 2.5));
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
void SmearHist(TH1F* h, float width1, float width2)
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
			if(j<i)
				m[i][j] = TMath::Gaus(j,i,width1);
			else
				m[i][j] = TMath::Gaus(j,i,width2);

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
			/*if(lep1.fl_lep != 10 && lep1.fl_lep!=21)
			{cont = FakeLep;}
			else {
			if(it_f -> Type == 'c')
				cont = continuum;
			else if(abs(event.lclass)==1 )
				cont = signal_decay;
			else if(lep1.fl_mother==3)
				cont = charm_lep;
			else 	cont = charm_s_lep;

			}*/
			/*string DDecay_t(sDDecay);
			DDecay_t = NoCharge(DDecay_t);
			//cout<<BDecay_t<<endl;
			string bad_part("_11_12");
			string bad_part2("_13_14");
*/
			if(it_f -> Type == 'c')
				cont = continuum;
			else if(event.lclass==0)
				cont = badB;
			/*else if(!(DDecay_t.find(bad_part) == string::npos && DDecay_t.find(bad_part2) == string::npos) )
			{
				if(abs(event.lclass)>4 && abs(event.lclass)<8)
				cont = rare;
				else
				cont = ulnu;
			}*/

			
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
					cont = charm_ss_lep;
				else
					cont = charm_ss_lep;
			
			}
			else if( it_f -> Type == 'u' || it_f->Type == 'r' )
				cont = charm_lep;
			else if(abs(event.lclass)==4)
				cont = double_charm;
			else if(abs(event.lclass)>=8)
				cont = otherB;
			else if(abs(event.lclass)>4)
				cont = charm_other;
			else if( lep1.fl_mother != 3  )
				cont = rare;
			
			if(abs(event.lclass) == 2){
			
				if(event.dclass == 5) cont = signal_decay;
				else if(event.dclass == 6) cont = charm_lep;
				else if(event.dclass == 7) cont = charm_ss_lep;
				else if(event.dclass == 8) cont = charm_s_lep;
				else if(event.dclass == 9) cont = FakeLep;
				else if(event.dclass == 10) cont = otherB;
				else if(event.dclass == 2) cont = rare;
				else if(event.dclass == 4) cont = charm_other;
			
			} else cont = continuum;
			
			/*	
			if(it_f -> Type == 'o'){ 
			//	if(gx.m2<3.25)
					cont = signal_decay;
			//	else cont = charm_lep;
			} else if(it_f -> Type == 's'){ 
				//if(btag.MCinfo>0)
					//cont = double_charm;
				 cont = FakeLep;
			}
			else cont = charm_lep;*/
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


double hadevtmcorr[] = {0, 2.46835, 0, 6.3939, 1.56408, 1.12852, 1.94338, 0, 1.11921, 0.659311, 1.08986, 0.823203, 1.78223, 0.84319, 0.769463, 1.11906, 0.845603, 1.00735, 1.0098, 1.01441, 0.946604, 1.11069, 0.963036, 1.05889, 1.01184, 1.06138, 1.01439, 1.12209, 1.04781, 1.08079, 1.03779, 1.0215, 1.00248, 1.00602, 0.959273, 0.948362, 0.927711, 0.949564, 0.968424, 1.00326, 1.04087, 1.06397, 1.07196, 1.05354, 1.03732, 1.00863, 0.985663, 0.980595, 0.962734, 0.985286, 1.00348, 0.940521, 0.949531, 0.933668, 0.908892, 0.909383, 0.986867, 0.908655, 0.932292, 0.930977, 0.912199, 0.914302, 0.925017, 0.950972, 0.869345, 0.900125, 0.928797, 0.880132, 0.980829, 0.933984, 0.904398, 1.01765, 1.40589, 0.919362, 0.304029, 1.08425, 1.68535, 0.667243, 0.717869, 2.28896};

//double hadevtmcorr[] = {0.976451, 1.00016, 1.02787, 1.02775, 1.02923, 1.17261, 1.02595, 1.11353, 1.18384, 1.02181, 1.20224, 1.1355, 1.13015, 1.15708, 1.11735, 1.19843, 1.11803, 1.19244, 1.15083, 1.08573, 1.21765, 1.25708, 1.28014, 1.20556, 1.17012, 1.24802, 1.36183, 1.20769, 1.22923, 1.27731, 1.19655, 1.3205, 1.17774, 1.1051, 1.53465, 1.24733, 1.84056, 0.881984, 1.15178, 0.804584, 0.846207, 1.06732, 0.846522, 2.81344, 4.03633, 1.82478, 3.81649, 1.84419, 13.4293, 0.225018, 1, 1.16686, 0.64684, 1, 1, 1, 0.927069, 1, 1.09958, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.91554, 1, 1, 1, 1, 1, 1, 1, 1};


float mm2[]={
0.837799	,
1.13376	,
1.18667	,
1.14855	,
0.92933	,
1.15564	,
1.03698	,
1.01785	,
1.09103	,
1.14161	,
1.13077	,
1.18325	,
1.10911	,
1.08299	,
1.04874	,
0.94439	,
0.954535	,
0.978097	,
0.958552	,
0.951883	,
0.985025	,
0.979419	,
0.958512	,
0.986949	,
0.988573	,
0.970326	,
0.944252	,
1.03253	,
0.992645	,
1.05505	,
0.951124	,
1.05568	,
0.992202	,
1.01805	,
1.0068	,
0.970205	,
1.0223	,
0.993993	,
0.998213	,
0.968549	,
0.985375	,
1.01729	,
1.01635	,
1.03719	,
1.12079	,
0.912672	,
1.05963	,
0.903222	,
1.03563	,
1.0159	,
0.988084	,
0.962257	,
1.17206	,
1.06457	,
1.13449	,
0.87553	,
1.04245	,
1.07889	,
0.741053	,
0.980703	
};

float EgammaCorr[] = {
1,
1,
1,
1.00834,
0.995654,
1.00019,
0.990906,
0.989635,
0.998663,
1.03341,
1.02301,
1.03025,
1.02524,
1.01925,
1.0293,
1.03999,
1.02202,
1.00991,
1.03195,
1.02079,
1.03262,
1.01989,
1.00048,
0.993257,
1.01549,
1.00912,
0.98267,
1.01492,
0.983652,
1.01485,
1.00555,
0.989415,
1.04824,
0.998664,
0.981177,
1.00024,
0.972118,
1.04896,
1.00154,
0.977615,
1.10514,
1.03276,
1.0641,
1.02965,
1.05636,
1.05275,
1.02918,
1.01447,
1.12833,
1.06277,
1.00304,
1.07065,
0.983975,
1.08735,
1.08479,
1.05551,
1.18538,
0.990338,
1.08392,
1.08712,
1.43477,
1.10126,
1.17703,
0.941079,
0.849435,
1.29239,
0.9678,
2.21791,
0.646688,
0.531839,
1.14199,
1.0108,
1.40323,
1.85771,
1.18201,
1.43822,
0.264603,
1.19087,
0.724259,
1.48191
};
float mm2_corr(float x)
{
double a               = 0.965744     ;
double b               = 0.0395465    ;   
double c               = -0.00930655   ;   
double d               = 0.000774376    ;  
double e               = -2.04358e-05   ;  
double z               = 1.03433        ;  
return a + b*(x-z) + c*pow(x-z,2) + d*pow(x-z,3) + e*pow(x-z,4);


}


float mm2_corr2 (float x)
{

double a               = 0.995916    ;//     // 5.398e+08    (5.42e+10%)
double b               = 0.0151975      ;//  // 3.917e+08    (2.577e+12%)
double c               = -0.00551708  ;//    // 6.175e+07    (1.119e+12%)
double d               = 0.000579994  ;//    // 2.099e+06    (3.618e+11%)
double e               = -1.47854e-05 ;//    // 1.65e-05     (111.6%)
double h               = -0.277708    ;//    // 2.089e+10    (7.522e+12%)
double i               = -1.28267     ;//    // 1.396e+10    (1.088e+12%)
double j               = -0.42851     ;//    // 2.972e+09    (6.935e+11%)
double k               = -0.0608041   ;//    // 2.1e+08      (3.454e+11%)
double l               = -0.00322239  ;//    // 0.004697     (145.8%)
double y               = 1.88136      ;//    // 1.63e+10     (8.662e+11%)
double z               = 2.51416      ;//    // 3.549e+10    (1.411e+12%)

double res = x>0 ? a + b*(x-z) + c*pow(x-z,2) + d*pow(x-z,3) + e*pow(x-z,4): h + i*(x-y) + j*pow(x-y,2) + k*pow(x-y,3) + l*pow(x-y,4);

return res;

}

float mm2_corr3 (float x)
{

double a               = 0.955128         ;// 2.953e+09    (3.091e+11%)
double b               = 0.060381         ;// 1.422e+09    (2.354e+12%)
double c               = -0.0145376       ;// 1.702e+08    (1.171e+12%)
double d               = 0.00116043       ;// 5.756e+06    (4.96e+11%)
double e               = -2.94386e-05     ;// 5.157e-06    (17.52%)
double f               = -0.00808814      ;// 7.637e+10    (9.443e+14%)
double g               = -0.850382        ;// 3.194e+10    (3.756e+12%)
double h               = -0.177768        ;// 2.675e+09    (1.505e+12%)
double i               = -0.00992677      ;// 0.001912     (19.26%)
double y               = 1.47704          ;// 8.984e+10    (6.082e+12%)
double z               = 1.79805          ;// 4.888e+10    (2.719e+12%)


double res = x>0 ? a + b*(x-z) + c*pow(x-z,2) + d*pow(x-z,3) + e*pow(x-z,4): f + g*(x-y) + h*pow(x-y,2) + i*pow(x-y,3) ;

return res;

}

float Umiss_corr(float x)
{

double a               = 0.986312         ;// 1.546e+08    (1.568e+10%)
double b               = -0.0277884       ;// 1.559e+08    (5.609e+11%)
double c               = 0.0140896        ;// 2.638e+08    (1.872e+12%)
double d               = 0.0158385        ;// 5.936e+07    (3.748e+11%)
double e               = 0.00267049       ;// 0.002146     (80.35%)
double f               = -1.89504         ;// 4.179e+10    (2.205e+12%)
double g               = -3.50531         ;// 3.127e+10    (8.92e+11%)
double h               = -1.31092         ;// 5.554e+09    (4.236e+11%)
double i               = -0.155187        ;// 0.06382      (41.12%)
double y               = 1.73211          ;// 1.193e+10    (6.887e+11%)
double z               = 2.76284          ;// 5.557e+09    (2.011e+11%)

double res = x>0 ? a + b*(x-z) + c*pow(x-z,2) + d*pow(x-z,3) + e*pow(x-z,4): f + g*(x-y) + h*pow(x-y,2) + i*pow(x-y,3) ;

return res;
}

void FillHist()
{


    int nCut = 0;
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
	LepEff ElecEff1 = LepEff(1,"elec eff1",0.9,0);


	KPiEff KPiEffFake = KPiEff(0.1,1);
	KPiEff PiKEffFake = KPiEff(0.9,1);

	CorrectBf CorrBf = CorrectBf();
	
	double countD1=0;

	int fCount =0;

	TH2F** h2d = new TH2F*[5];
	struct {int bin; float min; float max;} xa,ya;
	xa.bin = 30; xa.min = 0.; xa.max = 3;
	xa.bin = 30; xa.min = 0.3; xa.max = 2.5;
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
	ya.bin = 60; ya.min = 0; ya.max = 180;
	//ya.bin = 20; ya.min = -1.; ya.max = 1;
	h2d[0]= new TH2F("h2d0","h2d0",xa.bin,xa.min,xa.max,ya.bin,ya.min,ya.max);
	h2d[1]= new TH2F("h2d1","h2d1",xa.bin,xa.min,xa.max,ya.bin,ya.min,ya.max);
	h2d[2]= new TH2F("h2d2","h2d2",xa.bin,xa.min,xa.max,ya.bin,ya.min,ya.max);
	h2d[3]= new TH2F("h2d3","h2d3",xa.bin,xa.min,xa.max,ya.bin,ya.min,ya.max);
	h2d[4]= new TH2F("h2d4","h2d4",xa.bin,xa.min,xa.max,ya.bin,ya.min,ya.max);


	cout<<"open file..."<<endl;
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
	#ifdef smallFile
		it_f->Tree = "tree";
	#endif
        cout<<it_f->Path<<" "<<it_f->Name<<" "<<it_f->Tree<<" "<<it_f->FileName<<" "<<endl;

        TFile f( it_f -> FileName.c_str() ,"read");
        t = (TTree*)f.Get(it_f->Tree.c_str());
        cout<<t<<endl;

		//TTreePerfStats *ps =
		//	   new TTreePerfStats(it_f->Name.c_str(),t);


        if(_newTree)
        {
	   // string f_newName = string("/scratch/bctaunu/kakuno/new_0815/small/")+it_f->Name;	
            string f_newName = string("~/small/")+it_f->Name;	
            f_new = new TFile(f_newName.c_str(),"RECREATE","0");
            t_new = new TTree("tree","tree");
            t_new->Branch("run", &run, "etot/F:eher:eler:eth:exp/i:num/i");
            t_new->Branch("evtnr", &evtnr, "evtnr/i");

    	    t_new->Branch("btag",   &btag,   	"pcode_b/f:b_mode:sub1_mod:sub2_mod:sub3_mod:sub4_mod:MCinfo:NB:contNB:m_bc:delta_e:cos_thr:thr:p:p1:p2:p3:e:p1_reco:p2_reco:p3_reco:e_reco:NFS/I:fl_tag/I");
            t_new->Branch("lep1", &lep1,		"elid/f:muid:atckpi:atcpk:atcppi:e/f:p/f:th:dr:dz:es:ps:ths:ec:pc:thc:id_truth/I:q/I:mother_id/I:fl_lep:fl_mother");
            t_new->Branch("gmiss", &gmiss,	"e:p:th:m2:es:ps:ths:ec:pc:thc");
            t_new->Branch("gx", &gx,		"e:p:th:m2:es:ps:ths:ec:pc:thc");
            t_new->Branch("event", &event,	"Eecl:q2:cos_thrA:cos_thrA2:cos_thrAm:cos_thrB:cos_thrC:thrX:thrSigX:thrSigX2:thrSigXm:nch/I:nn/I:npi0/I:q/I:nk:lclass:dclass:BB");
            t_new->Branch("truth", &truth,	"q2/F:p_l:w:costh:m:th:thc:ths:p:pc:ps:miss_m2:miss_e:miss_p:miss_th:miss_es:miss_ps:miss_ths:nkl/I:nch/I");
            t_new->Branch("sBDecay",		&sBDecay, "sBDecay/C");
            t_new->Branch("sDDecay",		&sDDecay, "sDDecay/C");
            //t_new->Branch("thXMiss", &fThXMiss);

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
	//	t->SetBranchAddress("miss", &miss);
	//	t->SetBranchAddress("BDT1", &BDT1);
	//	t->SetBranchAddress("BDT2", &BDT2);
		t->SetBranchAddress("truth", &truth);
#ifndef smallFile
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
#endif
		/*
		t->SetBranchAddress("EFC_energy", &efc_e);
		t->SetBranchAddress("EFC_theta", &efc_th);
		t->SetBranchAddress("Ecl_energy", &cal_e);
		t->SetBranchAddress("Ecl_theta", &cal_th);
*/
		t->SetBranchAddress("sBDecay", &sBDecay);
		t->SetBranchAddress("sDDecay", &sDDecay);
/*
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
*/
		t->SetCacheSize(500000000);
		t->AddBranchToCache("*");

		int evtnr_old =-1,num_old=-1;

		cout<<"Read "<<it_f->FileName<<" with "<<t->GetEntries()<<" Events and weight "<<it_f->weight<<endl;
		for(int i =0; i<t->GetEntries(); ++i)
		{
			t->GetEntry(i);
		//	if(i>100000) break;


			int iLep = lep1.fl_lep%10;
			if(i%10000 == 0) cout<<"\r"<<i*100/t->GetEntries()<<flush;
		//	if(run.exp>37) continue;
			if( btag.m_bc<5.27) continue;
			if( btag.pcode_b*lep1.q>0) continue;
			if( log(btag.NB)<-4 ) continue;
			//if( log(btag.NB)<-0.15) continue;
			//if( run.etot>11.495) continue;
			//if(btag.p<0.3) continue;
			//if( abs(event.q)==0 ) continue;
			//if(lep1.ps < 0.5) continue;
			//if( gmiss.th>0.7 ) continue;	
			//if( gmiss.m2<0. ) continue;
			//if( gmiss.m2>2.5 ) continue;
			//if( iLep == 0 && (lep1.ps<0.5||lep1.ps>2.5) ) continue;
			//if( iLep == 1 && (lep1.ps<0.7||lep1.ps>2.5) ) continue;
			//if(iLep == 0 && lep1.th>0.787) continue;
			//if(event.nch != 4) continue;	
			if(event.cos_thrAm > 0.8) continue;
			//if(event.q2<1.7*1.7) continue;
			//if(it_f->Type == 'd' && lep1.ps<1.1) continue;
			//if(sqrt(sqr(lep1.es) - sqr(lep1.ps))>0.0009) continue;
			//if(event.nn<=3 || event.nn >5) continue;
			//if(event.nn !=4) continue;
			//event.nn+=event.nch;
           // if(1.776-lep1.ps>gmiss.ps) {nCut++;continue;}
			//if(!iLep) continue;
			//if(it_f->Type == 'd' && lep1.ps<1.2) continue;
			//if(btag.MCinfo>0) continue;
			//if(it_f->Type != 'd'){ 
			//	gmiss.m2 -= 0.04 ;
			//}			
			gx.es += gmiss.es;
			gx.m2 = sqrt(gx.m2);
			fHadEvtE = x.es/gx.es;
			fHadEvtP = x.ps/gx.ps;
			fHadEvtM = sqrt(x.m2)/gx.m2;
			//if(fHadEvtM <0.1) continue;
			
			//if(lep1.ps<0.5) continue;
			//if(!(lep1.ps<1.2  && gx.m2<2.6 && gmiss.m2>2.5)) continue;
			//if(gx.m2<1.8) continue;
			if(abs(btag.pcode_b) == 521) continue;	
			
			if(abs(event.lclass) != 2) continue; 
			if(event.dclass == 1 || event.dclass == 3) continue; 
			//	if((btag.pcode_b) == 511) continue;
			
			
			
			
			//if(gx.m2< 2.6) continue;
			//if(iLep != 0) continue;
			fUmiss = gmiss.es - gmiss.ps;
			abspmisspl = lep1.es + gmiss.es;
			pmisspl = 5.27 - gx.es - gmiss.es;
			//if(fUmiss<-1 || fUmiss>0.5) continue;

			
				//if(lep1.elid<0) continue;
			//if(miss.e < 0) continue;

			//if(event.nch + btag.NFS < 10 || event.nch + btag.NFS > 13) continue;
			#ifndef smallFile
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
			#endif

			
			//if(iLep == 0) continue;

			//if(event.cos_thrA>0.8) continue;

			//if(abs(btag.pcode_b) == 521) continue;
			//if(lep1.q*btag.pcode_b>0) continue;
			//cout<<lep1.q<< " "<<btag.pcode_b<<" "<<event.lclass<<" "<<btag.fl_tag<<endl;



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
			
			


			//if(event.nch !=4) continue;
			
			
			if(event.lclass == -2 && event.dclass == 11){
				string sBDec(sBDecay);
				CorrectBf::NoCharge(sBDec);
				if(sBDec.find("511_100411") == 0){	
					event.dclass = 9;		
				} else if (sBDec.find("511_100413") == 0){
					event.dclass = 10;
				}
				
			} else if(event.lclass == 2 && event.dclass == 9){
				string sBDec(sBDecay);
				CorrectBf::NoCharge(sBDec);
				if(sBDec.find("521_100421") == 0){
				    event.dclass = 9;       
				} else if (sBDec.find("521_100423") == 0){
				    event.dclass = 10;
				}
			}
			
		
			
			MC_Correction = 1;
			int iFake   = (lep1.fl_lep != 10) && (lep1.fl_lep != 21);
			
			///if( ((i+2)%4)!=0 && (it_f->Type == 's')){ continue; }
//			

			//if(abs(event.lclass) != 2) continue;
			//if(event.dclass != 5)  continue;

			//remove D**lnu from generic
			if( it_f->Type == 'o'){ if(  abs(event.lclass) == 2 && event.dclass > 4) continue; }
			// add only true D**lnu events
			else if( it_f->Type == 's'){ if(!(abs(event.lclass) == 2 && event.dclass > 4))continue;}
			





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
			//if(fnk==0) continue;

			int cont = GetCont(it_f);
			int cont_old =cont;
			
			fdclass = event.dclass;
			flclass = abs(event.lclass);
			

			//if (cont != charm_lep  && cont != charm_s_lep && cont != charm_ss_lep) continue;

			
            for(int i =0; i<8; i++)
				 Corrections[i] = -1;
			if(cont != data )
			{

				
				//MC_Correction *= hadevtmcorr[hists[0][hadevtmidx]->FindBin(gx.es)-1];
			

				//if(it_f->Type != 'c')
					Corrections[0] = tagcorr(btag.b_mode, (double)btag.NB);
					if(it_f->Type == 'o' || it_f->Type == 'c')
						Corrections[0] *= genMCCorr(run.exp);
				//else Corrections[0]=0.7;
				/*
				float f_p0 = 1.056;
				if(event.lclass>0) //f+0
					Corrections[0] *= f_p0/(1+f_p0)/0.5;
				else if(event.lclass<0)
					Corrections[0] *= 1./(1+f_p0)/0.5;
			*/



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
			//if( it_f->Type == 's' )
			{
				if(abs(event.lclass) == 2 )
				{
						Corrections[1] = CorrBf.weightB(event.lclass, event.dclass);
				//cout<<event.dclass<<" "<<CorrBf.weightB(event.lclass, event.dclass)<<endl;
				}
				else
				{
					Corrections[1] = CorrBf.weightB(sBDecay);

				}	
			}//else Corrections[1] =1;
			Corrections[1] *= CorrBf.weightD(sDDecay);
			//cout<<CorrBf.weightD(sDDecay)<<endl;
			if(abs(event.lclass) == 2 && event.dclass >4){
				if(it_f->Type == 's'){
					Corrections[1] *= CorrBf.FixDss(sDDecay);
					//cout<<Corrections[1]<<endl;
					countD1+=CorrBf.FixDss(sDDecay)*0.39;
					//if(abs(it_f->weight-0.39)>0.001) cout<<"HiLfe!!!!  "<<it_f->weight<<endl;
					//cout<<CorrBf.FixDss(sDDecay)*it_f->weight<<" "<<CorrBf.FixDss(sDDecay)*0.39<<endl;
				}
			}
			
			
			//if(event.dclass == 4 && abs(event.lclass) == 2) cout<<event.lclass<<" "<<Corrections[1]<<endl;
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




				for(int i = 0; i<8; i++) //8
				{
					//if(i ==1 ) continue;
					if(Corrections[i]>=0)
						MC_Correction*=Corrections[i];
				}
				
			
			//MC_Correction *= mm2[hists[0][nVarPs1Elec+5]->FindBin(gmiss.m2)];
			
		//	if(cont>charm_ss_lep)	
				//MC_Correction*=Umiss_corr(fUmiss);
			//if(cont == charm_lep ) MC_Correction *= 0.9;

			} //end of Mc corection
			nEvents[cont]+=MC_Correction*it_f->weight;

			std::vector<MyParticle> vParticles;
			#ifndef smallFile

			float theta = 0.022;
//			
			TLorentzVector lvBeam( run.eher*sin( theta ), 0.0, run.eher*cos(theta ) - run.eler, run.eher + run.eler );
			fMy4s = lvBeam.M();
			TLorentzVector lvBreco(btag.p1,btag.p2,btag.p3,btag.e);
			TLorentzVector lvBsig = lvBeam - lvBreco;
			//TLorentzVector lvBreco(btag.p1_reco,btag.p2_reco,btag.p3_reco,btag.e_reco);
//		

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
            		
			vParticles.clear();
			/*std::vector<float>::iterator itp1_gen = vtrk_gen_p1->begin();
			std::vector<float>::iterator itp2_gen = vtrk_gen_p2->begin();
			std::vector<float>::iterator itp3_gen = vtrk_gen_p3->begin();
			std::vector<float>::iterator ite_gen = vtrk_gen_e->begin();
			std::vector<int>::iterator itid_gen = vtrk_gen_id->begin();
			std::vector<int>::iterator itis_gen = vtrk_gen_ist->begin();
			std::vector<int>::iterator it_moid_gen = vtrk_gen_mo_id->begin();
			*/int nTracks = 0;
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
			TLorentzVector p_true, p_true_reco;
			
			for(; itp1!=vtrk_p1->end();++itp1, ++itp2, ++itp3, ++ite, ++itsep,++itdr, ++itdz, ++itq/*,++itp1_gen, ++itp2_gen, ++itp3_gen, ++ite_gen, ++itid_gen, ++itis_gen, ++it_moid_gen,++itkpid*/)
			{

				//if(cont!=data && gRand.Uniform(1.)<0.37e-2) continue;

				MyParticle mp_tmp;
				mp_tmp.pid = *itsep%10;
				mp_tmp.pid_true = int(*itsep/10);
				mp_tmp.charge = *itq;
				mp_tmp.dr = *itdr;
				mp_tmp.dz = *itdz;
				mp_tmp.P = TLorentzVector(*itp1,*itp2,*itp3,*ite);

				/*mp_tmp.gen_id = *itid_gen;
				mp_tmp.gen_ist = *itis_gen;
				mp_tmp.gen_mo_id = *it_moid_gen;
				mp_tmp.genP = TLorentzVector(*itp1_gen, *itp2_gen,*itp3_gen, *ite_gen);
				
				p_true += mp_tmp.genP;
				*/

				
				if(mp_tmp.P.Perp()<0.1) continue;
				if(abs(mp_tmp.dr)>.5) continue;
				if(abs(mp_tmp.dz)>1.5) continue;

				if( (mp_tmp.P.CosTheta()) > 0.956 || mp_tmp.P.CosTheta()<-0.866 ) continue;

				if(mp_tmp.pid ==5) continue;
				p_true_reco += mp_tmp.genP;
				if(mp_tmp.pid >=4) continue;
				
				
					TLorentzVector tmp = mp_tmp.P;
					//tmp.Boost(-lvBeam.BoostVector());
					SumPsig += tmp.Rho();
				
				
				
				//if(mp_tmp.pid == 3 && mp_tmp.pid_true != 3 && mp_tmp.P.Rho()<.3 ) nFake++;
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
					//if(mp_tmp.pid>1) mp_tmp.P = TLorentzVector(*itp1,*itp2,*itp3,sqrt(pow(mp_tmp.P.Rho(),2) + pow(0.138,2)));
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
				    		

			}
			if(!SignalLep) continue;
		
		/*	for(auto it = vParticles.begin(); it!=vParticles.end();++it) {
				for(auto jt = it+1; jt!=vParticles.end();++jt) {
					if(it->pid_true == jt->pid_true && it->P == jt->P)
					{
						if(it->flag&2==false) it->flag+=2;
						if(jt->flag&2==false) jt->flag+=2;
						cout<<"equal "<< it->flag<<endl;
					}
				}
			}
*/
			TLorentzVector lvGamma;
			itp1 = gam_p1->begin();
			itp2 = gam_p2->begin();
			itp3 = gam_p3->begin();
			int nGamma = 0;
			TLorentzVector lvMaxGamma;
			for(; itp1!=gam_p1->end();++itp1, ++itp2, ++itp3)
			{
				//if(cont!=data && gRand.Uniform(1.)<0.02) continue;				
				TLorentzVector tmp  = TLorentzVector(*itp1, *itp2, *itp3,
						sqrt(pow(*itp1,2)+ pow(*itp2,2)+ pow(*itp3,2)) );
				if(tmp.E()<0.15) continue;
		           
				lvGamma += tmp;
				
				tmp.Boost(-lvBeam.BoostVector());
				if(tmp.E() > lvMaxGamma.E()) lvMaxGamma = tmp;
				nGamma++;
			}
			
			
            event.nn = nGamma;
			
			lvPmisswTag -= lvTracks;
			TLorentzVector lvPmissWoGamma = lvPmisswTag;
			lvPmisswTag -= lvGamma;
			



		//SumPsig -= lvBreco.Rho();
t_track.costh = lvTracks.CosTheta();

		TLorentzVector lvVis = lvTracks + lvGamma-lvSignal;
		TLorentzVector lvVis2 = lvTracks + lvGamma;
        event.q2 = (lvBeam-lvBreco-lvVis).M2();
        sort(vParticles.begin(),vParticles.end(),[](const MyParticle &a, const MyParticle &b){return a.P.Rho() < b.P.Rho();});
        
        
        	TLorentzVector lvPmissLepSys = lvPmisswTag;
        	lvPmissLepSys.Boost(-lvSignal.BoostVector());
        	EmissLepSys = lvPmissLepSys.E();
        
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
		cout<<t_PmisswTag.m2<<endl;
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
				t_gamma.e = lvMaxGamma.E();// lvGamma.E();
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

			#endif

			btag.contNB = log(btag.contNB);
			btag.NB = log(btag.NB);
			fpi0 = (float) event.npi0;
			fnch = (float) event.nch;
			fnn = (float) event.nn;
			fq = (float) abs(event.q);
			fNFS = (float) btag.NFS;
			fb_mode = 2*(fabs(btag.pcode_b) - 511) + fabs(getnum(btag.b_mode));
			
			fl_lep1 = (float)lep1.fl_lep;
			fl_lep2 = (float)lep2.fl_lep;
			fl_mother1 = (float)lep1.fl_mother;
			fl_mother2 = (float)lep2.fl_mother;
			fl_mother12 = fl_mother1+fl_mother2;
			

			fexp = (int)run.exp;
			frunno = (float)run.num;

		



			if(lep1.fl_lep%10 == 0) lep1.muid = 2;
			else if(lep1.fl_lep%10 == 1) lep1.elid = 2;
			if(lep2.fl_lep%10 == 0) lep2.muid = 2;
			else if(lep2.fl_lep%10 == 1) lep2.elid = 2;

			lep1.atckpi = log(lep1.atckpi);
			lep2.atckpi = log(lep2.atckpi);



			//gx.m2 = mysqrt(gx.m2);
			

            MyParticle p[2];
           #ifndef smallFile

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
                    hists[cont][ihistKK+2]->Fill((it->P + jt->P).M(),MC_Correction*it_f->weight);
                    hists_uw[cont][ihistKK+2]->Fill((it->P + jt->P).M());
                }

            }

			

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
				
			}
		    }
			//if(nSlow>1) continue;

		#endif
			
			
			
			if(cont != data)
			{
					h2d[1]->Fill(lep1.p,acos(lep1.th)/M_PI*180.,MC_Correction*it_f->weight);
				}else{
					h2d[0]->Fill(lep1.ps,acos(lep1.th)/M_PI*180.,MC_Correction*it_f->weight);
					//if(cont == charm_lep || cont == charm_s_lep )
						h2d[3]->Fill(lep1.ps,acos(lep1.th)/M_PI*180.,MC_Correction*it_f->weight);
				//	else 
						h2d[4]->Fill(lep1.ps,acos(lep1.th)/M_PI*180.,MC_Correction*it_f->weight);
						
				}
			if(cont == data) {
				MC_Correction=1;
            }
            		if(cont <= charm_ss_lep && cont!=signal_decay)
            			hOverlap[0]->Fill(lep1.ps, MC_Correction*it_f->weight);
            		else
            			hOverlap[1]->Fill(lep1.ps,MC_Correction*it_f->weight);

			for(int j = 0; j<nVariables; ++j)
			{

				if((j==nVarPs1Elec-1) && lep1.fl_lep%10 != 0) continue;
				if((j==nVarPs1Elec) && lep1.fl_lep%10 != 1) continue;
				if((j==nVarPs1Elec+1) && lep1.fl_lep%10 != 0) continue;
				if((j==nVarPs1Elec+2) && lep1.fl_lep%10 != 1) continue;
				if(j<(nVariables -nExtraVariables))
				{
					


					if(j<12){
						cont = GetCont(it_f);
						if(event.nch == j+1)
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
				#ifndef smallFile
				if(j == nVariables -nExtraVariables)
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
				#endif//smallFile
			}

			cont = GetCont(it_f);
			nEventsSplit[iLep][fabs(btag.pcode_b)==511?0:1][cont]+=MC_Correction*it_f->weight ;
			//+=MC_Correction*it_f->weight ;
			nEventsPerFile[0]+=MC_Correction*it_f->weight;
			if(lep1.ps<1.2 && xlep.m2 <2.3)
				nEventsSigReg[cont]+=MC_Correction*it_f->weight;
			if(1)
			{
				if(BDecay[0] != '0')
				{
					BDecays->Fill(NoCharge(BDecay).c_str(),1);
					if(DDecay[0] != '0')
					{
						string stmp = NoCharge(DDecay);
						DDecays->Fill(NoSubDecay(stmp).c_str(),1);
					}
				}
			}
			
			
			btag.contNB = exp(btag.contNB);
			btag.NB = exp(btag.NB);

			if(_newTree) t_new->Fill();

		}

		cout<<"\rDone" <<endl;

		//ps->Print();
		//char sNam[100];
		//sprintf(sNam,"ReadPerf_%s.root",it_f->Name.c_str());
		//ps->SaveAs();
		//delete ps;
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
		
		for(int i=0;i<nContributions;++i)
			cout<<nEventsPerFile[i]<<"\t";
		


	}

	
	cout<<endl;
	for(int i =1; i<=h2d[0]->GetNbinsX();i++)
	{
		for(int j =1; j<=h2d[0]->GetNbinsY();j++)
		{
			float dat = h2d[1]->GetBinContent(i,j);
			float mc = h2d[0]->GetBinContent(i,j);
			float mcc = h2d[3]->GetBinContent(i,j); //const mc
			float mcf = h2d[4]->GetBinContent(i,j); //scale mc

			float tmp=0,dtmp =0;
			if( mc )
			{
				//tmp = (dat-mc)/sqrt(dat);
				tmp = dat/mc;
				dtmp = sqrt(dat)/mc;
			}
			else tmp = 1;
			h2d[2]->SetBinContent(i,j,tmp);
			cout<<h2d[0]->GetXaxis()->GetBinCenter(i)<<" "<<h2d[0]->GetYaxis()->GetBinCenter(j)<<" "<<tmp<<" "<<dtmp<<endl;
			
		} 
		cout<<endl;
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
	
	c2->Clear();
	
	
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



cout<<"cutted "<<nCut<<endl;


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
			//SmearHist(hists[i][j],0.3,0.5);
			//SmearHist(hists_uw[i][j],0.3,0.5);
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
		if(vVar[j].Name == string("mm2gtrans"))
		{
			FILE *fp = fopen("mmiss.txt","w");

			for(int i = 1; i<= sumMC[j]->GetNbinsX(); i++)
            {
                fprintf(fp,"%i\t%f\t",i,hists[data][j]->GetBinCenter(i));
                for(int u = 0; u<nContributions; u++)
                    fprintf(fp,"%f\t",hists[u][j]->GetBinContent(i));
                fprintf(fp,"%f\t%f\n", sumMC[j]->GetBinContent(i), hists[data][j]->GetBinContent(i));
            }
            fclose(fp);
		}
		if(Scale_MC)
		{
			MC2DataScale = hists[data][j]->Integral()/sumMC[j]->Integral();
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

		cout<<i<<" "<<sumMC[j]->GetBinCenter(i)<<" "<<dat/val<<" "<<1./val*sqrt(dat + dat*dat/val)<<"\n";

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
		pad -> Print(sPadName);

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

            if(hists[i][j]->Integral() >0.0001)
			    hists[i][j] -> Scale(1./hists[i][j]->Integral());


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
		
		if(_drawData){
			hists[data][j]->Draw("E1 X0 same");
			sumMC[j]->Draw("HIST same");
			hists[data][j]->GetXaxis()->SetTitle(vVar[j].xTitle.c_str());
			hists[data][j]->GetYaxis()->SetTitle("Events");
			Ymax = hists[data][j]->GetMaximum();
			Ymax +=sqrt(Ymax);
			if(Ymax < sumMC[j]->GetMaximum()) Ymax = sumMC[j]->GetMaximum();
			StackMC[j]->SetMaximum(Ymax*1.05);
		}
		
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
