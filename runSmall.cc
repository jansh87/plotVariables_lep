#include"TROOT.h"
void runSmall()
{
	gROOT->ProcessLine(".L ../shared/SemiLepWeights2D.cc+");
	gROOT->ProcessLine(".L ../shared/SemiLepDssWeights2D.cc+");
	gROOT->ProcessLine(".L ../shared/SemiLepWeightsXc2D.cc+");
	gROOT->ProcessLine(".L ../shared/pid/LepFake.cc+");
	gROOT->ProcessLine(".L ../shared/pid/LepEff.cc+");
	gROOT->ProcessLine(".L ../shared/pid/KPiEff.cc+");
    gROOT->ProcessLine(".L ../shared/CFile.cc+");
	//gROOT->ProcessLine(".L ../shared/B2pilnuFF/formfactor.C+O");
	gROOT->ProcessLine(".L CreateSmall.c+");
	CreateSmall();
}
