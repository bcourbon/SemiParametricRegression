{

//{

TString libstr(Form("%s/lib/%s/%s",gSystem->Getenv("CMSSW_BASE"),gSystem->Getenv("SCRAM_ARCH"),"libHiggsAnalysisGBRLikelihood.so")); 
gSystem->Load(libstr);

//}

gInterpreter->AddIncludePath((TString(":")+TString(gSystem->Getenv("ROOFITSYS"))+TString("/include")).Data());

gSystem->AddIncludePath("-I$CMSSW_BASE/src/HiggsAnalysis/GBRLikelihood/interface");
gSystem->AddIncludePath("-I$ROOFITSYS/include");

}

