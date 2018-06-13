/*
	for the generation of photon-jet observables through Pythia pp collisions 
	by Francesco Vassalli
*/
/* stop trying to make one file do so much.  Make test files */
#include <sstream>
#include <queue>
#include <fstream>
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h" //added plugin for HepMC, think we will need some new library in pythia for this
using namespace Pythia8;
using namespace std;
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TMath.h"
#include "Utils.C" 
#include "TROOT.h"

float deltaPhi(Photon p, Jet j);
float deltaR(Parton,Jet);

template<class T>
void swapPointer(T* a, T* b){
	T* t=a;
	a=b;
	b=t;
}

//inclusive
queue<Jet> getSignificantJet(SlowJet* antikT, float minGeV, float rad){
	queue<Jet> r;
	int i=0;
	while (antikT->pT(i)>=minGeV)
	{
		r.push(Jet(antikT->pT(i),antikT->phi(i),antikT->y(i),rad));
		i++;
	}
	return r;
}

template<class T>
T positivePhi(T in){
	if (in<0)
	{
		in = in+2*TMath::Pi();
	}
	return in;
}

inline bool quickPhotonCheck(Particle p, float gammaCut){
	return p.id()==22&&p.isFinal()&&p.pT()>gammaCut&&TMath::Abs(p.eta())<1.1;
}

queue<myParticle> EventToQueue(Event e){
	myParticle temp;
	queue<myParticle> r;
	for (int i = 0; i < e.size(); ++i)
	{
		temp = myParticle(e[i].id(),e[i].eT(),e[i].phi(),e[i].y());
		r.push(temp);
	}
	return r;
}

void makeData(std::string filename, long nEvents, string pTHat, float gammaCut, bool genHEP){
	using namespace HepMC;
	string hepName = filename+".dat";
    HepMC::IO_GenEvent ascii_io(hepName, std::ios::out); //file where HepMC events will be stored.
	HepMC::Pythia8ToHepMC ToHepMC;    // Interface for conversion from Pythia8::Event to HepMC event.
	filename+=".root";
	TFile* f = new TFile(filename.c_str(),"RECREATE");
  	TTree* interest = new TTree("interest","interest");
  	interest->SetAutoSave(30000);
  	
	/*pythia set up*/
    Pythia pythiaengine;
    pythiaengine.readString("Beams:eCM = 50200."); //LHC VS RHIC
  	pythiaengine.readString("promptphoton:all = on");
 	pythiaengine.readString("HardQCD:all = on");
  	pythiaengine.readString("Random::setSeed = on");
  	pythiaengine.readString("Random::seed =0");
  	pTHat = "PhaseSpace:pTHatMin = "+pTHat+".";
  	pythiaengine.readString(pTHat);
  	pythiaengine.init();
  	/* Tbranching  */
  	SlowJet *antikT2 = new SlowJet(-1,.4,10,2,2,1); 
  	/*init for the TTree*/
  	float asymmetry;
  	float deltaPhi;
  	float e1,e2,pldeltaPhi,psdeltaPhi;
  	float photonpT;
  	/* setting up the branches*/
  	interest->Branch("asymmetry",&asymmetry);
  	interest->Branch("deltaPhi",&deltaPhi);
  	interest->Branch("pldeltaPhi",&pldeltaPhi);
  	interest->Branch("psdeltaPhi",&psdeltaPhi);
  	interest->Branch("e1",&e1);
  	interest->Branch("e2",&e2);
  	interest->Branch("photonpT",&photonpT);

  	/* generation loop*/
    for (int iEvent = 0; iEvent < nEvents; ++iEvent)
  	{
  		if (!pythiaengine.next()){
      		cout<<"pythia.next() failed"<<"\n";
      		iEvent--;
      		continue;
    	}

    	for (int i = 0; i < pythiaengine.event.size(); ++i)
    	{
    		if (isDirect(pythiaengine.info.code())&&quickPhotonCheck(pythiaengine.event[i],gammaCut)) //eta, pT, and photon cut
    		{
    			/*make the event*/	
    			antikT2->analyze(pythiaengine.event);
    			DiJet dJTemp(antikT2,.2,pythiaengine.event[i].phi(),TMath::Pi()/2.0);
    			/*fill the tree*/ 
    			asymmetry=dJTemp.getR2J2();
    			e1=dJTemp.getleading().getpT().value;
    			e2=dJTemp.getsubleading().getpT().value;
    			deltaPhi=dJTemp.getDeltaPhi();
    			pldeltaPhi=dJTemp.getleading().deltaPhi(pythiaengine.event[i].phi());
    			psdeltaPhi=dJTemp.getsubleading().deltaPhi(pythiaengine.event[i].phi());
    			photonpT=pythiaengine.event[i].pT();
    			interest->Fill();
  				if (genHEP)
    			{
    				HepMC::GenEvent* hepmcevt = new HepMC::GenEvent(); //create HepMC "event"
          			ToHepMC.fill_next_event( pythiaengine, hepmcevt ); //convert event from pythia to HepMC
          			ascii_io << hepmcevt;//write event to file
          			delete hepmcevt; //delete event so it can be redeclared next time
    			}
     			break;
    		}
    	}
  	}
  	delete antikT2; //clear the mem
  	interest->Write();
  	f->Close();
  	delete f;
  	f=NULL;
}

int main(int argc, char const *argv[] )
{
	string fileOut = string(argv[1]);
	string pTHat = string(argv[2]);
	float gammaCut= strtod(argv[3],NULL);
	long nEvents =strtol(argv[4],NULL,10);  // 5000000;
	bool genHEP=false;
	makeData(fileOut,nEvents, pTHat, gammaCut,genHEP);
	return 0;
}

inline float deltaPhi(Photon p, Jet j){
	Scalar r= Scalar(TMath::Abs((p.getphi()-j.getphi()).value));
	if (r>TMath::Pi())
	{
		r= r*(-1)+2*TMath::Pi();
	}
	return r.value;
}
inline float deltaR(Parton p, Jet j){
	return TMath::Power(TMath::Power(TMath::Abs(p.getphi()-j.getphi().value),2)+TMath::Power(TMath::Abs(p.gety()-j.gety().value),2),.5);
}