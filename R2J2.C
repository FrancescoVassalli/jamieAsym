using namespace std;

#include <iostream>

int plotCount=0;
inline string getNextPlotName(){
	return "plot"+to_string(plotCount++);
}

queue<Jet> makeJets(float radius,float photonPhi,float* jetphi,float* jety, float* jetpT, float* jetR, float* pz, float* jetm,int SIZE){
	queue<Jet> r;
	
	for (int i = 0; i < SIZE; ++i)
	{
		//phi and R cuts
		if(jetR[i]==radius&&deltaPhi(photonPhi,jetphi[i])>=TMath::Pi()/2.0){ //make the jets that are the radius I want and far enough in phi space from the photon
			r.push(Jet(jetpT[i],jetphi[i],jety[i],jetR[i],pz[i],jetm[i]));
		}
	}

	return r;
}

DiJet makeDiJet(queue<Jet> jQ){
	if (jQ.size()<2)
	{
		return DiJet(false);
	}
	else{
		Jet j1, j2;
		while(!jQ.empty()){ //get the two highest energy jets 
			if (jQ.front().getEnergy().value>j1.getEnergy().value)
			{
				j2=j1;
				j1=jQ.front();
			}
			else if (jQ.front().getEnergy().value>j2.getEnergy().value)
			{
				j2=jQ.front();
			}
			jQ.pop();
		}
		DiJet r =DiJet(j1,j2,true);
		//cout<<j1.getEnergy().value<<'\n';
		return r;
	}
}

inline float deltaPhi(Photon p, Jet j){
	float r= TMath::Abs((p.getphi()-j.getphi()).value);
	if (r>TMath::Pi())
	{
		r= r*(-1)+2*TMath::Pi();
	}
	return r;
}

void plot(TH2F *plot){
	TCanvas *tc = new TCanvas();
	tc->SetRightMargin(.15);
	axisTitles(plot,"#Delta#phi Jet1-Jet2","(E1-E2)/(E1+E2)");
	gPad->SetLogz();
	plot->Scale(1/plot->Integral());
	plot->Draw("colz");
}

void pickR2J2(TChain* interest){
	float asymmetry;
  	float deltaPhi;
  	float deltaR;
  	float deltaEta;
  	float photonpT;
	interest->SetBranchAddress("asymmetry",&asymmetry);
    interest->Branch("deltaEta",&deltaEta);
  	interest->SetBranchAddress("deltaPhi",&deltaPhi);
  	interest->SetBranchAddress("photonpT",&photonpT);
  	interest->Branch("deltaR",&deltaR);
	
	TH2F *p_r2j2 = new TH2F(getNextPlotName().c_str(),"",15,0,.6,15,0,1); 
	for (int i = 0; i < interest->GetEntries(); ++i)
	{
		interest->GetEntry(i);
		p_r2j2->Fill(deltaEta,asymmetry);
	}
	//cout<<"Entries:"<<p_r2j2->GetEntries()<<'\n';
	cout<<interest->GetEntries()<<endl;
	plot(p_r2j2);
}

void handleFile(string name, string extension, int filecount){
	TChain *all = new TChain("interest");
	string temp;
	for (int i = 0; i < filecount; ++i)
	{
		temp = name+to_string(i)+extension;
		all->Add(temp.c_str());
	}
	pickR2J2(all);
}

void R2J2(){
	string fileLocation = "/home/user/Droptemp/LHCj/";
	string filename = "LHCGamma";
	string extension = ".root";
	string temp = fileLocation+filename;
	handleFile(temp,extension,100);

}