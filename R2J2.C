using namespace std;

#include <iostream>
namespace {
	int plotCount=0;
}
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

std::vector<TH1F*> makeTH1Farray(float min, float max, float plotwidth, float Nbins){
	int Nplots = (max-min)/plotwidth;
	vector<TH1F*> rarry; //= new TH1F*[Nplots];
	for (int i = 0; i < Nplots; ++i)
	{
		rarry[i]= new TH1F(getNextPlotName(&plotCount).c_str(),"",Nbins,min+i*plotwidth,min+(i+1)*plotwidth);
	}
	return rarry;
}

void plot(TH2F *plot, TProfile* prof){
	TCanvas *tc = new TCanvas();
	tc->SetRightMargin(.15);
	axisTitles(plot,"#DeltaR Jet1-Jet2","(E1-E2)/(E1+E2)");
	gPad->SetLogz();
	plot->Scale(1/plot->Integral());
	plot->Draw("colz");
	prof->Draw("same");
}
void plotAve(TH2F *plot, TProfile* prof){
	TCanvas *tc = new TCanvas();
	tc->SetRightMargin(.15);
	axisTitles(plot,"Energy Sum","(E1-E2)/(E1+E2)");
	gPad->SetLogz();
	plot->Scale(1/plot->Integral());
	plot->Draw("colz");
	prof->Draw("same");
}
void plot1d(TH1F *plot,string xTitle, string yTitle){
	TCanvas *tc = new TCanvas();
	//tc->SetRightMargin(.15);
	axisTitles(plot,xTitle.c_str(),yTitle.c_str());
	//gPad->SetLogz();
	plot->Scale(1/plot->Integral());
	plot->Draw("P0");
}

void pickR2J2(TChain* interest){
	float asymmetry;
  	float deltaPhi;
  	float deltaR;
  	float deltaEta;
  	float photonpT;
	float e1;
	float e2;  	
	interest->SetBranchAddress("asymmetry",&asymmetry);
    interest->SetBranchAddress("deltaEta",&deltaEta);
  	interest->SetBranchAddress("deltaPhi",&deltaPhi);
  	interest->SetBranchAddress("photonpT",&photonpT);
  	interest->SetBranchAddress("deltaR",&deltaR);
  	interest->SetBranchAddress("e1",&e1);
  	interest->SetBranchAddress("e2",&e2);
	
	TH2F *p_r2j2 = new TH2F(getNextPlotName(&plotCount).c_str(),"",100,0,4,100,0,1); 
	TProfile *profp2j2 = new TProfile(getNextPlotName(&plotCount).c_str(),"",100,0,4,0,1);
	TH2F *ave = new TH2F(getNextPlotName(&plotCount).c_str(),"",35,20,350,20,0,1);
	TProfile *aveProf = new TProfile(getNextPlotName(&plotCount).c_str(),"",35,20,350,0,1);
	TH1F *delR1 = new TH1F(getNextPlotName(&plotCount).c_str(),"",500,0,4);
	TH1F *asym1 = new TH1F(getNextPlotName(&plotCount).c_str(),"",200,0,1);
	//vector<TH1F*> splits = makeTH1Farray(0,1,.2,20);
	for (int i = 0; i < interest->GetEntries(); ++i)
	{
		interest->GetEntry(i);
		p_r2j2->Fill(deltaR,asymmetry);
		profp2j2->Fill(deltaR,asymmetry);
		ave->Fill(e1+e2,asymmetry);
		aveProf->Fill(e1+e2,asymmetry);
		delR1->Fill(deltaR);
		asym1->Fill(asymmetry);
		//cout<<deltaR<<"\n";
	}
	//cout<<"Entries:"<<p_r2j2->GetEntries()<<'\n';
	//cout<<interest->GetEntries()<<endl;
	plot(p_r2j2, profp2j2);
	//plotAve(ave,aveProf);
	//plot1d(delR1,"#DeltaR","count");
	//plot1d(asym1,"asymmetry","count");
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