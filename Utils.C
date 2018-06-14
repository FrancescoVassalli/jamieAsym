#include "TMath.h"
#include <queue>
#include <iostream>

template<class T>
void clear( std::queue<T> &q )
{
   std::queue<T> empty;
   std::swap( q, empty );
}

template<class T>
T* vectorToArray(std::vector<T> v){
	T *rarry= new T[v.size()];
	for (int i=0; i<v.size();i++)
	{
		rarry[i]=v[i];
	}
	return rarry;
}

float deltaPhi(float i1, float i2);
#ifndef Scalar_h
#define Scalar_h
using namespace std;
class Scalar
{
public:
	Scalar(){
		value=uncertainty=0;
	}
	Scalar(float value){
		this->value=value;
		uncertainty=0;
	}
	Scalar(float value, float uncertainty){
		this->value=value;
		this->uncertainty=uncertainty;
	}
	Scalar(const Scalar &s){
		this->value=s.value;
		uncertainty=s.uncertainty;
	}
	~Scalar(){}; 
	//operators might need to return class types but I'm not sure
	Scalar operator+(const Scalar &s){
		//cout<<"Addition "<<value<<" + "<<s.value<<"\n";;
		Scalar next = Scalar();
		next.value = this->value + s.value;
		next.uncertainty = quadrature(this->uncertainty, s.uncertainty);
		return next;
	}
	Scalar operator/(const Scalar &s){
		Scalar next = Scalar();
		next.uncertainty = (value/s.value)*quadrature(this->uncertainty/value, s.uncertainty/s.value);
		next.value = value/s.value;
		if (value==0||s.value==0)
		{
			next.uncertainty=0;
		}
		return next;
	}
	void operator+=(Scalar s){
		value+=s.value;
		uncertainty = quadrature(this->uncertainty, s.uncertainty);
	}
	void operator*=(Scalar s){
		*this = *this * s;
	}
	void operator/=(Scalar s){
		*this = *this / s;
	}
	void operator-=(Scalar s){
		*this = *this - s;
	}
	Scalar operator+(float s){
		Scalar next;
		next.value = this->value + s;
		return next;
	}
	Scalar operator-(Scalar s){
		Scalar next;
		next.value = this->value - s.value;
		next.uncertainty = quadrature(this->uncertainty, s.uncertainty);
		return next;
	}
	Scalar operator-(float s){
		Scalar next;
		next.value = this->value - s;
		return next;
	}
	void operator=(Scalar s){
		value=s.value;
		uncertainty=s.uncertainty;
	}
	bool operator==(Scalar s){
		return value==s.value;
	}
	bool operator>(Scalar s){
		return value>s.value;
	}
	bool operator<(Scalar s){
		return value<s.value;
	}
	bool operator>=(Scalar s){
		return value>=s.value;
	}
	bool operator<=(Scalar s){
		return value<=s.value;
	}
	bool operator!=(Scalar s){
		return value!=s.value;
	}
	Scalar operator*(Scalar s){
		Scalar next;
		next.uncertainty = (value*s.value)*quadrature(this->uncertainty/value, s.uncertainty/s.value);
		if (value==0||s.value==0)
		{
			next.uncertainty=0;
		}
		next.value = value*s.value;
		return next;
	}
	Scalar operator*(float s){
		Scalar next;
		next.value = value*s;
		next.uncertainty = (value*s)*uncertainty/value;
		if (value==0||s==0)
		{
			next.uncertainty=0;
		}
		return next;
	}
	Scalar operator/(float s){ //progate uncertainty
		Scalar next;
		next.uncertainty=(value/s)*uncertainty/value;
		if (value==0)
		{
			next.uncertainty=0;
		}
		next.value = value/s;
		return next;
	}
	Scalar pow(double n){
		Scalar r;
		if (n==0)
		{
			r.value=1;
			r.uncertainty=0;
			return r;
		}
		r.value = TMath::Power((double)value,n);
		r.uncertainty=(r.value*n*uncertainty/value);
		return r;
	}
	Scalar log(float base){
		Scalar r;
		r.value = TMath::Log(value)/TMath::Log(base);
		r.uncertainty = uncertainty/(value*TMath::Log(base));
		return r;
	}
	friend std::ostream& operator<<(std::ostream& os, Scalar const & tc) {
       return os <<"Scalar:" << tc.value <<char(241)<<tc.uncertainty<<'\n';
    }

	float value;
	float uncertainty;
protected:
	template<class T>
T quadrature(T d1, T d2){
	return TMath::Sqrt((double)d1*d1+d2*d2);
}

template<class T>
T quadrature(T* a, int SIZE){
	T* b = clone(a);
	arrayMultiply(b,b,SIZE);
	T r = sum(b,SIZE);
	return TMath::Sqrt(r);
}
	
};
#endif
inline bool inRange(float in, float low, float high){ // inclusive-exclusive
	return in>=low&&in<high;
}
#ifndef Angle_h
#define Angle_h
//not sure of all the inheritence rules on virtualization/abstracttion want it to do all the same things as a Scalar but check for <PI b4 all the returns
class Angle : public Scalar //all the members of Scalar become members of Angle with equal scope
{ //cannot be used yet it is not correct 
public:
	Angle(){
		value=uncertainty=0;
	}
	~Angle(){}
	Angle(float value){
		this->value=value;
		uncertainty=0;
	}
	Angle(float value, float uncertainty){
		this->value=value;
		this->uncertainty=uncertainty;
	}
	Angle(const Scalar &s){
		this->value=s.value;
		uncertainty=s.uncertainty;
	}
	Angle operator+(float s){
		Angle next;
		next.value = this->value + s;
		return next;
	}
	Angle operator-(Scalar s){
		Angle next;
		next.value = this->value - s.value;
		next.uncertainty = quadrature(this->uncertainty, s.uncertainty);
		return next;
	}
	Angle operator-(float s){
		Angle next;
		next.value = this->value - s;
		return next;
	}
	Angle operator*(Scalar s){
		Angle next;
		next.uncertainty = (value*s.value)*quadrature(this->uncertainty/value, s.uncertainty/s.value);
		if (value==0||s.value==0)
		{
			next.uncertainty=0;
		}
		next.value = value*s.value;
		makeIn2Pi();
		return next;
	}
	Angle operator*(float s){
		Angle next;
		next.value = value*s;
		next.uncertainty = (value*s)*uncertainty/value;
		if (value==0||s==0)
		{
			next.uncertainty=0;
		}
		makeIn2Pi();
		return next;
	}
	Angle operator/(float s){ //progate uncertainty
		Angle next;
		next.uncertainty=(value/s)*uncertainty/value;
		if (value==0)
		{
			next.uncertainty=0;
		}
		next.value = value/s;
		makeIn2Pi();
		return next;
	}
	Angle operator+(const Scalar &s){
		//cout<<"Addition "<<value<<" + "<<s.value<<"\n";;
		Angle next = Scalar();
		next.value = this->value + s.value;
		next.uncertainty = quadrature(this->uncertainty, s.uncertainty);
		makeIn2Pi();
		return next;
	}
	Angle operator/(const Scalar &s){
		Angle next = Scalar();
		next.uncertainty = (value/s.value)*quadrature(this->uncertainty/value, s.uncertainty/s.value);
		next.value = value/s.value;
		if (value==0||s.value==0)
		{
			next.uncertainty=0;
		}
		makeIn2Pi();
		return next;
	}
	void operator+=(Scalar s){
		value+=s.value;
		uncertainty = quadrature(this->uncertainty, s.uncertainty);
	}
	void operator*=(Scalar s){
		*this = *this * s;
	}
	void operator/=(Scalar s){
		*this = *this / s;
	}
	void operator-=(Scalar s){
		*this = *this - s;
	}
	void operator=(Scalar s){
		value=s.value;
		uncertainty=s.uncertainty;
	}
	bool operator==(Scalar s){
		return value==s.value;
	}
	bool operator>(Scalar s){
		return value>s.value;
	}
	bool operator<(Scalar s){
		return value<s.value;
	}
	bool operator>=(Scalar s){
		return value>=s.value;
	}
	bool operator<=(Scalar s){
		return value<=s.value;
	}
	bool operator!=(Scalar s){
		return value!=s.value;
	}
	Angle pow(double n){
		Angle r;
		if (n==0)
		{
			r.value=1;
			r.uncertainty=0;
			return r;
		}
		r.value = TMath::Power((double)value,n);
		r.uncertainty=(r.value*n*uncertainty/value);
		makeIn2Pi();
		return r;
	}
	Angle log(float base){
		Angle r;
		r.value = TMath::Log(value)/TMath::Log(base);
		r.uncertainty = uncertainty/(value*TMath::Log(base));
		makeIn2Pi();
		return r;
	}
	friend std::ostream& operator<<(std::ostream& os, Angle const & tc) {
       return os <<"Radian:" << tc.value <<char(241)<<tc.uncertainty<<'\n';
    }
private:
	inline void makeIn2Pi(){
		if(value<0) value+=2*TMath::Pi();
		if (value>2*TMath::Pi()) value-=2*TMath::Pi();
	}
};
#endif

#ifndef Point_h
#define Point_h
struct Point
{
	Scalar x;
	Scalar y;
};
#endif
#ifndef Parton_h
#define Parton_h
#include "TLorentzVector.h"
class Parton
{
public:
	Parton(){}
	~Parton(){}
	Parton(int ID, float phi, float y){
		quark=isQuark(ID);
		this->phi=phi;
		this->y=y;
	}
	Parton(int ID, float phi, float eta, float eT, float e){
		quark=isQuark(ID);
		TLorentzVector tlv;
		tlv.SetPtEtaPhiE(eTTopT(eT,ID),eta,phi,e);
		this->phi=phi;
		y=tlv.Rapidity();
		this->eta = eta;
	}
	float getphi(){
		return phi;
	}
	float gety(){
		return y;
	}
	float geteta(){
		return eta;
	}
	bool getQuark(){
		return quark;
	}
protected:
	bool quark;
	float phi;
	float y;	
	float eta;
	bool isQuark(int ID){
		return TMath::Abs(ID)>0&&TMath::Abs(ID)<9;
	}
	float eTTopT(float eT, int ID){
		return TMath::Power(eT*eT-idToMass(ID)*idToMass(ID),.5);
	}
	float idToMass(int ID){
		ID=TMath::Abs(ID);
		switch(ID){
			case 2: return 0.0022;
				break;
			case 1: return 0.0047;
				break;
			case 3: return 0.096;
				break;
			case 4: return 1.28;
				break;
			case 5: return 4.18;
				break;
			case 6: return 173;
				break;
			case 21: return 0;
				break;
			default: return -1;
				break;
		}
	}
};
#endif

template<class T>
T bigger(T x, T y){
	if (x>y)
	{
		return x;
	}
	else{
		return y;
	}
}
template<class T>
T smaller(T x, T y){
	if (x<y)
	{
		return x;
	}
	else{
		return y;
	}
}

bool isDirect(int i)
{
  if(i > 200 and i < 267)
  {
    return true;
  }
  else
  {
    return false;
  }
}
#ifndef myParticle_h
#define myParticle_h
class myParticle
{
public:
	myParticle(){}
	~myParticle(){}
	myParticle(int id, float eT, float phi, float y){
		this->id=id;
		this->eT=eT;
		this->phi=phi;
		this->y=y;
	}
	int getId(){
		return id;
	}
	float geteT(){
		return eT;
	}
	float getphi(){
		return phi;
	}
	float gety(){
		return y;
	}

private:
	int id;
	float eT;
	float phi;
	float y;
};
#endif
#ifndef Photon_h
#define Photon_h
class Photon
{
public:
	Photon(){}
	Photon(float _pT){
		this->pT = Scalar(_pT);
	}
	Photon(float _pT,float _phi){
		this->pT = Scalar(_pT);
		this->phi= Scalar(_phi);
	}
	Photon(double _pT,double _phi){
		this->pT = Scalar((float)_pT);
		this->phi= Scalar((float)_phi);
	}
	Photon(double _pT,double _phi,bool process){
		this->pT = Scalar((float)_pT);
		this->phi= Scalar((float)_phi);
		direct=process;
	}
	Photon(float _pT,float _phi, float _eta){
		this->pT = Scalar(_pT);
		this->phi= Scalar(_phi);
		this->eta= Scalar(_eta);
	}
	Photon(double _pT,double _phi, double _eta){
		this->pT = Scalar((float)_pT);
		this->phi= Scalar((float)_phi);
		this->eta= Scalar((float)_eta);
	}
	Photon(double _pT,double _phi, double _eta, bool process){
		this->pT = Scalar((float)_pT);
		this->phi= Scalar((float)_phi);
		this->eta= Scalar((float)_eta);
		direct=process;
	}
	Photon(double _pT,double _phi, double _eta, bool process, std::queue<myParticle> all){
		this->pT = Scalar((float)_pT);
		this->phi= Scalar((float)_phi);
		this->eta= Scalar((float)_eta);
		direct=process;
		findIsoEt(all);
	}
	Photon(int index,double _pT,double _phi, double _eta, bool process, queue<myParticle> all){
		this->pT = Scalar((float)_pT);
		this->phi= Scalar((float)_phi);
		this->eta= Scalar((float)_eta);
		direct=process;
		position=index;
		findIsoEt(all);
	}
	Photon(double _pT,double _phi, double _eta,queue<myParticle> all){
		this->pT = Scalar((float)_pT);
		this->phi= Scalar((float)_phi);
		this->eta= Scalar((float)_eta);
		findIsoEt(all);
	}
	Photon(int SIZE, int position, float* eT, float* phi, float* eta, bool direct,float eTCone){
		etCone=eTCone;
		this->position=position;
		pT=Scalar(eT[position]);
		this->phi=Scalar(phi[position]);
		this->eta=Scalar(eta[position]);
		this->direct=direct;
		findIsoEt(phi,eta,eT,SIZE);
	}
	Photon(int SIZE, int* id, float* eT, float* phi, float* eta, bool direct,float eTCone,float eTCut){
		etCone=eTCone;
		findPosition(SIZE,id,eT,eTCut);
		pT=Scalar(eT[position]);
		this->phi=Scalar(phi[position]);
		this->eta=Scalar(eta[position]);
		this->direct=direct;
		findIsoEt(phi,eta,eT,SIZE);
	}
	~Photon(){}
	Scalar getpT(){
		return pT;
	}
	void setpT(float _pT){
		this->pT = Scalar(_pT);
	}
	Scalar getphi(){
		return phi;
	}
	void setphi(float _phi){
		this->phi = Scalar(_phi);
	}
	Scalar geteta(){
		return eta;
	}
	void seteta(float _eta){
		this->eta = Scalar(_eta);
	}
	void setParton(Parton p){
		parton=p;
	}
	bool isDirect(){
		return direct;
	}
	int getPosition(){
		return position;
	}
	float findIsoEt(queue<myParticle> all){
		isoEt=0;
		while(!all.empty()){
			if (inCone(all.front().gety(),all.front().getphi()))
			{
				isoEt+=all.front().geteT();
			}
			all.pop();
		}
		return isoEt;
	}

	float getIsoEt(){
		return isoEt;
	}
	
	Point getAngle(){
		Point r;
		r.x=phi;
		r.y=eta;
		return r;
	}

private:
	Scalar pT;
	Scalar phi;
	Scalar eta;
	bool direct;
	float isoEt;
	Parton parton;
	float etCone = 0.3;
	int position;
	bool inCone(float geta, float gphi)
	{
	  if( sqrt(TMath::Power(TMath::Abs(geta-eta.value),2)+TMath::Power(TMath::Abs(gphi-phi.value),2)) < etCone )
	  {
	    return true;
	  }
	  else
	  {
	    return false;
	  }
	}
	float findIsoEt(float* phi, float* eta, float* eT, int SIZE){
		isoEt=0;
		for(int i=0;i<SIZE;i++){
			if (inCone(eta[i],phi[i]))
			{
				isoEt+=eT[i];
			}
		}
		isoEt-=pT.value; //take the photon out
		return isoEt;
	}
	inline bool isPhoton(int id){
		return id==22;
	}
	int findPosition(int SIZE, int* id, float* et,float eTCut){
		for (int i = 0; i < SIZE; ++i)
		{
			if (isPhoton(id[i])&&et[i]>eTCut)
			{
				position=i;
				return position;
			}
		}
		return -1;
	}
};
#endif
#ifndef Jet_h
#define Jet_h
class Jet
{
public:
	Jet(){
		pT=0;
		energy=0;
	}
	Jet(float _pT, float _phi, float _y, float _r){
		this->pT =Scalar(_pT);
		this->phi = Scalar(_phi);
		this->y = Scalar(_y);
		this->r = Scalar(_r);
	}
	Jet(double _pT, double _phi, double _y){
		pT=Scalar(_pT);
		phi=Scalar(_phi);
		y=Scalar(_y);
		r=0;
	}
	Jet(float _pT, float _phi, float _y, float _r, float pz){ // calculate eta
		this->pT =Scalar(_pT);
		this->phi = Scalar(_phi);
		this->y = Scalar(_y);
		this->r = Scalar(_r);
		eta= Scalar(calculateEta(_pT,pz));
	}
	Jet(float _pT, float _phi, float _y, float _r, float pz, float mass){ // calculate eta and e
		this->pT =Scalar(_pT);
		this->phi = Scalar(_phi);
		this->y = Scalar(_y);
		this->r = Scalar(_r);
		this->mass = Scalar(mass);
		energy = Scalar(calculateEnergy(pz));
		eta= Scalar(calculateEta(_pT,pz));
	}
	Jet(float _pT, float _phi, float _y, float _r, float pz, float mass,float energy){ // calculate eta
		this->pT =Scalar(_pT);
		this->phi = Scalar(_phi);
		this->y = Scalar(_y);
		this->r = Scalar(_r);
		this->mass = Scalar(mass);
		this->energy = Scalar(energy);
		eta= Scalar(calculateEta(_pT,pz));
	}
	Jet(float _pT, float _phi, float _y, float _r, float pz, float mass,float energy,std::vector<int> constiutuentIndices,Event e){ // calculate eta set constiutuent particles
		this->pT =Scalar(_pT);
		this->phi = Scalar(_phi);
		this->y = Scalar(_y);
		this->r = Scalar(_r);
		this->mass = Scalar(mass);
		this->energy = Scalar(energy);
		eta= Scalar(calculateEta(_pT,pz));
		mult=constiutuentIndices.size();
		for (std::vector<int>::iterator i = constiutuentIndices.begin(); i != constiutuentIndices.end(); ++i)
		{
			pTs.push_back(e[*i].pT());
			etas.push_back(e[*i].eta());
			phis.push_back(e[*i].phi());
		}
	}
	Jet(float _pT, float _phi, float _y, float _r, float pz, float mass,float energy,int mult){ // calculate eta set constiutuent particles
		this->pT =Scalar(_pT);
		this->phi = Scalar(_phi);
		this->y = Scalar(_y);
		this->r = Scalar(_r);
		this->mass = Scalar(mass);
		this->energy = Scalar(energy);
		eta= Scalar(calculateEta(_pT,pz));
		this->mult=mult;
	}
	~Jet(){	}
	void setMult(int m){
		mult=m;
	}
	void setConstituents(std::vector<int> indices, Event e){
		mult=indices.size();
		for (std::vector<int>::iterator i = indices.begin(); i != indices.end(); ++i)
		{
			pTs.push_back(e[*i].pT());
			etas.push_back(e[*i].eta());
			phis.push_back(e[*i].phi());
		}
		cout<<"here:"<<etas[0]<<"\n";
	}
	Parton setParton(Parton p1, Parton p2){  
		if (deltaR(p1)<deltaR(p2))//comparision
		{
			parton=p1;	
		}
		else{
			parton=p2;
		}
		return parton;
	}
	inline float deltaPhi(float in){
		float r = TMath::Abs(in-phi.value);
		if (r>TMath::Pi())
		{
			r= 2*TMath::Pi()-r;
		}
		return r;
	}
	void setParton(Parton p){
		parton=p;
	}
	bool isJetQuark(){
		return parton.getQuark();
	}
	int getmult(){
		return mult;
	}
	void fill(int* size, float* phi, float* eta, float* pT){
		*size=mult;
		phi= vectorToArray(phis);
		eta = vectorToArray(etas);
		pT= vectorToArray(pTs);
		cout<<"there"<<etas[0]<<"\n";
	}
	Scalar getpT(){
		return pT;
	}
	Scalar getphi(){
		return phi;
	}
	Scalar gety(){
		return y;
	}
	Scalar getr(){
		return r;
	}
	Scalar getEnergy(){
		return energy;
	}
	Scalar geteta(){
		return eta;
	}
	Scalar operator/(float s){ 
		return pT/s;
	}
	Scalar operator/(Jet s){ 
		return pT/s.getpT();
	}
	bool operator==(Jet s){
		return pT==s.pT&&phi==s.getphi();
	}
	bool operator>(Jet s){
		return pT>s.getpT();
	}
	bool operator<(Jet s){
		return pT<s.getpT();
	}
	bool operator>=(Jet s){
		return pT>=s.getpT();
	}
	bool operator<=(Jet s){
		return pT<=s.getpT();
	}
	bool operator!=(Jet s){
		return pT!=s.getpT();
	}
	friend std::ostream& operator<<(std::ostream& os, Jet const & tc) {
       return os <<"Jet energy:" << tc.energy.value <<" phi:"<<tc.phi.value<<" eta:"<<tc.eta.value<<'\n';
    }
private:
	Scalar pT=-1;
	Scalar phi =7;
	Scalar y=0;
	Scalar r=-1;
	Scalar eta;
	int mult=0;
	Scalar mass;
	Scalar energy;
	Parton parton;
	std::vector<float> pTs;
	std::vector<float> etas;
	std::vector<float> phis;

	inline float deltaPhi(float i1, float i2){
		float r = TMath::Abs(i1-i2);
		if (r>TMath::Pi())
		{
			r= 2*TMath::Pi()-r;
		}
		return r;
	}
	float deltaR(Parton p){
	  return TMath::Power((TMath::Power(TMath::Abs(p.geteta()-eta.value),2)+TMath::Power(deltaPhi(p.getphi(),phi.value),2)),.5);
	}
	float calculateEta(float pt, float pz){
		return .5* TMath::Log((TMath::Power(pt*pt+pz*pz,.5)+pt))/((TMath::Power(pt*pt+pz*pz,.5)-pt));
	}
	float calculateEnergy(float pz){
		return TMath::Power((pT.value)*(pT.value)+pz*pz+mass.value*mass.value,.5);
	}
	
};
#endif
#ifndef XjPhi_h
#define XjPhi_h
class XjPhi
{
public:
	XjPhi(float pt1, float phi1, float pt2, float phi2){
		//calculate scalar Xjs from parameters 
		xj=Scalar(smaller(pt1,pt2)/bigger(pt1,pt2));
		phi=Scalar(deltaPhi(phi1,phi2));
	}
	XjPhi(Jet leading, Jet subleading){
		xj=Scalar(subleading/leading);
		phi = Scalar(deltaPhi(leading.getphi().value,subleading.getphi().value));
	}
	XjPhi(Photon p, Jet j){
		xj=j.getpT()/p.getpT();
		phi = Scalar(deltaPhi(p.getphi().value,j.getphi().value));
		photon=p;
		jet=j;
	}
	XjPhi(Scalar xj, Scalar phi){
		this->xj=xj;
		this->phi=phi;
	}
	XjPhi(){}
	~XjPhi(){}
	Scalar getXj(){
		return xj;
	}
	Scalar getphi(){
		return phi;
	}
	Photon getPhoton(){
		return photon;
	}
	bool isQuark(){
		return jet.isJetQuark();
	}
	Jet getJet(){
		return jet;
	}
	friend ostream& operator<<(ostream& os, XjPhi const & tc) {
        return os << tc.xj.value <<','<<tc.phi.value<<'\n';
    }
	
private:
	Scalar xj;
	Scalar phi;
	Photon photon;
	Jet jet;
};
#endif
#ifndef DiJet_h
#define DiJet_h
class DiJet
{
public:
	DiJet(Jet j1, Jet j2){
		leading = bigger(j1,j2);
		subleading=smaller(j1,j2);
		calculateR2J2();
	}
	DiJet(Jet j1, Jet j2,bool t){ //for sorted jets 
		if (t)
		{
			isDijet=true;
			leading=j1;
			subleading=j2;
		}
		else{
			leading = bigger(j1,j2);
			subleading=smaller(j1,j2);
		}
		calculateR2J2();
	}
	DiJet(double pt1, double phi1, double pt2,double phi2){
		if (pt1>pt2)
		{
			leading=Jet(pt1,phi1,0,0);
			subleading=Jet(pt2,phi2,0,0);
		}
		else{
			subleading=Jet(pt1,phi1,0,0);
			leading=Jet(pt2,phi2,0,0);
		}
		calculateR2J2();
		makeJetDeltaPhi();
	}
	//makes a Dijet event given a SlowJet that has already been filled it will exclude jets with deltaphi less than phiRange from phi0
	DiJet(SlowJet* antikT, float radius,float phi0, float phiRange){
		int sizeJet = antikT->sizeJet();
		int count=0;
		isDijet=false;
	if (sizeJet>2)
		{
			for (int i = 0; i < antikT->sizeJet(); ++i)
			{
				if (deltaPhi(antikT->phi(i),phi0)>phiRange)
				{
					if(count==0){
						leading = Jet(antikT->pT(i),antikT->phi(i),antikT->y(i),radius,(antikT->p(i)).pz(),antikT->m(i),(antikT->p(i)).e());
						jet1Consit=antikT->constituents(i);
						count++;
					}
					else if (count==1)
					{
						subleading = Jet(antikT->pT(i),antikT->phi(i),antikT->y(i),radius,(antikT->p(i)).pz(),antikT->m(i),(antikT->p(i)).e());
						jet2Consit=antikT->constituents(i);
						count++;
						break;
					}
				}
			}
			isDijet= (count==2);
			calculateR2J2();
			makeJetDeltaPhi();
			makeJetDeltaEta();
			makeJetDeltaR();
			
		}
	}
	DiJet(bool f){
		isDijet=f;
	}
	DiJet(){}

	~DiJet(){}
	void setConstituents(Event e){
		leading.setConstituents(jet1Consit,e);
		subleading.setConstituents(jet2Consit,e);
	}
	Jet getleading(){
		return leading;
	}
	Jet getsubleading(){
		return subleading;
	}
	void fill(bool leading, int* size, float* phi, float* eta, float* pT){
		if (leading)
		{
			this->leading.fill(size, phi,eta,pT);
		}
		else{
			subleading.fill(size, phi,eta,pT);
		}
	}
	float getR2J2(){
		return r2j2;
	}
	float getDeltaPhi(){
		return jetDeltaPhi;
	}
	float getDeltaR(){
		return jetDeltaR;
	}
	float getDeltaEta(){
		return jetDeltaEta;
	}
	void operator=(DiJet d2){
		isDijet=(bool)d2;
		leading=d2.getleading();
		subleading=d2.getsubleading();
		jetDeltaPhi=d2.getDeltaPhi();
	}
	operator bool(){
		return isDijet;
	}
	friend std::ostream& operator<<(std::ostream& os, DiJet const & tc) {
       return os <<"Leading " << tc.leading <<"Subleading "<<tc.subleading;
    }
private:
	Jet leading;
	Jet subleading;
	float jetDeltaPhi;
	float jetDeltaR;
	float jetDeltaEta;
	float photonDeltaPhi;
	float r2j2;
	bool isDijet;
	std::vector<int> jet1Consit;
	std::vector<int> jet2Consit;

	inline float deltaPhi(float i1, float i2){
		float r = TMath::Abs(i1-i2);
		if (r>TMath::Pi())
		{
			r= 2*TMath::Pi()-r;
		}
		return r;
	}

	inline float makeJetDeltaPhi(){
		jetDeltaPhi = TMath::Abs(leading.getphi().value-subleading.getphi().value);
		if (jetDeltaPhi>TMath::Pi())
		{
			jetDeltaPhi= 2*TMath::Pi()-jetDeltaPhi;
		}
		return jetDeltaPhi;
	}

	inline void calculateR2J2(){
		r2j2 = (leading.getEnergy().value-subleading.getEnergy().value)/(leading.getEnergy().value+subleading.getEnergy().value);
	}
	inline void makeJetDeltaR(){
	  jetDeltaR=TMath::Power((TMath::Power(jetDeltaEta,2)+TMath::Power(jetDeltaPhi,2)),.5);
	}
	inline void makeJetDeltaEta(){
		jetDeltaEta=TMath::Abs(leading.geteta().value-subleading.geteta().value);
	}
	
};
#endif


float** qXjPhiTo2DArray(queue<XjPhi> in){
	float** out = new float*[2];
	float *o1 = new float[in.size()];
	float *o2 = new float[in.size()];
	int i=0;
	while(!in.empty()){
		o1[i]=in.front().getXj().value;
		o2[i]=in.front().getphi().value;
		i++;
		in.pop();
	}
	out[0]=o1;
	out[1]=o2;
	return out;
}

//takes two angles and returns their difference in Phi must be less than Pi
inline float deltaPhi(float i1, float i2){
	float r = TMath::Abs(i1-i2);
	if (r>TMath::Pi())
	{
		r= 2*TMath::Pi()-r;
	}
	return r;
}

