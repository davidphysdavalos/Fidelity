#include <iostream>
#include "cpp/dev_random.cpp"
#include <tclap/CmdLine.h>
#include <itpp/itbase.h>
//#include <itpp/base/algebra/cholesky.h>
#include <itpp/stat/histogram.h>
#include "cpp/itpp_ext_math.cpp"
#include "cpp/spinchain.cpp"
#include <itpp/stat/misc_stat.h>
#include <fstream>
#include "cpp/RMT.cpp"

using namespace std; 
using namespace itpp;
using namespace itppextmath;
using namespace cfpmath;
using namespace spinchain;
using namespace RMT;

cmat Isingterm2(int i,int N){

//cmat term=zeros_c(pow(2,N),pow(2,N));

cmat term;

    if(N < 2+i)
        printf("indice de termino muy grande, procure que sea i<=N-2");

    if(N>2+i){
        cmat C=TensorPow(eye_c(2),N-2-i);
        term=TensorProduct(TensorPow(sigma(3), 2), C);
           }
        
    if(N==2+i)
        term=TensorPow(sigma(3),2);
    
	if (i>0){
        cmat term2=TensorProduct(TensorPow(eye_c(2),i),term);
        return term2;
	}
	else{
		return term;
		}
        
}

cmat Isingterm(int i,int j,int N){
	
cmat lista[N];

for(int h=0;h < N;h++){
	
	if(h==i || h==j){
	 lista[h]=sigma(3);
 }
 else{
	 lista[h]=eye_c(2);
	 }
 
 }

cmat termino=kron(lista[0],lista[1]);

for(int h=2;h<N;h++) termino=kron(termino,lista[h]);

return termino;

}

cmat sigmaddress(int i,int qubit,int N){
	
	cmat lista[N];
	
	for(int h=0;h < N;h++){
	
	if(h==qubit){
	 lista[h]=sigma(i);
 }
 else{
	 lista[h]=eye_c(2);
	 }
	 
 }
	 
cmat termino=kron(lista[0],lista[1]);

for(int h=2;h<N;h++) termino=kron(termino,lista[h]);

return termino;

}

cmat HamiltonianChainU(double J, vec b, int Nqubits, double delta){
	
	cmat HI=zeros_c(pow(2,Nqubits),pow(2,Nqubits));
	
	cmat HK=HI;
	
	for(int i=0; i<Nqubits-1; i++){
		
		HI=Isingterm(i,i+1,Nqubits)+HI;
		}
		
		HI=J*HI+J*Isingterm(Nqubits-1,0,Nqubits);
		
	for(int i=0;i<Nqubits;i++){
		HK=b(0)*sigmaddress(1,i,Nqubits)+b(2)*sigmaddress(3,i,Nqubits)+HK;
	}
	
	return exponentiate_nonsym(-complex <double>(0,1)*(HK+delta*RandomGUE(pow(2,Nqubits))))*exponentiate_nonsym(-complex <double>(0,1)*HI);
	
}

cmat HamiltonianChainUbothSites(double J, vec b, int Nqubits, cmat V){
	
	cmat HI=zeros_c(pow(2,Nqubits),pow(2,Nqubits));
	
	cmat HK=HI;
	
	for(int i=0; i<Nqubits-1; i++){
		
		HI=Isingterm(i,i+1,Nqubits)+HI;
		}
		
		HI=J*HI+J*Isingterm(Nqubits-1,0,Nqubits);
		
	for(int i=0;i<Nqubits;i++){
		HK=b(0)*sigmaddress(1,i,Nqubits)+b(2)*sigmaddress(3,i,Nqubits)+HK;
	}
	
	return exponentiate_nonsym(-complex <double>(0,1)*HK)*exponentiate_nonsym(-complex <double>(0,1)*(HI+V));
	
}

int main(int argc, char* argv[])
{


TCLAP::CmdLine cmd("No reputas estupidas mames",' ', "0.1");
TCLAP::ValueArg<string> optionArg("o","option", "Option" ,false,"normalito", "string",cmd);
TCLAP::ValueArg<string> medida("m","medida", "Medida" ,false,"BLP", "string",cmd);
TCLAP::ValueArg<string> optionArg2("","option2", "onesite or bothsites Hamiltonian RMT perturbation" ,false,"onesite", "string",cmd);
TCLAP::ValueArg<unsigned int> seed("s","seed", "Random seed [0 for urandom]",false, 243243,"unsigned int",cmd);
TCLAP::ValueArg<int> qubits("q","qubits", "number of qubits",false, 4,"int",cmd);
TCLAP::ValueArg<double> J("J","ising_coupling", "Ising interaction in the z-direction",false, 1.0,"double",cmd);
TCLAP::ValueArg<double> bx("","bx", "Magnetic field in x direction",false, 1.4,"double",cmd);
TCLAP::ValueArg<double> by("","by", "Magnetic field in y direction",false, 0.,"double",cmd);
TCLAP::ValueArg<double> bz("","bz", "Magnetic field in z direction",false, 1.4,"double",cmd);
TCLAP::ValueArg<double> meshtheta("","meshtheta", "polar angle mesh size",false, 0.1,"double",cmd);
TCLAP::ValueArg<double> meshphi("","meshphi", "azimultal angle mesh size",false, 0.1,"double",cmd);
TCLAP::ValueArg<double> deltapert("","deltapert", "perturbation on RMT",false, 0.1,"double",cmd);
TCLAP::ValueArg<int> steps("","steps","steps",false, 100,"int",cmd);
TCLAP::ValueArg<double> Jpert("","Jpert","Perturbation on Ising",false, 0.0,"double",cmd);
//TCLAP::ValueArg<int> i("i","addressi","Direccion del qubit para el Ising term i",false, 1,"int",cmd);
//TCLAP::ValueArg<int> j("j","addressj","Direccion del qubit para el Ising term j",false, 2,"int",cmd);

cmd.parse( argc, argv );
cout.precision(12);

vec b(3); 
b(0)=bx.getValue(); 
b(1)=by.getValue();
b(2)=bz.getValue();

// {{{ Set seed for random
unsigned int semilla=seed.getValue();
if (semilla == 0){
  Random semilla_uran; semilla=semilla_uran.strong();
} 
RNG_reset(semilla);
// }}}


cvec state, staterev, qustate;

string option2=optionArg2.getValue();

cmat U, Udelta;

if(option2=="onesite"){
	
U=HamiltonianChainU(J.getValue(),b , qubits.getValue(), 0.0);

Udelta=HamiltonianChainU(J.getValue(),b , qubits.getValue(), deltapert.getValue());
}

if(option2=="bothsites"){
	
cmat V=RandomGUE(pow(2, qubits.getValue()));

U=HamiltonianChainUbothSites(J.getValue(),b , qubits.getValue(), deltapert.getValue()*V);

Udelta=HamiltonianChainUbothSites(J.getValue(),b , qubits.getValue(), -1.0*deltapert.getValue()*V);
}

double theta=0.0, phi=0.0;

while(theta<3.14159+meshtheta.getValue()){
	
	phi=0.0;
	
while(phi<6.28319+meshphi.getValue()){

qustate=BlochToQubit(theta,phi);

state=TensorPow(qustate,qubits.getValue());

staterev=state;

for(int i=0;i<steps.getValue()+1;i++){
	
	cout<<pow( abs( dot( conj(staterev),state)),2)<<endl;
	//cout<<theta<<" "<<phi<<endl;
	
	state=U*state;
	//apply_chain(state, J.getValue(), b);
	staterev=Udelta*staterev;
	
	}
	
	phi=phi+meshphi.getValue();
}

	theta=theta+meshtheta.getValue();
}


//cout<<U<<endl;

//cout<< HamiltonianChainU(J.getValue(), b, qubits.getValue()) <<endl;

//cout << HamiltonianChainU(J.getValue(), b, qubits.getValue()) <<endl;

//cout<< Isingterm(i.getValue(),j.getValue(),qubits.getValue())<<endl;

return 0;
}


