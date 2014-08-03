#include <iostream>
#include <cpp/dev_random.cpp>
#include <tclap/CmdLine.h>
#include <itpp/itbase.h>
//#include <itpp/base/algebra/cholesky.h>
#include <itpp/stat/histogram.h>
#include <cpp/itpp_ext_math.cpp>
#include <cpp/spinchain.cpp>
#include <itpp/stat/misc_stat.h>
#include <fstream>
#include <cpp/RMT.cpp>

using namespace std; 
using namespace itpp;
using namespace itppextmath;
using namespace cfpmath;
using namespace spinchain;
using namespace RMT;

TCLAP::CmdLine cmd("No reputas estupidas mames",' ', "0.1");
TCLAP::ValueArg<string> optionArg("o","option", "Option" ,false,"normalito", "string",cmd);
TCLAP::ValueArg<string> optionArg2("","option2", "Option2" ,false,"fidelity", "string",cmd);
TCLAP::ValueArg<unsigned int> seed("s","seed", "Random seed [0 for urandom]",false, 243243,"unsigned int",cmd);
TCLAP::ValueArg<int> qubits("q","qubits", "number of qubits",false, 4,"int",cmd);
TCLAP::ValueArg<double> J("J","ising_coupling", "Ising interaction in the z-direction",false, 1.0,"double",cmd);
TCLAP::ValueArg<double> bx("","bx", "Magnetic field in x direction",false, 1.4,"double",cmd);
TCLAP::ValueArg<double> by("","by", "Magnetic field in y direction",false, 0.,"double",cmd);
TCLAP::ValueArg<double> bz("","bz", "Magnetic field in z direction",false, 1.4,"double",cmd);
TCLAP::ValueArg<double> theta("","theta", "polar angle",false, 1.0,"double",cmd);
TCLAP::ValueArg<double> phi("","phi", "azimultal angle",false, 1.0,"double",cmd);
TCLAP::ValueArg<double> deltabx("","deltabx", "perturbation",false, 0.1,"double",cmd);
TCLAP::ValueArg<int> steps("","steps","steps",false, 100,"int",cmd);
TCLAP::ValueArg<double> Jpert("","Jpert","Perturbation on Ising",false, 0.0,"double",cmd);
TCLAP::ValueArg<int> i("i","addressi","Direccion del qubit para el Ising term i",false, 1,"int",cmd);
TCLAP::ValueArg<int> j("j","addressj","Direccion del qubit para el Ising term j",false, 2,"int",cmd);


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

cmat HamiltonianChainU(double J, cvec b(3), int qubits){
	
	cmat HI=zeros_c(pow(2,qubits),pow(2,qubits));
	
	HK=HI;
	
	for(int i=0; i<qubits-1; i++){
		
		HI=Isingterm(i,i+1,qubits)+HI;
		}
		
		HI=J*HI+J*Isingterm(qubits-1,0,qubits);
		
	for(int i=0;i<qubits;i++){
		HK=b(1)*sigmaddress(1,i,qubits)+b(3)*sigmaddress(2,i,qubits)+HK;
	}
	
	
	
}

int main(int argc, char* argv[])
{

cmd.parse( argc, argv );
cout.precision(12);


//cout<< Isingterm(i.getValue(),j.getValue(),qubits.getValue()) <<endl;




}


