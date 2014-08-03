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
TCLAP::ValueArg<int> address("A","address","Direccion del qubit para el Ising term",false, 0,"int",cmd);


cmat Isingterm(int i,int N){

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

int main(int argc, char* argv[])
{

cmd.parse( argc, argv );
cout.precision(12);

//cout<< Isingterm(address.getValue(),qubits.getValue()) <<endl;

cmat lista[10];

lista[1]=sigma(1);
lista[2]=sigma(1);
lista[3]=sigma(1);
lista[4]=sigma(1);
lista[5]=sigma(1);

for(int h=1;h<6,h++) TensorProduct(lista[i],lista[i+1])

cout<< lista[1] <<endl;
}


