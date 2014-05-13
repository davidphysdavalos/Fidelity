#include <iostream>
#include <cpp/dev_random.cpp>
#include <tclap/CmdLine.h>
#include <itpp/itbase.h>
#include <itpp/stat/histogram.h>
#include <cpp/itpp_ext_math.cpp>
#include <cpp/spinchain.cpp>
#include <itpp/stat/misc_stat.h>
#include <fstream>

using namespace std; 
using namespace itpp;
using namespace itppextmath;
using namespace cfpmath;
using namespace spinchain;


TCLAP::CmdLine cmd("Command description message", ' ', "0.1");
TCLAP::ValueArg<int> qubits("q","qubits", "number of qubits",false, 4,"int",cmd);
TCLAP::ValueArg<double> J("J","ising_coupling", "Ising interaction in the z-direction",false, 1.0,"double",cmd);
TCLAP::ValueArg<double> bx("","bx", "Magnetic field in x direction",false, 1.4,"double",cmd);
TCLAP::ValueArg<double> by("","by", "Magnetic field in y direction",false, 0.,"double",cmd);
TCLAP::ValueArg<double> bz("","bz", "Magnetic field in z direction",false, 1.4,"double",cmd);
TCLAP::ValueArg<double> theta("","theta", "polar angle",false, 1.0,"double",cmd);
TCLAP::ValueArg<double> phi("","phi", "azimultal angle",false, 1.0,"double",cmd);
TCLAP::ValueArg<int> steps("","steps","steps",false, 100,"int",cmd);

cvec TensorPow(cvec state, int qub){

cvec newstate;

newstate=state;

for(int i=0;i<qub-1;i++){
newstate=TensorProduct(newstate,state);
}
return newstate;

}


int main(int argc, char* argv[])
{

cmd.parse( argc, argv );
cout.precision(12);

vec b(3); 
b(0)=bx.getValue(); 
b(1)=by.getValue();
b(2)=bz.getValue();

cvec state, staterev, qustate;

cmat ro;

//ofstream fidelity;
//fidelity.open("fidelity.dat");

//qustate=RandomState(64);

//int dim=pow_2(qubits.getValue());

qustate=BlochToQubit(theta.getValue(),phi.getValue());

//qustate=RandomState(2);

//for(int i=0; i<qubits.getValue()+1;i++){

//list(i)=qustate;

//}

state=TensorPow(qustate,qubits.getValue());

//cout<< qustate ;

complex<double> sigmax, sigmay, sigmaz;

int flag;

for(int i=0;i<steps.getValue();i++){

apply_chain(state, J.getValue(), b); 

sigmax=0;

sigmay=0;

sigmaz=0;

for(int j=0; j<qubits.getValue(); j++){

flag=pow(2,j);

ro=partial_trace_qubits(state, flag);

sigmax=trace(sigma(1)*ro)+sigmax;

sigmay=trace(sigma(2)*ro)+sigmay;

sigmaz=trace(sigma(3)*ro)+sigmaz;

}
cout<< real(sigmax) <<" "<< real(sigmay) <<" "<< real(sigmaz) <<endl; 
}


}
