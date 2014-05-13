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
TCLAP::ValueArg<unsigned int> seed("s","seed", "Random seed [0 for urandom]",false, 243243,"unsigned int",cmd);
TCLAP::ValueArg<int> qubits("q","qubits", "number of qubits",false, 4,"int",cmd);
TCLAP::ValueArg<double> J("J","ising_coupling", "Ising interaction in the z-direction",false, 1.0,"double",cmd);
TCLAP::ValueArg<double> bx("","bx", "Magnetic field in x direction",false, 1.4,"double",cmd);
TCLAP::ValueArg<double> by("","by", "Magnetic field in y direction",false, 0.,"double",cmd);
TCLAP::ValueArg<double> bz("","bz", "Magnetic field in z direction",false, 1.4,"double",cmd);
TCLAP::ValueArg<double> theta("","theta", "polar angle",false, 1.0,"double",cmd);
TCLAP::ValueArg<double> phi("","phi", "azimultal angle",false, 1.0,"double",cmd);
TCLAP::ValueArg<int> steps("","steps","steps",false, 100,"int",cmd);

int main(int argc, char* argv[])
{

cmd.parse( argc, argv );
cout.precision(12);

vec b(3); 
b(0)=bx.getValue(); 
b(1)=by.getValue();
b(2)=bz.getValue();

cvec state, qustate;

qustate=BlochToQubit(theta.getValue(),phi.getValue());

state=TensorPow(qustate,qubits.getValue());

for(int i=0;i<steps.getValue();i++){

apply_chain(state, J.getValue(), b);


}

for(int i=0; i<state.size();i++)
	cout << real(state(i)) << imag(state(i)) <<endl;

}
