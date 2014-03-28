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
TCLAP::ValueArg<double> Jout("J","Jout", "Outer Ising interactions",false, 1.0,"double",cmd);
TCLAP::ValueArg<double> Jint("","Jint", "Internal Ising Interactions",false, 1.0,"double",cmd);
TCLAP::ValueArg<double> bx("","bx", "Magnetic field in x direction",false, 1.4,"double",cmd);
TCLAP::ValueArg<double> by("","by", "Magnetic field in y direction",false, 0.,"double",cmd);
TCLAP::ValueArg<double> bz("","bz", "Magnetic field in z direction",false, 1.4,"double",cmd);
TCLAP::ValueArg<double> theta("","theta", "polar angle",false, 1.0,"double",cmd);
TCLAP::ValueArg<double> phi("","phi", "azimultal angle",false, 1.0,"double",cmd);
TCLAP::ValueArg<int> steps("","steps","steps",false, 100,"int",cmd);
TCLAP::ValueArg<double> deltabx("","deltabx", "perturbation on bx",false, 0.0,"double",cmd);

int main(int argc, char* argv[])
{

cmd.parse( argc, argv );
cout.precision(12);

vec b(3), bpert(3); 
b(0)=bx.getValue(); 
b(1)=by.getValue();
b(2)=bz.getValue();
bpert=b;
bpert(0)=b(0)+deltabx.getValue();

//Ising parameters

double Jeps=Jout.getValue()-Jint.getValue();

//Construccion de estado coherentre

cvec state, staterev, qustate;

qustate=BlochToQubit(theta.getValue(),phi.getValue());

state=TensorPow(qustate,qubits.getValue());

staterev=state;

//Lista de la fidelidad

vec list(steps.getValue());

for(int i=0;i<steps.getValue();i++){

list(i)=pow( abs( dot( conj(staterev),state)),2);

//cout<< pow( abs( dot( conj(staterev),state)),2) <<endl;

cout << list(i) <<endl;
// cout<< i<< " " << list(i) <<endl;

list(i)=sqrt(list(i));

apply_ising_allvsall(state, Jint.getValue());
//Symmetry Breaking

apply_chain(state, Jeps, b);

//reversing

apply_ising_allvsall(staterev, Jint.getValue());
//Symmetry Breaking

apply_chain(state, Jeps, bpert);

}
 
//fidelity.close();

//cout << staterev;

cout<< sum_positive_derivatives(list)<< endl;


}
