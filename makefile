INCLUDE = ~/libs

helloworld :: helloworld.cpp makefile
	g++ helloworld.cpp -o helloworld -I $(INCLUDE) -litpp

hello :: hello.cpp makefile
	g++ hello.cpp -o hello -I $(INCLUDE) -litpp

Fidelity :: Fidelity.cpp makefile
	g++ Fidelity.cpp -o Fidelity -I $(INCLUDE) -litpp

FidelityKlimov :: FidelityKlimov.cpp makefile
	g++ FidelityKlimov.cpp -o FidelityKlimov -I $(INCLUDE) -litpp

Averages :: Averages.cpp makefile
	g++ Averages.cpp -o Averages -I $(INCLUDE) -litpp

Dephasing :: Dephasing.cpp makefile
	g++ Dephasing.cpp -o Dephasing -I $(INCLUDE) -litpp

StateEvolving :: StateEvolving.cpp makefile
	g++ StateEvolving.cpp -o StateEvolving -I $(INCLUDE) -litpp

FidelityInhom :: FidelityInhom.cpp makefile
	g++ FidelityInhom.cpp -o FidelityInhom -I $(INCLUDE) -litpp

FidelityBxBroken :: FidelityBxBroken.cpp makefile
	g++ FidelityBxBroken.cpp -o FidelityBxBroken -I $(INCLUDE) -litpp

FidelityBxInhom :: FidelityBxInhom.cpp makefile
	g++ FidelityBxInhom.cpp -o FidelityBxInhom -I $(INCLUDE) -litpp

FidelityAllvsAll :: FidelityAllvsAll.cpp makefile
	g++ FidelityAllvsAll.cpp -o FidelityAllvsAll -I $(INCLUDE) -litpp 
	
FidelityAllvsAll_Special :: FidelityAllvsAll_Special.cpp makefile
	g++ FidelityAllvsAll_Special.cpp -o FidelityAllvsAll_Special -I $(INCLUDE) -litpp
	
FidelityAllvsAll_Special_second :: FidelityAllvsAll_Special_second.cpp makefile
	g++ FidelityAllvsAll_Special_second.cpp -o FidelityAllvsAll_Special_second -I $(INCLUDE) -litpp
	
FidelityAllvsAll_Special_BrokenPermutations :: FidelityAllvsAll_Special_BrokenPermutations.cpp makefile
	g++ FidelityAllvsAll_Special_BrokenPermutations.cpp -o FidelityAllvsAll_Special_BrokenPermutations -I $(INCLUDE) -litpp
