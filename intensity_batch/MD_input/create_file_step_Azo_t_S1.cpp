#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int main() {
	ofstream fout;

	string molecule_name = "Azo_t_S1";
	string input_name = "Azo_t_S1.in";

	string file_name[16] = {"pre_therm", "therm1", "therm2", "therm3", "therm4", "50ps", "100ps", "200ps", "300ps", "400ps", "500ps", "600ps", "700ps", "800ps", "900ps", "1ns"};

	for(int i=0; i<16; i++) {
		fout.open((file_name[i]+"_"+molecule_name).c_str());

		fout << "# Input system files" << endl
			<< "sys-spec-file=" << input_name << endl;

		if(i!=0) {
			fout << "restart-file=" << (file_name[i-1]+"_"+molecule_name) << ".mds" << endl;
		}

		fout << endl;

		fout << "# Characteristics" << endl
			<< "density=" << "0.779" << endl
			<< "temperature=" << "300" << endl
			<< endl;

		fout << "# timesteps" << endl;
		if(i==0) {
			// pre_therm : 0.0000001*5000=0.0005ps
			fout << "nsteps=" << "50000" << endl
				<< "step=" << "0.00000001" << endl;
		}
		else if(i>=1 && i<=2) {
			// therm1~2 : 0.000001*10000=0.01ps
			fout << "nsteps=" << "10000" << endl
				<< "step=" << "0.000001" << endl;
		}
		else if(i>=3 && i<=4) {
			// therm3~4 : 0.0005*20000=10ps
			fout << "nsteps=" << "20000" << endl
				<< "step=" << "0.0005" << endl;
		}
		else if(i>=5 && i<=6) {
			// 50ps, 100ps : 0.0005*100000=50ps
			fout << "const-temp=1" << endl;
			fout << "nsteps=" << "100000" << endl
				<< "step=" << "0.0005" << endl;
		}
		else {
			// 200ps~1ns : 0.0005*200000=100ps
			fout << "const-temp=1" << endl;
			fout << "nsteps=" << "200000" << endl
				<< "step=" << "0.0005" << endl;
		}
		fout << endl;

		if(i<=4) {
			// pre_therm, therm1~4 : scaling
			fout << "# scaling" << endl
			<< "scale-options=" << "6" << endl
			<< "scale-interval=" << "100" << endl
			<< "scale-end=" << "6000" << endl
			<< endl;
		}
	
		fout << "# when to calculate averages and print" << endl
			<< "print-interval=" << "100" << endl
			<< "roll-interval=" << "100" << endl
			<< "average-interval=" << "100" << endl
			<< "begin-average=" << "0" << endl
			<< endl;

		fout << "# values" << endl
			<< "#alpha=" << "-1" << endl
			<< "subcell=" << "2.5" << endl
			<< "strict-cutoff=" << "1" << endl
			<< "cutoff=" << "11.86" << endl
			<< "surface-dipole=" << "0" << endl
			<< endl;
		
		fout << "# rdf" << endl
			<< "begin-rdf=" << "0" << endl
			<< "rdf-limit=" << "20" << endl
			<< "rdf-out=" << "400" << endl
			<< "nbins=" << "800" << endl
			<< endl;

		fout << "# dump files" << endl
			<< "dump-file=" << (file_name[i]+"_"+molecule_name) << ".dump.\%d" << endl
			<< "dump-level=" << "1" << endl
			<< "dump-interval=" << "20" << endl
			<< "ndumps=" << "1000" << endl
			<< endl;

		fout << "save-file=" << (file_name[i]+"_"+molecule_name) << ".mds";

		fout.close();
	}

	// Saves execution record on ah_debug
	fout.open("execution_script");
	for(int i=0; i<16; i++) {
		fout << "date >> ah_debug" << endl
			<< "echo Start : " << file_name[i] << " on `hostname` >> ah_debug" << endl
			<< "moldy < " << (file_name[i]+"_"+molecule_name) << " > " << (file_name[i]+"_"+molecule_name+".out") << endl
			<< "chmod 777 " << (file_name[i]+"_"+molecule_name+".out") << endl
			<< "echo Finish : " << file_name[i] << " >> ah_debug" << endl
			<< "date >> ah_debug" << endl
			<< "echo =================================== >> ah_debug" << endl;
	}
	fout.close();
	
	return 0;
}
