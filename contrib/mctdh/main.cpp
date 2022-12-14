#include "Core/stdafx.h"
#include "mctdh.h"

int main(int argc, char* argv[])
{

	cout << "=================================================\n";
	cout << "=====     mctdh++     ===========================\n";
	cout << "=================================================\n";

	if (argc != 2){
		cerr << "Please provide a mctdh input as argument." << endl;
		exit(1);
	}

	string filename(argv[1]);
	ifstream ifs {filename};

	if (ifs.fail()){
		cerr << "Error opening the input file. Please make sure the Spelling\n"
		     << "is correct." << endl;
		exit(1);
	}

	mctdh_state state = parser::run(filename);

	return 0;
}

