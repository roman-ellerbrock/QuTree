#pragma once
#include "Vector.h"
#include "Matrix.h"

class MNISTdata
{
public:
	MNISTdata(int _nsets);
	MNISTdata(int _nsets, int blocksize, int x=0, int y=0);
	~MNISTdata() = default;

	const Vectord& GetPicture(size_t i)const {
		assert(i < pictures.size());
		return pictures[i];
	}

	void print(size_t p)const;
	void print()const;
	void read(const string& filename);
	void readlabels(string filename);
	size_t size()const { return pictures.size(); }

	vector<Vectord>::const_iterator begin()const { return pictures.begin(); }
	vector<Vectord>::const_iterator end()const { return pictures.end(); }

private:
	void readheader(ifstream & file);
	void readlabelheader(ifstream & file);
	Vectord readpicture(ifstream& file);
	Vectord readlabel(ifstream & file);

	int nsets;
	int nin;
	int nout;
	int blocksize;
	int x, y;

	vector <Vectord> pictures;
	vector <Vectord> labels;
	uint32_t magic;
	uint32_t labelmagic;
	uint32_t nimages, nlabels;
	uint32_t nrows;
	uint32_t ncols;
};

typedef Vectord Picture;


