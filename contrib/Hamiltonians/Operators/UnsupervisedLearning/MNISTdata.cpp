#include "stdafx.h"
#include "MNISTdata.h"
//#include <Winsock2.h>
#include <arpa/inet.h>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <stdint.h>
//#include "Targa.h"

MNISTdata::MNISTdata(int _nsets)
	: nsets(_nsets), nin(28 * 28), nout(28 * 28), blocksize(28), x(0), y(0) {
	// Constructor to create an object that reads the whole picture
}

MNISTdata::MNISTdata(int _nsets, int _blocksize, int _x, int _y)
	: nsets(_nsets), nin(28 * 28), nout(28 * 28), blocksize(_blocksize), x(_x), y(_y) {
	// Constructor to create an object that reads a block of the picture
}

void MNISTdata::read(const string& filename) {
	cout << "Reading MNIST-data set..." << endl;

	// open file
	ifstream file;
	file.open(filename.c_str(), ios::in | ios::ate | ios::binary);
	if (!file.is_open()) {
		cout << "Unable to open file." << endl;
		return;
	}

	// get filesize
	size_t end = (size_t) file.tellg();
	file.seekg(0, ios::beg);
	size_t begin = (size_t) file.tellg();
	size_t filesize;
	filesize = end - begin;
	cout << "File-size of training set:\t" << filesize / 1024 << " Kb" << endl;

	// read header of MNIST data file
	readheader(file);

	nsets = min((int) nsets, (int) nimages);
	//  read images
	for (int i = 0; i < nsets; i++) {
		pictures.push_back(readpicture(file));
	}

/*	// write to hard drive as .tga
	for (Vectord picture : pictures)
	{
		string name = "test.tga";
		Targa tga(name, 28, 28, picture);
//		tga.getimage(picture);
		tga.write(name);
	}*/

	cout << "... done reading." << endl;
}

void MNISTdata::readlabels(string filename) {
	cout << "Reading MNIST-data labels..." << endl;

	// open file
	ifstream file;
	file.open(filename.c_str(), ios::in | ios::ate | ios::binary);
	if (!file.is_open()) {
		cout << "Unable to open file." << endl;
		return;
	}

	// get filesize
	size_t end = (size_t) file.tellg();
	file.seekg(0, ios::beg);
	size_t begin = (size_t) file.tellg();
	size_t filesize;
	filesize = end - begin;
	cout << "File-size of label set:\t" << filesize / 1024 << " Kb" << endl;

	// read header of MNIST data file
	readlabelheader(file);

	nsets = min((int) nsets, (int) nlabels);
	//  read labels
	for (int i = 0; i < nsets; i++) {
		labels.push_back(readlabel(file));
	}
	cout << "... done reading." << endl;
}

Vectord MNISTdata::readlabel(ifstream& file) {

	// magic reading
	unsigned char charlabel;
	file.read((char *) &charlabel, sizeof(char));
	uint32_t ilabel = (uint32_t) charlabel;

//	cout << "Label: " << ilabel << endl;
	// just to make sure..
	ilabel = max((uint32_t) 0, ilabel);
	ilabel = min((uint32_t) 9, ilabel);

	Vectord label(10);

	for (int i = 0; i < 10; i++) {
		//label(i) = 0.3;
		//label(i) = -1.0;
		//label(i) = 0.001;
		label(i) = 0;
	}
	//label(ilabel) = 0.7;
	//label(ilabel) = 0.999;
	label(ilabel) = 1.0;

	return label;
}

Vectord MNISTdata::readpicture(ifstream& file) {
	int dim = 28;

	Matrixd image(dim, dim);

	// read an image
	for (unsigned int j = 0; j < ncols; j++) {
		for (unsigned int i = 0; i < nrows; i++) {
			// magic reading
			unsigned char charpixel;
			file.read((char *) &charpixel, sizeof(char));
			uint32_t pixel = (uint32_t) charpixel;
			image(j, i) = (double) pixel / 255;
			//image(j, i) = 2*((double) pixel/255)-1;
		}
	}

	// read block from image
	int xstart = x * blocksize;
	int ystart = y * blocksize;
	int xend = (x + 1) * blocksize;
	int yend = (y + 1) * blocksize;
	Vectord picture(blocksize * blocksize);
	if (xend > (int) ncols) {
		cout << "ERROR." << endl;
		return picture;
	}
	if (yend > (int) nrows) {
		cout << "ERROR." << endl;
		return picture;
	}

	int k = 0;
//	for (int j = ystart; j < yend; ++j) {
//		for (int i = xstart; i < xend; i++) {
	for (int j = yend - 1; j >= ystart; j--) {
		for (int i = xend - 1; i >= xstart; i--) {
			picture(k) = image(j, i);
			k++;
		}
	}
	return picture;
}

void MNISTdata::readheader(ifstream& file) {

	uint32_t magic2;
	uint32_t nimages2;
	uint32_t nrows2;
	uint32_t ncols2;
	file.read((char *) &magic2, sizeof(uint32_t));
	file.read((char *) &nimages2, sizeof(uint32_t));
	file.read((char *) &nrows2, sizeof(uint32_t));
	file.read((char *) &ncols2, sizeof(uint32_t));

	// Machine: Little-Endian, File: Big-Endian --> Convert
	magic = htonl(magic2);
	if (magic != 2051) {
		cout << "Magic number for imagesis wrong." << endl;
	}
	nimages = htonl(nimages2);
	nrows = htonl(nrows2);
	ncols = htonl(ncols2);
	cout << "nimages=" << nimages << endl;
}

void MNISTdata::readlabelheader(ifstream& file) {

	uint32_t magic2;
	uint32_t nimages2;
	file.read((char *) &magic2, sizeof(uint32_t));
	file.read((char *) &nimages2, sizeof(uint32_t));

	// Machine: Little-Endian, File: Big-Endian --> Convert
	labelmagic = htonl(magic2);
	if (labelmagic != 2049) {
		cout << "Magic number for labels is wrong!" << endl;
	}
	nlabels = htonl(nimages2);
	if (nlabels != nimages) {
		cout << "Number of labels doesnt fit number of images." << endl;
	}
}

void MNISTdata::print()const {
	for (size_t i = 0; i < size(); ++i) {
		print(i);
		cout << endl;
	}
}

void MNISTdata::print(size_t p)const {
	const Vectord& pic = GetPicture(p);
	for (size_t i = 0; i < pic.Dim(); ++i) {
		if (pic(pic.Dim()-1-i) > 0) {
			cout << "1";
		} else {
			cout << "0";
		}
		if ((i + 1) % 28 == 0) {
			cout << endl;
		}
	}
}

