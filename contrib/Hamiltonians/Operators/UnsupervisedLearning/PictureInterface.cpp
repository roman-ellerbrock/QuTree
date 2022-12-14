//
// Created by Hayley Weir on 16/10/2019.
//

#include <Operators/SumOfProductsOperator.h>
#include "PictureInterface.h"
#include "MNISTdata.h"

namespace PictureInterface {

	SOP DatasetOperator(const MNISTdata& dataset, float threshold, bool adjungate) {

		SOP DatasetOps;

		for (const Picture& Pic : dataset) {
			// Get picture in operator form and add to dataset of operators
			MPO PicOp = PictureOperator(Pic, threshold, adjungate);
			DatasetOps.push_back(PicOp, 1.);
		}
		return DatasetOps;
	}

	MPO PictureOperator(const Picture& pic, float threshold, bool adjungate) {
		MPO PicOp;
		using namespace GateOperators;

		for (size_t i = 0; i < pic.Dim(); ++i) {
			int id = GateOperators::ConvertIdx2D(i, 32, 32);
			if (pic(i) > threshold) {
//				PicOp.push_back(GateOperators::X, GateOperators::ConvertIdx2D(i, 32, 32));
//				PicOp.push_back(GateOperators::SetBlack, GateOperators::ConvertIdx2D(i, 32, 32));
				shared_ptr<BottomLayerSPO> black = make_shared<SetBlack>(0.90, adjungate);
//				shared_ptr<BottomLayerSPO> R = make_shared<Rk>(k, conj_phase);
				PicOp.push_back(black, GateOperators::ConvertIdx2D(i, 32, 32));
			} else {
			}
/*			cout << id << " ";
			if (id < 10) {
				cout << "  ";
			}
			if (id < 100) {
				cout << " ";
			}
			if ((i + 1) % 32 == 0) {
				cout << endl;
			}*/
		}
//		cout << endl << endl;
		return PicOp;
	}
}

