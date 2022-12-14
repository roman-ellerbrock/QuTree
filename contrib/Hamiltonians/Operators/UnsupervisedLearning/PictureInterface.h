//
// Created by Hayley Weir on 16/10/2019.
//

#ifndef MCTDH_PICTUREINTERFACE_H
#define MCTDH_PICTUREINTERFACE_H

#include "MNISTdata.h"
#include "GateOperators.h"

namespace PictureInterface {
    SOP DatasetOperator(const MNISTdata& dataset, float threshold, bool adjungate);
    MPO PictureOperator(const Picture& pic, float threshold, bool adjungate);
}

#endif //MCTDH_PICTUREINTERFACE_H
