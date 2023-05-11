/*
 * idl.cpp
 *
 *  Created on: 2023/05/11
 *      Author: ando
 */
#include <iostream>
#include <cmath>
#include "fdtd3d_mur.h"

double Idl(double t){
  constexpr double tDx { 2.0*DX/C0 };
  constexpr double sig { 5.0 * tDx };
  constexpr double t0 { 6.0 * sig };

//  std::cout << t0 << " , " << Tmax << std::endl;
//  exit(0);

  return (t - t0) / sig * std::exp( - (t - t0) * (t - t0) / 2.0 / sig / sig );

}
