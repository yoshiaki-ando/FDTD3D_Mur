/*
 * fdtd3d_mur.h
 *
 *  Created on: 2023/05/10
 *      Author: ando
 */

#ifndef FDTD3D_MUR_H_
#define FDTD3D_MUR_H_

#include <cmath>

constexpr double C0 { 3.0e8 };
constexpr double MU0 { 4.0*M_PI };
constexpr double EPS0 { 1.0/C0/C0/MU0 };

constexpr int NX { 200 };
constexpr int NY { 200 };
constexpr int NZ { 200 };

constexpr double DX { 1.0e0 };
constexpr double DY { 1.0e0 };
constexpr double DZ { 1.0e0 };

constexpr double XI { 0.99e0 };
constexpr double DT { XI / std::sqrt( 1.0/DX/DX + 1.0/DY/DY + 1.0/DZ/DZ ) / C0 };

constexpr double Tmax { 1000.0 * DT };

constexpr int NT { int(Tmax / DT) };

/*** プロトタイプ ***/

void update_ex(double ***Ex, double ***Hy, double ***Hz);
void update_ey(double ***Ey, double ***Hz, double ***Hx);
void update_ez(double ***Ez, double ***Hx, double ***Hy);
void update_hx(double ***Hx, double ***Ey, double ***Ez);
void update_hy(double ***Hy, double ***Ez, double ***Ex);
void update_hz(double ***Hz, double ***Ex, double ***Ey);

void mur_ex(double ***Ex, double ***Exy, double ***Exz);
void mur_ey(double ***Ey, double ***Eyz, double ***Eyx);
void mur_ez(double ***Ez, double ***Ezx, double ***Ezy);
void store_ex(double ***Ex, double ***Exy, double ***Exz);
void store_ey(double ***Ey, double ***Eyz, double ***Eyx);
void store_ez(double ***Ez, double ***Ezx, double ***Ezy);

double Idl(double t); /* 電流モーメント */

#endif /* FDTD3D_MUR_H_ */
