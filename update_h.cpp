/*
 * update_h.cpp
 *
 *  Created on: 2023/05/11
 *      Author: ando
 */
#include "fdtd3d_mur.h"

void update_hx(double ***Hx, double ***Ey, double ***Ez){

  for(int i = 1; i < NX; i++){
    for(int j = 0; j < NY; j++){
      for(int k = 0; k < NZ; k++){
        Hx[i][j][k] = Hx[i][j][k] - DT/MU0/DY * ( Ez[i][j+1][k] - Ez[i][j][k] )
            + DT/MU0/DZ * ( Ey[i][j][k+1] - Ey[i][j][k] );
      }
    }
  }
}

void update_hy(double ***Hy, double ***Ez, double ***Ex){

  for(int i = 0; i < NX; i++){
    for(int j = 1; j < NY; j++){
      for(int k = 0; k < NZ; k++){
        Hy[i][j][k] = Hy[i][j][k] - DT/MU0/DZ * ( Ex[i][j][k+1] - Ex[i][j][k] )
            + DT/MU0/DX * ( Ez[i+1][j][k] - Ez[i][j][k] );
      }
    }
  }
}

void update_hz(double ***Hz, double ***Ex, double ***Ey){

  for(int i = 0; i < NX; i++){
    for(int j = 0; j < NY; j++){
      for(int k = 1; k < NZ; k++){
        Hz[i][j][k] = Hz[i][j][k] - DT/MU0/DX * ( Ey[i+1][j][k] - Ey[i][j][k] )
            + DT/MU0/DY * ( Ex[i][j+1][k] - Ex[i][j][k] );
      }
    }
  }
}
