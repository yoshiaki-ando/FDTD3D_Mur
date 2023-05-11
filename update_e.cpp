/*
 * update_e.cpp
 *
 *  Created on: 2023/05/10
 *      Author: ando
 */
#include "fdtd3d_mur.h"

void update_ex(double ***Ex, double ***Hy, double ***Hz){

  for(int i = 0; i < NX; i++){
    for(int j = 1; j < NY; j++){
      for(int k = 1; k < NZ; k++){
        Ex[i][j][k] = Ex[i][j][k] + DT/EPS0/DY * ( Hz[i][j][k] - Hz[i][j-1][k] )
            - DT/EPS0/DZ * ( Hy[i][j][k] - Hy[i][j][k-1] );
      }
    }
  }
}


void update_ey(double ***Ey, double ***Hz, double ***Hx){

  for(int i = 1; i < NX; i++){
    for(int j = 0; j < NY; j++){
      for(int k = 1; k < NZ; k++){
        Ey[i][j][k] = Ey[i][j][k] + DT/EPS0/DZ * ( Hx[i][j][k] - Hx[i][j][k-1] )
            - DT/EPS0/DX * ( Hz[i][j][k] - Hz[i-1][j][k] );
      }
    }
  }
}

void update_ez(double ***Ez, double ***Hx, double ***Hy){

  for(int i = 1; i < NX; i++){
    for(int j = 1; j < NY; j++){
      for(int k = 0; k < NZ; k++){
        Ez[i][j][k] = Ez[i][j][k] + DT/EPS0/DX * ( Hy[i][j][k] - Hy[i-1][j][k] )
            - DT/EPS0/DY * ( Hx[i][j][k] - Hx[i][j-1][k] );
      }
    }
  }
}
