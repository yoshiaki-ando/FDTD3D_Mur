/*
 * mur.cpp
 *
 *  Created on: 2023/05/11
 *      Author: ando
 */
#include "fdtd3d_mur.h"

constexpr double CX { (C0*DT - DX) / (C0*DT + DX) };
constexpr double CY { (C0*DT - DY) / (C0*DT + DY) };
constexpr double CZ { (C0*DT - DZ) / (C0*DT + DZ) };

void mur_ex(double ***Ex, double ***Exy, double ***Exz){

  for(int i = 0; i < NX; i++){
    for(int k = 1; k < NZ; k++){
      Ex[i][0][k]  = Exy[i][0][k] + CY * ( Ex[i][1][k]    - Ex[i][0][k]  );
      Ex[i][NY][k] = Exy[i][1][k] + CY * ( Ex[i][NY-1][k] - Ex[i][NY][k] );
    }
  }

  for(int i = 0; i < NX; i++){
    for(int j = 1; j < NY; j++){
      Ex[i][j][0]  = Exz[i][j][0] + CZ * ( Ex[i][j][1]    - Ex[i][j][0] );
      Ex[i][j][NZ] = Exz[i][j][1] + CZ * ( Ex[i][j][NZ-1] - Ex[i][j][NZ] );
    }
  }
}

void mur_ey(double ***Ey, double ***Eyz, double ***Eyx){

  for(int i = 1; i < NX; i++){
    for(int j = 0; j < NY; j++){
      Ey[i][j][0]  = Eyz[i][j][0] + CZ * ( Ey[i][j][1]    - Ey[i][j][0] );
      Ey[i][j][NZ] = Eyz[i][j][1] + CZ * ( Ey[i][j][NZ-1] - Ey[i][j][NZ] );
    }
  }

  for(int j = 0; j < NY; j++){
    for(int k = 1; k < NZ; k++){
      Ey[0][j][k]  = Eyx[0][j][k] + CX * ( Ey[1][j][k]    - Ey[0][j][k] );
      Ey[NX][j][k] = Eyx[1][j][k] + CX * ( Ey[NX-1][j][k] - Ey[NX][j][k] );
    }
  }
}

void mur_ez(double ***Ez, double ***Ezx, double ***Ezy){

  for(int j = 1; j < NY; j++){
    for(int k = 0; k < NZ; k++){
      Ez[0][j][k]  = Ezx[0][j][k] + CX * ( Ez[1][j][k]    - Ez[0][j][k] );
      Ez[NX][j][k] = Ezx[1][j][k] + CX * ( Ez[NX-1][j][k] - Ez[NX][j][k] );
    }
  }

  for(int i = 1; i < NX; i++){
    for(int k = 0; k < NZ; k++){
      Ez[i][0][k]  = Ezy[i][0][k] + CY * ( Ez[i][1][k]    - Ez[i][0][k] );
      Ez[i][NY][k] = Ezy[i][1][k] + CY * ( Ez[i][NY-1][k] - Ez[i][NY][k] );
    }
  }
}

void store_ex(double ***Ex, double ***Exy, double ***Exz){

  for(int i = 0; i < NX; i++){
    for(int k = 1; k < NZ; k++){
      Exy[i][0][k] = Ex[i][1][k];
      Exy[i][1][k] = Ex[i][NY-1][k];
    }
  }

  for(int i = 0; i < NX; i++){
    for(int j = 1; j < NY; j++){
      Exz[i][j][0] = Ex[i][j][1];
      Exz[i][j][1] = Ex[i][j][NZ-1];
    }
  }
}

void store_ey(double ***Ey, double ***Eyz, double ***Eyx){

  for(int i = 1; i < NX; i++){
    for(int j = 0; j < NY; j++){
      Eyz[i][j][0] = Ey[i][j][1];
      Eyz[i][j][1] = Ey[i][j][NZ-1];
    }
  }

  for(int j = 0; j < NY; j++){
    for(int k = 1; k < NZ; k++){
      Eyx[0][j][k] = Ey[1][j][k];
      Eyx[1][j][k] = Ey[NX-1][j][k];
    }
  }
}

void store_ez(double ***Ez, double ***Ezx, double ***Ezy){

  for(int j = 1; j < NY; j++){
    for(int k = 0; k < NZ; k++){
      Ezx[0][j][k] = Ez[1][j][k];
      Ezx[1][j][k] = Ez[NX-1][j][k];
    }
  }

  for(int i = 1; i < NX; i++){
    for(int k = 0; k < NZ; k++){
      Ezy[i][0][k] = Ez[i][1][k];
      Ezy[i][1][k] = Ez[i][NY-1][k];
    }
  }
}



