/*
 * main.cpp
 *
 *  Created on: 2023/05/10
 *      Author: ando
 */

#include <iostream>
#include <fstream>
#include <string>
#include <memory_allocate.h>

#include "fdtd3d_mur.h"

/* 保存するディレクトリを指定する */
std::string work_dir("/home/ando/works/FDTD3D_Mur/data/");

int main(void){

  double ***Ex = AndoLab::allocate_memory3d(NX, NY+1, NZ+1, 0.0);
  double ***Ey = AndoLab::allocate_memory3d(NX+1, NY, NZ+1, 0.0);
  double ***Ez = AndoLab::allocate_memory3d(NX+1, NY+1, NZ, 0.0);

  double ***Exy = AndoLab::allocate_memory3d(NX, 2, NZ+1, 0.0);
  double ***Exz = AndoLab::allocate_memory3d(NX, NY+1, 2, 0.0);
  double ***Eyz = AndoLab::allocate_memory3d(NX+1, NY, 2, 0.0);
  double ***Eyx = AndoLab::allocate_memory3d(2, NY, NZ+1, 0.0);
  double ***Ezx = AndoLab::allocate_memory3d(2, NY+1, NZ, 0.0);
  double ***Ezy = AndoLab::allocate_memory3d(NX+1, 2, NZ, 0.0);

  double ***Hx = AndoLab::allocate_memory3d(NX+1, NY, NZ, 0.0);
  double ***Hy = AndoLab::allocate_memory3d(NX, NY+1, NZ, 0.0);
  double ***Hz = AndoLab::allocate_memory3d(NX, NY, NZ+1, 0.0);

  std::ofstream ofs_obs( work_dir + "at_qNX.dat" );

  for(int n = 1; n <= NT; n++){
    if ( n%5 == 0 ){
      std::cout << n << " / " << NT << std::endl;
    }

    store_ex(Ex, Exy, Exz);
    store_ey(Ey, Eyz, Eyx);
    store_ez(Ez, Ezx, Ezy);

    update_ex(Ex, Hy, Hz);
    update_ey(Ey, Hz, Hx);
    update_ez(Ez, Hx, Hy);

    mur_ex(Ex, Exy, Exz);
    mur_ey(Ey, Eyz, Eyx);
    mur_ez(Ez, Ezx, Ezy);

    double t = (n - 0.5) * DT;
    Ez[NX/2][NY/2][NZ/2] -= DT/EPS0 * Idl(t) / DX/DY/DZ;

//    if ( (n%5 == 0) && (n >= 50) ){
//      std::ofstream ofs( work_dir + std::to_string(n) + ".dat" );
//      for(int i = 0; i <= NX; i++){
//        for(int j = 0; j <= NY; j++){
//          ofs << i << " " << j << " " << Ez[i][j][NZ/2] << "\n";
//        }
//        ofs << "\n";
//      }
//    }

    ofs_obs << n*DT << " " << Ez[NX/2-25][NY/2-25][NZ/2] << std::endl;

    update_hx(Hx, Ey, Ez);
    update_hy(Hy, Ez, Ex);
    update_hz(Hz, Ex, Ey);
  }

  ofs_obs.close();

  AndoLab::deallocate_memory3d(Ex);
  AndoLab::deallocate_memory3d(Ey);
  AndoLab::deallocate_memory3d(Ez);
  AndoLab::deallocate_memory3d(Hx);
  AndoLab::deallocate_memory3d(Hy);
  AndoLab::deallocate_memory3d(Hz);

  AndoLab::deallocate_memory3d(Exy);
  AndoLab::deallocate_memory3d(Exz);
  AndoLab::deallocate_memory3d(Eyz);
  AndoLab::deallocate_memory3d(Eyx);
  AndoLab::deallocate_memory3d(Ezx);
  AndoLab::deallocate_memory3d(Ezy);

  return 0;
}
