#include <iostream>
#include <mpi/mpi.h>
#include <omp.h>
#include <chrono>
using namespace std;

/*
  An example of how to use nec2++ from within a C++ program.
  
  See the makefile target 'test_cpp' for how to compile this
  (as well as the nec2++) source code.

  Copyright (C) 2005,2015  Timothy C.A. Molteno (tim@molteno.net)
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
// #include "../src/c_geometry.h"
#include "c_geometry.h"
#include "nec_context.h"
#include "nec_exception.h"
#include "nec_radiation_pattern.h"

/*
atoi(string) --> convert string to integer
 it's better to use stro(primitive datatpe: i, l, d, f) --> convert string to i, l, d, or f 

char str[100] = "133";
i.e long int x = strtlong(str, NULL, 10);
*/
int main(int argc, char **argv)
{

  MPI_Init(&argc, &argv);
  using std::chrono::duration;
  using std::chrono::duration_cast;
  using std::chrono::high_resolution_clock;
  using std::chrono::milliseconds;

  auto t1 = high_resolution_clock::now();

  // MPI_CXX_DOUBLE_COMPLEX

  try
  {
    cout << "Nec2++ C++ example. Running (takes a few minutes...)" << endl;

    nec_context nec;
    nec.initialize();

    c_geometry *geo = nec.get_geometry();
    geo->wire(0, 512, -0.048, 0.021, -0.005, 0.035, 0.043, 0.014, 0.001, 1.0, 1.0);  // 70
    geo->wire(0, 512, 0.017, -0.015, 0.014, -0.027, 0.04, -0.031, 0.001, 1.0, 1.0);  //66
    geo->wire(0, 512, 0.046, -0.01, 0.028, -0.013, -0.005, 0.031, 0.001, 1.0, 1.0);  //47
    geo->wire(0, 512, -0.048, -0.038, -0.04, 0.049, -0.045, -0.04, 0.001, 1.0, 1.0); //77
    nec.geometry_complete(0);

    nec.gn_card(-1, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    nec.ld_card(5, 0, 0, 0, 3.72e7, 0.0, 0.0);
    nec.pt_card(-1, 0, 0, 0);
    nec.ex_card(EXCITATION_LINEAR, 1, 1, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    nec.fr_card(0, 2, 2400.0, 100.0);
    nec.rp_card(0, 10, 10, 0, 5, 0, 0, 0.0, 0.0, 9.0, 9.0, 0.0, 0.0);

    // now get the radiation pattern data. The result index is 0 since
    // this is the first (and only) radiation pattern.

    // nec_radiation_pattern* rp = nec.get_radiation_pattern(0);
    // int nth = rp->get_ntheta();
    // int nph = rp->get_nphi();

    // cout << endl << "Theta \tPhi \tHorizontal \tVertical \tTotal" << endl;
    // for (int j=0; j<nph; j++) {
    //   for (int i=0; i<nth; i++) {
    //     cout
    //       << rp->get_theta(i) << "  \t"
    //       << rp->get_phi(j) << "  \t"
    //       << rp->get_power_gain_horiz(i,j) << "  \t"
    //       << rp->get_power_gain_vert(i,j) << "  \t"
    //       << rp->get_power_gain(i,j) << "  \t"
    //       << rp->get_etheta_magnitude(i,j) << "  \t"
    //       << rp->get_etheta_phase(i,j) << "  \t"
    //       << rp->get_ephi_magnitude(i,j) << "  \t"
    //       << rp->get_ephi_phase(i,j)
    //       << endl;
    //   }
    // }

  auto t2 = high_resolution_clock::now();

  /* Getting number of milliseconds as an integer. */
  auto ms_int = duration_cast<duration<double>>(t2 - t1);

  /* Getting number of milliseconds as a double. */
  duration<double, std::milli> ms_double = t2 - t1;
  std::cout << ms_int.count() << "ms\n";
  }

  //  MPI_Barrier(MPI_COMM_WORLD);
  
  catch (nec_exception *e)
  {
    cout << e->get_message() << endl;
  }
  MPI_Finalize();
  return 0;
}
