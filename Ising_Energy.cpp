#include "RFIM.cpp"
#include "RBIM.cpp"

Vec_d Calc_Ising_Energy(Mat_d& New_H_value, Vec_d& New_Ising_Spin, Vec_d& New_Ising_Spin_type)
{
   Vec_d New_Ising_Atomistic_Energy; 
   New_Ising_Atomistic_Energy.resize(total_spins);
   New_Ising_Atomistic_Energy.setZero();

   #pragma omp parallel for schedule(static) shared (Neighbor_List, New_Ising_Atomistic_Energy, New_H_value, New_Ising_Spin, New_Ising_Spin_type)
      for(int spin = 0; spin < total_spins; spin++){
         int type = New_Ising_Spin_type[spin];  

	  for(int n = 0; n < nearest_neighbor_count; n++){
	     if( Neighbor_List(n, spin) > spin ){
 	       int nei_spin = New_Ising_Spin[Neighbor_List(n, spin)];
	       int nei_type = New_Ising_Spin_type[Neighbor_List(n, spin)];
	       double nei_dist = Distance_List(n, spin);

	       double J = J_matrix(type, nei_type) * ( 1.0 - ((nei_dist - 1.0)/grid_cut_off) );
               if(RBIM_)
                  J = RBIM(type, nei_type);
	       double new_energy = New_H_value(spin, Neighbor_List(n, spin)) * New_Ising_Spin[spin]
                                                   - 0.5 * (J * New_Ising_Spin[spin] * nei_spin);

               //Symmetric interaction...
	       #pragma omp atomic
               New_Ising_Atomistic_Energy[spin] += new_energy;
	       #pragma omp atomic
	       New_Ising_Atomistic_Energy[Neighbor_List(n, spin)] += new_energy;
	       }
	    }  
          }

    return New_Ising_Atomistic_Energy; 	 
} 
