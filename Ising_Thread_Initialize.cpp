
void Thread_Initialize()
{
   nearest_neighbor_count = 0;
   for(int spin = 1; spin < total_spins; spin++){
	double distance = 0;
	for(int d = 0; d < dimension; d++){
           int shift = 0;
  	   shift = (int)( 2 * (0 - Ising_Spin_Coordinates(d, spin)) / total_spins_per_dim[d]);
           double x = 0 - (Ising_Spin_Coordinates(d, spin) + shift * total_spins_per_dim[d]);
           distance += pow(x, 2);
	   }
	distance = sqrt(distance);

	if(distance <= grid_cut_off){
	   nearest_neighbor_count++;
	  }
   }

   Neighbor_List.resize(nearest_neighbor_count, total_spins);//Isotropic grid
   Distance_List.resize(nearest_neighbor_count, total_spins);

   for(int spin_i = 0; spin_i < total_spins; spin_i++){
      int nei_count = -1;
      for(int spin_j = 0; spin_j < total_spins; spin_j++){
         if( spin_i != spin_j ){
	   double distance = 0;
	   for(int d = 0; d < dimension; d++){
              int shift = 0;
  	      shift = (int)( 2 * (Ising_Spin_Coordinates(d, spin_i) - Ising_Spin_Coordinates(d, spin_j)) / total_spins_per_dim[d]);
              double x = Ising_Spin_Coordinates(d, spin_i) - (Ising_Spin_Coordinates(d, spin_j) + shift * total_spins_per_dim[d]);
              distance += pow(x, 2);
	      }
	   distance = sqrt(distance);
	   if(distance <= grid_cut_off){
	      nei_count++;
	      Neighbor_List(nei_count, spin_i) = spin_j;
	      Distance_List(nei_count, spin_i) = distance;
	      }
	   H_value(spin_i, spin_j) = Field;	   
	  }
        }    
     }
}
