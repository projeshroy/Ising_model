
void MasterThread_Initialize()
{
  EigenMultiDimArray<double> Grid;
  Grid.resize_DummyMD(total_spins_per_dim);
  for(int spin = 0; spin < total_spins; spin++){
      auto coord = Grid.gdim(spin);
      for(int d = 0; d < dimension; d++)	
      	  Ising_Spin_Coordinates(d, spin) = coord[d];
  }
  Ising_Spin.fill(type_spins[0]);
  Ising_Spin_type.fill(0);
  Ising_Spin_name.fill(type_names[0]);

  for(int i = 1; i < total_types; i++){
      for(int j = 0; j < type_count[i]; j++){
          int random_spin = 0;

          while(Ising_Spin[random_spin] != type_spins[0])
	  random_spin = getRand(0, total_spins-1);

          Ising_Spin[random_spin] = type_spins[i];
          Ising_Spin_type[random_spin] = i;
          Ising_Spin_name[random_spin] = type_names[i];
	  }  
     }
}
