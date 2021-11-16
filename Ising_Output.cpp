
void Write_Trajectory(std::ofstream& trajectory_file, Vec_s& New_Ising_Spin_name)
{

 trajectory_file << total_spins << std::endl;
 trajectory_file << std::endl;

 int d;
 for (unsigned int i = 0; i < total_spins; i++){
      trajectory_file << New_Ising_Spin_name[i] << "    ";

         for (d =0; d< dimension; d++)
             trajectory_file << Ising_Spin_Coordinates(d, i) << "   ";
      
          while(d <= 2){
	       trajectory_file << 0.0 << "   ";
	       d = d+1;}
          
	   trajectory_file << std::endl;
     }
 trajectory_file << std::endl;
}

/*
  int Plus_Spin_Counts = total_spins - (total_spins - New_Ising_Spin.sum())/2;

//printing the header
//  trajectory_file << "# STRUCTURE BLOCK" << std::endl;
//  trajectory_file << "# define the default atom" <<std::endl;
//  trajectory_file << "atom default radius   " << 0.8 << "  name   " << "H  " << std::endl;
//  trajectory_file << "# define other atoms, undefined atoms will be filled with the default" << std::endl;
  
  int Plus_Spin_Counts = total_spins - (total_spins - New_Ising_Spin.sum())/2;
 
  trajectory_file << "atom   " << 0
                  << ":" << (Plus_Spin_Counts - 1)
		  << "     radius    " << 0.5 << "   name   "
		  << std::string("P") << std::endl;

  trajectory_file << "atom   " << Plus_Spin_Counts
                  << ":" << (total_spins - 1)
		  << "     radius    " << 0.5 << "   name   "
		  << std::string("N") << std::endl;

//writing the coordinates 

// trajectory_file << "# TIMESTEP BLOCKS" <<std::endl;
// trajectory_file << "# start a new timestep (ordered by default)" << std::endl;
 trajectory_file << "timestep" <<std::endl;
// trajectory_file << "# set the unitcell" <<std::endl;

  trajectory_file << "pbc    ";
  unsigned int d;
  for (d = 0; d < dimension; d++)
      trajectory_file << total_spins_per_dim[d] << "   ";

  while(d <= 2){
     trajectory_file << 1.0 << "   ";
     d = d+1;}

     trajectory_file << std::endl;

//   trajectory_file << "# now define the coordinates" << std::endl;

  for (unsigned int i = 0; i < total_spins; i++){
      if( New_Ising_Spin[i] == 1){
         for (d =0; d< dimension; d++)
             trajectory_file << Ising_Spin_Coordinates(d, i) << "   ";
      
          while(d <= 2){
	       trajectory_file << 0.0 << "   ";
	       d = d+1;}
          
	   trajectory_file << std::endl;
	   }
      }

   for (unsigned int i = 0; i < total_spins; i++){   
      if( New_Ising_Spin[i] == -1){
         for (d = 0; d< dimension; d++)
             trajectory_file << Ising_Spin_Coordinates(d, i) << "   ";
      
          while(d <= 2){
	       trajectory_file << 0.0 << "   ";
	       d = d+1;}
          
	   trajectory_file << std::endl;
	   }
   
       }
*/       

