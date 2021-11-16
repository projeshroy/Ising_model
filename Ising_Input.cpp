
//Read Input file...............................................................................................   
void Read_Input(std::string& input_file_address)
{
   std::ifstream input_file(input_file_address.c_str());
   
   if (!input_file.is_open())
      std::cout << " Unable to find file " << std::endl;

   std::string read_string;
   input_file >> read_string;
   input_file >> read_string >> dimension;
   input_file >> read_string >> total_spins;

   Ising_Spin_Coordinates.resize(dimension, total_spins); Ising_Spin_Coordinates.setZero();
   Ising_Spin.resize(total_spins); Ising_Spin.setZero();
   Ising_Spin_type.resize(total_spins); Ising_Spin_type.setZero();
   Ising_Spin_name.resize(total_spins); 
   total_spins_per_dim.resize(dimension); total_spins_per_dim.setZero();
   Ising_Atomistic_Energy.resize(total_spins); Ising_Atomistic_Energy.setZero();
   J_value.resize(total_spins, total_spins); J_value.setZero();
   H_value.resize(total_spins, total_spins); H_value.setZero();

   input_file >> read_string;
   for(int d = 0; d < dimension; d++)
       input_file >> read_string >> total_spins_per_dim[d];

   
   input_file >> read_string >> grid_cut_off;//Nearest neighbor at dist 1.
   input_file >> read_string >> temperature;
   input_file >> read_string >> Field;

   input_file >> read_string >> total_types;
   type_names.resize(total_types); type_spins.resize(total_types); type_count.resize(total_types);
   input_file >> read_string >> read_string >> read_string;
   for(int i = 0; i < total_types; i++)
   input_file >> type_names[i] >> type_spins[i] >> type_count[i];

   input_file >> read_string >> FLIP_PROBABILITY;
   input_file >> read_string;
   flip_type.resize(total_types);
   for(int i = 0; i < total_types; i++)
   input_file >> read_string >> flip_type[i];
   input_file >> read_string >> SWAP_PROBABILITY;   

   input_file >> read_string >> read_string;
   if(read_string == std::string("yes"))
      RFIM_ = true;
   else RFIM_ = false;   
   if(RFIM_)
     input_file >> read_string >> Sigma_H;
   
   J_matrix.resize(total_types, total_types);
   Delta_J.resize(total_types, total_types);
   RBIM_p.resize(total_types, total_types);

   input_file >> read_string >> read_string;
   if(read_string == std::string("yes"))
      RBIM_ = true;
   else RBIM_ = false;

   input_file >> read_string;
   if(RBIM_)
     input_file >> read_string >> read_string;

   for(int i = 0; i < total_types; i++){
      for(int j = 0; j < total_types; j++){
         input_file >> read_string >> J_matrix(i, j);
         if(RBIM_)
         input_file >> Delta_J(i, j) >> RBIM_p(i, j);
         }
     }
   
   input_file >> read_string >> total_steps;
   input_file >> read_string >> equilib_steps;
   input_file >> read_string >> data_write_step;
   total_steps_count = std::floor((double)total_steps/(double)data_write_step) + 1;

   int thread_id = omp_get_thread_num();
   if(thread_id == 0){
   std::cout << " total_steps " << total_steps << std::endl;
   std::cout << " equilib_steps " << equilib_steps << std::endl; 
   std::cout << " data_write_step " << data_write_step << std::endl;
   std::cout << " total_steps_count " << total_steps_count << std::endl;
   std::cout << " dimension " << dimension << std::endl;
   std::cout << " total_spins " << total_spins << std::endl;
   std::cout << " total_spins_per_dim " << total_spins_per_dim << std::endl;
   std::cout << " grid_cut_off " << grid_cut_off << std::endl;
   std::cout << " temperature " << temperature << std::endl;
   std::cout << " Field " << Field << std::endl;
   std::cout << " FLIP_PROBABILITY " << FLIP_PROBABILITY << std::endl;
   std::cout << " SWAP_PROBABILITY " << SWAP_PROBABILITY << std::endl;
   std::cout << " RFIM_ " << RFIM_ << std::endl;
   if(RFIM_)
    std::cout << " Sigma_H " << Sigma_H << std::endl;
   std::cout << " RBIM_ " << RBIM_ << std::endl;
   if(RBIM_)
     std::cout << " Delta_J " << Delta_J << " RBIM_p " << RBIM_p << std::endl;
   std::cout << " total_types " << total_types << std::endl;
   std::cout << " type_names  type_spins  type_count " << std::endl;
   for(int i = 0; i < total_types; i++)
      std::cout << "   " << type_names[i] << "   " << type_spins[i] << "   " << type_count[i] << std::endl;
   std::cout << " FLIP_PROBABILITY " << FLIP_PROBABILITY << std::endl;
   std::cout << " flip_type " << flip_type.transpose() << std::endl;
   std::cout << " SWAP_PROBABILITY " << SWAP_PROBABILITY << std::endl;
   std::cout << " J_matrix \n" << J_matrix << std::endl;
    }
}

