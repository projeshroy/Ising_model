#include "../includeallglobals.h"
#include "../Ising_Input.cpp"
#include "../Ising_Energy.cpp"
#include "../Ising_Thread_Initialize.cpp"
#include "../Ising_MasterThread_Initialize.cpp"

int main (int argc, char** argv)
{
   std::string directory[2];
   directory[0]= argv[1];
   directory[1]= argv[2];
   std::string outdir = directory[1];
   std::string indir = directory[0];

   std::string input_file_address = directory[0].append(std::string("ising.in"));
   directory[0] = indir;
   Read_Input(input_file_address);

   int max_threads = 0;
#ifdef _OPENMP
   max_threads = omp_get_max_threads();
#endif   
   MasterThread_Initialize();
   Thread_Initialize();
#pragma omp parallel num_threads(max_threads)
 { 
   Read_Input(input_file_address);
   Thread_Initialize();
 }

   std::string triplet_input_file_address = directory[0].append(std::string("ising_triplet.in"));
   std::ifstream triplet_input_file(triplet_input_file_address.c_str());
   directory[0] = indir;
   
   if (!triplet_input_file.is_open())
      std::cout << " Unable to find file " << std::endl;

   std::string triplet_file_address = directory[1].append(std::string("Ising_Triplet.dat"));
   std::ofstream triplet_file(triplet_file_address.c_str());
   triplet_file.precision(14);
   directory[1] = outdir;

   std::string AAA_file_address = directory[1].append(std::string("Ising_AAA_Triplet_Energy_dist.dat"));
   std::ofstream AAA_file(AAA_file_address.c_str());
   AAA_file.precision(14);
   directory[1] = outdir;

   std::string AAB_file_address = directory[1].append(std::string("Ising_AAB_Triplet_Energy_dist.dat"));
   std::ofstream AAB_file(AAB_file_address.c_str());
   AAB_file.precision(14);
   directory[1] = outdir;

   std::string BAB_file_address = directory[1].append(std::string("Ising_BAB_Triplet_Energy_dist.dat"));
   std::ofstream BAB_file(BAB_file_address.c_str());
   BAB_file.precision(14);
   directory[1] = outdir;
 
   std::cout.unsetf ( std::ios::floatfield );// floatfield not set
   std::cout.precision(10);
   std::cout.setf( std::ios::fixed, std:: ios::floatfield );
//...............................................
//For One-dimensional triplet analysis only
//...............................................
   std::string read_string;

   std::string  spin_file_address, energy_file_address;
   triplet_input_file >> read_string;
   triplet_input_file >> read_string >> spin_file_address;
   triplet_input_file >> read_string >> energy_file_address;
   
   std::ifstream spin_file(spin_file_address.c_str());
   std::ifstream energy_file(energy_file_address.c_str());
   if ( (!spin_file.is_open()) || (!energy_file.is_open()) )
      std::cout << " Unable to find file " << std::endl;

   double min_E, max_E, bin_width;

   triplet_input_file >> read_string >> min_E;
   triplet_input_file >> read_string >> max_E;
   triplet_input_file >> read_string >> bin_width;
   Vec_d Energy_bin; 
   int bin_size = (int)((max_E - min_E)/bin_width) + 1;
   Energy_bin.resize(bin_size);
   for(int i = 0; i < bin_size; i++)
      Energy_bin[i] = min_E + i * bin_width;

//...........................   
   Vec_d temp_Triplet_AAA_Energy, temp_Triplet_AAB_Energy, temp_Triplet_BAB_Energy;
   temp_Triplet_AAA_Energy.resize(total_steps * total_spins); 
   temp_Triplet_AAB_Energy.resize(total_steps * total_spins); 
   temp_Triplet_BAB_Energy.resize(total_steps * total_spins);
        
   double AAA_count = -1;
   double AAB_count = -1;
   double BAB_count = -1;

   AAA_file << "#Energy   Probability" << std::endl;
   AAB_file << "#Energy   Probability" << std::endl;
   BAB_file << "#Energy   Probability" << std::endl;

   triplet_file << "#Triplet  permutation   Probability  Av_Energy" << std::endl;
   
   for(int step = 0; step < total_steps; step++){

      energy_file >> Ising_Energy;

      for(int spin = 0; spin < total_spins; spin++){
          spin_file >> Ising_Spin[spin];  
          energy_file >> Ising_Atomistic_Energy[spin];
	  }
      
      for(int spin = 0; spin < total_spins; spin++){
         double triplet_energy = Ising_Atomistic_Energy[spin] + 
	                      Ising_Atomistic_Energy[Neighbor_List(0, spin)] + 
			      Ising_Atomistic_Energy[Neighbor_List(1, spin)];

         if( (Ising_Spin[spin] == Ising_Spin[Neighbor_List(0, spin)]) && 
             (Ising_Spin[Neighbor_List(1, spin)] == Ising_Spin[Neighbor_List(0, spin)]) ){
	      AAA_count++;
	      temp_Triplet_AAA_Energy[(int)AAA_count] = triplet_energy;}

         if( (Ising_Spin[Neighbor_List(1, spin)] != Ising_Spin[Neighbor_List(0, spin)]) && 
             ((Ising_Spin[spin] == Ising_Spin[Neighbor_List(0, spin)]) || (Ising_Spin[spin] == Ising_Spin[Neighbor_List(1, spin)])) ){
	      AAB_count++;
              temp_Triplet_AAB_Energy[(int)AAB_count] = triplet_energy;}

         if( (Ising_Spin[Neighbor_List(1, spin)] == Ising_Spin[Neighbor_List(0, spin)]) && 
             (Ising_Spin[spin] != Ising_Spin[Neighbor_List(0, spin)]) ){
	     BAB_count++;
	     temp_Triplet_BAB_Energy[(int)BAB_count] = triplet_energy;}
         }
      print_progress_bar<int>(step, total_steps, 10);
      }

      AAA_count++;
      AAB_count++;
      BAB_count++;
      
      Vec_d Triplet_AAA_Energy, Triplet_AAB_Energy, Triplet_BAB_Energy;
      Triplet_AAA_Energy = temp_Triplet_AAA_Energy.topRows((int)AAA_count);
      Triplet_AAB_Energy = temp_Triplet_AAB_Energy.topRows((int)AAB_count);
      Triplet_BAB_Energy = temp_Triplet_BAB_Energy.topRows((int)BAB_count);

      triplet_file << 111 << "   " << 2 << "   " 
                   << AAA_count / (double)(total_steps * total_spins) << "   " 
		   << Triplet_AAA_Energy.sum() / AAA_count << std::endl;
      triplet_file << 112 << "   " << 4 << "   " 
                   << AAB_count / (double)(total_steps * total_spins) << "   " 
		   << Triplet_AAB_Energy.sum() / AAB_count << std::endl;
      triplet_file << 212 << "   " << 2 << "   " 
                   << BAB_count / (double)(total_steps * total_spins) << "   " 
		   << Triplet_BAB_Energy.sum() / BAB_count << std::endl;
      
      auto AAA_Energy_hist = hist(Triplet_AAA_Energy, Energy_bin);
      auto AAB_Energy_hist = hist(Triplet_AAB_Energy, Energy_bin);
      auto BAB_Energy_hist = hist(Triplet_BAB_Energy, Energy_bin);

      for(int i = 0; i < bin_size; i++){
         AAA_file << Energy_bin[i] << "   " << AAA_Energy_hist[i] / AAA_count << std::endl;
         AAB_file << Energy_bin[i] << "   " << AAB_Energy_hist[i] / AAB_count << std::endl;
         BAB_file << Energy_bin[i] << "   " << BAB_Energy_hist[i] / BAB_count << std::endl;
	 }

//..................................................................................................................
}
