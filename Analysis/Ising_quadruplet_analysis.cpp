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

   std::string quadruplet_input_file_address = directory[0].append(std::string("ising_quadruplet.in"));
   std::ifstream quadruplet_input_file(quadruplet_input_file_address.c_str());
   directory[0] = indir;
   
   if (!quadruplet_input_file.is_open())
      std::cout << " Unable to find file " << std::endl;

   std::string quadruplet_file_address = directory[1].append(std::string("Ising_Quadruplet.dat"));
   std::ofstream quadruplet_file(quadruplet_file_address.c_str());
   quadruplet_file.precision(14);
   directory[1] = outdir;

   std::string AAAA_file_address = directory[1].append(std::string("Ising_AAAA_Quadruplet_Energy_dist.dat"));
   std::ofstream AAAA_file(AAAA_file_address.c_str());
   AAAA_file.precision(14);
   directory[1] = outdir;

   std::string BAAA_file_address = directory[1].append(std::string("Ising_BAAA_Quadruplet_Energy_dist.dat"));
   std::ofstream BAAA_file(BAAA_file_address.c_str());
   BAAA_file.precision(14);
   directory[1] = outdir;

   std::string ABAA_file_address = directory[1].append(std::string("Ising_ABAA_Quadruplet_Energy_dist.dat"));
   std::ofstream ABAA_file(ABAA_file_address.c_str());
   ABAA_file.precision(14);
   directory[1] = outdir;

   std::string BBAA_file_address = directory[1].append(std::string("Ising_BBAA_Quadruplet_Energy_dist.dat"));
   std::ofstream BBAA_file(BBAA_file_address.c_str());
   BBAA_file.precision(14);
   directory[1] = outdir;

   std::string ABBA_file_address = directory[1].append(std::string("Ising_ABBA_Quadruplet_Energy_dist.dat"));
   std::ofstream ABBA_file(ABBA_file_address.c_str());
   ABBA_file.precision(14);
   directory[1] = outdir;

   std::string ABAB_file_address = directory[1].append(std::string("Ising_ABAB_Quadruplet_Energy_dist.dat"));
   std::ofstream ABAB_file(ABAB_file_address.c_str());
   ABAB_file.precision(14);
   directory[1] = outdir;
 
   std::cout.unsetf ( std::ios::floatfield );// floatfield not set
   std::cout.precision(10);
   std::cout.setf( std::ios::fixed, std:: ios::floatfield );
//...............................................
//For One-dimensional quadruplet analysis only
//...............................................
   std::string read_string;

   std::string  spin_file_address, energy_file_address;
   quadruplet_input_file >> read_string;
   quadruplet_input_file >> read_string >> spin_file_address;
   quadruplet_input_file >> read_string >> energy_file_address;
   
   std::ifstream spin_file(spin_file_address.c_str());
   std::ifstream energy_file(energy_file_address.c_str());
   if ( (!spin_file.is_open()) || (!energy_file.is_open()) )
      std::cout << " Unable to find file " << std::endl;

   double min_E, max_E, bin_width;

   quadruplet_input_file >> read_string >> min_E;
   quadruplet_input_file >> read_string >> max_E;
   quadruplet_input_file >> read_string >> bin_width;

   Vec_d Energy_bin; 
   int bin_size = (int)((max_E - min_E)/bin_width) + 1;
   Energy_bin.resize(bin_size);
   for(int i = 0; i < bin_size; i++)
      Energy_bin[i] = min_E + i * bin_width;

//...........................   
   Vec_d temp_Quadruplet_AAAA_Energy, temp_Quadruplet_BAAA_Energy, temp_Quadruplet_ABAA_Energy, 
   temp_Quadruplet_BBAA_Energy, temp_Quadruplet_ABBA_Energy, temp_Quadruplet_ABAB_Energy;

   temp_Quadruplet_AAAA_Energy.resize(total_steps * total_spins); 
   temp_Quadruplet_BAAA_Energy.resize(total_steps * total_spins); 
   temp_Quadruplet_ABAA_Energy.resize(total_steps * total_spins);
   temp_Quadruplet_BBAA_Energy.resize(total_steps * total_spins); 
   temp_Quadruplet_ABBA_Energy.resize(total_steps * total_spins); 
   temp_Quadruplet_ABAB_Energy.resize(total_steps * total_spins);
        
   double AAAA_count = -1;
   double BAAA_count = -1;
   double ABAA_count = -1;
   double BBAA_count = -1;
   double ABBA_count = -1;
   double ABAB_count = -1;

   AAAA_file << "#Energy   Probability" << std::endl;
   BAAA_file << "#Energy   Probability" << std::endl;
   ABAA_file << "#Energy   Probability" << std::endl;
   BBAA_file << "#Energy   Probability" << std::endl;
   ABBA_file << "#Energy   Probability" << std::endl;
   ABAB_file << "#Energy   Probability" << std::endl;

   quadruplet_file << "#Quadruplet  permutation   Probability  Av_Energy" << std::endl;
   
   for(int step = 0; step < total_steps; step++){

      energy_file >> Ising_Energy;

      for(int spin = 0; spin < total_spins; spin++){
          spin_file >> Ising_Spin[spin];  
          energy_file >> Ising_Atomistic_Energy[spin];
	  }
      
      for(int spin = 0; spin < total_spins; spin++){

         Vec_d quadruplet;
	 quadruplet.resize(4);
         int next_spin_index = spin; 
	 double quadruplet_energy = 0;

	 for(int i = 0; i < 4; i++){
	    quadruplet_energy += Ising_Atomistic_Energy[next_spin_index];
	    quadruplet[i] = Ising_Spin[next_spin_index];
	    next_spin_index++;
	    if(next_spin_index == total_spins)
	      next_spin_index = 0;
            }
//.............................

	if(sqrt(pow(quadruplet.sum(), 2)) == 4){
          AAAA_count++;
	  temp_Quadruplet_AAAA_Energy[(int)AAAA_count] = quadruplet_energy;}
        else if(sqrt(pow(quadruplet.sum(), 2)) == 2){
	  if(quadruplet[0] != quadruplet[3]){
            BAAA_count++;
	    temp_Quadruplet_BAAA_Energy[(int)BAAA_count] = quadruplet_energy;}
	  else {
	    ABAA_count++;
	    temp_Quadruplet_ABAA_Energy[(int)ABAA_count] = quadruplet_energy;}
	 }
        else if(sqrt(pow(quadruplet.sum(), 2)) == 0){
	  if(quadruplet[0] != quadruplet[3]){
	    if((quadruplet[0] != quadruplet[1]) && (quadruplet[1] != quadruplet[2]) && (quadruplet[2] != quadruplet[3])){
	       ABAB_count++;
               temp_Quadruplet_ABAB_Energy[(int)ABAB_count] = quadruplet_energy;}
	    else {
               BBAA_count++;
   	       temp_Quadruplet_BBAA_Energy[(int)BBAA_count] = quadruplet_energy;}
            }
	  else {
	    ABBA_count++;
	    temp_Quadruplet_ABBA_Energy[(int)ABBA_count] = quadruplet_energy;}
	 }
//........................	 
      }
      print_progress_bar<int>(step, total_steps, 10);
     }

      AAAA_count++;
      BAAA_count++;
      ABAA_count++;
      BBAA_count++;
      ABBA_count++;
      ABAB_count++;
      
      Vec_d Quadruplet_AAAA_Energy, Quadruplet_BAAA_Energy, Quadruplet_ABAA_Energy,
      Quadruplet_BBAA_Energy, Quadruplet_ABBA_Energy, Quadruplet_ABAB_Energy;

      Quadruplet_AAAA_Energy = temp_Quadruplet_AAAA_Energy.topRows((int)AAAA_count);
      Quadruplet_BAAA_Energy = temp_Quadruplet_BAAA_Energy.topRows((int)BAAA_count);
      Quadruplet_ABAA_Energy = temp_Quadruplet_ABAA_Energy.topRows((int)ABAA_count);
      Quadruplet_BBAA_Energy = temp_Quadruplet_BBAA_Energy.topRows((int)BBAA_count);
      Quadruplet_ABBA_Energy = temp_Quadruplet_ABBA_Energy.topRows((int)ABBA_count);
      Quadruplet_ABAB_Energy = temp_Quadruplet_ABAB_Energy.topRows((int)ABAB_count);
      
      quadruplet_file << 1111 << "   " << 2 << "   " 
                   << AAAA_count / (double)(total_steps * total_spins) << "   " 
		   << Quadruplet_AAAA_Energy.sum() / AAAA_count << std::endl;
      quadruplet_file << 2111 << "   " << 4 << "   " 
                   << BAAA_count / (double)(total_steps * total_spins) << "   " 
		   << Quadruplet_BAAA_Energy.sum() / BAAA_count << std::endl;
      quadruplet_file << 1211 << "   " << 4 << "   " 
                   << ABAA_count / (double)(total_steps * total_spins) << "   " 
		   << Quadruplet_ABAA_Energy.sum() / ABAA_count << std::endl;
      quadruplet_file << 2211 << "   " << 2 << "   " 
                   << BBAA_count / (double)(total_steps * total_spins) << "   " 
		   << Quadruplet_BBAA_Energy.sum() / BBAA_count << std::endl;
      quadruplet_file << 1221 << "   " << 2 << "   " 
                   << ABBA_count / (double)(total_steps * total_spins) << "   " 
		   << Quadruplet_ABBA_Energy.sum() / ABBA_count << std::endl;
      quadruplet_file << 1212 << "   " << 2 << "   " 
                   << ABAB_count / (double)(total_steps * total_spins) << "   " 
		   << Quadruplet_ABAB_Energy.sum() / ABAB_count << std::endl;

      
      auto AAAA_Energy_hist = hist(Quadruplet_AAAA_Energy, Energy_bin);
      auto BAAA_Energy_hist = hist(Quadruplet_BAAA_Energy, Energy_bin);
      auto ABAA_Energy_hist = hist(Quadruplet_ABAA_Energy, Energy_bin);
      auto BBAA_Energy_hist = hist(Quadruplet_BBAA_Energy, Energy_bin);
      auto ABBA_Energy_hist = hist(Quadruplet_ABBA_Energy, Energy_bin);
      auto ABAB_Energy_hist = hist(Quadruplet_ABAB_Energy, Energy_bin);

      for(int i = 0; i < bin_size; i++)
         {
         AAAA_file << Energy_bin[i] << "   " << AAAA_Energy_hist[i] / AAAA_count << std::endl;
         BAAA_file << Energy_bin[i] << "   " << BAAA_Energy_hist[i] / BAAA_count << std::endl;
         ABAA_file << Energy_bin[i] << "   " << ABAA_Energy_hist[i] / ABAA_count << std::endl;
         BBAA_file << Energy_bin[i] << "   " << BBAA_Energy_hist[i] / BBAA_count << std::endl;
         ABBA_file << Energy_bin[i] << "   " << ABBA_Energy_hist[i] / ABBA_count << std::endl;
         ABAB_file << Energy_bin[i] << "   " << ABAB_Energy_hist[i] / ABAB_count << std::endl;
	 }
//..................................................................................................................
}
