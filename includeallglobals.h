//all global functions ...........................................
#include "declarations.h"
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

 void omp_set_num_threads(int num_threads);

 typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vec_d;
 typedef Eigen::Matrix<int, Eigen::Dynamic, 1> Vec_i;
 typedef Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> Vec_ui;
 typedef Eigen::Matrix<std::string, Eigen::Dynamic, 1> Vec_s;
 typedef Eigen::Matrix<bool, Eigen::Dynamic, 1> Vec_b;
 typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> Vec_Compl_d;

 typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> Mat_i;
 typedef Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic> Mat_ui;
 typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Mat_d;
 typedef Eigen::Matrix<std::string, Eigen::Dynamic, Eigen::Dynamic> Mat_s;
 typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> Mat_b;
 typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> Mat_Compl_d;
 std::mt19937 generator = std::mt19937(std::random_device()());

template <class T>
T oddpow(T& x, int power)
{
   T x_power = 1;
   T x_cube = x * x * x;
   
   int mult_power = std::floor((T)power/3.0);       
   int residual_power = power - mult_power * 3;
   if(residual_power == 2) 
      x_power *= x * x;
   else x_power *= x;
   for(int i = 0; i < mult_power; i++)
       x_power *= x_cube;

   return x_power;
}

template <class T>
T evenpow(T& x, int power)
{
   T x_power = 1;
   T x_sq = x * x;
   
   int mult_power = std::floor((T)power/2.0);       
   for(int i = 0; i < mult_power; i++)
       x_power *= x_sq;

   return x_power;
}
 
template <class T>
T tothepower(T& x, int power)
{
  if( std::floor(power/2.0) == std::ceil(power/2.0) )
     return evenpow<T>(x ,power);
  else return oddpow<T>(x ,power);
}

  
double getRand(double range_ini, double range_fin)
{
   std::mt19937 generator = std::mt19937(std::random_device()());
   std::uniform_real_distribution<double> distribution(range_ini, range_fin);
   double _rand = distribution(generator);

   return (_rand);
}

double getRand(double range_ini, double range_fin, double range_gap)
{
   std::mt19937 generator = std::mt19937(std::random_device()());
   std::uniform_real_distribution<double> distribution(range_ini, range_fin);
   double _rand = distribution(generator);
   double new_rand = std::round( (_rand - range_ini)/range_gap )*range_gap + range_ini;

   return (new_rand);
}

Vec_d getUniqueRandVector(int& rand_vector_size, double range_ini, double range_fin, double range_gap)
{
   Vec_d rand_vector; rand_vector.resize(rand_vector_size);
   Vec_d linear_vector; linear_vector.resize(rand_vector_size); linear_vector.fill(1);

   while(linear_vector.sum() != 0){
   int rand_index = (int)getRand(0, rand_vector_size+1);
   if(rand_index < rand_vector_size){
   if(linear_vector[rand_index] == 1){
   rand_vector[rand_index] = rand_index*range_gap + range_ini;
   linear_vector[rand_index] = 0;
   }}}

   return rand_vector;    
}

Vec_d getUniqueRandVector(int& rand_vector_size, double range_ini, double range_fin)
{
   Vec_d rand_vector; rand_vector.resize(rand_vector_size);
   double large_number = getRand(1000, 10000);

   #pragma omp parallel for schedule(static) shared(rand_vector)
   for(int i = 0; i < rand_vector.size(); i++){
    int new_range_ini = range_ini * large_number;
    int new_range_fin = range_fin * large_number;
    rand_vector[i] = getRand(new_range_ini, new_range_fin);
    rand_vector[i] /= large_number;
   }
   return rand_vector;
}

double getMax(Vec_d& array, std::string getindex)
{
   auto max_array = array[0];
   double max_index = 0;
   auto array_size = array.size();
   for(int i = 0 ; i < (array_size - 1); i++)
{
   if ((array[i+1] > array[i]) && (array[i+1] > max_array))
      {
       max_array = array[i+1];
       max_index = i+1;
      }
}
   if (getindex == std::string("getindex"))
      return max_index;
   else return max_array;
}

double getMin(Vec_d& array, std::string getindex)
{
   auto min_array = array[0];
   double min_index = 0;
   auto array_size = array.size();
   for(int i = 0 ; i < (array_size - 1); i++)
{
   if ((array[i+1] < array[i]) && (array[i+1] < min_array))
      {
       min_array = array[i+1];
       min_index = i+1;
      }
}
   if (getindex == std::string("getindex"))
      return min_index;
   else return min_array;
}

double determ(Mat_d& matrix, int n) {

  int det=0;
  int p, h, k, i, j;
  if(matrix.rows() != matrix.cols())
    return 0;

  Mat_d temp;
  temp.resize(matrix.rows(), matrix.cols());

  if(n==1) {
    return matrix(0, 0);
  } else if(n==2) {
    det=(matrix(0, 0)*matrix(1, 1) - matrix(0, 1)*matrix(1, 0));
    return det;
  } else {
    for(p=0;p<n;p++) {
      h = 0;
      k = 0;
      for(i=1;i<n;i++) {
        for( j=0;j<n;j++) {
          if(j==p) {
            continue;
          }
          temp(h, k) = matrix(i, j);
          k++;
          if(k==n-1) {
            h++;
            k = 0;
          }
        }
      }
      det=det + matrix(0, p)*pow(-1,p)*determ(temp, n-1);
    }
    return det;
  }
}

double dotProduct(Vec_d& Vec_A, Vec_d& Vec_B)
{
       auto shape1 = Vec_A.size();
       //auto shape2 = Vec_B.size();
       //if(shape1 != shape2 )
       //  PC_Obj->_Terminate->Error_signal(std::string("size mismatch while dot product product calculation"));

       double dot = 0;
       for(int i = 0; i < shape1 ; i++)
           dot = Vec_A[i] * Vec_B[i] + dot;

       return dot;
}

double dotProduct(Mat_d& Mat_A, Mat_d& Mat_B)
{
       auto shape1 = Mat_A.rows();
       //auto shape2 = Vec_B.size();
       //if(shape1 != shape2 )
       //  PC_Obj->_Terminate->Error_signal(std::string("size mismatch while dot product product calculation"));

       double dot = 0;
       for(int i = 0; i < shape1 ; i++)
           dot = Mat_A(i, 0) * Mat_B(i, 0) + dot;

       return dot;
}

template <class T>
T hist(T& data_vector, T& bins, T& bins_index, bool roundoff)
{
//make sure bin are in low to high order
//bins_index is the same size as bins
      T count_vector;
      count_vector.resize(bins.size());
      count_vector.setZero();
      bins_index.resize(data_vector.size());
      bins_index.setZero();
      
      for(int i = 0; i < data_vector.size(); i++)
         {
         if(data_vector[i] < bins[0]){
            count_vector[0]++;
            bins_index[i] = 0;}
         else if (data_vector[i] >= bins[bins.size()-1]){
             count_vector[bins.size()-1]++;
             bins_index[i] = bins.size()-1;}

         for (int j = 1; j < bins.size(); j++)
             {
              if(roundoff){
              if((data_vector[i] < (bins[j]-((bins[j]-bins[j-1])/2)) )&&(data_vector[i] >= bins[j-1])){
                 count_vector[j-1]++;
                 bins_index[i] =j-1;
                 break;}
              else if((data_vector[i] >= (bins[j]-((bins[j]-bins[j-1])/2)) )&&(data_vector[i] < bins[j])){
                 count_vector[j]++;
                 bins_index[i] = j;
		 break;}
		 }
               else {
               if((data_vector[i] < bins[j])&&(data_vector[i] >= bins[j-1]))
                 {count_vector[j-1]++;
                 bins_index[i] = j-1;
		 break;}
                 }
                }
	     }

       return count_vector;
}

template <class T>
T hist(T& data_vector, T& bins)
{
  bool roundoff = false;
  T bins_index; bins_index.resize(data_vector.size());
  T count_vector = hist<T>(data_vector, bins, bins_index, roundoff);
  return count_vector;
}

template <class T>
T lowtohigh(T& Vec_A)
{
T reshaped_a_array = Vec_A;
T flatted_array;
flatted_array.resize(Vec_A.size());

double smallest_number = getMin(Vec_A, std::string(" "));
double highest_number = getMax(Vec_A, std::string(" "));

int count = 0;

flatted_array[Vec_A.size() - 1] = highest_number;
flatted_array[0] = smallest_number;

int array_size = Vec_A.size();

if(smallest_number == highest_number){
    for (int i = 0; i < array_size; i++)
        flatted_array[i] = smallest_number;
 }
else
{
 while(count < array_size){
      if(smallest_number != highest_number){

         for (int i = 0; i < array_size; i++){
             if(reshaped_a_array[i] == smallest_number){
 	        count++;
	        flatted_array[count-1] = smallest_number;
	        reshaped_a_array[i] = highest_number;
	       }}
		
          smallest_number = getMin(reshaped_a_array, std::string(" "));
	  highest_number = getMax(reshaped_a_array, std::string(" "));
	 }
      else 
      {
         for (int i = count; i < array_size; i++)
	     flatted_array[i] = highest_number;
         count = array_size;
      }					     
}}
 return flatted_array;
}

//print on screen .........   
bool global_progress_bar = false;
template <class T>
T print_progress_bar(T cur_val, T tot_val, int percent, bool& progress_bar)
{
   if( (cur_val > 0)&&(((int)(cur_val * 100/tot_val) % percent) > 1e-5) )
     progress_bar = true;
   if( (progress_bar) && (((int)(cur_val * 100/tot_val) % percent) < 1e-5) ){
     std::cout << (int)(cur_val * 100/tot_val) << " % ..... Completed " << std::endl;     
     progress_bar = false;
     }
}

template <class T>
T print_progress_bar(T cur_val, T tot_val, int percent)
{
print_progress_bar<T>(cur_val, tot_val, percent, global_progress_bar);
}

std::string getFileAddress(std::string& directory, std::string file_name)
{
std::string temp_directory = directory;
temp_directory.append(file_name);
return (temp_directory);
}

template <class T>
int search_array(T& value, Eigen::Matrix<T, Eigen::Dynamic, 1>& ref_array, int division)
{
	int value_id = -1;
	int ref_size = ref_array.size();

//Assuming ref_array is monotonously increasing...................................................................................
	if((value >= ref_array[0]) && (value <= ref_array[ref_size-1])){
	for(int i = 0; i < division; i++){
	int ini_ref_id = i*std::ceil((ref_size/division));
	int fin_ref_id = (i+1)*std::ceil((ref_size/division))-1;

	if(ini_ref_id >= ref_size)
	break;
	if(fin_ref_id >= ref_size)
	fin_ref_id = ref_size-1;
	if((value >= ref_array[ini_ref_id]) && (value <= ref_array[fin_ref_id])){
	for(int j = ini_ref_id; j <= fin_ref_id; j++){
	if(value == ref_array[j]){
	value_id = j;
	break;
	}}}}}

	return value_id;
}

double round(double& x, int n){
	int d = 0;
    	if( (x * pow(10, n + 1)) - (10*(std::floor(x * pow(10, n))) ) > 4) 
          d = 1;
    	x = (floor(x * pow(10, n)) + d) / pow(10, n);
    	return x;
}
// all variables ...........................................................
   Mat_d Ising_Spin_Coordinates, J_value, H_value,  J_matrix, Delta_J, RBIM_p, Neighbor_List, Distance_List;
   Vec_d Ising_Spin, Ising_Spin_type, Ising_Atomistic_Energy, type_spins, flip_type, total_spins_per_dim, type_count;
   Vec_s type_names, Ising_Spin_name;
   int total_spins, total_steps, total_steps_count, equilib_steps, data_write_step, dimension, nearest_neighbor_count, total_types, equil_random_spin_size, prod_random_spin_size;
   double Field, Ising_Energy, temperature, grid_cut_off, Sigma_H, FLIP_PROBABILITY, SWAP_PROBABILITY;
   bool RFIM_, RBIM_;
//....................................................................

 int getType(int spin)
 {
 for(int i = 0; i < total_types; i++){
 if(spin == type_spins[i]){
 return i;
 break;
 }
 } 
}


