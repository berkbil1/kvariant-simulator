#include <iostream> // allows program to perform input and output
#include <string>
#include <cstdlib>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include <windows.h>
#include <random>

using namespace std;
#define PI 3.14159
#define MEMORY_SIZE 50000

//Global variables

int N;								// Memory
int n;								// Vulnerable memory space
int M;								// Number of attempts
int K;								// Number of variants
int a;                              // Vulnerable memory starting index
int b;								// Size of the dummy buffer
int variantFunct;					// Function id that is used when generating variants

double iteration  = 100000;        //# of iteration for simuator
double success = 0;                 //# of successful attack
float range = 1000;                 //range of distribution
int x = 0;                          //index of attack
int pos = 0;                        //index of first attack

int i = 0;                          //dummy buffer 1 
int j = 0;                          //dummy buffer 2 
int attackType = 0;                 //Attack type 1 is for uniform dist. && Attack type 2 is for normal dist.  

int k_variant[20][MEMORY_SIZE];     // First dimension repesents variant and the second dimension represents the memory, [i][0] reserved for checking bit
int currentVariant[MEMORY_SIZE];    // To track current variants for efficency

float probs_array[MEMORY_SIZE];     //To keep probability of each index for normal distribution
float array_sums[MEMORY_SIZE];      //Cummulative function of probs_array
int x_array[MEMORY_SIZE];           //index array for probs_array and array_sums

int spaces[10];                     //It holds space index of each variant
int indexes[10];                    //It holds vulnerable memory starting index of each variant

//Global Variables End


//Merhods Start

//Method Description:Round function
double round(double x)
{
	return ceil(x);
}
//Method Description: Produce random number between min and max
int RandomNumber(int min, int max)
{
	int r = rand() %(max - min + 1);
	return min + r;
}

//Method Description: Normal distribution function
float f_of_x (float x, float sigma) {

	float first_term = 0;
	float second_term = 0;
	float result = 0;

	first_term = 1.0 / (sigma*sqrt(2.0*PI));
	second_term = exp(-(x*x)/(2.0*sigma*sigma));
	result = first_term*second_term;
	return result;
}


//Initialices the array of probabilities and the array of positions
void initialize_arrays (int range_prima) {
	int i;
	for (i = 0; i < range_prima; i++) {
		probs_array[i] = 0;
		x_array[i] = 0;
		array_sums [i] = 0;
	}
}

//Method Description: Probability arrays are calculated according to normal distribution
void calculate_Probs_array_and_array_sums()
{
	float pm=0.001;
	float sigma = 1;
	float increment=0.01;
	float border;
	float p = 0;
	float delta = 0;
	float sum = 0;
	int x = 0;
	float aux = 0;
	int i = 0;

	for (border = 0;; border = border + increment) {		// Find Border
		p = f_of_x (border, sigma);
		if (p < pm) break;
	}

	delta = border / (range/2);

	for (x = -range/2; x < range/2; x++) {
		aux = x*delta;
		sum = sum + p*delta;
		p = f_of_x (aux,sigma);
		probs_array[i] = p;
		x_array [i] = x;
		array_sums[i] = sum;
		i++;
	}
}

//Method Description:Return the index x for the attack according to probability arrays
int calculate_x (int position) {

	int x = 0;
	int x_aux = 0;
	float gen = 1+rand()%9999;			// Random gen: 1 <= pos <= 100
	float sum = 0;

	gen = gen/100;

	while (sum < gen) {
		sum = array_sums [x_aux];
		x_aux++;
	}

	x = x_array[x_aux-1] + position;
	return x;
}

//Method Description: Check that given x is in the determined range N --> 1 <= x <= N
int isoutofrange (int x, int N) {
	if ((x < 1) || (x > N)) return -1;
	else return 1;
}

//Method Description: Multiply probs multiply each probability stored in probs_array *100
void multiply_probs (int range_prima) {
	int i;
	for (i = 0; i < range_prima; i++) {
		array_sums[i] = array_sums[i]*100;
	}
}

//Method Description: Initialize k_variant array with 0 
void initialize_K_variant_array()
{
	int i = 0;
	int j = 0;

	for(i = 0; i < 20; i++)//first dimension
	{
		for (j = 0; j < N; j++)//second dimension
		{
			k_variant[i][j] = 0;
		}
	}
}

//Method Description: Initialize checking bit (second dimension first index) --> k_variant[i][0] = 0
void initialize_check_bit_of_K_variant_array()
{
	int i = 0;

	for(i = 0; i < K; i++)//first dimension
	{	
		k_variant[i][0] = 0;
	}
}

//Method Description: Variants are created according to index of "a" and buffer size "b". Vulnerable memory is reprenseted as 1. 
void produce_variants()
{
	int i = 0;
	int j = 0;
	int index = a;

	for(i = 0; i < K; i++)//first dimension (variant)
	{	
		for (j = 0; j < n; j++)//second dimension
		{
			k_variant[i][(j + index) % N] = 1; //starting from "a" (index) to a+ n- 1, make them 1 to show the vulnerabile memory
		}
		index = index + b; //shifted as buffers size
	}
}
void produce_random_variants(){
	int v; // variants  first address
	int x;
	int variantList[20];
	int variantIndex = 0;
	//initialize the variant list 
	for (int m = 0; m<20; m++)
	{
		variantList[m] = -1;
	}

	//generate variants that allows overlapping
	for (int i = 0; i<K; i++)
	{
		std::random_device rd;     // only used once to initialise (seed) engine
		std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
		std::uniform_int_distribution<int> uni(0, 1000); // guaranteed unbiased
		v = (uni(rng) % (N - n - 2)) + 1; //first variant is the master and other variants generated from the first variant
			
		variantIndex++;
		for (int j = 0; j < n; j++) //assign the first variant
		{
			k_variant[i][(v + j)] = 1;
		}
	}
}

void produce_variants_with_overlapping(){
	int v1; // variants  first address
	int x;
	int variantList[20];
	int variantIndex = 0;
	//initialize the variant list 
	for(int m = 0; m<20; m++)
	{
		variantList[m] = -1;
	}

	v1 = ( rand() % (N-n-2)) + 1; //first variant is the master and other variants generated from the first variant
	variantList[variantIndex] = v1;
	variantIndex++;
	//set the first variant
	for(int j = 0; j < n; j++) //assign the first variant
	{
		k_variant[0][(v1 + j) ] = 1;	
	}

	//generate other overlapping variants
	for(int i = 1; i<K; i++)
	{
GENERATENEWVARIANT:
		int shifting = ( rand() % n ) + 1; //generate random shifting for overlapping variant 
		int sign = rand() %2; //left or right shifting 
		if( sign == 1 ){ // 1 is left shifting
			shifting = shifting * (-1); //assign shifting way
		}

		if(v1 + shifting >= N)  //if new variant is out of boundary
		{
			goto GENERATENEWVARIANT;
		}
		x = v1 + shifting;

		//check that this variant exist or not 
		for(int j=0; j<i; j++)
		{
			if(variantList[j] == x)
				goto GENERATENEWVARIANT;
		}

		variantList[variantIndex] = x;
		variantIndex++;
		for(int j = 0; j < n; j++)
		{
			k_variant[i][(x + j) ] = 1;	
		}
	}
} 

//It produces random gap "s" for each variant according to left memory, then find corresponding indexes and produce variants 
void generateVariantsByUsingGapMethod()
{
	int s = 0;	               //total space left
	int i = 0;                 //dummy variable
	int j = 0;                 //dummy variable
	int v = 0;                 //random variant
	int leftVariants[10];      //track left variants. It is used because of random selection of gaps for each variant
	int cummulativeSpaces = 0; //summation of previos spaces

	if(  (N - K * n) >= 0 ) // control that producing K variants is possible with current N and n 
	{
		//find spaces
		s = N - K * n; 

		for(i = 0; i < K; i++) //produce space for each variant
		{
			if( s != 0)
			{
				spaces[i] =  rand() % s + 1;  //generate random space between 0 to S
				s = s - spaces[i];        //substract the generated space to left space 

			}
			else
			{
				spaces[i] = 0;
			}
		}

		//find indexes according to spaces
		cummulativeSpaces = 1;
		for(j = 0; j < K; j++)
		{
			indexes[j] =  spaces[j] + (n * j) + cummulativeSpaces;
			cummulativeSpaces = cummulativeSpaces + spaces[j];
		}

		//generate variants
		//initialize_K_variant_array(); //clean array for new variant
		for(i = 0; i < K; i++)        //Initialize leftVariants array which will be used to track generated random variants
		{
			leftVariants[i] = 0;
		}
		j = 0;

		//variants are generated randomly based on random spaces
		while(j < K)
		{
			v =  rand() % K; //select random variant number 

			while(leftVariants[v] !=0) // if it not already seleceted 
				v =  rand() % K; //reselect variant again 

			leftVariants[v] = 1; //now, this k is selected

			for(i = 0; i < n; i++)
			{
				k_variant[j][(indexes[v] + i) % N] = 1;	
			}
			j++;
		}
	}
}

//checks two variants are non-overlap
bool variantOK(int v1, int v2)
{
	if( v2 >=v1 && v2 <= v1+n-1) //overlapping
		return false;
	else
	{
		if(v2 < v1)
		{
			if(v2 + n <= v1)
				return true;
			else
				return false;
		}
		else 
			if(v2>=v1+n )
				return true;
	}

	return false;
}

void generateTwoVariantWithSecondVariantVulnerableAtTheEnd()
{

	int i; 
	int v1;
	int v2;
	v1 = ( rand() % (N-n+1) ) + 1;
	for(i = 0; i < n; i++)
	{
		k_variant[0][(v1 + i) ] = 1;	
	}

	v2 = ( rand() % (2 * N / 10 ) ) + (8 * N / 10);
	while(variantOK(v1,v2) == false)
	{
		v2 = ( rand() % (2 * N / 10 ) ) + (8 * N / 10);

	}

	for(i = 0; i < n; i++)
	{
		k_variant[1][(v2 + i) ] = 1;	
	}

}

void VariantGeneratorWithVulnerableEnd()
{
	int i; 
	int v1;
	int v2;
	int v3;
	v1 = ( rand() % (N-n-2)) + 1; //first variant is totaly random from whole memory
	for(i = 0; i < n; i++)
	{
		k_variant[0][(v1 + i) ] = 1;	
	}
VARIANT2:
	v2 = ( rand() % (N-n-2)) + 1; //second variant is totaly random from whole memory
	if(v2 >= v1 && v2<v1+n) //if v2 overlap with v1 
		goto VARIANT2;
	if(v1 >=v2 && v1< v2+n) //if v1 overlap with v2
		goto VARIANT2;

	for(i = 0; i < n; i++)
	{
		k_variant[1][(v2 + i) ] = 1;	
	}

VARIANT3: 
	v3 = ( rand() % ( (2 * N / 10 ) + 1 )  ) + (8 * N / 10) - n; //right 80% of memory

	if(v2 >= v3 && v2<v3+n) //if v2 overlap with v3
		goto VARIANT3;
	if(v3 >=v2 && v3< v2+n) //if v3 overlap with v2
		goto VARIANT3;

	if(v3 >= v1 && v3<v1+n) //if v3 overlap with v1 
		goto VARIANT3;
	if(v1 >=v3 && v1< v3+n) //if v1 overlap with v3
		goto VARIANT3;

	
	for(i = 0; i < n; i++)
	{
		k_variant[2][(v3 + i) ] = 1;	
	}	
}

//Method Description: This function is used in randomDistributingK(). x is random generated "a" and k is variant number
bool isSuitable(int x, int k)
{
	int i;
	int j;

	for(i = 0;  i < k; i++)
	{
		if(k_variant[i][x] == 1 || k_variant[i][x+n] == 1)
			return false;
	}

	return true;
}

//Method Description: This function is used in randomDistributingK()
bool isItInValidRange(int x)
{
	if( x + n > N)
		return false;
	else
		return true;
}


//Method Description: Check wheather all variants are compromised or not. If all variants comromised, it returns "true". This method will be used inside the attack function
bool is_all_variants_compromised()
{
	int i = 0;

	for(i = 0; i < K; i++)//check all variant 
	{
		if(k_variant[i][0] == 0) // if it is not compromised, return false
			return false;
	}

	//if all variants are comprised (k_variant[i][0] == 1 for all i), return true 
	return true; 
}

//Method Description: Takes the attack index x, returns "true", if attack is successful.
bool attack(int x)
{
	int i = 0;
	int j = 0;

	for(i = 0; i < K; i++) //(variant)
	{
		if(k_variant[i][0] == 0) //if this variant is not compromised, check the variant. If it is 1, it is already compromised.
		{
			if(k_variant[i][x] == 1 ) //If attacked index for variant 'i' is vulnerable 
			{
				k_variant[i][0] = 1; //Now, this 'i'th variant is compromised
			}
		}
	}
	//after this line all variants are cheched....

	if(is_all_variants_compromised() == true) // if all variants are compromised, return true
	{ 
		return true; 
	}
	else // all variants are not compromised, return false 
	{
		return false;
	}
}


//Method Description: Retrun delta which will be used for a step
int deltaStepFunction(int N, int M)
{
	double k = 1;

	while(M > (pow(2,k) + 1) )
	{
		k++;
	}

	return N/ pow(2,k);
}

//Methods End


void main (int argc, char *argv[]) {
	
	
	N = atoi(argv[1]);
	n = atoi(argv[2]);
	M = atoi(argv[3]);
	K = atoi(argv[4]);
	range = atoi(argv[5]);
	attackType = atoi(argv[6]);
	variantFunct = atoi(argv[7]);

	//cout<<"N :" <<N <<endl; 
	//cout<<"n :" <<n <<endl; 
	//cout<<"M :" <<M <<endl; 
	//cout<<"K :" <<K <<endl; 
	//cout<<"range :" <<range <<endl; 
	//cout<<"attackType :" <<attackType <<endl;
	//cout<<"variantFunct :" <<variantFunct <<endl;
	
	
	/*
	//Get Input values
	cout << "Please, type memory space                     (N): ";
	cin >> N;
	cout << "Please, type vulnerable memory space          (n): ";
	cin >> n;
	cout << "Please, type the number of attempts           (M): ";
	cin >> M;
	cout << "Please, type the number of variants           (K): ";
	cin >> K;  
	cout << "Please, type the range                    (Range): ";
	cin >> range;
	//cout << "Please, type vulnerable memory starting index (a): ";
	//cin >> a;
	//cout << "Please, type the buffer size                  (b): ";
	//cin >> b;
	//cout << "Please, type the atack type(1-Uniform 2- Normal 3-Delta 4-Mixed 5-Binary 6-Delta-step): ";
	cout << "Please, type the atack type(1-Uniform 2- Normal 5-Binary 6-Delta-step): ";
	cin >> attackType;
	cout << "Please select variant Gen. Function (1-Non-overlap 2- Overlap 3- End): ";
	cin >>variantFunct;
	//Get Input values End
	*/


	if(attackType == 1)//if uniform distribution is selected.
	{
		//START success counter for the graph 
		int successArray[500]; //keep number of success for each M
		for(int i=0; i<M; i++)
		{
			successArray[i] = 0;
		}
		float cumPu[500]; //keep cummulative Pu for each M
		//END success counter for the graph 

		//Start Simulation
		success = 0; //intialize the #of successful attack
		int overlapCount = 0;
		for(i = 0; i < iteration; i++) //for each iteration
		{

			initialize_K_variant_array();

			switch(variantFunct) //Select variant generator function
			{
			case 1:	generateVariantsByUsingGapMethod(); break; //non-overlapping 
			case 2:	produce_variants_with_overlapping();break; //overlapping 
			case 3: VariantGeneratorWithVulnerableEnd(); break;	//vulnerable at the end 
			case 4: produce_random_variants(); break; //total random allowing overlapping and it can also be non-overlapping
			}

			//Start Attack
			for(j = 0; j < M; j++)  //The attacker can attack maximum M times
			{
				x =  rand() % N + 1 ; //genetate random attack index x

				if(attack(x) == true) //if attack is successful (all variants are compromised) 
				{
					successArray[j] = successArray[j] + 1;  //success counter for the current M  
					success++; //#of successful attack is incremented.
					break; //terminate this attack, because it is successful.
				}
			}
			//End Attack
		}
		//End Simulation

		//Calculate Pu
		cout<< "Pu :" << 1- success/iteration   <<endl; 
		
	}
	if(attackType == 2)//if normal distribution is selected.
	{
		//START success counter for the graph 
		int successArray[500]; //keep number of success for each M
		for(int i=0; i<M; i++)
		{
			successArray[i] = 0;
		}
		float cumPu[500]; //keep cummulative Pu for each M
		//END success counter for the graph 

		//Start Simulation
		range = N;
		initialize_arrays(int(range));
		calculate_Probs_array_and_array_sums();
		multiply_probs (int(range));

		success = 0; //intialize the #of successful attack

		for(i = 0; i < iteration; i++) //for each iteration
		{

			initialize_K_variant_array();

			switch(variantFunct) //Select variant generator function
			{
			case 1:	generateVariantsByUsingGapMethod(); break; //non-overlapping 
			case 2:	produce_variants_with_overlapping();break; //overlapping 
			case 3: VariantGeneratorWithVulnerableEnd(); break;	//vulnerable at the end 
			case 4: produce_random_variants(); break; //total random allowing overlapping and it can also be non-overlapping

			}

			pos =  rand() % N + 1 ; //genetate random first attack spos

			//checking weather first attemp is successful or not
			if(attack(pos) == true) //if attack is successful (all variants are compromised) 
			{
				successArray[0] = successArray[0] + 1; //increment the first attack success
				success++; //#of successful attack is incremented.
				continue; //terminate this attack, because it is successful.
			}
			else // first attempt is unsuccessful
			{
				//from 2nd  to Mth attacks

				j = 1;
				while(j < M ) //The attacker can attack maximum M times
				{

					x = calculate_x(pos);
					if (isoutofrange (x,N) == -1) continue; //check that x is in the given interval
					j++;

					if(attack(x) == true) //if attack is successful (all variants are compromised) 
					{
						successArray[j] = successArray[j] + 1; 
						success++; //#of successful attack is incremented.
						break; //terminate this attack, because it is successful.
					}
				}
			}
		}
		//Calculate Pu
		cout<< "Pu :" << 1- success/iteration   <<endl;

		
		//START Print Pu for each M until Mmax
		for(int j = 0 ; j<M; j++)
		{
			int temp = 0;
			for(int i = 0; i<=j; i++)
			{
				temp = temp +  successArray[i];
			} 
			cumPu[j] = 1 - (temp/iteration);
		}
		
	}
	///
	if(attackType == 3)//
	{
		//Start Simulation

		float delta, N1, M1;
		int tempDelta = 0; // step location, incemented in each attempt
		float tempDeltafloat=0; success = 0; //intialize the #of successful attack
		N1=N;
		M1=M;
		delta = N1/(M1-1);
		success = 0; //intialize the #of successful attack

		for(i = 0; i < iteration; i++) //for each iteration
		{
			tempDeltafloat=0;

			tempDelta = 1; // start with first slot
			initialize_K_variant_array();

			switch(variantFunct) //Select variant generator function
			{
			case 1:	generateVariantsByUsingGapMethod(); break; //non-overlapping 
			case 2:	produce_variants_with_overlapping();break; //overlapping 
			case 3: VariantGeneratorWithVulnerableEnd(); break;	//vulnerable at the end 
			case 4: produce_random_variants(); break; //total random allowing overlapping and it can also be non-overlapping

			}

			//Start Attack
			for(j = 0; j < M; j++)  //The attacker can attack maximum M times
			{
				if(attack(tempDelta) == true) //if attack is successful (all variants are compromised) 
				{
					success++; //#of successful attack is incremented.
					break; //terminate this attack, because it is successful.
				} 

				tempDeltafloat = tempDeltafloat + delta; //increment by delta
				tempDelta= round(tempDeltafloat);
				if (tempDelta>N) tempDelta=N;
			}

			//End Attack
		}
		//End Simulation

		//Calculate Pu
		cout<< "Pu :" << 1- success/iteration   <<endl; 
	}
	///
	if(attackType == 4)//
	{
		//Start Simulation

		float A = 3; //# of jumps 
		float numberOfAttempsPerPos = M/A;
		int flag = 0; //to show successful attack

		//initialization functions
		initialize_arrays(int(range));
		calculate_Probs_array_and_array_sums();
		multiply_probs (int(range));

		success = 0; //intialize the #of successful attack

		for(i = 0; i < iteration; i++) //for each iteration
		{
			int tempA = A; //number of distribution 

			initialize_K_variant_array();

			switch(variantFunct) //Select variant generator function
			{
			case 1:	generateVariantsByUsingGapMethod(); break; //non-overlapping 
			case 2:	produce_variants_with_overlapping();break; //overlapping 
			case 3: VariantGeneratorWithVulnerableEnd(); break;	//vulnerable at the end 
			case 4: produce_random_variants(); break; //total random allowing overlapping and it can also be non-overlapping

			};

			pos =  RandomNumber(1,N); //generate random pos between 1 and N, 0 is reserved slot

NEXT: //next distribution 

			if(tempA > 0)
			{
				if(attack(pos) == true) //if attack is successful (all variants are compromised) 
				{
					success++; //#of successful attack is incremented.
					continue; //terminate this attack, because it is successful.
				}
				else // first attempt is unsuccessful
				{
					j = 1;
					flag = 0;
					while(j < numberOfAttempsPerPos) //The attacker can attack maximum M times
					{
						x = calculate_x(pos);
						if (isoutofrange (x,N) == -1) continue; //check that x is in the given interval
						j++;
						if(attack(x) == true) //if attack is successful (all variants are compromised) 
						{
							success++; //#of successful attack is incremented.
							flag = 1; //set flag to prevent next iteration 
							break; //terminate this attack, because it is successful.

						}
					}

					if(flag != 1)//if attack is unsuccessful
					{
						pos =  RandomNumber(1,N); //generate random pos between 1 and N, 0 is reserved slot
						tempA--; // decrement the number of attack
						goto NEXT;
					}					
				}
			}
		}
		//Calculate Pu
		cout<< "Pu :" << 1- success/iteration   <<endl;
	}
	///
	if(attackType == 5)//if binary attack is selected.
	{	
		//START success counter for the graph 
		int successArray[500]; //keep number of success for each M
		for(int i=0; i<M; i++)
		{
			successArray[i] = 0;
		}
		float cumPu[500]; //keep cummulative Pu for each M
		//END success counter for the graph

		//Start Simulation
		for(i = 0; i < iteration; i++) //for each iteration
		{

			initialize_K_variant_array();

			switch(variantFunct) //Select variant generator function
			{
			case 1:	generateVariantsByUsingGapMethod(); break; //non-overlapping 
			case 2:	produce_variants_with_overlapping();break; //overlapping 
			case 3: VariantGeneratorWithVulnerableEnd(); break;	//vulnerable at the end 
			case 4: produce_random_variants(); break; //total random allowing overlapping and it can also be non-overlapping

			}

			int flag = 0;
			int m = 0;
			int denominator  = 2;
			//Start Attack
			while(m < M)  //The attacker can attack maximum M times
			{

				if(flag  == 1) //to terminate attack
					break;

				j = 1;
				while(j < denominator && m <M)
				{
					x = ceil (double((N*j))/denominator); //calculate attack location 

					if(attack(x) == true) //if attack is successful (all variants are compromised) 
					{
						successArray[m] = successArray[m] + 1;
						success++; //#of successful attack is incremented.
						flag = 1; //to terminate the attack
						break;
					}
					m++; // incremenet the number of attempts

					if(m >= M)
					{
						flag = 1;
						break;
					}	
					j = j + 2;
				}
				denominator = denominator *2;
			}
			//End Attack
		}
		//End Simulation
		cout<< "Pu :" << 1- success/iteration   <<endl;
		
	}
	if(attackType == 6)
	{
		//START success counter for the graph 
		int successArray[500]; //keep number of success for each M
		for(int i=0; i<M; i++)
		{
			successArray[i] = 0;
		}
		float cumPu[500]; //keep cummulative Pu for each M
		//END success counter for the graph

		success = 0; //intialize the #of successful attack
		int delta;
		float df;
		df = ((float)N) / (M - 1); //generate random shift constant .
		delta = floor(df);

		for(i = 0; i < iteration; i++) //for each iteration
		{

			initialize_K_variant_array();

			switch(variantFunct) //Select variant generator function
			{
			case 1:	generateVariantsByUsingGapMethod(); break; //non-overlapping 
			case 2:	produce_variants_with_overlapping();break; //overlapping 
			case 3: VariantGeneratorWithVulnerableEnd(); break;	//vulnerable at the end 
			case 4: produce_random_variants(); break; //total random allowing overlapping and it can also be non-overlapping

			}

			x =  1; //start from the beginning
			//Start Attack
			for(j = 0; j < M; j++)  //The attacker can attack maximum M times
			{

				if(attack(x) == true) //if attack is successful (all variants are compromised) 
				{		
					successArray[j] = successArray[j] + 1;
					success++; //#of successful attack is incremented.
					break; //terminate this attack, because it is successful.
				}
				x = ( x + delta) % (N + 1); //shift by delta			
			}
			//End Attack
		}
		//End Simulation

		//Calculate Pu
		cout<< "Pu :" << 1- success/iteration   <<endl;
		
	}
	if(attackType == 7)
	{
		//new attack will be here
	}

}