#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// Compare Memory Blocks
#include <string.h>

//====================
#include <cstdlib>
#include <ctime>


void init_random() {
    std::srand(static_cast<unsigned int>(std::time(NULL)));
}

// 生成 0 到 359 之間的隨機整數
int random_theta() {
    return std::rand() % 360;
}


//====================



using namespace std;
int const n = 90; /* n-qubit state */
static uint64_t s[4];
static uint64_t x = 0xbf3749f5b97cd3b9; /* The state can be seeded with any value. */
double r; // random number

int correct(int (*Error_vector)[n], int a, int B);

int SQerror(int (*Record_error)[n], int a, double probability); // Single qubit error
int TQerror(int (*Record_error)[n], int a, int b, double probability); // Two qubit error
int CNOT(int (*Error_vector)[n], int a, int b);

int SQerror2(int (*Record_error2)[n], int a, double probability); // Single qubit error
int TQerror2(int (*Record_error2)[n], int a, int b, double probability); // Two qubit error
int CNOT2(int (*Error_vector2)[n], int a, int b);

int TQerror3(int (*Record_error)[n], int (*Record_error2)[n], int a, int b, double probability); // Two qubit error
int CNOT3(int (*Error_vector)[n], int (*Error_vector2)[n], int a, int b);


int CZ(int (*Error_vector)[n], int a, int b);
int perfect(int (*Error_vector)[n]);
int Pauli_correction(int (*Error_vector)[n]);
int measuement(int (*Error_vector)[n], int a, int b);
int pseudo_CE(int (*Error_vector)[n], int a, int theta);

int main_circuit_z(int (*Error_vector)[n], double probability, double alpha, double beta, double gamma, double delta, double theta, double theta2);
int main_circuit_x(int (*Error_vector)[n], double probability, double alpha, double beta, double gamma, double delta, double theta, double theta2);




uint64_t next64() {  // random number
    uint64_t z = (x += 0x9e3779b97f4a7c15);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
    z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
    return z ^ (z >> 31);
}

void rnd256_init() { // random number
    s[0] = next64();  s[1] = next64();  s[2] = next64();  s[3] = next64();
}



int const table_length_FT = 6927;
int const table_width_FT = 22;

int syndrome_FT[table_length_FT][table_width_FT]= { 0 };
int syndrome_result_1[table_width_FT] = {0};
int test1_FT[table_width_FT] = { 0 };  //
int recovery_FT[table_length_FT][24]= { 0 };


int const table_length_x = 16;
int const table_length_z = 84;

int syndrome_x[table_length_x][4]= { 0 };
int syndrome_z[table_length_z][7]= { 0 };
int recovery_x[table_length_x][12]= { 0 };
int recovery_z[table_length_z][12]= { 0 };
int symdorme_result_x[4];
int symdorme_result_z[7];
int test_x[4]; // syndrome talbe for function used
int test_z[7]; // syndrome talbe for function used
int syndrome_result[11] = {0};


int compare_result;
int table_index;



int main(){
    // 數秒counter
    int input;
    clock_t start, end;
    start = clock();
    // print time
    time_t current_time;
    char* c_time_string;
    current_time = time(NULL);
    c_time_string = ctime(&current_time);
    printf("%% Current time is %s", c_time_string);
    
    
    //讀取table
//============================================================================================================================
    
    
    ifstream fin1a("table\\syndore_FT.txt");
    if(!fin1a) {
        cout << "無法讀入檔案\n";
        system("pause");
        return 0;
    }
    for(int i=0;i<table_length_FT;i++)
        for(int k=0;k<table_width_FT;k++)
            fin1a >> syndrome_FT[i][k];
    fin1a.close();
    
    
    
    ifstream fin1b("table\\recovery_FT.txt");
    if(!fin1b) {
        cout << "無法讀入檔案\n";
        system("pause");
        return 0;
    }
    for(int i=0;i<table_length_FT;i++)
        for(int k=0;k<24;k++)
            fin1b >> recovery_FT[i][k];
    fin1b.close();
    
    
//============================================================================================================================
    
    
    ifstream fin_1x("table\\syndore_x.txt");
    if(!fin_1x) {
        cout << "無法讀入檔案\n";
        system("pause");
        return 0;
    }
    for(int i=0;i<table_length_x;i++)
        for(int k=0;k<4;k++)
            fin_1x >> syndrome_x[i][k];
    fin_1x.close();
    

    
    ifstream fin_1z("table\\syndore_z.txt");
    if(!fin_1z) {
        cout << "無法讀入檔案\n";
        system("pause");
        return 0;
    }
    for(int i=0;i<table_length_x;i++)
        for(int k=0;k<7;k++)
            fin_1z >> syndrome_z[i][k];
    fin_1z.close();
    

    ifstream fin_2x("table\\recovery_x.txt");
    if(!fin_2x) {
        cout << "無法讀入檔案\n";
        system("pause");
        return 0;
    }
    for(int i=0;i<table_length_x;i++)
        for(int k=0;k<12;k++)
            fin_2x >> recovery_x[i][k];
    fin_2x.close();
    
    ifstream fin_2z("table\\recovery_z.txt");
    if(!fin_2z) {
        cout << "無法讀入檔案\n";
        system("pause");
        return 0;
    }
    for(int i=0;i<table_length_z;i++)
        for(int k=0;k<12;k++)
            fin_2z >> recovery_z[i][k];
    fin_2z.close();
    

    
    
//============================================================================================================================
    
    
    
    double Total_number, Error_number, detected_number, probability;
    int s1, s2, s3, s4, s5, s6, s7, s8, m1, m2, m3, m4,m5 ,m6 ,q1 ,q2 ,q3 ,ga1, ga2;
    int s9,s10,s11,s12,s13,s14,s15,s16, counter;
    int N = 15; // 切 N 個格子
    int Error_vector[5][n]={0};
    
    double Display_p[N], Display_e[N], Display_d[N];
    
    
    double alpha = 1; // CNOT error
    double beta = 1;  // measurement error
    double gamma = 0.01; // idle error
    double delta = 1;  // ancilla state preparation error
    
    
    cout << "% [[12,1,3]] CE code" << " ,loop=" << N-1 <<" ,gamma = " << gamma <<  ", alpha = " << alpha << ", beta = "  << beta << endl << endl << endl;
    std::ofstream outputFile("1213_FT_Coherent.txt");
    
    for(int i=1; i<N; i++){
        
  
        probability = 0.0001+0.0001*(i-1)*(i); 
        Total_number = 0; //
        Error_number = 0; //
        while(Total_number<10000000||Error_number<600){
        //while(Total_number<1000000){    
        int test_no = -1;    
            

            
            
            //int theta = random_theta();
            //int theta2 = random_theta();
            int theta = 0;
            int theta2 = 0;
            
            s1 = 0;            s2 = 0;            s3 = 0;            s4 = 0;
            s5 = 0;            s6 = 0;            s7 = 0;            s8 = 0;
            s9 = 0;            s10= 0;            s11= 0;            s12= 0;
            s13= 0;            s14= 0;            s15= 0;            s16= 0;
            m1 = 0;            m2 = 0;            m3 = 0;            m4 = 0;
            m5 = 0;            m6 = 0;            q1 = 0;            q2 = 0;
            q3 = 0;            ga1= 0;            ga2= 0;
            
            //**get time***************************
            time_t current_time;
            char* c_time_string;
            current_time = time(NULL);
            c_time_string = ctime(&current_time);
            //*************************************
            
            
            for(int i=0; i<3; i++){
                for(int j=0; j<n; j++)	Error_vector[i][j] = 0; // 在進入新的電路之前, 將error vector歸零 (第一列: 紀錄 X error, 第二列: 紀錄 Z error)
            }
            
            
             SQerror(Error_vector, 1, probability);
             SQerror(Error_vector, 2, probability);
             SQerror(Error_vector, 3, probability);
             SQerror(Error_vector, 4, probability);
             SQerror(Error_vector, 5, probability);
             SQerror(Error_vector, 6, probability);
             SQerror(Error_vector, 7, probability);
             SQerror(Error_vector, 8, probability);
             SQerror(Error_vector, 9, probability);
             SQerror(Error_vector, 10, probability);
             SQerror(Error_vector, 11, probability);
             SQerror(Error_vector, 12, probability);
            
            
//------------------------------------------------------------------------------------------------------------------------------
            
            
            main_circuit_z(Error_vector, probability, alpha, beta, gamma, delta, theta, theta2);
            

            
            main_circuit_x(Error_vector, probability, alpha, beta, gamma, delta, theta, theta2);

				//============================================		                
                          if ( Total_number == test_no ){
                       
				cout << "No."<<Total_number; 
                cout << " & " ;
                for (int i=0; i<12; i++){
                    cout << Error_vector[0][i];
                }
                cout << " & " ;
                for (int i=0; i<12; i++){
                    cout << Error_vector[1][i];
                }
     
                cout << " & " ;
                for (int i=30; i<37; i++){
                    cout << Error_vector[1][i];
                }
                cout << " & " ;
                for (int i=20; i<24; i++){
                    cout << Error_vector[1][i];
                }
                cout <<  " 1st ";		               
                cout << endl;
                        
                              }                            
                //============================================         
            
            syndrome_result_1[0] = Error_vector[1][20];  //m1
            syndrome_result_1[1] = Error_vector[1][21];  //m2
            syndrome_result_1[2] = Error_vector[1][22];  //m3
            syndrome_result_1[3] = Error_vector[1][23];  //m4
            
            syndrome_result_1[8] = Error_vector[1][30];  //m5
            syndrome_result_1[9] = (Error_vector[1][31]+1)%2;  //m6
            syndrome_result_1[10] = (Error_vector[1][32]+1)%2;  //m7
            syndrome_result_1[11] = (Error_vector[1][33]+1)%2;  //m8
            syndrome_result_1[12] = (Error_vector[1][34]+1)%2;  //m9
            syndrome_result_1[13] = (Error_vector[1][35]+1)%2;  //m10
            syndrome_result_1[14] = (Error_vector[1][36]+1)%2;  //m11
            
            
            for (int i=0;i<4;i++){
               ga1 = (syndrome_result_1[i] + ga1);
            }
            for (int i=8;i<15;i++){
               ga1 = (syndrome_result_1[i] + ga1)%2;
            } 
            if (ga1==0){
            	
            	 if ( Total_number == test_no ){
            	cout << endl <<" goto" << endl;
            	}           	
            	goto ckp2;
			}
            
            
            for(int i=0; i<3; i++){
                for(int j=12; j<n; j++)	Error_vector[i][j] = 0; // 在進入新的電路之前, 將error vector歸零 (第一列: 紀錄 X error, 第二列: 紀錄 Z error)
            }
            
            

            
            main_circuit_z(Error_vector, probability, alpha, beta, gamma, delta, theta, theta2);
            

            
            main_circuit_x(Error_vector, probability, alpha, beta, gamma, delta, theta, theta2);
            
            
				//============================================		                
                          if ( Total_number == test_no ){
                       
				cout << "No."<<Total_number; 
                cout << " & " ;
                for (int i=0; i<12; i++){
                    cout << Error_vector[0][i];
                }
                cout << " & " ;
                for (int i=0; i<12; i++){
                    cout << Error_vector[1][i];
                }
                
                cout << " & " ;
                for (int i=30; i<37; i++){
                    cout << Error_vector[1][i];
                }
                cout << " & " ;
                for (int i=20; i<24; i++){
                    cout << Error_vector[1][i];
                }    
                cout <<  " 2nd ";
                cout << endl;
                        
                              }                            
                //============================================ 

            
            // For standard table
            //X-type stabilizer
            syndrome_result[0] = Error_vector[1][20];  //m1
            syndrome_result[1] = Error_vector[1][21];  //m2
            syndrome_result[2] = Error_vector[1][22];  //m3
            syndrome_result[3] = Error_vector[1][23];  //m4
            
            //Z-type stabilizer
            syndrome_result[4] = Error_vector[1][30];  //m5
            syndrome_result[5] = (Error_vector[1][31]+1)%2;  //m6
            syndrome_result[6] = (Error_vector[1][32]+1)%2;  //m7
            syndrome_result[7] = (Error_vector[1][33]+1)%2;  //m8
            syndrome_result[8] = (Error_vector[1][34]+1)%2;  //m9
            syndrome_result[9] = (Error_vector[1][35]+1)%2;  //m10
            syndrome_result[10] = (Error_vector[1][36]+1)%2;  //m11
            
            // For addition  table
            syndrome_result_1[4] = Error_vector[1][20];  //m1
            syndrome_result_1[5] = Error_vector[1][21];  //m2
            syndrome_result_1[6] = Error_vector[1][22];  //m3
            syndrome_result_1[7] = Error_vector[1][23];  //m4
            
            syndrome_result_1[15] = Error_vector[1][30];  //m5
            syndrome_result_1[16] = (Error_vector[1][31]+1)%2;  //m6
            syndrome_result_1[17] = (Error_vector[1][32]+1)%2;  //m7
            syndrome_result_1[18] = (Error_vector[1][33]+1)%2;  //m8
            syndrome_result_1[19] = (Error_vector[1][34]+1)%2;  //m9
            syndrome_result_1[20] = (Error_vector[1][35]+1)%2;  //m10
            syndrome_result_1[21] = (Error_vector[1][36]+1)%2;  //m11
            
            
            
//------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------
            
          
            //standard table
            
            //loading syndorme result to array
            for (int i=0; i<4; i++){
                symdorme_result_x[i] = syndrome_result[i];
            }
            for (int i=0; i<7; i++){
                symdorme_result_z[i] = syndrome_result[i+4];
            }
            
            
            
            
            //compare test1 and symdorme_result_z array and get index j
            for (int j=0; j<table_length_x; j++){
                
                for (int i=0; i<4; i++){
                    test_x[i]=syndrome_x[j][i];
                }
                // 排除syndrome = 0000的情快
                
                if (symdorme_result_x[0] == 0 &&  symdorme_result_x[1] == 0 && symdorme_result_x[2] == 0 && symdorme_result_x[3] == 0   ){
                    break;
                }
                
                compare_result = memcmp(test_x, symdorme_result_x,sizeof(symdorme_result_x));
                
                if (compare_result == 0) {
                    //   table_index = j;
                    //使用recovery進行 x error
                    
                    for (int k=0; k<12; k++){
                        Error_vector[1][k]= (Error_vector[1][k]+recovery_x[j][k])%2;
                    }
                    break;
                    
                }
            }
            
            
            //compare test1 and symdorme_result_x array and get index j
            
            
            for (int j=0; j<table_length_z; j++){
                
                for (int i=0; i<7; i++){
                    test_z[i]=syndrome_z[j][i];
                }
                
                
                if (symdorme_result_z[0] == 0 &&  symdorme_result_z[1] == 1 && symdorme_result_z[2] == 1 && symdorme_result_z[3] == 1 && symdorme_result_z[4] == 1 &&  symdorme_result_z[5] == 1 && symdorme_result_z[6] == 1  ){
                    break;
                }
                
                compare_result = memcmp(test_z, symdorme_result_z,sizeof(symdorme_result_z));
                
                if (compare_result == 0) {
                    //table_index = j;
                    
                    //使用table_r進行 Z error
                    for (int k=0; k<12; k++){
                        Error_vector[0][k]= (Error_vector[0][k]+recovery_z[j][k])%2;
                        
                    }
                //    goto  ckp2;
                    break;
                    
                }
            }
            
            
				//============================================		                
                          if ( Total_number == test_no ){
                       
 				cout << "No."<<Total_number; 
                cout << " & " ;
                for (int i=0; i<12; i++){
                    cout << Error_vector[0][i];
                }
                cout << " & " ;
                for (int i=0; i<12; i++){
                    cout << Error_vector[1][i];
                }
                
                cout << " & " ;
                for (int i=30; i<37; i++){
                    cout << Error_vector[1][i];
                }
                cout << " & " ;
                for (int i=20; i<24; i++){
                    cout << Error_vector[1][i];
                }  
                cout  << " after recovery ";
                cout << endl;
                                      
                }                            
                //============================================   
            

            ckp2:
                
                
//------------------------------------------------------------------------------------------------------------------------------
                
                
                
//---Perfect circuit------------------------------------------------------------------------------------------------------------
                
                
                for(int i=0; i<3; i++){
                    for(int j=12; j<n; j++)	Error_vector[i][j] = 0; // 在進入新的電路之前, 將error vector歸零 (第一列: 紀錄 X error, 第二列: 紀錄 Z error)
                }
                
                perfect(Error_vector );
                
                //X-type stabilizer
                syndrome_result[0] = Error_vector[1][20];  //m1
                syndrome_result[1] = Error_vector[1][21];  //m2
                syndrome_result[2] = Error_vector[1][22];  //m3
                syndrome_result[3] = Error_vector[1][23];  //m4
                
                //Z-type stabilizer
                syndrome_result[4] = Error_vector[1][30];  //m5
                syndrome_result[5] = (Error_vector[1][31]+1)%2;  //m6
                syndrome_result[6] = (Error_vector[1][32]+1)%2;  //m7
                syndrome_result[7] = (Error_vector[1][33]+1)%2;  //m8
                syndrome_result[8] = (Error_vector[1][34]+1)%2;  //m9
                syndrome_result[9] = (Error_vector[1][35]+1)%2;  //m10
                syndrome_result[10] = (Error_vector[1][36]+1)%2;  //m11
                
                
                //loading syndorme result to array
                for (int i=0; i<4; i++){
                    symdorme_result_x[i] = syndrome_result[i];
                }
                for (int i=0; i<7; i++){
                    symdorme_result_z[i] = syndrome_result[i+4];
                }
                
                
                
                
                //compare test1 and symdorme_result_z array and get index j
                for (int j=0; j<table_length_x; j++){
                    
                    for (int i=0; i<4; i++){
                        test_x[i]=syndrome_x[j][i];
                    }
                    // 排除syndrome = 0000的情快
                    
                    if (symdorme_result_x[0] == 0 &&  symdorme_result_x[1] == 0 && symdorme_result_x[2] == 0 && symdorme_result_x[3] == 0   ){
                        break;
                    }
                    
                    compare_result = memcmp(test_x, symdorme_result_x,sizeof(symdorme_result_x));
                    
                    if (compare_result == 0) {
                        //   table_index = j;
                        //使用recovery進行 x error
                        
                        for (int k=0; k<12; k++){
                            Error_vector[1][k]= (Error_vector[1][k]+recovery_x[j][k])%2;
                        }
                        break;
                        
                    }
                }
                
                
                //compare test1 and symdorme_result_x array and get index j
                
                
                for (int j=0; j<table_length_z; j++){
                    
                    for (int i=0; i<7; i++){
                        test_z[i]=syndrome_z[j][i];
                    }
                    
                    
                    if (symdorme_result_z[0] == 0 &&  symdorme_result_z[1] == 1 && symdorme_result_z[2] == 1 && symdorme_result_z[3] == 1 && symdorme_result_z[4] == 1 &&  symdorme_result_z[5] == 1 && symdorme_result_z[6] == 1  ){
                        break;
                    }
                    
                    compare_result = memcmp(test_z, symdorme_result_z,sizeof(symdorme_result_z));
                    
                    if (compare_result == 0) {
                        //table_index = j;
                        
                        //使用table_r進行 Z error
                        for (int k=0; k<12; k++){
                            Error_vector[0][k]= (Error_vector[0][k]+recovery_z[j][k])%2;
                            
                        }
                        break;
                        
                    }
                }
                
                
				//============================================		                
                          if ( Total_number == test_no ){
                       
				cout << "No."<<Total_number; 
                cout << " & " ;
                for (int i=0; i<12; i++){
                    cout << Error_vector[0][i];
                }
                cout << " & " ;
                for (int i=0; i<12; i++){
                    cout << Error_vector[1][i];
                }
                
                cout << " & " ;
                for (int i=30; i<37; i++){
                    cout << Error_vector[1][i];
                }
                cout << " & " ;
                for (int i=20; i<24; i++){
                    cout << Error_vector[1][i];
                }
                cout << " after perfect recovery ";
                cout << endl;
                        
                              }                            
                //============================================ 
                
                
                
                
//------------------------------------------------------------------------------------------------------------------------------
                
                
                //排除logci X, logic Z
                                   
            //logical X operator
            s1 =  (Error_vector[1][4] + Error_vector[1][5] + Error_vector[1][10] + Error_vector[1][11])%2 ;
            //logical X operator
            s2 =  (Error_vector[0][6] + Error_vector[0][8] + Error_vector[0][11])%2 ;
     
            if (s1+s2 ==0){
                    
                    
                    Total_number = Total_number + 1;
                }
                else {
                          
				//============================================		                
                         if ( Total_number == test_no ){
                       
				cout << "No."<<Total_number; 
                cout << " & " ;
                for (int i=0; i<12; i++){
                    cout << Error_vector[0][i];
                }
                cout << " & " ;
                for (int i=0; i<12; i++){
                    cout << Error_vector[1][i];
                }
                
                cout << " & " ;
                for (int i=30; i<37; i++){
                    cout << Error_vector[1][i];
                }
                cout << " & " ;
                for (int i=20; i<24; i++){
                    cout << Error_vector[1][i];
                }
                cout << " after decision ";
                cout << endl;
                        
                              }                            
                //============================================                        
                    Error_number = Error_number + 1;
                    Total_number = Total_number + 1;
                    
                }
                
                
                
        } // while
        Display_p[i-1] = probability;
        Display_e[i-1] = Error_number/(Total_number-detected_number);
        Display_d[i-1] = detected_number/Total_number;
        cout  <<  "% No." << i << " Physical error rate : " << probability << "; Error number = " << Error_number << ", total_number = " << Total_number << ", logical error rate = " <<  Error_number/(Total_number) << "    |" << c_time_string ;//<< endl;
        outputFile  <<  "% No." << i << " Physical error rate : " << probability << "; Error number = " << Error_number << ", total_number = " << Total_number << ", logical error rate = " <<  Error_number/(Total_number) << "    |" << c_time_string ;//<< endl;
        
        
    } // i loop
    
    cout << "Physical_error_rate = [ ";
    for(int i=0; i<N-1; i++){
        if(i<N-2)	cout <<  Display_p[i] << ", ";
        if(i==N-2) cout << Display_p[i];
    } cout << " ]" << endl;
    
    cout << "Logical_error_rate = [ ";
    for(int i=0; i<N-1; i++){
        if(i<N-2)	cout <<  Display_e[i] << ", ";
        if(i==N-2) cout << Display_e[i];
    } cout << " ]" << endl;
    
    
    outputFile << "Physical_error_rate = [ ";
    for(int i=0; i<N-1; i++){
        if(i<N-2)	outputFile <<  Display_p[i] << ", ";
        if(i==N-2) outputFile << Display_p[i];
    } outputFile << " ]" << endl;
    
    outputFile << "Logical_error_rate = [ ";
    for(int i=0; i<N-1; i++){
        if(i<N-2)	outputFile <<  Display_e[i] << ", ";
        if(i==N-2) outputFile << Display_e[i];
    } outputFile << " ]" << endl;
    
    
    end = clock();
    double diff = end - start; // ms
    printf("%% %f  sec", diff / CLOCKS_PER_SEC );
    scanf("%d", &input);
    
    outputFile.close();
    
    return 0;
}

int correct(int (*Error_vector)[n], int a, int B){
    a = a - 1;
    if( B==2 ){
        *(*(Error_vector+0)+a) = (*(*(Error_vector+0)+a) + 1)%2;   // X error時 ,於該qubit補上一個error count
        //	cout << "X error on " << a+1 << endl;
    }
    if( B==3 ){
        *(*(Error_vector+1)+a) = (*(*(Error_vector+1)+a) + 1)%2;    // Z error時,於該qubit補上一個error count
        //	cout << "Z error on " << a+1 << endl;
    }
    if( B==4 ){
        *(*(Error_vector+0)+a) = (*(*(Error_vector+0)+a) + 1)%2;    // Y error時,於該qubit補上一個X & Y error count
        *(*(Error_vector+1)+a) = (*(*(Error_vector+1)+a) + 1)%2;
        //	cout << "Y error on " << a+1 << endl;
    }
    
}

int SQerror(int (*Error_vector)[n], int a, double probability){
    a = a - 1;
    r = next64()/(pow(2,64));
    if( r < (probability/3) ){ // depolarizing channel, r < (1/3)*p 發生 X error
        *(*(Error_vector+0)+a) = (*(*(Error_vector+0)+a) + 1)%2;
    }
    if( r > (probability/3) && r < (2*probability/3) ){ // (1/3)*p < r < (2/3)*p 發生 Z error
        *(*(Error_vector+1)+a) = (*(*(Error_vector+1)+a) + 1)%2;
    }
    if( r > (2*probability/3) && r < (probability) ){ // (2/3)*p < r < p 發生 Y error
        *(*(Error_vector+0)+a) = (*(*(Error_vector+0)+a) + 1)%2;
        *(*(Error_vector+1)+a) = (*(*(Error_vector+1)+a) + 1)%2;
    }
    
}

int TQerror(int (*Error_vector)[n], int a, int b, double probability){
    SQerror(Error_vector, a, probability);
    SQerror(Error_vector, b, probability);
}

int CNOT(int (*Error_vector)[n], int a, int b){
    a = a - 1;
    b = b - 1;
    if(Error_vector[0][a]==1){                                     // X error 發生
        Error_vector[0][b] = (Error_vector[0][b] + 1)%2;
    }
    if(Error_vector[1][b]==1){    								   // Z error 發生
        Error_vector[1][a] = (Error_vector[1][a] + 1)%2;
    }
    
}

int CZ(int (*Error_vector)[n], int a, int b){
    a = a - 1;
    b = b - 1;
    if(Error_vector[0][a]==1){
        Error_vector[1][b] = (Error_vector[1][b] + 1)%2;
    }
    if(Error_vector[0][b]==1){
        Error_vector[1][a] = (Error_vector[1][a] + 1)%2;
    }
}


int CNOT_rev(int (*Error_vector)[n], int b, int a){
    a = a - 1;
    b = b - 1;
    if(Error_vector[0][a]==1){                                     // X error 發生
        Error_vector[0][b] = (Error_vector[0][b] + 1)%2;
    }
    if(Error_vector[1][b]==1){    								   // Z error 發生
        Error_vector[1][a] = (Error_vector[1][a] + 1)%2;
    }
    
}


int CNOT_error(int (*Error_vector)[n], int a, int b, double probability){
    a = a - 1;
    b = b - 1;
    if(Error_vector[0][a]==1){                                     // X error 發生
        Error_vector[0][b] = (Error_vector[0][b] + 1)%2;
    }
    if(Error_vector[1][b]==1){    								   // Z error 發生
        Error_vector[1][a] = (Error_vector[1][a] + 1)%2;
    }
    
    SQerror(Error_vector, a, probability);
    SQerror(Error_vector, b, probability);
    
}

int CZ_error(int (*Error_vector)[n], int a, int b, double probability){
    a = a - 1;
    b = b - 1;
    if(Error_vector[0][a]==1){
        Error_vector[1][b] = (Error_vector[1][b] + 1)%2;
    }
    if(Error_vector[0][b]==1){
        Error_vector[1][a] = (Error_vector[1][a] + 1)%2;
    }
    
    SQerror(Error_vector, a, probability);
    SQerror(Error_vector, b, probability);
    
}


int CNOT_rev_error(int (*Error_vector)[n], int b, int a, double probability){
    a = a - 1;
    b = b - 1;
    if(Error_vector[0][a]==1){                                     // X error 發生
        Error_vector[0][b] = (Error_vector[0][b] + 1)%2;
    }
    if(Error_vector[1][b]==1){    								   // Z error 發生
        Error_vector[1][a] = (Error_vector[1][a] + 1)%2;
	    // phare error rule
	    if( Error_vector[2][b] !=0 ){
	        Error_vector[2][a] = Error_vector[2][b];
	    }  
    }
    
    SQerror(Error_vector, a, probability);
    SQerror(Error_vector, b, probability);
    
}


int pseudo_CE(int (*Error_vector)[n], int a, int theta){
    a = a-1;
    if (Error_vector[0][a] !=0 ){
        Error_vector[2][a] = Error_vector[2][a] + 2*theta;
    }
    
}


int measuement(int (*Error_vector)[n], int a, int b){

	a = a-1; // ancilla qubit
	b = b-1; // data qubit
	
	    double threshold1, threshold2;
		int theta_temp;
        //cout << "start decesion " << endl;
        theta_temp =  Error_vector[2][a];
        
        if ( theta_temp !=0){
	        
	        threshold1 = std::pow(std::cos(theta_temp*M_PI/180), 2); //將theta轉換為degree並且取cos後平方
	        threshold2 = std::pow(std::sin(theta_temp*M_PI/180), 2); //將theta轉換為degree並且取sin後平方
	        
	        r = next64()/(pow(2,64));
	        
	        if ( r < threshold2) { //如果選到sin^2 \theta
				Error_vector[1][a] = 1;
				Error_vector[1][b] = 1;
	            
	        }  
			else { //如果選到cos^2 \theta
	           Error_vector[2][a] = 0;
	           Error_vector[2][b] = 0;
	              
	        }
        
	}	
	
}



int perfect(int (*Error_vector)[n]){
    
    
// gz1=[1,3,5,7,9,11];
    CZ(Error_vector, 1, 31);
    CZ(Error_vector, 3, 31);
    CZ(Error_vector, 5, 31);
    CZ(Error_vector, 7, 31);
    CZ(Error_vector, 9, 31);
    CZ(Error_vector, 11, 31);
    
    
// gz2=[1,2 ];
    CZ(Error_vector, 1, 32);
    CZ(Error_vector, 2, 32);
    
// gz3=[3,4 ];
    CZ(Error_vector, 3, 33);
    CZ(Error_vector, 4, 33);
    
    
    
// gz4=[5,6 ];
    CZ(Error_vector, 5, 34);
    CZ(Error_vector, 6, 34);
    
    
// gz5=[7,8 ];
    CZ(Error_vector, 7, 35);
    CZ(Error_vector, 8, 35);
    
    
    
// gz6=[9,10 ];
    CZ(Error_vector, 9, 36);
    CZ(Error_vector, 10, 36);
    
    
    
// gz7=[11,12 ];
    CZ(Error_vector, 11, 37);
    CZ(Error_vector, 12, 37);
    
    
//=====================================================================
    
    
// gx1=[1,2,3,4 ];
    
    CNOT_rev(Error_vector, 1, 21);
    CNOT_rev(Error_vector, 2, 21);
    CNOT_rev(Error_vector, 3, 21);
    CNOT_rev(Error_vector, 4, 21);
    
    
// gx2=[3,4,5,6 ];
    
    CNOT_rev(Error_vector, 3, 22);
    CNOT_rev(Error_vector, 4, 22);
    CNOT_rev(Error_vector, 5, 22);
    CNOT_rev(Error_vector, 6, 22);
    
    
// gx3=[7,8,9,10 ];
    
    CNOT_rev(Error_vector, 7, 23);
    CNOT_rev(Error_vector, 8, 23);
    CNOT_rev(Error_vector, 9, 23);
    CNOT_rev(Error_vector, 10, 23);
    
    
// gx4=[9,10,11,12 ];
    
    CNOT_rev(Error_vector, 9, 24);
    CNOT_rev(Error_vector, 10, 24);
    CNOT_rev(Error_vector, 11, 24);
    CNOT_rev(Error_vector, 12, 24);
    
    
    
    
}


int Pauli_correction(int (*Error_vector)[n]){
    
    //Z-type
    //==Pauli correction (perfect)========================================
    correct(Error_vector, 3, 3); // qubit 3, apply Z correction
    correct(Error_vector, 7, 3); // qubit 3, apply Z correction
    correct(Error_vector, 11, 3); // qubit 3, apply Z correction
    //==Pauli correction (perfect)========================================
    correct(Error_vector, 2, 3); // qubit 3, apply Z correction
    //==Pauli correction (perfect)========================================
    correct(Error_vector, 4, 3); // qubit 3, apply Z correction
    //==Pauli correction (perfect)========================================
    correct(Error_vector, 6, 3); // qubit 3, apply Z correction
    //==Pauli correction (perfect)========================================
    correct(Error_vector, 8, 3); // qubit 3, apply Z correction
    //==Pauli correction (perfect)========================================
    correct(Error_vector, 10, 3); // qubit 3, apply Z correction
    //==Pauli correction (perfect)========================================
    correct(Error_vector, 11, 3); // qubit 3, apply Z correction
    
    //X-type
    //==Pauli correction (perfect)========================================
    correct(Error_vector, 2, 2); // qubit 3, apply X correction
    correct(Error_vector, 4, 2); // qubit 3, apply X correction
    //==Pauli correction (perfect)========================================
    correct(Error_vector, 4, 2); // qubit 3, apply X correction
    correct(Error_vector, 6, 2); // qubit 3, apply X correction
    //==Pauli correction (perfect)========================================
    correct(Error_vector, 8, 2); // qubit 3, apply X correction
    correct(Error_vector, 10, 2); // qubit 3, apply X correction
    //==Pauli correction (perfect)========================================
    correct(Error_vector, 10, 2); // qubit 3, apply X correction
    correct(Error_vector, 12, 2); // qubit 3, apply X correction
    
}


int main_circuit_z(int (*Error_vector)[n], double probability, double alpha, double beta, double gamma, double delta,double theta, double theta2){
    
     
//*********************************************************************************************************************************
    
//++time slot++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

    //==Ancilla preparation error ========================================
    SQerror(Error_vector, 57, probability*delta);
    SQerror(Error_vector, 58, probability*delta);
    SQerror(Error_vector, 59, probability*delta);
    SQerror(Error_vector, 60, probability*delta);
    SQerror(Error_vector, 61, probability*delta);
    SQerror(Error_vector, 62, probability*delta);

	// add coherent error
	pseudo_CE(Error_vector, 1, theta2);
	pseudo_CE(Error_vector, 2, theta2);
	pseudo_CE(Error_vector, 3, theta2);
	pseudo_CE(Error_vector, 4, theta2);
	pseudo_CE(Error_vector, 5, theta2);
	pseudo_CE(Error_vector, 6, theta2);
	pseudo_CE(Error_vector, 7, theta2);
	pseudo_CE(Error_vector, 8, theta2);
	pseudo_CE(Error_vector, 9, theta2);
	pseudo_CE(Error_vector, 10, theta2);
	pseudo_CE(Error_vector, 11, theta2);
	pseudo_CE(Error_vector, 12, theta2);
	pseudo_CE(Error_vector, 13, theta2);

	// gz1=[1,3,5,7,9,11];
    CZ_error(Error_vector, 1, 57,probability*alpha);
    CZ_error(Error_vector, 3, 58,probability*alpha);
    CZ_error(Error_vector, 5, 59,probability*alpha);
    CZ_error(Error_vector, 7, 60,probability*alpha);
    CZ_error(Error_vector, 9, 61,probability*alpha);
    CZ_error(Error_vector, 11, 62,probability*alpha);
    
	//==idle error ========================================
    SQerror(Error_vector, 2, probability*gamma);
    SQerror(Error_vector, 4, probability*gamma);
    SQerror(Error_vector, 6, probability*gamma);
    SQerror(Error_vector, 8, probability*gamma);   
    SQerror(Error_vector, 10, probability*gamma);
    SQerror(Error_vector, 12, probability*gamma);

	// add coherent error
	pseudo_CE(Error_vector, 1, theta);
	pseudo_CE(Error_vector, 2, theta);
	pseudo_CE(Error_vector, 3, theta);
	pseudo_CE(Error_vector, 4, theta);
	pseudo_CE(Error_vector, 5, theta);
	pseudo_CE(Error_vector, 6, theta);
	pseudo_CE(Error_vector, 7, theta);
	pseudo_CE(Error_vector, 8, theta);
	pseudo_CE(Error_vector, 9, theta);
	pseudo_CE(Error_vector, 10, theta);
	pseudo_CE(Error_vector, 11, theta);
	pseudo_CE(Error_vector, 12, theta);
	pseudo_CE(Error_vector, 13, theta);
    pseudo_CE(Error_vector, 57, theta);
    pseudo_CE(Error_vector, 58, theta);
    pseudo_CE(Error_vector, 59, theta);
    pseudo_CE(Error_vector, 60, theta);
    pseudo_CE(Error_vector, 61, theta);
    pseudo_CE(Error_vector, 62, theta);
		    
	//==idle error ========================================
    SQerror(Error_vector, 1, probability*gamma);
    SQerror(Error_vector, 2, probability*gamma);
    SQerror(Error_vector, 3, probability*gamma);
    SQerror(Error_vector, 4, probability*gamma); 
    SQerror(Error_vector, 5, probability*gamma);
    SQerror(Error_vector, 6, probability*gamma);
    SQerror(Error_vector, 7, probability*gamma);
    SQerror(Error_vector, 8, probability*gamma); 
    SQerror(Error_vector, 9, probability*gamma);
    SQerror(Error_vector, 10, probability*gamma);
    SQerror(Error_vector, 11, probability*gamma);
    SQerror(Error_vector, 12, probability*gamma); 
    
	measuement(Error_vector,57,1);
    measuement(Error_vector,58,3);
    measuement(Error_vector,59,5);
    measuement(Error_vector,60,7);
    measuement(Error_vector,61,9);
    measuement(Error_vector,62,11);
        
		    
    //==measurement error ========================================
    SQerror(Error_vector, 57, probability*beta);
    SQerror(Error_vector, 58, probability*beta);
    SQerror(Error_vector, 59, probability*beta);
    SQerror(Error_vector, 60, probability*beta);
    SQerror(Error_vector, 61, probability*beta);
    SQerror(Error_vector, 62, probability*beta);
    
    for (int i = 56; i <= 61; ++i) {
        Error_vector[1][30] += Error_vector[1][i];
    }
    Error_vector[1][30] %= 2;

	    
//++time slot++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
    
    
//==Ancilla preparation error ========================================
    SQerror(Error_vector, 63, probability*delta);
    SQerror(Error_vector, 64, probability*delta);
    SQerror(Error_vector, 65, probability*delta);
    SQerror(Error_vector, 66, probability*delta);
    SQerror(Error_vector, 67, probability*delta);
    SQerror(Error_vector, 68, probability*delta);
    SQerror(Error_vector, 69, probability*delta);
    SQerror(Error_vector, 70, probability*delta);
    SQerror(Error_vector, 71, probability*delta);
    SQerror(Error_vector, 72, probability*delta);
    SQerror(Error_vector, 73, probability*delta);
    SQerror(Error_vector, 74, probability*delta);
    
    
	// add coherent error
	pseudo_CE(Error_vector, 1, theta2);
	pseudo_CE(Error_vector, 2, theta2);
	pseudo_CE(Error_vector, 3, theta2);
	pseudo_CE(Error_vector, 4, theta2);
	pseudo_CE(Error_vector, 5, theta2);
	pseudo_CE(Error_vector, 6, theta2);
	pseudo_CE(Error_vector, 7, theta2);
	pseudo_CE(Error_vector, 8, theta2);
	pseudo_CE(Error_vector, 9, theta2);
	pseudo_CE(Error_vector, 10, theta2);
	pseudo_CE(Error_vector, 11, theta2);
	pseudo_CE(Error_vector, 12, theta2);
	pseudo_CE(Error_vector, 13, theta2);
    
    
// gz2=[1,2 ];
    CZ_error(Error_vector, 1, 63,probability*alpha);
    CZ_error(Error_vector, 2, 64,probability*alpha);

// gz3=[3,4 ];
    CZ_error(Error_vector, 3, 65,probability*alpha);
    CZ_error(Error_vector, 4, 66,probability*alpha);
  
// gz4=[5,6 ];
    CZ_error(Error_vector, 5, 67,probability*alpha);
    CZ_error(Error_vector, 6, 68,probability*alpha);
  
// gz5=[7,8 ];
    CZ_error(Error_vector, 7, 69,probability*alpha);
    CZ_error(Error_vector, 8, 70,probability*alpha);

// gz6=[9,10 ];
    CZ_error(Error_vector, 9, 71,probability*alpha);
    CZ_error(Error_vector, 10, 72,probability*alpha);

// gz7=[11,12 ];
    CZ_error(Error_vector, 11, 73,probability*alpha);
    CZ_error(Error_vector, 12, 74,probability*alpha);

  
            // add coherent error
            pseudo_CE(Error_vector, 1, theta);
            pseudo_CE(Error_vector, 2, theta);
            pseudo_CE(Error_vector, 3, theta);
            pseudo_CE(Error_vector, 4, theta);
            pseudo_CE(Error_vector, 5, theta);
            pseudo_CE(Error_vector, 6, theta);
            pseudo_CE(Error_vector, 7, theta);
            pseudo_CE(Error_vector, 8, theta);
            pseudo_CE(Error_vector, 9, theta);
            pseudo_CE(Error_vector, 10, theta);
            pseudo_CE(Error_vector, 11, theta);
            pseudo_CE(Error_vector, 12, theta);
            pseudo_CE(Error_vector, 13, theta);
		    pseudo_CE(Error_vector, 63, theta);
		    pseudo_CE(Error_vector, 64, theta);
		    pseudo_CE(Error_vector, 65, theta);
		    pseudo_CE(Error_vector, 66, theta);
		    pseudo_CE(Error_vector, 67, theta);
		    pseudo_CE(Error_vector, 68, theta);
		    pseudo_CE(Error_vector, 69, theta);
		    pseudo_CE(Error_vector, 70, theta);
		    pseudo_CE(Error_vector, 71, theta);
		    pseudo_CE(Error_vector, 72, theta);
		    pseudo_CE(Error_vector, 73, theta);
		    pseudo_CE(Error_vector, 74, theta);

		    
	//==idle error ========================================
    SQerror(Error_vector, 1, probability*gamma);
    SQerror(Error_vector, 2, probability*gamma);
    SQerror(Error_vector, 3, probability*gamma);
    SQerror(Error_vector, 4, probability*gamma); 
    SQerror(Error_vector, 5, probability*gamma);
    SQerror(Error_vector, 6, probability*gamma);
    SQerror(Error_vector, 7, probability*gamma);
    SQerror(Error_vector, 8, probability*gamma); 
    SQerror(Error_vector, 9, probability*gamma);
    SQerror(Error_vector, 10, probability*gamma);
    SQerror(Error_vector, 11, probability*gamma);
    SQerror(Error_vector, 12, probability*gamma); 
    



	measuement(Error_vector,63,1);
    measuement(Error_vector,64,2);
	    
    //==measurement error ========================================
    SQerror(Error_vector, 63, probability*beta);
    SQerror(Error_vector, 64, probability*beta);
    //============================================================
    Error_vector[1][31] = (Error_vector[1][62] + Error_vector[1][63])%2;
    
    
	measuement(Error_vector,65,3);
    measuement(Error_vector,66,4);
	    
    
    //==measurement error ========================================
    SQerror(Error_vector, 65, probability*beta);
    SQerror(Error_vector, 66, probability*beta);
    //============================================================
    Error_vector[1][32] = (Error_vector[1][64] + Error_vector[1][65])%2;
    
    

	measuement(Error_vector,67,5);
    measuement(Error_vector,68,6);

    
    //==measurement error ========================================
    SQerror(Error_vector, 67, probability*beta);
    SQerror(Error_vector, 68, probability*beta);
    //============================================================
    Error_vector[1][33] = (Error_vector[1][66] + Error_vector[1][67])%2;
    

	measuement(Error_vector,69,7);
    measuement(Error_vector,70,8);
	    

    
    
    //==measurement error ========================================
    SQerror(Error_vector, 69, probability*beta);
    SQerror(Error_vector, 70, probability*beta);
    //============================================================
    Error_vector[1][34] = (Error_vector[1][68] + Error_vector[1][69])%2;
    

	measuement(Error_vector,71,9);
    measuement(Error_vector,72,10);

    
    //==measurement error ========================================
    SQerror(Error_vector, 71, probability*beta);
    SQerror(Error_vector, 72, probability*beta);
    //============================================================
    Error_vector[1][35] = (Error_vector[1][70] + Error_vector[1][71])%2;
    

	measuement(Error_vector,73,11);
    measuement(Error_vector,74,12);

    
    //==measurement error ========================================
    SQerror(Error_vector, 73, probability*beta);
    SQerror(Error_vector, 74, probability*beta);
    //============================================================
    Error_vector[1][36] = (Error_vector[1][72] + Error_vector[1][73])%2;
    
//++time slot++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
    

    
    
}

int main_circuit_x(int (*Error_vector)[n], double probability, double alpha, double beta, double gamma, double delta, double theta, double theta2){
    

  
//++time slot++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
    
    
//==Ancilla preparation error ========================================
    
    SQerror(Error_vector, 41, probability*delta);
    SQerror(Error_vector, 42, probability*delta);
    SQerror(Error_vector, 43, probability*delta);
    SQerror(Error_vector, 44, probability*delta);
    SQerror(Error_vector, 49, probability*delta);
    SQerror(Error_vector, 50, probability*delta);
    SQerror(Error_vector, 51, probability*delta);
    SQerror(Error_vector, 52, probability*delta);
    
	// add coherent error
	pseudo_CE(Error_vector, 1, theta2);
	pseudo_CE(Error_vector, 2, theta2);
	pseudo_CE(Error_vector, 3, theta2);
	pseudo_CE(Error_vector, 4, theta2);
	pseudo_CE(Error_vector, 5, theta2);
	pseudo_CE(Error_vector, 6, theta2);
	pseudo_CE(Error_vector, 7, theta2);
	pseudo_CE(Error_vector, 8, theta2);
	pseudo_CE(Error_vector, 9, theta2);
	pseudo_CE(Error_vector, 10, theta2);
	pseudo_CE(Error_vector, 11, theta2);
	pseudo_CE(Error_vector, 12, theta2);
	pseudo_CE(Error_vector, 13, theta2);
    
// gx1=[1,2,3,4 ];
    
    CNOT_rev_error(Error_vector, 1, 41,probability*alpha);
    CNOT_rev_error(Error_vector, 2, 42,probability*alpha);
    CNOT_rev_error(Error_vector, 3, 43,probability*alpha);
    CNOT_rev_error(Error_vector, 4, 44,probability*alpha);
    

// gx3=[7,8,9,10 ];
    
    CNOT_rev_error(Error_vector, 7, 49,probability*alpha);
    CNOT_rev_error(Error_vector, 8, 50,probability*alpha);
    CNOT_rev_error(Error_vector, 9, 51,probability*alpha);
    CNOT_rev_error(Error_vector, 10, 52,probability*alpha);
    

	//==idle error ========================================
    SQerror(Error_vector, 5, probability*gamma);
    SQerror(Error_vector, 6, probability*gamma);
    SQerror(Error_vector, 11, probability*gamma);
    SQerror(Error_vector, 12, probability*gamma);   

	 
    // add coherent error
    pseudo_CE(Error_vector, 1, theta);
    pseudo_CE(Error_vector, 2, theta);
    pseudo_CE(Error_vector, 3, theta);
    pseudo_CE(Error_vector, 4, theta);
    pseudo_CE(Error_vector, 5, theta);
    pseudo_CE(Error_vector, 6, theta);
    pseudo_CE(Error_vector, 7, theta);
    pseudo_CE(Error_vector, 8, theta);
    pseudo_CE(Error_vector, 9, theta);
    pseudo_CE(Error_vector, 10, theta);
    pseudo_CE(Error_vector, 11, theta);
    pseudo_CE(Error_vector, 12, theta);
    pseudo_CE(Error_vector, 13, theta);
    pseudo_CE(Error_vector, 41, theta);
    pseudo_CE(Error_vector, 42, theta);
    pseudo_CE(Error_vector, 43, theta);
    pseudo_CE(Error_vector, 44, theta);
    pseudo_CE(Error_vector, 49, theta);
    pseudo_CE(Error_vector, 50, theta);
    pseudo_CE(Error_vector, 51, theta);
    pseudo_CE(Error_vector, 52, theta);
		    
    
	//==idle error ========================================
    SQerror(Error_vector, 1, probability*gamma);
    SQerror(Error_vector, 2, probability*gamma);
    SQerror(Error_vector, 3, probability*gamma);
    SQerror(Error_vector, 4, probability*gamma); 
    SQerror(Error_vector, 5, probability*gamma);
    SQerror(Error_vector, 6, probability*gamma);
    SQerror(Error_vector, 7, probability*gamma);
    SQerror(Error_vector, 8, probability*gamma); 
    SQerror(Error_vector, 9, probability*gamma);
    SQerror(Error_vector, 10, probability*gamma);
    SQerror(Error_vector, 11, probability*gamma);
    SQerror(Error_vector, 12, probability*gamma); 
    
    
    //measurement    
	measuement(Error_vector,41,1);
    measuement(Error_vector,42,2);
    measuement(Error_vector,43,3);
    measuement(Error_vector,44,4);
    

    //==measurement error ========================================
    SQerror(Error_vector, 41, probability*beta);
    SQerror(Error_vector, 42, probability*beta);
    SQerror(Error_vector, 43, probability*beta);
    SQerror(Error_vector, 44, probability*beta);
    //============================================================
    
    for (int i = 40; i <= 43; ++i) {
        Error_vector[1][20] += Error_vector[1][i];
    }
    Error_vector[1][20] %= 2;
    
    
    //measurement
    //measurtment_4(Error_vector,49,50,51,52);
  

    
	measuement(Error_vector,49,7);
    measuement(Error_vector,50,8);
    measuement(Error_vector,51,9);
    measuement(Error_vector,52,10);
    

    
    //==measurement error ========================================
    SQerror(Error_vector, 49, probability*beta);
    SQerror(Error_vector, 50, probability*beta);
    SQerror(Error_vector, 51, probability*beta);
    SQerror(Error_vector, 52, probability*beta);
    //============================================================
    
    for (int i = 48; i <= 51; ++i) {
        Error_vector[1][22] += Error_vector[1][i];
    }
    Error_vector[1][22] %= 2;
    
    

//++time slot++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

//==Ancilla preparation error ========================================
    SQerror(Error_vector, 45, probability*delta);
    SQerror(Error_vector, 46, probability*delta);
    SQerror(Error_vector, 47, probability*delta);
    SQerror(Error_vector, 48, probability*delta);
    SQerror(Error_vector, 53, probability*delta);
    SQerror(Error_vector, 54, probability*delta);
    SQerror(Error_vector, 55, probability*delta);

	// add coherent error
	pseudo_CE(Error_vector, 1, theta2);
	pseudo_CE(Error_vector, 2, theta2);
	pseudo_CE(Error_vector, 3, theta2);
	pseudo_CE(Error_vector, 4, theta2);
	pseudo_CE(Error_vector, 5, theta2);
	pseudo_CE(Error_vector, 6, theta2);
	pseudo_CE(Error_vector, 7, theta2);
	pseudo_CE(Error_vector, 8, theta2);
	pseudo_CE(Error_vector, 9, theta2);
	pseudo_CE(Error_vector, 10, theta2);
	pseudo_CE(Error_vector, 11, theta2);
	pseudo_CE(Error_vector, 12, theta2);
	pseudo_CE(Error_vector, 13, theta2);
	
// gx2=[3,4,5,6 ];
    
    CNOT_rev_error(Error_vector, 3, 45,probability*alpha);
    CNOT_rev_error(Error_vector, 4, 46,probability*alpha);
    CNOT_rev_error(Error_vector, 5, 47,probability*alpha);
    CNOT_rev_error(Error_vector, 6, 48,probability*alpha);
    

    
// gx4=[9,10,11,12 ];
    
    CNOT_rev_error(Error_vector, 9, 53,probability*alpha);
    CNOT_rev_error(Error_vector, 10, 54,probability*alpha);
    CNOT_rev_error(Error_vector, 11, 55,probability*alpha);
    CNOT_rev_error(Error_vector, 12, 56,probability*alpha);
    
	//==idle error ========================================
    SQerror(Error_vector, 1, probability*gamma);
    SQerror(Error_vector, 2, probability*gamma);
    SQerror(Error_vector, 7, probability*gamma);
    SQerror(Error_vector, 8, probability*gamma);   
    
    
    // add coherent error
    pseudo_CE(Error_vector, 1, theta);
    pseudo_CE(Error_vector, 2, theta);
    pseudo_CE(Error_vector, 3, theta);
    pseudo_CE(Error_vector, 4, theta);
    pseudo_CE(Error_vector, 5, theta);
    pseudo_CE(Error_vector, 6, theta);
    pseudo_CE(Error_vector, 7, theta);
    pseudo_CE(Error_vector, 8, theta);
    pseudo_CE(Error_vector, 9, theta);
    pseudo_CE(Error_vector, 10, theta);
    pseudo_CE(Error_vector, 11, theta);
    pseudo_CE(Error_vector, 12, theta);
    pseudo_CE(Error_vector, 13, theta);
    pseudo_CE(Error_vector, 45, theta);
    pseudo_CE(Error_vector, 46, theta);
    pseudo_CE(Error_vector, 47, theta);
    pseudo_CE(Error_vector, 48, theta);
    pseudo_CE(Error_vector, 53, theta);
    pseudo_CE(Error_vector, 54, theta);
    pseudo_CE(Error_vector, 55, theta);
    
    
	//==idle error ========================================
    SQerror(Error_vector, 1, probability*gamma);
    SQerror(Error_vector, 2, probability*gamma);
    SQerror(Error_vector, 3, probability*gamma);
    SQerror(Error_vector, 4, probability*gamma); 
    SQerror(Error_vector, 5, probability*gamma);
    SQerror(Error_vector, 6, probability*gamma);
    SQerror(Error_vector, 7, probability*gamma);
    SQerror(Error_vector, 8, probability*gamma); 
    SQerror(Error_vector, 9, probability*gamma);
    SQerror(Error_vector, 10, probability*gamma);
    SQerror(Error_vector, 11, probability*gamma);
    SQerror(Error_vector, 12, probability*gamma); 
    
    
	measuement(Error_vector,45,3);
    measuement(Error_vector,56,4);
    measuement(Error_vector,47,5);
    measuement(Error_vector,48,6);
    

    //==measurement error ========================================
    SQerror(Error_vector, 45, probability*beta);
    SQerror(Error_vector, 46, probability*beta);
    SQerror(Error_vector, 47, probability*beta);
    SQerror(Error_vector, 48, probability*beta);
    //============================================================
    
    for (int i = 44; i <= 47; ++i) {
        Error_vector[1][21] += Error_vector[1][i];
    }
    Error_vector[1][21] %= 2;
    
    

    
	measuement(Error_vector,53,9);
    measuement(Error_vector,54,10);
    measuement(Error_vector,55,11);
    measuement(Error_vector,56,12);

    
    //==measurement error ========================================
    SQerror(Error_vector, 53, probability*beta);
    SQerror(Error_vector, 54, probability*beta);
    SQerror(Error_vector, 55, probability*beta);
    SQerror(Error_vector, 56, probability*beta);
    //============================================================
    
    for (int i = 52; i <= 55; ++i) {
        Error_vector[1][23] += Error_vector[1][i];
    }
    Error_vector[1][23] %= 2;
    
//++time slot++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
    
}
