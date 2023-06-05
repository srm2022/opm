
# This is a generic C++implementation of IGV. Please customise it as per 
# your own requirements, by changing its termination condition, reported 
# performance indicators, etc. We customised it differently for different 
# sections of the paper. Please also make sure the path points to your data 
# folder. To keep it simple, it only runs once, for 10 seconds, on each 
# instance of the large dataset and prints the best objective value and the 
# duration time (which is around 10 seconds, excluding the time spent on 
# reading input data). 
# Please note that this code is meant to be simple so there is room for 
# optimising it.
# Comments and questions are welcome to seyed_r_mousavi@yahoo.com and new 
# versions (if any) will be made available at https://github.com/srm2022/opm, 
# as long as these accounts remain active. Thank you.

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <time.h>
#include <string>
#include <vector>
#include <math.h>

using namespace std;

const int MAX_N = 900; //NOTE: please use a tight upper bound on the number of vertices (MAX_N >= the maximum number of nodes in all instances).

static const std::string WHITESPACE = " \n\r\t\f\v";

static int MY_INF = 1000000; //must be larger than the sum of the weights of all the edges in the given graph but no more than half of the maximum int value (to avoid oberflow when adding two of them)
static int n, p, q;
static int* D[MAX_N];
static int* N[MAX_N];
static int* Ninverse[MAX_N];
static bool x[MAX_N]; //list of indicator variables xj which is true iff facility j is in P

static double gamma = 0.6;
static double tau = 6;




string ltrim(const string& s) {
    size_t start = s.find_first_not_of(WHITESPACE);

    return (start == std::string::npos) ? "" : s.substr(start);
}

string rtrim(const string& s) {
    size_t end = s.find_last_not_of(WHITESPACE);

    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

string trim(const string& s) {

    return rtrim(ltrim(s));
}
//---------------

vector<string> split(string s, char d) {
    vector<string> r;
    int j = 0;
    for (int i = 0; i < s.length(); i++) {
        if (s[i] == d) {
            string cur = s.substr(j, i - j);
            if (cur.length()) {
                r.push_back(cur);
            }
            j = i + 1;
        }
    }
    if (j < s.length()) {
        r.push_back(s.substr(j));
    }
    return r;
}//split
//-------------



    static void ee(string message) {
        cout << message <<endl;
        exit(1);
    }//ee
    //----------


    static int nextpos(int c, int curpos)
    {
        int pos = 1 + curpos;
        while(!x[N[c][pos]])
            pos++;
        return pos;
    }//nextpos
    //----------


    static double g(int radius){
        return gamma * pow(2, tau*(1-radius)) + (1-gamma) * (rand()*1.0) / RAND_MAX;
    }//g
    //----------
    
    static int calcf(int pP[], int pp)
    {
        int i, c, tmpmin;
        int ret_val;
               
        ret_val = 0;
        for(c=0; c<n; c++)
        {
            tmpmin = MY_INF;
            for(i=0; i<pp; i++)
            {
                if(D[pP[i]][c] < tmpmin)
                    tmpmin = D[pP[i]][c];
            }
            ret_val += tmpmin;
        }
        return ret_val;
    }//calcf
    //----------

    static int initRandSol(int P[], int U[], bool x[], int F[], int* C[MAX_N], int Cinx[], int nC[])
    {
        int k, fa, tmpf, i, c, minpos, minfa, mind, tmpinx;
        double tmpPbyN = p * 1.0 / n;
        for (fa = 0; fa < n; fa++) {
            x[fa] = false;
        }

        k = 0;
        fa = rand() % n; //random starting point
        while (k < p) {
            if (!x[fa]) {
                if ( (rand()*1.0)/RAND_MAX < tmpPbyN) {
                    P[k] = fa;
                    x[fa] = true;
                    k++;
                }
            }
            fa = (fa + 1) % n;
        }

        //copy all unselected ones into U
        k = 0;
        for (fa = 0; fa < n; fa++) {
            if (!x[fa]) {
                U[k] = fa; 
                k++;
            }
        }
        
        for (i = 0; i < n; i++) { //n not p (for all potential fa)
            nC[i] = 0;
        }
        tmpf = 0;
        for (c = 0; c < n; c++) {
            minpos = nextpos(c, -1);
            minfa = N[c][minpos];                
            mind = D[minfa][c];

            F[c] = minfa;
            tmpinx = nC[minfa];
            C[minfa][tmpinx] = c;
            Cinx[c] = tmpinx;
            nC[minfa]++;
            tmpf += mind;
        }

        return tmpf;
    }//initRandSol
    //----------
    
    static long openfa(int P[], int pp, int U[], bool x[], int F[], int* C[MAX_N], int Cinx[], int nC[], double alpha)
    {
        int j, j_best_fa, fa, c, best_fa, oldfa, oldinx, tmpVlast, newinx;
        long mindf;
        long tmpabsdf[n];
        
        int tmpq = n - pp;
        long ret_val;
               
        mindf = MY_INF;
        j_best_fa = -1;
        for (j = 0; j < tmpq; j++) {
            if ((rand()*1.0)/RAND_MAX > alpha) {
                continue;
            }
            fa = U[j];
            tmpabsdf[j] = 0;
            for (c = 0; c < n; c++) {
                //calc contributation of c to tmpdf
                if (D[fa][c] < D[F[c]][c]) {
                    tmpabsdf[j] += D[F[c]][c] - D[fa][c];
                }
            }

            if (tmpabsdf[j] < mindf) {
                mindf = tmpabsdf[j];
                j_best_fa = j;
            }
        }

        if (j_best_fa == -1) {
            j_best_fa = rand() % tmpq;

            fa = U[j_best_fa];
            tmpabsdf[j_best_fa] = 0;
            for (c = 0; c < n; c++) {
                //calc contributation of c to tmpdf
                if (D[fa][c] < D[F[c]][c]) {
                    tmpabsdf[j_best_fa] += D[F[c]][c] - D[fa][c];
                }
            }

            mindf = tmpabsdf[j_best_fa];
        }

        //update
        ret_val = mindf;

        best_fa = U[j_best_fa];
        P[pp] = best_fa;
        pp++;
        U[j_best_fa] = U[tmpq - 1];
        tmpq--;
        x[best_fa] = true;

        nC[best_fa] = 0; //CRUCIAL

        for (c = 0; c < n; c++) {
            if (Ninverse[c][best_fa] < Ninverse[c][F[c]]) { //NOTE thant we CANNOT use if (D[best_fa][c] < D[F[c]][c])
                //remove from old list
                oldfa = F[c];
                oldinx = Cinx[c];
                tmpVlast = C[oldfa][nC[oldfa] - 1];
                C[oldfa][oldinx] = tmpVlast;
                Cinx[tmpVlast] = oldinx;
                nC[oldfa]--;

                //add to new list
                F[c] = best_fa;
                newinx = nC[best_fa];
                C[best_fa][newinx] = c;
                Cinx[c] = newinx;
                nC[best_fa]++;
            }
        }

        return ret_val;
    }//openfa
    //----------

    static long closefa(int P[], int pp, int U[], bool x[], int F[], int* C[MAX_N], int Cinx[], int nC[], double alpha)
    {
        int i, k, i_best_fa, fa, c, best_fa, tmpnc, mind, tmpminfa, nextminfa, tmpnextpos, tmpinx;
        long maxdf;
        long tmpabsdf[n];
        
        int tmpq = n - pp;
        long ret_val;
        
        i_best_fa = -1;
        maxdf = -1;
        for (i = 0; i < pp; i++) {
            if ((rand()*1.0)/RAND_MAX > alpha) {
                continue;
            }
            fa = P[i];
            tmpnc = nC[fa];
            //calc the amount of change in f if fa dropped
            tmpabsdf[i] = 0;
            for (k = 0; k < tmpnc; k++) {
                c = C[fa][k];
                
                tmpnextpos = nextpos(c, Ninverse[c][F[c]]);
                nextminfa = N[c][tmpnextpos];                
                mind = D[nextminfa][c];
                
                tmpabsdf[i] += mind - D[fa][c];
            }
            if (tmpabsdf[i] > maxdf) {
                maxdf = tmpabsdf[i];
                i_best_fa = i;
            }
        }

        if (i_best_fa == -1) {            
            i_best_fa = rand() % pp;
            fa = P[i_best_fa];
            tmpnc = nC[fa];
            //calculate the amount of change in f if fa dropped
            tmpabsdf[i_best_fa] = 0;
            for (k = 0; k < tmpnc; k++) {
                c = C[fa][k];

                tmpnextpos = nextpos(c, Ninverse[c][F[c]]);
                nextminfa = N[c][tmpnextpos];                
                mind = D[nextminfa][c];

                tmpabsdf[i_best_fa] += mind - D[fa][c];
            }           
            
            maxdf = tmpabsdf[i_best_fa];
        }

        //update info
        ret_val = maxdf;

        best_fa = P[i_best_fa];
        U[tmpq] = best_fa;
        tmpq++;

        P[i_best_fa] = P[pp-1];
        pp--;
        x[best_fa] = false;

        //update F and C as well
        tmpnc = nC[best_fa];
        for (k = 0; k < tmpnc; k++) {
            c = C[best_fa][k];

            tmpnextpos = nextpos(c, Ninverse[c][F[c]]);
            tmpminfa = N[c][tmpnextpos];
            
            F[c] = tmpminfa;
            tmpinx = nC[tmpminfa];
            C[tmpminfa][tmpinx] = c;
            Cinx[c] = tmpinx;
            nC[tmpminfa]++;
        }
        nC[best_fa] = 0; //not required, just for more reliability
        
        return ret_val;
    }//closefa
    //----------
    
    
    static void quicksort(int c, int from, int to)
    {
        int i, j, tmp, pivotitem;
        
        if(from < to)
        {
            pivotitem = D[N[c][from]][c];
            j = from;
            for(i = from + 1; i <= to; i++)
            {
                if(D[N[c][i]][c] < pivotitem)
                {
                    j++;
                    tmp = N[c][i];
                    N[c][i] = N[c][j];
                    N[c][j] = tmp;
                }
            }

            tmp = N[c][j];
            N[c][j] = N[c][from];
            N[c][from] = tmp;
            
            quicksort(c, from, j - 1);
            quicksort(c, j + 1, to);
        }
    }//quicksort
    //----------


    static void fillin_N_and_Ninverse()
    {
        int j, c;
        
        for(c=0; c < n; c++)
        {
            for(j=0; j < n; j++)
                N[c][j] = j;
            quicksort(c, 0, n-1);
        }
        
        //now Ninverse
        for(c=0; c < n; c++)
            for(j = 0; j < n; j++)
                Ninverse[c][N[c][j]] = j;
    }//fillin_N_and_Ninverse
    //----------
    

    
static int readParams(string filename, int pmeds[], int ps[], string abs[]) //popuates pmed, ps, and abs
{
    string tmprow;
    int t = 0;

    try {
        ifstream infile(filename);
        getline(infile, tmprow);

        vector<string> paramnp;
        while (getline(infile, tmprow)) {
            paramnp = split(trim(tmprow), ' ');
            pmeds[t] = stoi(trim(paramnp[0]));
            ps[t] = stoi(trim(paramnp[1]));
            abs[t] = trim(paramnp[2]);
            t++;
        }
        infile.close();
    }    catch (exception e) {
        ee("could not read the parameters table");
    }

    return t;
}//readParams
//----------

    static void read_n_and_D(string theRawFilename) {
        int i, j, posFrom, posEnd, tmpM, tmpP, c, f;
        string tmpS, tmpRow;
        vector<string> tmpVals;
        float tmpfloat;

        try {
        ifstream input_file(theRawFilename);
        if (!input_file.is_open()) {
            ee("Could not open the input file: " + theRawFilename);
            exit(EXIT_FAILURE);
        }
        tmpS = string((std::istreambuf_iterator<char>(input_file)), std::istreambuf_iterator<char>());
        
            
            posFrom = tmpS.find("n=") + 2;
            posEnd = tmpS.find("m=") - 1;
            n = stoi(trim(tmpS.substr(posFrom, posEnd - posFrom)));

            posFrom = tmpS.find("m=") + 2;
            posEnd = tmpS.find("p=") - 1;
            tmpM = stoi(trim(tmpS.substr(posFrom, posEnd - posFrom)));
            if (tmpM != n) {
                ee("error: inconsistent n");
            }

            posFrom = tmpS.find("p=") + 2;
            posEnd = tmpS.find("client") - 1;
            tmpP = stoi(trim(tmpS.substr(posFrom, posEnd - posFrom)));
            if (tmpP != p) {
                ee("error: inconsistent p");
            }

            for (i = 0; i < n; i++) {
                for (j = 0; j < n; j++) {
                    D[i][j] = -1;
                }
            }

            posFrom = tmpS.find("table = ");
            posEnd = tmpS.length() - 1;
            while (tmpS[posEnd] != '}') {
                posEnd--;
            }
            tmpS = trim(tmpS.substr(posFrom + 10, posEnd-posFrom));
            for (f = 0; f < n; f++) {
                posFrom = tmpS.find("{");
                posEnd = tmpS.find("}");
                tmpRow = trim(tmpS.substr(posFrom + 1, posEnd-posFrom));
                tmpVals = split(tmpRow, ',');
                if (f < n - 1) {
                    tmpS = trim(tmpS.substr(posEnd + 2, tmpS.length() - (posEnd+1)));
                }
                for (c = 0; c < n; c++) {
                    tmpfloat = stof ( trim(tmpVals[c]) );
                    D[f][c] = (int) tmpfloat;
                    if (tmpfloat != D[f][c]) {
                        ee("error: inconsistent weight");
                    }
                }
            }
        } catch (exception e) {
            ee("error reading input file");
        }
    }//read_n_and_D
    //----------

    static void prepareDandNandP(string theRawFilename)
    {
        int i, j, k, w, tmp_inf;
        int tmpN = 0, tmpM = 0, tmpP = 0; //dummy init
        //read the graph      
        try {
            //first read the number of nodes and edges:
            ifstream infile(theRawFilename);
            string tmpFirstLine;
            getline(infile, tmpFirstLine);
            tmpFirstLine = trim(tmpFirstLine);
            vector<string> tmpParameters;
            tmpParameters = split(tmpFirstLine, ' ');
            tmpN = stoi(tmpParameters.at(0));
            tmpM = stoi(tmpParameters.at(1));
            tmpP = stoi(tmpParameters.at(2));

            tmp_inf = MY_INF; //NOTE this must be greater than the distance between the farthest nodes in the graph
            for (i = 0; i < tmpN; i++) {
                D[i][i] = 0;
                for (j = i + 1; j < tmpN; j++) {
                    D[i][j] = D[j][i] = tmp_inf;
                }
            }
            
            int tmpCntM = 0;
            string e;
            vector<string> ijw;
            while (getline(infile, e)) {
                ijw = split(trim(e), ' ');
                i = stoi(trim(ijw.at(0)));
                j = stoi(trim(ijw.at(1)));
                w = stoi(trim(ijw.at(2)));
                if ((i == j) || (w <= 0) || (i > tmpN) || (j > tmpN)) {
                    ee("error: inconsistency in the input file!");
                } else {
                    D[i - 1][j - 1] = D[j - 1][i - 1] = w;
                    //cout<< D[i - 1][j - 1]<<endl;
                    tmpCntM++;
                }
            }
            infile.close();
            if(tmpCntM != tmpM)
                ee("tmpCntM != tmpM");
        } catch (exception e) {
            ee("error reading input file!");
        }

        //run Floyd
        for (k = 0; k < tmpN; k++) {
            for (i = 0; i < tmpN; i++) {
                for (j = 0; j < tmpN; j++) {
                    if (D[i][k] + D[k][j] < D[i][j]) {
                        D[i][j] = D[i][k] + D[k][j];
                    }
                }
            }
        }

        n = tmpN;
        p = tmpP;
    }//prepareDandNandP
    //----------    
    
    

int main(int argc, char** argv) {
        string path = "C:\\Users\\ad0204\\Desktop\\OpM_LIB_2016\\pmeds\\"; //NOTE: Please change this path so it points to the folder (directory) on your machine where the data files are located.
        string filename;
        int t; 
        int F[MAX_N];//this is the list fi in the paper (for client's facility)
        int* C[MAX_N];//this is the list of cj in the paper (for facility's clients)
        int Cinx[MAX_N]; //to indicate the index of each client in C
        int nC[MAX_N]; 
        for(int i=0; i < MAX_N; i++)
        {
            D[i] = (int*)malloc(MAX_N*sizeof(int));
            N[i] = (int*)malloc(MAX_N*sizeof(int));
            Ninverse[i] = (int*)malloc(MAX_N*sizeof(int));
            C[i] = (int*)malloc(MAX_N*sizeof(int));
        }
        
        
        int P[MAX_N], U[MAX_N], bestP[MAX_N]; //lists of centres and non-centres vertices

        int i, d, radius, f, worstf, bestf, deltaf;
        double alpha, gap;

        bool flgImproved;
        long start_time;
        srand(time(NULL));
        double g1 = g(1);
        
        for (t = 20; t < 40; t++){
            filename = path + "pmed" + to_string(t + 1) + ".txt"; // e.g. pmed1.txt
            prepareDandNandP(filename);

            for (int tmpk = 0; tmpk <= 1; tmpk++)
            {   // overwrite p with n/4 or n/3
                if (tmpk == 0)
                {
                    p = n / 4;
                }
                else
                {
                    p = n / 3;
                }
                cout << "pmed" << t+1 <<", new p = " << p << "\t" << flush;

                start_time = clock();
                q = n - p;

                fillin_N_and_Ninverse();

                f = initRandSol(P, U, x, F, C, Cinx, nC);

                for(i=0; i < p; i++)
                    bestP[i] = P[i];
                worstf = bestf = f;
                radius = 1;
                alpha = g1;
                int duration = 10000; //for 10 seconds per instance

                while ( clock() - start_time < duration ){ //please adjust the termination condition based on your own requirements
                    //RLS1     
                    flgImproved = true;
                    while (flgImproved && clock() - start_time < duration){
                        flgImproved = false;
                        
                        deltaf = 0;
                        for (d = 0; d < radius; d++) {
                            deltaf += closefa(P, p, U, x, F, C, Cinx, nC, alpha);
                            p--;
                            q++;
                        }//for d

                        for (d = 0; d < radius; d++) {
                            deltaf -= openfa(P, p, U, x, F, C, Cinx, nC, alpha);
                            p++;
                            q--;
                        }//for d

                        f += deltaf;

                        if (deltaf > 0) {
                            flgImproved = true;
                            radius = 1;
                            alpha = g1;
                            if (f > bestf) {
                                for(i= 0; i < p; i++)
                                    bestP[i] = P[i];
                                bestf = f;
                            }
                        }
                        if (f < worstf) {
                            worstf = f;
                        }
                    }//while RLS1
                    
                    //RLS2
                    flgImproved = true;
                    while (flgImproved && clock() - start_time < duration){
                        flgImproved = false;
                        
                        deltaf = 0;
                        for (d = 0; d < radius; d++) {
                            deltaf -= openfa(P, p, U, x, F, C, Cinx, nC, alpha);
                            p++;
                            q--;
                        }
                        for (d = 0; d < radius; d++) {
                            deltaf += closefa(P, p, U, x, F, C, Cinx, nC, alpha);
                            p--;
                            q++;
                        }

                        f += deltaf;

                        if (deltaf > 0) {
                            flgImproved = true;
                            radius = 1;
                            alpha = g1;
                            if (f > bestf) {
                                for(i=0; i < p; i++)
                                    bestP[i] = P[i];
                                bestf = f;
                            }
                        }
                        if (f < worstf) {
                            worstf = f;
                        }
                    }//while RLS2
                    
                    if (worstf == bestf) {
                        gap = 0;
                    } else {
                        gap = ((double) bestf - f) / (bestf - worstf);
                    }
                    radius += (int) ((p - 2) * gap + 1);
                    if (radius > p - 1) {
                        radius = p - 1;
                    }             
                    alpha = g(radius);
                }//while termination condition not satisfied

                if (bestf != calcf(bestP, p)) {// to double-check the correctness of the algorithm
                    ee("The final answer is incorrect!");
                }

                cout << bestf << "\t" << (clock() - start_time)/1000.0 << endl << flush; 
                
            }//for p 
        }//for pmed file

    return 0;
}

