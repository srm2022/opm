
# This is a generic Java implementation of IGV. Please customise it as per 
# your own requirements, by changing its termination condition, reported 
# performance indicators, etc. We customised it differently for different 
# sections of the paper. Please also make sure the path points to your data 
# folder. To keep it simple, it only runs once, for 10 seconds, on each 
# instance of the large dataset and prints the best objective value and the 
# duration time (which is around 10 seconds, excluding the time spent on 
# reading input data). 
# Comments and questions are welcome to seyed_r_mousavi@yahoo.com and new 
# versions (if any) will be made available at https://github.com/srm2022/opm, 
# as long as these accounts remain active. Thank you.

import java.io.File;
import java.io.FileNotFoundException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Scanner;

public class IGV_generic {

    static int MY_INF = Integer.MAX_VALUE;
    static int n, p, q;
    static int[][] D;
    static int[][] N;
    static int[][] Ninverse;
    static boolean[] x; //list of indicator variables xj which is true iff facility j is in P

    static double gamma = 0.6;
    static double tau = 6;

    public static void main(String[] args) {
        String path = "C:\\Users\\ad0204\\Desktop\\OpM_LIB_2016\\pmeds\\"; //NOTE: Please change this path so it points to the folder (directory) on your machine where the data files are located.
        String filename;
        int t;
        int[] F;//this is the list fi in the paper (for client's facility)
        int[][] C;//this is the list of cj in the paper (for facility's clients)
        int[] Cinx; //to indicate the index of each client in C
        int[] nC;

        int[] P, U, bestP; //lists of centres and non-centres vertices

        int i, d, radius, f, worstf, bestf, deltaf;
        double alpha, gap;

        boolean flgImproved;
        long start_time;
        double g1 = g(1);

        for (t = 20; t < 40; t++) { // from 0 to 40, for the whole pmed benchmark
            filename = path + "pmed" + (t + 1) + ".txt"; //e.g. pmed1.txt
            prepareDandNandP(filename);

            for (int tmpk = 0; tmpk <=1; tmpk++) { // overwrite p with n/4 or n/3
                if (tmpk == 0) {
                    p = n/4 ;
                }else{
                    p =n/3;
                }
                System.out.print("pmed" + (t+1) + ", new p = " + p + "\t");
                
                start_time = System.currentTimeMillis();
                q = n - p;
                P = new int[n];
                bestP = new int[n];
                U = new int[n];

                N = new int[n][n];
                Ninverse = new int[n][n];
                fillin_N_and_Ninverse();

                x = new boolean[n];
                F = new int[n];
                C = new int[n][n];
                Cinx = new int[n];
                nC = new int[n];

                f = initRandSol(P, U, x, F, C, Cinx, nC);

                for (i = 0; i < p; i++) {
                    bestP[i] = P[i];
                }
                worstf = bestf = f;
                radius = 1;
                alpha = g1;
                int duration = 10000; //for 10 seconds per instance- please adjust the termination condition as you require
                
                while (System.currentTimeMillis() - start_time < duration) { //please adjust the termination condition based on your own requirements
                    //RLS1     
                    flgImproved = true;
                    while (flgImproved && System.currentTimeMillis() - start_time < duration) {
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
                                for (i = 0; i < p; i++) {
                                    bestP[i] = P[i];
                                }
                                bestf = f;
                            }
                        }
                        if (f < worstf) {
                            worstf = f;
                        }
                    }//while RLS1

                    //RLS2
                    flgImproved = true;
                    while (flgImproved && System.currentTimeMillis() - start_time < duration) {
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
                                for (i = 0; i < p; i++) {
                                    bestP[i] = P[i];
                                }
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

                if (bestf != calcf(bestP, p)) {// to double-ccheck the correctness of the algorithm
                    ee("The final answer is incorrect!");
                }

                System.out.println(bestf + "\t" + (System.currentTimeMillis() - start_time) / 1000.0); 
            }//for p
        }//for pmed file
    }//main
    //----------

    static double g(int radius) {
        return gamma * Math.pow(2, tau * (1 - radius)) + (1 - gamma) * Math.random();
    }//g
    //----------

    static int calcf(int[] pP, int pp) {
        int i, c, tmpmin;
        int ret_val;

        ret_val = 0;
        for (c = 0; c < n; c++) {
            tmpmin = MY_INF;
            for (i = 0; i < pp; i++) {
                if (D[pP[i]][c] < tmpmin) {
                    tmpmin = D[pP[i]][c];
                }
            }
            ret_val += tmpmin;
        }
        return ret_val;
    }//calcf
    //----------

    static int initRandSol(int[] P, int[] U, boolean[] x, int[] F, int[][] C, int[] Cinx, int[] nC) {
        int k, fa, tmpf, i, c, minpos, minfa, mind, tmpinx;
        double tmpPbyN = p * 1.0 / n;
        for (fa = 0; fa < n; fa++) {
            x[fa] = false;
        }

        k = 0;
        fa = (int) (Math.random() * n); //random starting point
        while (k < p) {
            if (!x[fa]) {
                if (Math.random() < tmpPbyN) {
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

    static long openfa(int[] P, int pp, int[] U, boolean[] x, int[] F, int[][] C, int[] Cinx, int[] nC, double alpha) {
        int j, j_best_fa, fa, c, best_fa, oldfa, oldinx, tmpVlast, newinx;
        long mindf;
        long[] tmpabsdf = new long[n];

        int tmpq = n - pp;
        long ret_val;

        mindf = MY_INF;
        j_best_fa = -1;
        for (j = 0; j < tmpq; j++) {
            if (Math.random() > alpha) {
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
            j_best_fa = (int) (Math.random() * tmpq);

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

    static long closefa(int[] P, int pp, int[] U, boolean[] x, int[] F, int[][] C, int[] Cinx, int[] nC, double alpha) {
        int i, k, i_best_fa, fa, c, best_fa, tmpnc, mind, tmpminfa, nextminfa, tmpnextpos, tmpinx;
        long maxdf;
        long[] tmpabsdf = new long[n];

        int tmpq = n - pp;
        long ret_val;

        i_best_fa = -1;
        maxdf = -1;
        for (i = 0; i < pp; i++) {
            if (Math.random() > alpha) {
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
            i_best_fa = (int) (Math.random() * pp);
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

        P[i_best_fa] = P[pp - 1];
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

    static int nextpos(int c, int curpos) {
        int pos = 1 + curpos;
        while (!x[N[c][pos]]) {
            pos++;
        }
        return pos;
    }//nextpos
    //----------

    static void fillin_N_and_Ninverse() {
        int j, c;

        for (c = 0; c < n; c++) {
            for (j = 0; j < n; j++) {
                N[c][j] = j;
            }
            quicksort(c, 0, n - 1);
        }

        //now Ninverse
        for (c = 0; c < n; c++) {
            for (j = 0; j < n; j++) {
                Ninverse[c][N[c][j]] = j;
            }
        }
    }//fillin_N_and_Ninverse
    //----------

    static void quicksort(int c, int from, int to) {
        int i, j, tmp, pivotitem;

        if (from < to) {
            pivotitem = D[N[c][from]][c];
            j = from;
            for (i = from + 1; i <= to; i++) {
                if (D[N[c][i]][c] < pivotitem) {
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

    static void prepareDandNandP(String theRawFilename) {
        int i, j, k, w, tmp_inf;
        int tmpN = 0, tmpM = 0, tmpP = 0;//dummy init
        //read in the graph      
        try {
            File infile = new File(theRawFilename);
            Scanner myReader = new Scanner(infile);
            String tmpFirstLine = myReader.nextLine().trim();
            String[] tmpParameters = new String[3];
            tmpParameters = tmpFirstLine.split(" ");
            tmpN = Integer.valueOf(tmpParameters[0]);
            tmpM = Integer.valueOf(tmpParameters[1]);
            tmpP = Integer.valueOf(tmpParameters[2]);

            tmp_inf = Integer.MAX_VALUE / 2; //for infinity; devided by 2 to avoid overflow when adding two of them in Floyd alg.
            D = new int[tmpN][tmpN];
            for (i = 0; i < tmpN; i++) {
                D[i][i] = 0;
                for (j = i + 1; j < tmpN; j++) {
                    D[i][j] = D[j][i] = tmp_inf;
                }
            }

            int tmpCntM = 0;
            String e;
            String[] ijw;
            while (myReader.hasNextLine()) {
                e = myReader.nextLine().trim();
                //System.out.println(e);
                ijw = e.split(" ");
                i = Integer.valueOf(ijw[0]);
                j = Integer.valueOf(ijw[1]);
                w = Integer.valueOf(ijw[2]);
                if ((i == j) || (w <= 0) || (i > tmpN) || (j > tmpN)) {
                    ee("error: inconsistency in the input file");
                } else {
                    D[i - 1][j - 1] = D[j - 1][i - 1] = w;
                    tmpCntM++;
                }
            }
            myReader.close();
            if (tmpM != tmpCntM) {
                ee("error: inconsistent m");
            }
        } catch (FileNotFoundException e) {
            ee("error reading input file");
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

    static void ee(String message) {
        System.out.println(message);
        System.exit(1);
    }//ee
    //----------

}//the class
