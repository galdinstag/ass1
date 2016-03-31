import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;


/**
 * Created by gal on 3/27/2016.
 */
public class SequenceAlignment {
    HashMap<String,String> sequences1;
    HashMap<String,String> sequences2;
    TransitionMatrix matrix;
    int a;
    int b;



    public SequenceAlignment(String matrixFile,String fasteFile1,String fasteFile2){
        sequences1 = new HashMap<>();
        FastaSequence reader = new FastaSequence(fasteFile1);
        for (int i=0; i< reader.size(); i++)
        {
            sequences1.put(reader.getDescription(i), reader.getSequence(i));
        }
        sequences2 = new HashMap<>();
        reader = new FastaSequence(fasteFile2);
        for (int i=0; i< reader.size(); i++)
        {
            sequences2.put(reader.getDescription(i), reader.getSequence(i));
        }
        matrix = new TransitionMatrix(fixMatrix(matrixFile));
        new File("mid.txt").delete();
    }

    public File fixMatrix(String orgMatrix){
        try {
            BufferedReader bf = new BufferedReader(new FileReader(orgMatrix));
            BufferedWriter writer = new BufferedWriter(new FileWriter("mid.txt"));
            String line = bf.readLine();
            while((line != null) && line.contains("#")){    //dont read lines with #
                line = bf.readLine();
            }
            while((line != null) && (!line.contains("#"))){ //read matrix
                writer.write(line);
                writer.write("\n");
                writer.flush();
                line = bf.readLine();
            }
            while((line != null) && line.contains("#")){ //dont read lines with #
                line = bf.readLine();
            }
            //read A for: W(x) = Ax + B
            String num = line.substring(2);
            a = Integer.parseInt(num.toString());

            //read B for: W(x) = Ax + B
            line = bf.readLine();
            num = line.substring(2);
            b = Integer.parseInt(num.toString());

            writer.close();
            return new File("mid.txt");
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }

    //run global alignment on two strings
    //return the score of the highest rank + the path
    public void globalAlignment(String sequenceA, String sequenceB) {
        int maxi = 0;
        int maxj = 0;
        cellMatrix[][] M = new cellMatrix[sequenceA.length() + 1][sequenceB.length() + 1];
        //init new cells in matrix
        for (int i = 0; i <= sequenceA.length(); i++)
            for (int j = 0; j <= sequenceB.length(); j++)
                M[i][j] = new cellMatrix();

        //init zero in the left and up
        for (int i = 0; i <= sequenceA.length(); i++)
            M[i][0].setScore(0);
        for (int j = 0; j <= sequenceB.length(); j++)
            M[0][j].setScore(0);

        //mainloop:
        //full the matrix with scores and pies
        for (int i = 1; i <= sequenceA.length(); i++) {
            for (int j = 1; j <= sequenceB.length(); j++) {
                int maxScore = M[i - 1][j - 1].getScore() + matrix.score(sequenceA.charAt(i - 1), sequenceB.charAt(j - 1));
                M[i][j].setScore(maxScore);
                M[i][j].setPi(M[i - 1][j - 1]);

                if (M[i - 1][j].getScore() + matrix.score(sequenceA.charAt(i - 1), '*') > maxScore) {
                    maxScore = M[i - 1][j].getScore() + matrix.score(sequenceA.charAt(i - 1), '*');
                    M[i][j].setScore(maxScore);
                    M[i][j].setPi(M[i - 1][j]);
                }
                if (M[i][j - 1].getScore() + matrix.score(sequenceB.charAt(j - 1), '*') > maxScore) {
                    maxScore =M[i][j - 1].getScore() + matrix.score(sequenceB.charAt(j - 1), '*');
                    M[i][j].setScore(maxScore);
                    M[i][j].setPi(M[i][j - 1]);
                }
            }
        }
        //get the end of path cell
        int maxScore = M[sequenceA.length()][sequenceB.length()].getScore();

        for (int i = 1; i <= sequenceA.length(); i++)
            if (M[i][sequenceB.length()].getScore() >= maxScore){
                maxScore = M[i][sequenceB.length()].getScore();
                maxi = i;
                maxj = sequenceB.length();
            }
        for (int j = 1; j <= sequenceB.length(); j++)
            if (M[sequenceA.length()][j].getScore() >= maxScore){
                maxScore = M[sequenceA.length()][j].getScore();
                maxi = sequenceA.length();
                maxj = j;
            }

        // get best score path
        findPath(maxi, maxj, M, sequenceA, sequenceB);

    }

    //run local alignment on two strings
    //return the score of the two highest ranks + the path
    public void localAlignment(String sequenceA, String sequenceB) {
        cellMatrix[][] M = new cellMatrix[sequenceA.length()+1][sequenceB.length()+1];
        for (int i = 0; i <= sequenceA.length(); i++) {
            for (int j = 0; j <= sequenceB.length(); j++) {
                M[i][j] = new cellMatrix();
            }
        }
        int highestScore;

        //main loop
        for (int i = 1; i <= sequenceA.length(); i++) {
            for (int j = 1; j <= sequenceB.length(); j++) {
                highestScore = 0;
                if (M[i - 1][j - 1].getScore() + matrix.score(sequenceA.charAt(i - 1), sequenceB.charAt(j - 1)) > highestScore) {
                    M[i][j].setScore(M[i - 1][j - 1].getScore() + matrix.score(sequenceA.charAt(i - 1), sequenceB.charAt(j - 1)));
                    M[i][j].setPi(M[i - 1][j - 1]);
                    highestScore = M[i - 1][j - 1].getScore() + matrix.score(sequenceA.charAt(i - 1), sequenceB.charAt(j - 1));
                }
                if (M[i - 1][j].getScore() + matrix.score(sequenceA.charAt(i - 1), '*') > highestScore) {
                    M[i][j].setScore(M[i - 1][j].getScore() + matrix.score(sequenceA.charAt(i - 1), '*'));
                    M[i][j].setPi(M[i - 1][j]);
                    highestScore = M[i - 1][j].getScore() + matrix.score(sequenceA.charAt(i - 1), '*');
                }
                if (M[i][j - 1].getScore() + matrix.score(sequenceB.charAt(j - 1), '*') > highestScore) {
                    M[i][j].setScore(M[i][j - 1].getScore() + matrix.score(sequenceB.charAt(j - 1), '*'));
                    M[i][j].setPi(M[i][j - 1]);
                }
            }
        }
        //search the matrix for the 2 largest scores, return them and the paths.
        int best = 0;
        int iBestCeller = 0;
        int jBestCeller = 0;
        int secondBest = 0;
        int iSecondBestCeller = 0;
        int jSecondBestCeller = 0;
        for(int i = 1; i <= sequenceA.length(); i++) {
            for (int j = 1; j <= sequenceB.length(); j++) {
                //second best score so far
                if (M[i][j].getScore() >= secondBest) {
                    secondBest = M[i][j].getScore();
                    iSecondBestCeller = i;
                    jSecondBestCeller = j;
                }
                //best so far, swap best and second best
                if (M[i][j].getScore() > best) {
                    secondBest = best;
                    iSecondBestCeller = iBestCeller;
                    jSecondBestCeller = jBestCeller;
                    best = M[i][j].getScore();
                    iBestCeller = i;
                    jBestCeller = j;
                }
            }
        }
        System.out.println("best local alignment:");
        findPath(iBestCeller,jBestCeller,M,sequenceA,sequenceB);
        System.out.println("second best local alignment:");
        findPath(iSecondBestCeller, jSecondBestCeller, M, sequenceA, sequenceB);
    }

    private void findPath(int i, int j, cellMatrix[][] M, String sequenceA, String sequenceB) {
        int score = M[i][j].getScore();
//        for(int t = 0; t < M.length; t++){
//            System.out.println();
//            for(int k = 0; k < M[0].length; k++){
//                System.out.print(M[t][k].getScore() + "  ");
//            }
//        }
//
//        System.out.println();
        //at least one of i,j should be not-zero otherwise M[0][0] =
        StringBuilder first = new StringBuilder();
        StringBuilder second = new StringBuilder();
        boolean found;
        if(i != 0 && j != 0){
            //initialize last cell
            cellMatrix currCell = M[i][j];
            while(currCell != null && i > 0 && j > 0) {
                found = false;
                //where did i came from?
                //replace
                if (M[i][j].getPI() == M[i - 1][j - 1]) {
                    first.append(sequenceA.charAt(i - 1));
                    second.append(sequenceB.charAt(j - 1));
                    i--;
                    j--;
                    found = true;
                }
                if (!found) {
                    //delete
                    if (M[i][j].getPI() == M[i - 1][j]) {
                        first.append(sequenceA.charAt(i - 1));
                        second.append("_");
                        i--;
                        found = true;
                    }
                }
                if (!found) {
                    //insert
                    {
                        if (M[i][j].getPI() == M[i][j - 1]) {
                            first.append("_");
                            second.append(sequenceB.charAt(j - 1));
                            j--;
                            found = true;
                        }
                    }
                }
                currCell = currCell.getPI();
            }
            System.out.println(first.reverse());
            System.out.println(second.reverse());
            System.out.println("Score: " + score);
        }
    }
    public void gapGlobalAlignment(String sequenceA, String sequenceB) {
        //MAX VALUE MATRIX
        cellMatrix[][] M = new cellMatrix[sequenceA.length() + 1][sequenceB.length() + 1];
        //init new cells in matrix M
        for (int i = 0; i <= sequenceA.length(); i++)
            for (int j = 0; j <= sequenceB.length(); j++)
                M[i][j] = new cellMatrix();
        //INSERSION MATRIX
        cellMatrix[][] IS = new cellMatrix[sequenceA.length() + 1][sequenceB.length() + 1];
        //init new cells in matrix S
        for (int i = 0; i <= sequenceA.length(); i++)
            for (int j = 0; j <= sequenceB.length(); j++)
                IS[i][j] = new cellMatrix();
        //DELETION MATRIX
        cellMatrix[][] IT = new cellMatrix[sequenceA.length() + 1][sequenceB.length() + 1];
        //init new cells in matrix T
        for (int i = 0; i <= sequenceA.length(); i++)
            for (int j = 0; j <= sequenceB.length(); j++)
                IT[i][j] = new cellMatrix();
        //mainloop: full all matrixes
        int max;
        for (int i = 1; i <= sequenceA.length(); i++)
            for (int j = 1; j <= sequenceB.length(); j++) {
                // find max M
                max = M[i - 1][j - 1].getScore();
                M[i][j].setScore(max);
                M[i][j].setPi(M[i - 1][j - 1]);

                if (IS[i - 1][j - 1].getScore() > max) {
                    max = IS[i - 1][j - 1].getScore();
                    M[i][j].setScore(max);
                    M[i][j].setPi(IS[i - 1][j - 1]);
                }
                if (IT[i - 1][j - 1].getScore() > max) {
                    max = IT[i - 1][j - 1].getScore();
                    M[i][j].setScore(max);
                    M[i][j].setPi(IT[i - 1][j - 1]);
                }
                M[i][j].setScore(max + matrix.score(sequenceA.charAt(i - 1), sequenceB.charAt(j - 1)));

                //find max IS
                max = IS[i - 1][j].getScore() - b;
                IS[i][j].setScore(max);
                IS[i][j].setPi(IS[i - 1][j]);

                if ((IT[i - 1][j].getScore() - a) > max) {
                    max = IT[i - 1][j - 1].getScore() - a;
                    IS[i][j].setScore(max);
                    IS[i][j].setPi(IT[i - 1][j]);
                }
                if ((M[i - 1][j].getScore() - a) > max) {
                    max = M[i - 1][j].getScore() - a;
                    IS[i][j].setScore(max);
                    IS[i][j].setPi(M[i - 1][j]);
                }
                //find max IT
                max = IT[i][j - 1].getScore() - b;
                IT[i][j].setScore(max);
                IT[i][j].setPi(IS[i - 1][j]);

                if ((IS[i][j - 1].getScore() - a) > max) {
                    max = IS[i][j - 1].getScore() - a;
                    IT[i][j].setScore(max);
                    IT[i][j].setPi(IS[i][j - 1]);
                }
                if ((M[i][j - 1].getScore() - a) > max) {
                    max = M[i][j - 1].getScore() - a;
                    IT[i][j].setScore(max);
                    IT[i][j].setPi(M[i][j - 1]);
                }
            }
        //get the end of path cell
        int maxScore = M[sequenceA.length()][sequenceB.length()].getScore();
        int maxi = sequenceA.length();
        int maxj = sequenceB.length();

        for (int i = 1; i <= sequenceA.length(); i++) {
            if (M[i][sequenceB.length()].getScore() >= maxScore) {
                maxScore = M[i][sequenceB.length()].getScore();
                maxi = i;
                maxj = sequenceB.length();
            }
            if (IS[i][sequenceB.length()].getScore() >= maxScore) {
                maxScore = IS[i][sequenceB.length()].getScore();
                maxi = i;
                maxj = sequenceB.length();
            }
            if (IT[i][sequenceB.length()].getScore() >= maxScore) {
                maxScore = IT[i][sequenceB.length()].getScore();
                maxi = i;
                maxj = sequenceB.length();
            }
        }
        for (int j = 1; j <= sequenceB.length(); j++){
            if (M[sequenceA.length()][j].getScore() >= maxScore) {
                maxScore = M[sequenceA.length()][j].getScore();
                maxi = sequenceA.length();
                maxj = j;
            }
            if (IS[sequenceA.length()][j].getScore() >= maxScore) {
                maxScore = IS[sequenceA.length()][j].getScore();
                maxi = sequenceA.length();
                maxj = j;
            }
            if (IT[sequenceA.length()][j].getScore() >= maxScore) {
                maxScore = IT[sequenceA.length()][j].getScore();
                maxi = sequenceA.length();
                maxj = j;
            }
        }
        findAffinePath(maxi,maxj,M,IS,IT,sequenceA,sequenceB);
    }
    private void findAffinePath(int i, int j, cellMatrix[][] M, cellMatrix[][] IS, cellMatrix[][] IT, String sequenceA, String sequenceB) {
        int score = M[i][j].getScore();
        cellMatrix currCell = M[i][j];
        if(IS[i][j].getScore() > score){
            currCell = IS[i][j];
        }
        if(IT[i][j].getScore() > score){
            currCell = IT[i][j];
        }
//        for(int t = 0; t < M.length; t++){
//            System.out.println();
//            for(int k = 0; k < M[0].length; k++){
//                System.out.print(M[t][k].getScore() + "  ");
//            }
//        }
//
//        System.out.println();
        //at least one of i,j should be not-zero otherwise M[0][0] =
        StringBuilder first = new StringBuilder();
        StringBuilder second = new StringBuilder();
        boolean found;
        if(i != 0 && j != 0){
            //initialize last cell
            while(currCell != null && i > 0 && j > 0) {
                found = false;
                //where did i came from?
                //replace
                if (currCell.getPI() == M[i - 1][j - 1] || currCell.getPI() == IS[i - 1][j - 1] || currCell.getPI() == IT[i - 1][j - 1]) {
                    first.append(sequenceA.charAt(i - 1));
                    second.append(sequenceB.charAt(j - 1));
                    i--;
                    j--;
                    found = true;
                }
                if (!found) {
                    //delete
                    if (currCell.getPI() == M[i - 1][j] || currCell.getPI() == IS[i - 1][j] || currCell.getPI() == IT[i - 1][j]) {
                        first.append(sequenceA.charAt(i - 1));
                        second.append("_");
                        i--;
                        found = true;
                    }
                }
                if (!found) {
                    //insert
                    {
                        if (currCell.getPI() == M[i][j - 1] || currCell.getPI() == IS[i][j - 1] || currCell.getPI() == IT[i][j - 1]) {
                            first.append("_");
                            second.append(sequenceB.charAt(j - 1));
                            j--;
                        }
                    }
                }
                currCell = currCell.getPI();
            }
            System.out.println(first.reverse());
            System.out.println(second.reverse());
            System.out.println("Score: " + score);
        }
    }
}