import java.io.*;
import java.util.ArrayList;


/**
 * Created by gal on 3/27/2016.
 */
public class SequenceAlignment {
    ArrayList<String> sequences;
    TransitionMatrix matrix;
    int a;
    int b;



    public SequenceAlignment(String fasteFile,String matrixFile){
        sequences = new ArrayList<>();
        FastaSequence reader = new FastaSequence(fasteFile);
        for (int i=0; i< reader.size(); i++)
        {
            sequences.add(reader.getSequence(i));
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
        //full matrix with scores and pies
        for (int i = 1; i <= sequenceA.length(); i++) {
            for (int j = 1; j <= sequenceB.length(); j++) {
                double maxScore = Math.max(Math.max(M[i-1][j-1].getScore(), M[i-1][j].getScore()) , M[i][j-1].getScore());
                M[i][j].setScore(maxScore);
                if(maxScore == M[i-1][j-1].getScore())
                    M[i][j].setPi(M[i-1][j-1]);
                else if(maxScore == M[i-1][j].getScore())
                    M[i][j].setPi(M[i-1][j]);
                else M[i][j].setPi(M[i][j-1]);
            }
        }

        //get the end of path cell
        double maxScore = 0;
        int maxi;
        int maxj;
        for (int i = 1; i <= sequenceA.length(); i++)
            if (i >= maxScore){
                maxScore = i;
                maxi = i;
                maxj = 0;
            }
        for (int j = 1; j <= sequenceB.length(); j++)
            if (j >= maxScore){
                maxScore = j;
                maxi = 0;
                maxj = j;
            }

        // get best score path

        //CHECK
        System.out.println(M);


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
        double highestScore;

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
        double best = 0;
        int iBestCeller = 0;
        int jBestCeller = 0;
        double secondBest = 0;
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
        findPath(iBestCeller,jBestCeller,M,sequenceA,sequenceB);
    }

    private void findPath(int i, int j, cellMatrix[][] M, String sequenceA, String sequenceB) {
        for(int t = 0; t < M.length; t++){
            System.out.println();
            for(int k = 0; k < M[0].length; k++ ){
                System.out.print(M[t][k].getScore() + " ");
            }
        }
        System.out.println();
        //at least one of i,j should be not-zero otherwise M[0][0] =
        StringBuilder first = new StringBuilder();
        StringBuilder second = new StringBuilder();
        boolean found;
        if(i != 0 && j != 0){
            //initialize last cell
            cellMatrix currCell = M[i][j];
            while(currCell != null){
                found = false;
                System.out.println(i + " " + j + " " + M[i][j].getScore());
                //where did i came from?
                if(i > 0 && j > 0) {
                    //replace
                    if (M[i][j].getPI() == M[i - 1][j - 1]) {
                        first.append(sequenceA.charAt(i - 1));
                        second.append(sequenceB.charAt(j - 1));
                        i--;
                        j--;
                        found = true;
                    }
                }
                if(i > 0 && !found) {
                    //delete
                    if (M[i][j].getScore() == (M[i - 1][j].getScore() + matrix.score(sequenceA.charAt(i - 1), '*'))) {
                        first.append(sequenceA.charAt(i - 1));
                        second.append("_");
                        i--;
                        found = true;
                    }
                }
                if(j > 0 && !found) {
                    //insert
                    {
                        first.append("_");
                        second.append(sequenceB.charAt(j - 1));
                        j--;
                        found = true;
                    }
                }
                currCell = currCell.getPI();
            }
            System.out.println(first.reverse());
            System.out.println(second.reverse());
        }
    }
}



