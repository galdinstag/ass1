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


    public SequenceAlignment(String fasteFile, String matrixFile) {
        sequences = new ArrayList<>();
        FastaSequence reader = new FastaSequence(fasteFile);
        for (int i = 0; i < reader.size(); i++) {
            sequences.add(reader.getSequence(i));
        }
        matrix = new TransitionMatrix(fixMatrix(matrixFile));
        new File("mid.txt").delete();
    }

    public void GlobalAlignment() {


    }

    public File fixMatrix(String orgMatrix) {
        try {
            BufferedReader bf = new BufferedReader(new FileReader(orgMatrix));
            BufferedWriter writer = new BufferedWriter(new FileWriter("mid.txt"));
            String line = bf.readLine();
            while ((line != null) && line.contains("#")) {    //dont read lines with #
                line = bf.readLine();
            }
            while ((line != null) && (!line.contains("#"))) { //read matrix
                writer.write(line);
                writer.write("\n");
                writer.flush();
                line = bf.readLine();
            }
            while ((line != null) && line.contains("#")) { //dont read lines with #
                line = bf.readLine();
            }
            //read A for: W(x) = Ax + B
            String num = line.substring(2);
            a = Integer.parseInt(num);

            //read B for: W(x) = Ax + B
            line = bf.readLine();
            num = line.substring(2);
            b = Integer.parseInt(num);
            
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
            for (int j = 0; i<sequenceB.length(); j++)
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



    }


    //run local alignment on two strings
    //return the score of the two highest ranks + the paths
    public void localAlignment(String sequenceA, String sequenceB) {
        cellMatrix[][] M = new cellMatrix[sequenceA.length() + 1][sequenceB.length() + 1];
        double highestScore;
        //initialization

        for (int i = 0; i <= sequenceA.length(); i++) {
            M[i][0].setScore(0);
        }
        for (int j = 0; j <= sequenceB.length(); j++) {
            M[0][j].setScore(0);
        }

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
        //TODO: 2 best with same score!
        double best = 0;
        cellMatrix bestCeller = null;
        double secondBest = 0;
        cellMatrix secondBestCeller = null;
        for (int i = 1; i <= sequenceA.length(); i++) {
            for (int j = 1; j <= sequenceB.length(); j++) {
                //second best score so far
                if (M[i][j].getScore() > secondBest) {
                    secondBest = M[i][j].getScore();
                    secondBestCeller = M[i][j];
                }
                //best so far, swap best and second best
                if (M[i][j].getScore() > best) {
                    secondBest = best;
                    secondBestCeller = bestCeller;
                    best = M[i][j].getScore();
                    bestCeller = M[i][j];
                }
            }
        }
        //print best path
        if (bestCeller != null) {
            cellMatrix currCell = bestCeller;
            while (currCell != null) {
                //current pi is M[i-1][j-1]
                //if(currCell.getPI() == M)
            }
        }
    }
}



