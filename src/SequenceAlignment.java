import java.io.*;
import java.nio.file.Files;
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
    }

    public void GlobalAlignment(){


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
}



