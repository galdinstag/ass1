import java.io.BufferedReader;
import java.io.InputStreamReader;

/**
 * Created by gal on 3/27/2016.
 */
public class main {
    public static void main(String [] args) throws Exception
    {
       SequenceAlignment seq = new SequenceAlignment(args[0],args[1]);
        seq.gapGlobalAlignment("AAACCAA","GGCCGGG");
        //seq.localAlignment("AAAAATGCTGCAAAAA","TTTTTTTTTTTTTTTTGCTGCTTT");

        //seq.globalAlignment(seq.sequences.get(0),seq.sequences.get(1));

    }
}
