import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.text.SimpleDateFormat;
import java.util.Date;

/**
 * Created by gal on 3/27/2016.
 */
public class Alignments {
    public static void main(String [] args) throws Exception
    {
        //parse score matrix and fasta files
        SequenceAlignment seq = new SequenceAlignment(args[1],args[2],args[3]);
        String timeStamp = new SimpleDateFormat("mm.ss").format(new Date());

        //compute all-against-all alignment according to the flag
        for(String sequenceA : seq.sequences1.keySet()){
            for(String sequenceB : seq.sequences2.keySet()){
                if (sequenceA.compareTo(sequenceB) != 0) {
                    System.out.println(sequenceA);
                    System.out.println(sequenceB);
                    switch (args[0]) {
                        case "-g":
                            seq.globalAlignment(seq.sequences1.get(sequenceA), seq.sequences2.get(sequenceB));
                            break;
                        case "-l":
                            seq.localAlignment(seq.sequences1.get(sequenceA), seq.sequences2.get(sequenceB));
                            break;
                        case "-a":
                            seq.gapGlobalAlignment(seq.sequences1.get(sequenceA), seq.sequences2.get(sequenceB));
                            break;
                    }

                }
            }
        }
        String timeStampEnd = new SimpleDateFormat("mm.ss").format(new Date());
        int minuts = (new Integer(timeStampEnd.substring(0,2)).intValue() - new Integer(timeStamp.substring(0,2)).intValue());
        int seconds = (new Integer(timeStampEnd.substring(3,5)).intValue() - new Integer(timeStamp.substring(3,5)).intValue());
        System.out.println("total time is: " + minuts + " minutes and " + seconds + " seconds");

    }
}
