/**
 * Created by gal on 3/28/2016.
 */
public class cellMatrix {
    private cellMatrix pi;
    private int score;

    public cellMatrix(){
        pi = null;
        score = 0;
    }

    public void setScore(int score) {
        this.score = score;
    }

    public int getScore() {
        return score;
    }

    public void setPi(cellMatrix pi) {
        this.pi = pi;
    }

    public cellMatrix getPI() {
        return pi;
    }
}
