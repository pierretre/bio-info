package bioinfo;

public class AlignmentNW {
    private String s1;
    private String s2;
    private int scoreMatch = 5;
    private int scoreMismatch = -4;
    private int scoreIndel = -3;
    private int[][] alignmentMatrix;

    public AlignmentNW(String s1, String s2) {
        this.s1 = s1;
        this.s2 = s2;
        this.alignmentMatrix = new int[s2.length()+1][s1.length()+1];
        
        // Fill in the matrix :
        for (int j = 1; j < this.alignmentMatrix[0].length; j++) 
            this.alignmentMatrix[0][j] = this.alignmentMatrix[0][j-1] + this.scoreIndel;

        for (int i = 1; i < this.alignmentMatrix.length; i++) {
            this.alignmentMatrix[i][0] = this.alignmentMatrix[i-1][0] + this.scoreIndel;

            for (int j = 1; j < alignmentMatrix[i].length; j++) {
                int top_left_diag = this.alignmentMatrix[i-1][j-1] + ((s1.charAt(j-1) == s2.charAt(i-1))? this.scoreMatch : this.scoreMismatch);
                int top = this.alignmentMatrix[i-1][j] + this.scoreIndel;
                int left = this.alignmentMatrix[i][j-1] + this.scoreIndel;
                this.alignmentMatrix[i][j] = Math.max(top_left_diag, Math.max(top, left));
            }
        }
    }

    public void printMatrix() {
        
        System.out.println();
        System.out.println("Alignment Matrix = ");
        System.out.println(" \t"+s1.replace("", "\t"));
        for (int i = 0; i < alignmentMatrix.length; i++) {
            if(i > 0) System.out.print(s2.charAt(i-1)+"\t");
            else System.out.print(" \t");
            for (int j = 0; j < alignmentMatrix[i].length; j++) {
                System.out.print(alignmentMatrix[i][j]+"\t");
            }
            System.out.println();
        }
    }

    private int getScoreMin(){
        return this.scoreIndel * (s1.length() + s2.length());
    }
    
    private int getScoreMax(){
        return this.scoreMatch * Math.min(s1.length(), s2.length());        
    }

    public double getDistance(){
        return (double)(this.getScoreMax() - this.alignmentMatrix[s2.length()][s1.length()]) / (this.getScoreMax() - this.getScoreMin());
    }

    public static void main(String[] args) {
		String seq1 = "ATTACG";
		String seq2 = "ATATCG";
		String seq3 = "ACCCCG";
		String seq4 = "GGGGAA";
		String seq5 = "TTTACG";

        AlignmentNW a1 = new AlignmentNW(seq1, seq2);
        a1.printMatrix();
        System.out.println("distance = "+a1.getDistance());

        AlignmentNW a2 = new AlignmentNW(seq1, seq3);
        a2.printMatrix();
        System.out.println("distance = "+a2.getDistance());

        AlignmentNW a3 = new AlignmentNW(seq1, seq3);
        a3.printMatrix();
        System.out.println("distance = "+a3.getDistance());

        AlignmentNW a4 = new AlignmentNW(seq1, "eeee");
        a4.printMatrix();
        System.out.println("distance = "+a4.getDistance());

        AlignmentNW a5 = new AlignmentNW(seq1, seq1);
        a5.printMatrix();
        System.out.println("distance = "+a5.getDistance());
	}
}
