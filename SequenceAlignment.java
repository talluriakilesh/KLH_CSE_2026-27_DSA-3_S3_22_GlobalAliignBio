import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class SequenceAlignment {

    // Scoring parameters
    private static final int MATCH = 2;
    private static final int MISMATCH = -1;
    private static final int GAP = -2;

    public static void align(String seq1, String seq2) {
        int n = seq1.length();
        int m = seq2.length();

        // 1. Initialize the Scoring Matrix (DP Table)
        int[][] dp = new int[n + 1][m + 1];

        for (int i = 0; i <= n; i++) {
            dp[i][0] = i * GAP;
        }
        for (int j = 0; j <= m; j++) {
            dp[0][j] = j * GAP;
        }

        // 2. Fill the DP Table
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= m; j++) {
                char c1 = seq1.charAt(i - 1);
                char c2 = seq2.charAt(j - 1);

                int scoreMatchMismatch = (c1 == c2) ? MATCH : MISMATCH;
                
                int diagonal = dp[i - 1][j - 1] + scoreMatchMismatch;
                int up = dp[i - 1][j] + GAP;
                int left = dp[i][j - 1] + GAP;

                dp[i][j] = Math.max(diagonal, Math.max(up, left));
            }
        }

        // 3. Traceback to find the optimal alignment
        int i = n;
        int j = m;
        List<Character> alignedSeq1 = new ArrayList<>();
        List<Character> alignedSeq2 = new ArrayList<>();

        while (i > 0 || j > 0) {
            char c1 = (i > 0) ? seq1.charAt(i - 1) : '\0';
            char c2 = (j > 0) ? seq2.charAt(j - 1) : '\0';

            int scoreMatchMismatch = (i > 0 && j > 0 && c1 == c2) ? MATCH : MISMATCH;

            if (i > 0 && j > 0 && dp[i][j] == dp[i - 1][j - 1] + scoreMatchMismatch) {
                alignedSeq1.add(c1);
                alignedSeq2.add(c2);
                i--;
                j--;
            } else if (i > 0 && dp[i][j] == dp[i - 1][j] + GAP) {
                alignedSeq1.add(c1);
                alignedSeq2.add('-');
                i--;
            } else {
                alignedSeq1.add('-');
                alignedSeq2.add(c2);
                j--;
            }
        }

        // Reverse the alignments since we traced back from end to start
        Collections.reverse(alignedSeq1);
        Collections.reverse(alignedSeq2);

        // 4. Output Results
        System.out.println("--- Needleman-Wunsch Alignment Results ---");
        System.out.println("Input Sequence 1: " + seq1);
        System.out.println("Input Sequence 2: " + seq2);
        System.out.println("Optimal Alignment Score: " + dp[n][m]);
        System.out.print("Aligned Sequence 1: ");
        for (char c : alignedSeq1) System.out.print(c + " ");
        System.out.println();
        System.out.print("Aligned Sequence 2: ");
        for (char c : alignedSeq2) System.out.print(c + " ");
        System.out.println("\n------------------------------------------");
    }

    public static void main(String[] args) {
        // Example DNA sequences
        String seq1 = "AGACTABDCBAFG";
        String seq2 = "GACTABDHEFG";

        align(seq1, seq2);
    }
}