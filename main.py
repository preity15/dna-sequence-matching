def needleman_wunsch_score(seq1, seq2):
    n = len(seq1)
    m = len(seq2)

    match = 1
    mismatch = -1
    gap = -2
    #create DP table
    dp =[[0 for _ in range(m + 1)] for _ in range(n + 1)]

    #initialiazing
    for i in range(n + 1):
        dp[i][0] = gap * i
    for j in range(m + 1):
        dp[0][j] = gap 

    #filling the table
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i - 1] == seq2[j - 1]:
              score = match 
            else:
                score = mismatch
            dp[i][j] = max(
                dp[i - 1][j - 1] + score,  # diagonal
                dp[i - 1][j] + gap,          # up
                dp[i][j - 1] + gap           # left
            )
    #traceback to get the score
    aligned1 = ""
    aligned2 = ""
    i, j = n, m

    while i > 0 and j > 0:
        current = dp[i][j]
        diagonal = dp[i - 1][j - 1]
        up = dp[i - 1][j]
        left = dp[i][j - 1]

        if seq1[i - 1] == seq2[j - 1]:
            score = match
        else:
            score = mismatch
        
        if current == diagonal + score:
            aligned1 = seq1[i - 1] + aligned1
            aligned2 = seq2[j - 1] + aligned2
            i -= 1
            j -= 1
        elif current == up + gap:
            aligned1 = '-' + aligned1
            aligned2 = seq2[j - 1] + aligned2
            i -= 1
        else:
            aligned1 = "-" + aligned1
            aligned2 = seq2[j-1] + aligned2
            j -= 1

    # Add any remaining characters from either sequence
    # remaining
    while i > 0:
        aligned1 = seq1[i-1] + aligned1
        aligned2 = "-" + aligned2
        i -= 1

    while j > 0:
        aligned1 = "-" + aligned1
        aligned2 = seq2[j-1] + aligned2
        j -= 1

    return dp[n][m], aligned1, aligned2


seq1 = "AGTACG"
seq2 = "ACATAG"

score, a1, a2 = needleman_wunsch_score(seq1, seq2)

print("Score:", score)
print("Alignment 1:", a1)
print("Alignment 2:", a2)