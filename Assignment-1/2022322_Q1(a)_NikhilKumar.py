import numpy as np
import pandas as pd

def nw(x, y, match=2, mismatch=-3, gap=-1):
    nx = len(x)
    ny = len(y)
    
    # Initialization of the matrix.
    F = np.zeros((nx + 1, ny + 1))

    """
    I changed the np.linspace(0, -gap * nx, nx + 1) to np.linspace(0, gap * nx, nx + 1), so as to have gaps 
    for the first row and first column. Initailly it was having positive values due to the -ve sign, now
    as we know that we are given gap as -ve then in this case we don't have to change it, whereas initially 
    it was getting positive due to multiplication of two negative numbers.
    """
    F[:, 0] = np.linspace(0, gap * nx, nx + 1)
    F[0, :] = np.linspace(0, gap * ny, ny + 1)

    # Pointers to trace through an optimal alignment.
    """
    Initially this P matrix was getting initialized later on i.e - after filling of the F matrix but as we
    know that this matrix will be used for tracing back the alignments, then in this case we have to fill 
    this matrix at the same time when we are filling the matrix F. Hence we have to declare it here.
    """
    P = np.zeros((nx + 1, ny + 1), dtype=int)
    P[:, 0] = 3
    P[0, :] = 4

    # Matrix filling.
    """
    Now we start filling the matrix F which is the scoring matrix. Nested loops are used to fill the matrix
    F, with row being referred by the variable i and column by j.
    Logic Behind Matrix filling:-
    1.  We start filling our matrix from the cell: (1,1) which represents 2nd row and 2nd column, as we have
        alredy filled the 1st row and 1st column by the gap penalty. Now the value at any cell of matrix F
        will either represent a match, mismatch, a gap in first sequence or a gap in second sequence. 
    2.  Now if a cell represents a match then in that case we have to take the score from the diagonal cell
        and add the match score to it, and if it represents a mismatch then we have to add the mismatch
        score to it. Now the diagonal score can be accessed by the (i-1,j-1) which will represent the 
        previous row and previous column and now if it's a match which can be checked by comparing the two
        sequences using the same i and j variable and subtracting 1 from that(as python have 0 based 
        indexing), if they have same character at that positon then we add the match score and if they are 
        not equal then we add mismatch.
    3.  Now if a cell represents a gap in the second sequence then in that case we have to take the score
        from the upper cell which represents taking the character of the first sequence and add the gap
        score to it. Now the upper score can be accessed by the (i-1,j) which will represent the previous
        row and same column.
    4.  Now if a cell represents a gap in the first sequence then in that case we have to take the score
        from the left cell which represents taking the character of the second sequence and add the gap 
        score to it. Now the left score can be accessed by the (i,j-1) which will represent the same row and
        previous column.
    5.  After obtaining these three scores we simply have to fill our cell of matrix F(i,j) as the maximum
        of these three. 
    """

    """
    While we are filling the matrix F, at the same time we have to fill our Matrix P so that we can 
    traceback the alignments.
    Logic Behind filling Matrix P:-
    1. Initially the first row of Matrix P was filled with 4 and first column with 3. This indicates that
    when we are taking the value of left cell(represents a gap in the first sequence) then we are adding 4 
    and when we are taking the value from upper cell(represents a gap in the second sequence) then we are adding 3. This tells us that when it's a match/mismatch (diagonal cell) then the value has to be 2,
    and the later part of traceback proves this.
    2. Now when we have got the value for F[i,j] then we can fill our P[i,j] depending upon that value
    filled at cell(i,j) of F. since the value at F[i,j] can be from more than one cell, hence we have
    to fill our matrix P so that while backtracking we backtrack both the paths. And hence we add the 2,3,4
    to our matrix P, inside the if conditions.
    3. If the value filled at F[i,j] represents the score from the diagonal, then we add 2 to the P[i,j].
    4. If the value filled at F[i,j] represents the score from the left cell, then we add 3 to the P[i,j].
    5. If the value filled at F[i,j] represents the score from the upper cell, then we add 4 to the P[i,j].
    """
    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            diagonal_score=0           # This represents the match/mismatch score
            vertical_score=0           # This represent a gap in the second sequence, score from upper cell
            horizontal_score=0         # This represent a gap in the first sequence, score from left cell
            if(x[i-1]==y[j-1]):
                diagonal_score=F[i-1][j-1]+match
            else:
                diagonal_score=F[i-1][j-1]+mismatch
            vertical_score=F[i-1][j]+gap
            horizontal_score=F[i][j-1]+gap
            max_score=max(diagonal_score,vertical_score,horizontal_score)
            F[i][j]=max_score
            if diagonal_score == max_score:
                P[i,j] = P[i,j] + 2
            if vertical_score == max_score:
                P[i,j] = P[i,j] + 3
            if horizontal_score == max_score:
                P[i,j] = P[i,j] + 4
            
    # Print scoring matrix using pandas DataFrame
    print("Scoring Matrix:")
    df = pd.DataFrame(F, index=['-'] + list(x), columns=['-'] + list(y))
    print(df)
    
    # Trace through an optimal alignment.
    """
    To traceback through an optimal alignment we follow these steps. First of all we start from the 
    rightmost bottom cell which in the global alignment represents the maximum score.
    1. Now if the value of a cell is from [2,5,6,9] then this represents a match/mismatch between
    the two sequences. This is because a value of 2 directly indicates a match, a value of 5
    indicates a match/mismatch + gap in y sequence, i.e- 2+3=5, similarly we can consider all the 
    combination of these like a value of 6 indicates a match/mismatch + gap in x sequence, a value of 9
    indicates a match/mismatch + gap in x sequence + gap in y sequence,i.e- 2+3+4=9. Now in this case 
    we append both the x and y sequence characters. Also we have taken both the characters in answer
    hence we subtract 1 from both i and j.
    2. Now if the value of a cell is from [3,5,7,9] then this represents a gap in y sequence and and a 
    character in x sequence. This is because a value of 3 directly indicates a gap in y, a value of 5
    indicates a match/mismatch + gap in y sequence, i.e- 2+3=5, similarly we can consider all the 
    combination of these, like a value of 7 indicates a gap in y sequence + gap in x sequence, a value of 9
    indicates a match/mismatch + gap in x sequence + gap in y sequence,i.e- 2+3+4=9. Now in this case 
    we append the x character and keep y as a gap. Also we have taken only the x character in answer
    hence we subtract 1 from i and keep j as it is.
    3. Now if the value of a cell is from [4,6,7,9] then this represents a gap in x sequence and a 
    character in y sequence. This is because a value of 4 directly indicates a gap in x, a value of 6
    indicates a match/mismatch + gap in x sequence, i.e- 2+4=6, similarly we can consider all the 
    combination of these, like a value of 7 indicates a gap in y sequence + gap in x sequence, a value of 9
    indicates a match/mismatch + gap in x sequence + gap in y sequence,i.e- 2+3+4=9. Now in this case 
    we append the x character and keep y as a gap. Also we have taken only the y character in answer
    hence we subtract 1 from j and keep i as it is.
    4. After this we return our answer by reversing the orders of both the strings.
    """
    i = nx
    j = ny
    rx = []
    ry = []
    x_character=''         #to store the character of x sequence(passed as an argument to the function nw())
    y_character=''         #to store the character of y sequence(passed as an argument to the function nw())
    while i > 0 or j > 0:
        if P[i, j] in [2, 5, 6, 9]:
            #complete the code
            x_character=x[i-1]
            y_character=y[j-1]
            rx.append(x_character)
            ry.append(y_character)
            i=i-1
            j=j-1
        elif P[i, j] in [3, 5, 7, 9]:
            #complete the code
            x_character=x[i-1]
            y_character='-'
            rx.append(x_character)
            ry.append(y_character)
            i=i-1
        elif P[i, j] in [4, 6, 7, 9]:
            #complete the code
            x_character='-'
            y_character=y[j-1]
            rx.append(x_character)
            ry.append(y_character)
            j=j-1

    # Reverse the strings.
    """
    The Score for the optimal alignment in case of global alignment is given by the rightmost bottom cell.
    i.e- F[nx,ny].
    """
    print()
    rx = ''.join(rx)[::-1]
    ry = ''.join(ry)[::-1]
    return '\n'.join([rx, ry])+"\n"+"Score: "+str(F[nx,ny])

seq1 = "GATGCGCAG"
seq2 = "GGCAGTA"
print(nw(seq1, seq2))