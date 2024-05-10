import numpy as np
import pandas as pd

"""
Globally declared lists to store the optimal alignments
"""
rx_final=[]
ry_final=[]


def sw(x, y, match=2, mismatch=-1, gap=-3):
    nx = len(x)
    ny = len(y)
    """
    The Initialization of the matrix F and P is same as the Q1(a).py. The only difference is that in case
    of the local alignment every cell has minimum value of 0. And hence we avoid filling first row and
    first column by gap penalty.
    """
    # Initialization of the matrix.
    F = np.zeros((nx + 1, ny + 1))
    # Pointers to trace through an optimal alignment.
    P = np.zeros((nx + 1, ny + 1), dtype=int)
    P[:, 0] = 3
    P[0, :] = 4
    # Matrix filling.
    """
    variable to store the max_score and the optimal alignment positions
    """
    current_max_score=-1 
    pointers=[]

    """
    Now we start filling the matrix F which is the scoring matrix. Nested loops are used to fill the matrix
    F, with row being referred by the variable i and column by j.
    Logic Behind Matrix filling:-
    1.  We start filling our matrix from the cell: (1,1) which represents 2nd row and 2nd column. Now the 
        value at any cell of matrix F will either represent a match, mismatch, a gap in first sequence or a gap in second sequence. 
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
        of these three while taking into account that in case of local alignment the minimum value of any
        cell is 0. Hence we take the maximum of all of these with 0 as well.
    """

    """
    Now we need to fill the Matrix P as well. Since this is local alignment then there can be no gaps
    in between the alignments, hence here we only need to check for the match/mismatch. We apply the same
    logic used in Q1(a).py here too. 
    Now there can be many optimal alignment with the same score. Hence we need to traceback all of these.
    Hence we store the current max score and if there's any cell with the same max score then we store its
    position so that we can later on traceback it. 
    And if we find a cell with score greater than the max score then we clear the stored positions and store
    this score as the maximum score along with the positions.
    """
    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            diagonal_score=0
            vertical_score=0
            horizontal_score=0
            default_score=0
            if(x[i-1]==y[j-1]):
                diagonal_score=F[i-1][j-1]+match
            else:
                diagonal_score=F[i-1][j-1]+mismatch
            vertical_score=F[i-1][j]+gap
            horizontal_score=F[i][j-1]+gap
            max_score=max(diagonal_score,vertical_score,horizontal_score,default_score)
            F[i][j]=max_score
            if diagonal_score == max_score:
                P[i,j] = P[i,j] + 2
            if (max_score==current_max_score):
                newPositions=[i,j]
                pointers.append(newPositions)
            elif(max_score>current_max_score):
                current_max_score=max_score
                newPositions=[i,j]
                pointers.clear()
                pointers.append(newPositions)

    # Print scoring matrix using pandas DataFrame
    print("Scoring Matrix:")
    df = pd.DataFrame(F, index=['-'] + list(x), columns=['-'] + list(y))
    print(df)

    """
    Now we start tracing back the optimal alignments. Since in case of local alignment the optimal alignment
    can begin from anywhere in the matrix, hence we use the pointers list which stores the positions
    of those cell which have the maximum score. 
    After getting the cell with the maximum score from the pointers list we start iterating over that. It is
    a 2d list which stores the i,j in a list. Now we simply use a while loop to traceback the sequence to that
    cell which has value as 0. And if we have reached a cell with value as 0 then we have completed the 
    traceback for that path. And hence we store this sequence into the globally declared lists.
    """
    for pointer in pointers:
        i=pointer[0]
        j=pointer[1]
        x_character=''
        y_character=''
        rx = []
        ry = []
        while P[i,j] !=0:
            if P[i, j] in [2]:
                x_character=x[i-1]
                y_character=y[j-1]
                rx.append(x_character)
                ry.append(y_character)
                i=i-1
                j=j-1
        rx_final.append(rx)
        ry_final.append(ry)



"""
This is a function to calculate the score for the alignment. It uses a simple logic, it iterates over the
two sequences and finds whether it's a match,mismatch or a gap in any of the sequence, if it's match then
it adds the match score which is passed as an argument to this function and if it's a gap in any of the
sequence and it adds the gap score which again is passed as an argument to this function and finally if 
there's no gap or no match then it would be a mismatch, hence it adds the mismatch score to the final score.
And after iterating over the sequences it returns the score.
"""
def score(rx,ry,match,mismatch,gap):
    score=0
    for i in range(0,len(rx)):
        if(rx[i]==ry[i]):
            score=score+match
        elif rx[i]=='-' or ry[i]=='-':
             score=score+gap
        else:
             score=score+mismatch
    return score


seq1 = "GATGCGCAG"
seq2 = "GGCAGTA"
sw(seq1, seq2)


"""
This is a simple for loop to print all the optimal alignment along with the scores. 
This uses the list rx_final and ry_final which stores the optimal alignment as a 2D List.
And it prints them along with their score. 
To find the score we use the score function and pass the arguments appropriately.
"""
print("\nTotal Optimal Alignments: ",len(rx_final))
for i in range(0,len(rx_final)):
    rx=rx_final[i]
    ry=ry_final[i]
    rx = ''.join(rx)[::-1]
    ry = ''.join(ry)[::-1]
    print()
    print('\n'.join([rx, ry]))
    print("Score:",score(rx,ry,match=2,mismatch=-1,gap=-3))