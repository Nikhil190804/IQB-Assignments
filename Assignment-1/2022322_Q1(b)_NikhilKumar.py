import numpy as np
import pandas as pd

"""
Globally declared lists to store the optimal alignments
"""
rx_final = []
ry_final = []

"""
This function finds all the optimal alignments using recursion.
The main idea of tracing back the alignment is same as in the Q1(a).py.
The difference is that when there's a cell whose value has came from more than one cell, then we need to 
traceback both the paths here. So to achieve this we need to check every condition hence we use if blocks
rather than if,elif. 
Now since this is recursion hence we need to have a Base Case, here the base case occurs when the value
of both i and j have became less than or equal to zero, which represents that we have iterated over both the
x and y sequence (at matrix level). Hence now we just store the optimal alignment in a list which is declared
globally and then we return from the function. Also in the base case we create new copies of the lists so
as to ensure that list are newly created and not being passed by refrence. Also after the function call 
returns we need to pop the last element from both the rx and ry, since these are passed by refrence in python.
"""
def traceback(i, j, P, rx, ry, x, y):
    if (i <= 0 and j <= 0):
        newx = list.copy(rx)
        newy = list.copy(ry)
        rx_final.append(newx)
        ry_final.append(newy)
        return
    if P[i, j] in [2, 5, 6, 9]:
        # complete the code
        x_character = x[i-1]
        y_character = y[j-1]
        rx.append(x_character)
        ry.append(y_character)
        traceback(i-1, j-1, P, rx, ry, x, y)
        rx.pop()
        ry.pop()
    if P[i, j] in [3, 5, 7, 9]:
        # complete the code
        x_character = x[i-1]
        y_character = '-'
        rx.append(x_character)
        ry.append(y_character)
        traceback(i-1, j, P, rx, ry, x, y)
        rx.pop()
        ry.pop()
    if P[i, j] in [4, 6, 7, 9]:
        # complete the code
        x_character = '-'
        y_character = y[j-1]
        rx.append(x_character)
        ry.append(y_character)
        traceback(i, j-1, P, rx, ry, x, y)
        rx.pop()
        ry.pop()


"""
This is a function to calculate the score for the alignment. It uses a simple logic, it iterates over the
two sequences and finds whether it's a match,mismatch or a gap in any of the sequence, if it's match then
it adds the match score which is passed as an argument to this function and if it's a gap in any of the
sequence and it adds the gap score which again is passed as an argument to this function and finally if 
there's no gap or no match then it would be a mismatch, hence it adds the mismatch score to the final score.
And after iterating over the sequences it returns the score.
"""
def score(rx, ry,match,mismatch,gap):
    score = 0
    for i in range(0, len(rx)):
        if (rx[i] == ry[i]):
            score = score+ match
        elif rx[i] == '-' or ry[i] == '-':
            score =score+ gap
        else:
            score =score+ mismatch
    return score


"""
This function is same as the function nw() in Q1(a).py. This simply fills the scoring matrix using DP.
The logic behind this is also same. The only difference here is that after filling the scoring matrix
we call the traceback() function to trace all the optimal alignments.
"""
def nw(x, y, match=2, mismatch=-3, gap=-1):
    nx = len(x)
    ny = len(y)

    # Initialization of the matrix.
    F = np.zeros((nx + 1, ny + 1))
    F[:, 0] = np.linspace(0, gap * nx, nx + 1)
    F[0, :] = np.linspace(0, gap * ny, ny + 1)

    # Pointers to trace through an optimal alignment.
    P = np.zeros((nx + 1, ny + 1), dtype=int)
    P[:, 0] = 3
    P[0, :] = 4
    # Matrix filling.
    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            diagonal_score = 0
            vertical_score = 0
            horizontal_score = 0
            if (x[i-1] == y[j-1]):
                diagonal_score = F[i-1][j-1]+match
            else:
                diagonal_score = F[i-1][j-1]+mismatch
            vertical_score = F[i-1][j]+gap
            horizontal_score = F[i][j-1]+gap
            max_score = max(diagonal_score, vertical_score, horizontal_score)
            F[i][j] = max_score
            if diagonal_score == max_score:
                P[i, j] = P[i, j] + 2
            if vertical_score == max_score:
                P[i, j] = P[i, j] + 3
            if horizontal_score == max_score:
                P[i, j] = P[i, j] + 4

    # Print scoring matrix using pandas DataFrame
    print("Scoring Matrix:")
    df = pd.DataFrame(F, index=['-'] + list(x), columns=['-'] + list(y))
    print(df)
    
    # Trace through an optimal alignment.
    traceback(nx, ny, P, [], [], x, y)

seq1 = "GATGCGCAG"
seq2 = "GGCAGTA"
nw(seq1, seq2)


"""
This is a simple for loop to print all the optimal alignment along with the scores. 
This uses the list rx_final and ry_final which stores the optimal alignment as a 2D List.
And it prints them along with their score. 
The score would be same for all and it can be accessed by the rightmost bottom cell of the scoring matrix
but here i have wrote a function to find the score as well.
"""
print("\nTotal Optimal Alignments: ",len(rx_final))
for i in range(0, len(rx_final)):
    rx = rx_final[i]
    ry = ry_final[i]
    rx = ''.join(rx)[::-1]
    ry = ''.join(ry)[::-1]
    print()
    print('\n'.join([rx, ry]))
    print("Score: ", score(rx, ry,match=2,mismatch=-3,gap=-1))