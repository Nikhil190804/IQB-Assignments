-----------------------------:The Main Idea of Algorithm used:-------------------------------

1. first we scan for a valid window of appropriate size (6 for helix and 5 for strands)
2. After obtaining a valid window we store its start and end position in a list.
3. Now for every valid window stored in that list we extend this to both left and right.
4. After extension of that window we update our start and end postions respectively.
5. Now after getting the start and end postions for a secondary structure we now make a new list which
is filled of 0 or 1 or 2. 
6. If the ith position has 0 then it means that ith amino acid in original protein sequence has been assigned 
to no secondary structure. If it is 1 then it means it is assigned to helix and if its 2 then it is assigned to strand.
7. Now we simply iterate over this zero_one list and use appropriate if and else statemnets to print the 
desired secondary structure.
8. Now for resolving conflicts if the ith postion has 1 and 2 then this represents a conflicting region and 
hence we solve this considering which score is greater. And assign a final number to that ith postion.
9. Finally we print the desired output as done in step 7.



-----------------------------:NOTE:----------------------------------------
1. While calculating the score and adding them i have made them round off to 2 decimal places after 0, as 
with floating point numbers the comparison is not that much precise.
2. I have printed 80 characters in one line i.e - 80 amino acids will be printed in one line and below 
them would be there respectively secondary structures.
To change this you can edit the following function calls:-

Line number : 411 printAnswer(helix_zero_one_list,1,80)
Line number : 419 printAnswer(sheet_zero_one_list,2,80)
Line number : 426 displayConflicts(helix_zero_one_list,sheet_zero_one_list,80)
Line number : 434 resolveConflicts(helix_zero_one_list,sheet_zero_one_list,80) 

Here the last argument passed i.e- 80 represents that 80 characters would be printed in one line