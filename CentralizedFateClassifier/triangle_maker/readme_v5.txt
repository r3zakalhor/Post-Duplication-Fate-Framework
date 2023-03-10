The v5 refers to the prob calculation version (the current version).  The one we used in Nov/Dec was v4.  
csv_1_v5, csv_2_v5, csv_3_v5 contain the triangles from the runs Celine did with the large triangles.
testset_v5 contains triangles from a set of runs in November.  The script avoids repetition, that is, sets of m,w,h already encountered are ignored, so there may be gaps in the numbers.  
other/testset_v4 contains the same triangles, but with the v4 calculations 
other/testset_v5_different contains only the triangles that do not have the same probs as v5
other/testset_v5_with_repeats contains the same triangles, but includes repeated triangles.

Celine had comments on testset_v5_with_repeats, as follows.  
- what happens here 67 ?
- 103 (strange fate) 104 (why neo? This happens a lot in this example ) 125 (a bit sub?) 227 (cons is too low wrt 231 ?) 316 (really neo?) 418 (should we remove these g?) 443 (why a bit of spec?)
Some have changed, some haven't.  IN the cases that haven't changed, I'm fine with the current probs.





