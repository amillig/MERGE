# MERGE

MERGE represents a method that combines direct coupling analysis and machine learning techniques to predict a protein's fitness from sequence. It requires a binary parameter file outputted by [plmc](https://github.com/debbiemarkslab/plmc/tree/master) and sequence-fitness pairs.

![image](https://github.com/amillig/MERGE/assets/58852023/c11ff5d1-0c75-4307-a930-8ceb1a2edea2)


# Usage
See 'example/example.py' 
  1) Generate sequence representations
  2) Construct model
  3) Scrape fitness landscape

# References
“Combining evolutionary probability and machine learning enables data-driven protein engineering with minimized experimental effort” by Alexander-Maurice Illig, Niklas E. Siedhoff, Mehdi D. Davari*, and Ulrich Schwaneberg*

# Author
MERGE was developed and written by Alexander-Maurice Illig at RWTH Aachen University.

# Credits
MERGE uses binary parameter files that are generated with [plmc](https://github.com/debbiemarkslab/plmc/tree/master) written by [John Ingraham](https://github.com/jingraham).
