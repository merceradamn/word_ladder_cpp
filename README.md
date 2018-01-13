# Word Ladder in C++

The Word Ladder Problem: Given a start word and an end word, find a shortest sequence of words (shortest word ladder) from the start to the end so that:
  1. Each word is a legal word from a given dictionary
  2. The next word is obtained from the previous word by substituting
				exactly one letter.
			
If there is no such sequence, then say so. Note: Assume all words are in lowercase.
Example: Given a start word zero and an end word five, here is the ladder:

  > zero => hero => here => hire => fire => five
  
  We assume that all the words are legal.

This program has two parts. The first finds the unweighted word ladder from the starting word to the end word. Then the second part finds the weighted word ladder using the accompanied letter weight file.
