# Shrew Foraging

This repo contains scripts and code associated with a behavioural experiment conducted in common shrews (*Sorex araneus*) between 20xx-20-xx.

The set-up consists of a squared arena (110x110 cm), base covered in sand and closed on top by a transparent plastic layer. 
Inside the arena, 4 hexagonal feeders are disposed symmetrically at the same distance to each corner of the apparatus. Each of the 4 feeders (named A, B, C, D) features 6 openings covered by a little sliding door. 

*Maybe add a picture*

During the experiment, two of the feeders (A and D) are filled with one mealworm for each opening (6 in total), while the other two (feeders B & C) are empty.
The experiment consists of 4 sessions over two consecutive days. Each trial was composed by two session (S1 and S2), so in total every shrew performed the experiment in four sessions among two different days (T1_S1, T1_S2 - first day and T2_S1, T2_S2 second day
Each session lasts 30 minutes, but the shrew is free to roam inside the arena, interact with the feeders and open the feeders' doors. 

During the first session of the first day (T1S1), the feeders A and D are filled with food (12 mealworms). The shrew can explore all the feeders and get any of the food items. After 30 minutes, the shrew can leave the arena.
During the second session of the first day (T1S2), the feeders are *not* replenished and the arena is in the same "conditions" the shrew left it in (e.g. trail marks, urine marks). The shrew must remember the location they already got food from, or keep exploring to find more.
The sessions on the second day are identical: in T2S1, the arena is clean, the feeders are filled and the shrew can explore/find food in the same spots as the first time. In T2S2, the only food available is the one not consumed during the previous trial.

All shrews were recorded in the arena with a security camera. The recordings have been tracked with the tracking software Trex. 
The results of the first tracking gave the coordinates of the position of the animal at each frame of the recording. Afterwards, each video has been visually analysed to assess the effective visit of the islands: every visit was manually registered with the letter of the island visited, the number corresponding to the door, the presence or absence of food and if the door was opened towards left or right.

## Scripts in this repo

1. `foraging.R` enable to merge the tracking data from Trex and the island visit data manually recorded.
2. `trajectory_debug.R` solves some problem on the trajectory of the shrew movement.
3. `mismatch_solving.R` solves the mismatches between the tracking observation from Trex and the island visits manually registered.
