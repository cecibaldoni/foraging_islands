# Shrew Foraging
## Methods

Shrews for the behavioural experiments were captured in Möggingen, Germany (47°46'04.70"N, 8°59'47.11"E).
The individuals were sampled in three different seasons among two years: 26 individuals in spring 2020, 35 in summer 2021 and 44 in winter 2020 for a total of 105 individuals.

Trapping was carried out with wooden live-traps baited with mealworms (*Tenebrio molitor*) and checked at <2 h intervals. 
Each individual will be housed in two connected cages (measures). The front cage is arranged with a running wheel, a pot and dishes with food and water. The back cage is composed of several layers of hay and soil to imitate their natural environment. 
Food and water will be placed in ceramic dishes that will be changed every day. Shrews are fed with a mix of mealworms and raw meat, the amount of which changes during the year according to their body size. 
Cages will be held in an outdoor aviary subject to natural ambient temperature but not affected by weather conditions (i.e. rain, snow, wind).  

The set-up consists of a squared arena (110x110 cm), base covered in sand and closed on top by a transparent plastic layer to avoid escapes of the animals. 
Inside the arena 4 hexagonal feeders are disposed symmetrically at the same distance to each corner of the apparatus. Each of the 4 feeders (named A, B, C, D) feature 6 openings covered by a little door. 
Two of the feeders are filled with one mealworm for each opening (feeders A & D) while the other two (feeders B & C) are empty.

The tested shrew is able to enter the arena through a plastic pipe connected to its enclosure. 
After the entrance in the arena the animal could not go back since the opening will be automatically closed by an Arduino gate activated after the passage. 
Subsequently, the shrew is free to wonder inside the arena, interact with the feeders and open the feeders doors. 

Each individual was tested in two trials, one per day, during two days (T1 and T2). Each trial was composed by two session (S1 and S2), so in total every shrew performed the experiment in four sessions among two different days (T1_S1, T1_S2 - first day and T2_S1, T2_S2 second day)  
Each shrew has been tested and recorded for the entire duration of the experiment in the arena. Each video can last from 5 to 30 minutes, because the experiment ended when the shrew stopped to move or to interact with any feeder.
Subsequently, the recordings have been tracked with the tracking software Trex. 
The results of the first tracking gave the coordinates of the position of the animal at each frame of the recording. Afterwards, the video have been visually analysed to assess the effective visit of the islands.
Every visit was manually registered with the letter of the island visited, the number corresponding to the door, presence or absence of food and if the door was opened towards left or right.

## Scripts in this repo

1. foraging.R enable to merge the tracking data from Trex and the island visit data manually recorded.
2. trajectory_debug.R solve some problem on the trajectory of the shrew movement.
3. Mismatch solving.R solve the mismatches between the tracking observation from Trex and the island visits manually registered.
