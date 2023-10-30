# metaSPPstu: solving the Set Packing Problem - ElizabethG - AwenJ
Implementation in Julia (compliant Julia v1.9) of tools related to the Set Packing Problem (SPP) for pedagogical purposes.

This implementation is the base of an exercice of the course "metaheuristics".

------

Fonctionnalités des codes :

Eléments de soutien à l'implémentation en Julia v1.x des devoirs maison à réaliser dans le cadre du cours "métaheuristiques" en master 1 informatique". Révision pour l'année académique 2023-2024.

- `loadSPP.jl`       : lecture d'une instance de SPP au format OR-library
- `setSPP.jl`        : construction d'un modèle JuMP de SPP
- `getfname.jl`      : collecte les noms de fichiers non cachés présents dans un répertoire donné
- `experiment.jl`    : protocole pour mener une expérimentation numérique avec sorties graphiques
- `amelioration.jl`  : contient la fonction du glouton amelioration, celui de la plus profonde descente
- `construction.jl`  : contient la fonction du glouton contruction
- `destruction.jl`   : contient la fonction qui sert à détruire notre solution, elle est utilisée dans notre destroy and repair
- `grasp.jl`         : contient la fonction grasp et la fonction grasp_DR, qui est notre solution opérationnelle: le destroy and repair
- `reconstruction.jl`: contient la fonction qui reconstruite notre solution à partir de celle détruite. Elle se base sur l'algorithme de VND.
------

Le répertoire `Data` contient une sélection d'instances numériques de SPP au format OR-library :

- didactic.dat
- pb_100rnd0100.dat
- pb_100rnd0200.dat
- pb_100rnd0300.dat
- pb_200rnd0100.dat
- pb_200rnd0200.dat 
- pb_500rnd0100.dat
- pb_500rnd0300.dat
- pb_500rnd0500.dat
- pb_1000rnd0100.dat
- pb_2000rnd0100.dat

------

Exemple d'utilisation (`main.jl`) avec chemins d'accès correspondant à une configuration standard sur macOS :
- chargement de l'instance `didactic.dat` de SPP
- résolution d'une instance de SPP avec la solution de contruction gloutonne, puis d'amélioration gloutonne, de Grasp et enfin de notre solution opérationnelle: destroy and repair.
- collecte des noms d'instances présentes dans le répertoire `Data`
